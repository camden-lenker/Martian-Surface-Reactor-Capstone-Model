"""
Martian surface fission reactor -- OpenMC model builder (v0)

All materials, geometry, and simulation settings live in reactor_config.json.
This script only reads that file and turns it into an OpenMC model.

Usage (from your openmc-env, inside WSL):
    python reactor_model.py --mode plot              # geometry PNGs
    python reactor_model.py --mode run               # one k-eigenvalue run
    python reactor_model.py --mode sweep             # control drum worth
    python reactor_model.py --mode kinf              # no-leakage k diagnostic

    # any config value can be overridden from the command line:
    python reactor_model.py --mode run \
        --set geometry.lattice.n_rings=10 \
        --set geometry.control_drums.angle_deg=0

    python reactor_model.py --config other_design.json --mode plot

Requires OpenMC >= 0.14 and an ENDF/B-VIII.0 HDF5 library
(OPENMC_CROSS_SECTIONS must be set).

Units: cm, g/cm3, K.

Known simplifications (see also the _notes in the JSON):
  * Gaps and pin interstitials are void; heat pipes stop at the core
  * Drums are BeO cylinders with a B4C arc (no sleeve, bearings, actuators)
  * Temperatures are fixed inputs, not coupled to a thermal-hydraulics solution (Will change in future with T-H Calculations)
  * Shield is a flat disk stack, not a shadow cone sized to a dose limit
"""
import argparse
import json
import math
import os
from copy import deepcopy

import numpy as np
import openmc

#Set the path to read the config.json file
DEFAULT_CONFIG = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "reactor_config.json")


# --------------------------------------------------------------------------
# 1. CONFIG LOADING, OVERRIDES, VALIDATION
# --------------------------------------------------------------------------
def load_config(path):
    with open(path) as f:
        return json.load(f)


def apply_overrides(cfg, overrides):
    """Apply 'a.b.c=value' strings; value is parsed as JSON if possible."""
    cfg = deepcopy(cfg)
    for item in overrides or []:
        if "=" not in item:
            raise ValueError(f"--set expects path=value, got '{item}'")
        path, raw = item.split("=", 1)
        try:
            value = json.loads(raw)
        except json.JSONDecodeError:
            value = raw
        node = cfg
        keys = path.split(".")
        for k in keys[:-1]:
            if k not in node:
                raise KeyError(f"'{path}': no key '{k}' in config")
            node = node[k]
        if keys[-1] not in node:
            raise KeyError(f"'{path}': no key '{keys[-1]}' in config")
        node[keys[-1]] = value
    return cfg


def validate(cfg):
    g, mats = cfg["geometry"], cfg["materials"]

    # every referenced material must exist
    refs = [g["fuel_pin"]["fuel_material"], g["fuel_pin"]["clad_material"],
            g["heat_pipe"]["wall_material"], g["heat_pipe"]["wick_material"],
            g["lattice"]["filler_material"], g["core"]["barrel_material"],
            g["reflector"]["material"], g["control_drums"]["body_material"],
            g["control_drums"]["absorber_material"]]
    if g["lattice"].get("interstitial_material"):
        refs.append(g["lattice"]["interstitial_material"])
    if g["shield"]["enabled"]:
        refs += [layer["material"] for layer in g["shield"]["layers"]]
    for key in refs:
        if key not in mats:
            raise KeyError(f"Material '{key}' is referenced in geometry but "
                           f"not defined under 'materials'")

    # pin and heat pipe radii must nest
    fp, hp = g["fuel_pin"], g["heat_pipe"]
    if not fp["r_fuel"] + fp["gap"] < fp["r_clad"]:
        raise ValueError("fuel_pin: need r_fuel + gap < r_clad")
    if not hp["r_vapor"] < hp["r_wick"] < hp["r_outer"]:
        raise ValueError("heat_pipe: need r_vapor < r_wick < r_outer")
    if g["lattice"]["pitch"] < 2 * max(fp["r_clad"], hp["r_outer"]):
        raise ValueError("lattice.pitch must be >= 2 x pin outer radius "
                         "(pins would overlap)")

    # drums must not overlap each other
    dims = _radial_dims(cfg)
    d = g["control_drums"]
    chord = 2 * dims["R_drum_c"] * math.sin(math.pi / d["count"])
    if chord < 2 * d["radius"]:
        raise ValueError(
            f"Control drums overlap: center spacing {chord:.2f} cm < "
            f"diameter {2 * d['radius']:.2f} cm. Reduce count/radius.")
    if not d["absorber_thickness"] < d["radius"]:
        raise ValueError("drum absorber_thickness must be < drum radius")
    if not 0 < d["absorber_half_angle_deg"] < 90:
        raise ValueError("absorber_half_angle_deg must be between 0 and 90")


def _radial_dims(cfg):
    g = cfg["geometry"]
    R_core = ((g["lattice"]["n_rings"] - 1) * g["lattice"]["pitch"]
              + g["fuel_pin"]["r_clad"] + g["core"]["margin_to_barrel"])
    R_barrel = R_core + g["core"]["barrel_thickness"]
    d = g["control_drums"]
    R_drum_c = R_barrel + d["gap_to_barrel"] + d["radius"]
    R_ref = R_drum_c + d["radius"] + g["reflector"]["outboard_margin"]
    return dict(R_core=R_core, R_barrel=R_barrel,
                R_drum_c=R_drum_c, R_ref=R_ref)


# --------------------------------------------------------------------------
# 2. MATERIALS (built lazily from the JSON, only if referenced)
# --------------------------------------------------------------------------
class MaterialLibrary:
    def __init__(self, material_specs):
        self.specs = material_specs
        self.cache = {}

    def get(self, key):
        if key not in self.cache:
            self.cache[key] = self._build(key, self.specs[key])
        return self.cache[key]

    @staticmethod
    def _build(key, spec):
        mat = openmc.Material(name=spec.get("name", key))
        for comp in spec["components"]:
            ptype = comp.get("percent_type", "ao")
            if "nuclide" in comp:
                mat.add_nuclide(comp["nuclide"], comp["percent"], ptype)
            else:
                kw = {k: comp[k] for k in
                      ("enrichment", "enrichment_target", "enrichment_type")
                      if k in comp}
                mat.add_element(comp["element"], comp["percent"],
                                percent_type=ptype, **kw)
        density = spec["density_g_cm3"] * spec.get("fraction_of_theoretical",
                                                   1.0)
        mat.set_density("g/cm3", density)
        for sab in spec.get("s_alpha_beta", []):
            mat.add_s_alpha_beta(sab)
        return mat

    def colors(self, cfg_materials):
        return {mat: cfg_materials[key].get("plot_color", "gray")
                for key, mat in self.cache.items()}


# --------------------------------------------------------------------------
# 3. PIN / HEAT PIPE UNIVERSES AND THE HEX LATTICE
# --------------------------------------------------------------------------
def _interstitial(g, lib):
    """Material between pins (None = void) and its temperature."""
    lat = g["lattice"]
    key = lat.get("interstitial_material")
    if not key:
        return None, None
    return lib.get(key), lat["interstitial_temperature_K"]


def fuel_pin_universe(g, lib):
    fp = g["fuel_pin"]
    s_fuel = openmc.ZCylinder(r=fp["r_fuel"])
    s_gap = openmc.ZCylinder(r=fp["r_fuel"] + fp["gap"])
    s_clad = openmc.ZCylinder(r=fp["r_clad"])

    c_fuel = openmc.Cell(name="fuel", fill=lib.get(fp["fuel_material"]),
                         region=-s_fuel)
    c_fuel.temperature = fp["fuel_temperature_K"]
    c_gap = openmc.Cell(name="gap (void)", region=+s_fuel & -s_gap)
    c_clad = openmc.Cell(name="pin clad", fill=lib.get(fp["clad_material"]),
                         region=+s_gap & -s_clad)
    c_clad.temperature = fp["clad_temperature_K"]
    inter, inter_T = _interstitial(g, lib)
    c_out = openmc.Cell(name="interstitial", fill=inter, region=+s_clad)
    if inter is not None:
        c_out.temperature = inter_T
    return openmc.Universe(name="fuel pin",
                           cells=[c_fuel, c_gap, c_clad, c_out])


def heat_pipe_universe(g, lib):
    hp = g["heat_pipe"]
    s_vap = openmc.ZCylinder(r=hp["r_vapor"])
    s_wick = openmc.ZCylinder(r=hp["r_wick"])
    s_wall = openmc.ZCylinder(r=hp["r_outer"])

    c_vap = openmc.Cell(name="HP vapor (void)", region=-s_vap)
    c_wick = openmc.Cell(name="HP wick", fill=lib.get(hp["wick_material"]),
                         region=+s_vap & -s_wick)
    c_wall = openmc.Cell(name="HP wall", fill=lib.get(hp["wall_material"]),
                         region=+s_wick & -s_wall)
    inter, inter_T = _interstitial(g, lib)
    c_out = openmc.Cell(name="interstitial", fill=inter, region=+s_wall)
    if inter is not None:
        c_out.temperature = inter_T
    for c in (c_wick, c_wall):
        c.temperature = hp["temperature_K"]
    return openmc.Universe(name="heat pipe",
                           cells=[c_vap, c_wick, c_wall, c_out])


def build_lattice(g, lib, fuel_u, hp_u):
    lat_cfg = g["lattice"]
    filler = openmc.Cell(name="lattice filler",
                         fill=lib.get(lat_cfg["filler_material"]))
    filler.temperature = lat_cfg["filler_temperature_K"]
    filler_u = openmc.Universe(name="filler", cells=[filler])

    n, every = lat_cfg["n_rings"], lat_cfg["heat_pipe_every"]
    rings = []
    for k in range(n - 1, 0, -1):          # outermost ring first
        rings.append([hp_u if (i % every == 0) else fuel_u
                      for i in range(6 * k)])
    rings.append([fuel_u])                 # center pin

    lat = openmc.HexLattice(name="core lattice")
    lat.center = (0.0, 0.0)
    lat.pitch = (lat_cfg["pitch"],)
    lat.orientation = lat_cfg["orientation"]
    lat.universes = rings
    lat.outer = filler_u

    n_fuel = sum(1 for ring in rings for u in ring if u is fuel_u)
    n_hp = sum(1 for ring in rings for u in ring if u is hp_u)
    return lat, n_fuel, n_hp


# --------------------------------------------------------------------------
# 4. CONTROL DRUM HELPER
# --------------------------------------------------------------------------
def plane_through(center, theta_deg):
    """Vertical plane through `center` containing the ray at angle theta.
    The + side is the counter-clockwise side of that ray."""
    th = math.radians(theta_deg)
    a, b = -math.sin(th), math.cos(th)
    return openmc.Plane(a=a, b=b, c=0.0, d=a * center[0] + b * center[1])


# --------------------------------------------------------------------------
# 5. FULL MODEL
# --------------------------------------------------------------------------
def build_model(cfg):
    validate(cfg)
    g = cfg["geometry"]
    lib = MaterialLibrary(cfg["materials"])
    L = g["core"]["active_length"]
    dims = _radial_dims(cfg)
    R_core, R_barrel = dims["R_core"], dims["R_barrel"]
    R_drum_c, R_ref = dims["R_drum_c"], dims["R_ref"]

    # ---- lattice core ----
    fuel_u = fuel_pin_universe(g, lib)
    hp_u = heat_pipe_universe(g, lib)
    lattice, n_fuel, n_hp = build_lattice(g, lib, fuel_u, hp_u)

    # ---- axial planes ----
    z_fb, z_ft = -L / 2, L / 2
    z_rb = z_fb - g["reflector"]["axial_thickness"]
    z_rt = z_ft + g["reflector"]["axial_thickness"]
    shield_layers = g["shield"]["layers"] if g["shield"]["enabled"] else []

    s_fb = openmc.ZPlane(z0=z_fb)
    s_ft = openmc.ZPlane(z0=z_ft)
    s_rt = openmc.ZPlane(z0=z_rt, boundary_type="vacuum")
    s_rb = openmc.ZPlane(z0=z_rb,
                         boundary_type="transmission" if shield_layers
                         else "vacuum")

    # ---- radial cylinders ----
    s_core = openmc.ZCylinder(r=R_core)
    s_barrel = openmc.ZCylinder(r=R_barrel)
    s_ref = openmc.ZCylinder(r=R_ref, boundary_type="vacuum")

    cells = []

    # core
    core = openmc.Cell(name="core", fill=lattice,
                       region=-s_core & +s_fb & -s_ft)
    cells.append(core)

    # barrel
    barrel = openmc.Cell(name="barrel",
                         fill=lib.get(g["core"]["barrel_material"]),
                         region=+s_core & -s_barrel & +s_fb & -s_ft)
    barrel.temperature = g["core"]["barrel_temperature_K"]
    cells.append(barrel)

    # control drums: absorber arc defined in global coordinates with planes
    d = g["control_drums"]
    h = d["absorber_half_angle_deg"]
    drum_cyls = []
    for k in range(d["count"]):
        phi = 360.0 * k / d["count"]
        cx = R_drum_c * math.cos(math.radians(phi))
        cy = R_drum_c * math.sin(math.radians(phi))
        gamma = phi + 180.0 + d["angle_deg"]   # arc centre direction

        c_out = openmc.ZCylinder(x0=cx, y0=cy, r=d["radius"])
        c_in = openmc.ZCylinder(x0=cx, y0=cy,
                                r=d["radius"] - d["absorber_thickness"])
        lo = plane_through((cx, cy), gamma - h)
        hi = plane_through((cx, cy), gamma + h)
        drum_cyls.append(c_out)

        zr = +s_fb & -s_ft
        absorber = openmc.Cell(
            name=f"drum {k} absorber", fill=lib.get(d["absorber_material"]),
            region=-c_out & +c_in & +lo & -hi & zr)
        absorber.temperature = d["absorber_temperature_K"]
        body = openmc.Cell(
            name=f"drum {k} body", fill=lib.get(d["body_material"]),
            region=-c_out & (-c_in | -lo | +hi) & zr)
        body.temperature = d["body_temperature_K"]
        cells += [absorber, body]

    # radial reflector: annulus minus the drums
    refl_mat = lib.get(g["reflector"]["material"])
    refl_t = g["reflector"]["temperature_K"]
    region = +s_barrel & -s_ref & +s_fb & -s_ft
    for c_out in drum_cyls:
        region = region & +c_out
    radial = openmc.Cell(name="radial reflector", fill=refl_mat,
                         region=region)
    bot_refl = openmc.Cell(name="bottom axial reflector", fill=refl_mat,
                           region=-s_ref & +s_rb & -s_fb)
    top_refl = openmc.Cell(name="top axial reflector", fill=refl_mat,
                           region=-s_ref & +s_ft & -s_rt)
    for c in (radial, bot_refl, top_refl):
        c.temperature = refl_t
        cells.append(c)

    # shield layers stacked downward from the bottom reflector
    z_cursor, upper = z_rb, s_rb
    for i, layer in enumerate(shield_layers):
        z_next = z_cursor - layer["thickness"]
        last = i == len(shield_layers) - 1
        s_next = openmc.ZPlane(z0=z_next,
                               boundary_type="vacuum" if last
                               else "transmission")
        cell = openmc.Cell(name=layer["name"],
                           fill=lib.get(layer["material"]),
                           region=-s_ref & +s_next & -upper)
        cell.temperature = layer["temperature_K"]
        cells.append(cell)
        z_cursor, upper = z_next, s_next
    z_bot = z_cursor

    geometry = openmc.Geometry(openmc.Universe(cells=cells))

    # ---- fuel volume (needed for masses and, later, depletion) ----
    fuel_mat = lib.get(g["fuel_pin"]["fuel_material"])
    fuel_mat.volume = n_fuel * math.pi * g["fuel_pin"]["r_fuel"] ** 2 * L

    # ---- settings ----
    s = cfg["settings"]
    st = openmc.Settings()
    st.run_mode = s["run_mode"]
    st.batches = s["batches"]
    st.inactive = s["inactive"]
    st.particles = s["particles"]
    if s.get("seed") is not None:
        st.seed = s["seed"]
    st.temperature = {"method": s["temperature"]["method"],
                      "range": tuple(s["temperature"]["range"])}
    box = openmc.stats.Box((-R_core, -R_core, z_fb), (R_core, R_core, z_ft))
    st.source = openmc.IndependentSource(space=box,
                                         constraints={"fissionable": True})
    ent = openmc.RegularMesh()
    ent.dimension = tuple(s["entropy_mesh_dimension"])
    ent.lower_left = (-R_core, -R_core, z_fb)
    ent.upper_right = (R_core, R_core, z_ft)
    st.entropy_mesh = ent

    # ---- tallies ----
    t = cfg["tallies"]
    ll, ur = (-R_core, -R_core, z_fb), (R_core, R_core, z_ft)

    mesh_xy = openmc.RegularMesh()
    mesh_xy.dimension = tuple(t["fission_xy_dimension"])
    mesh_xy.lower_left, mesh_xy.upper_right = ll, ur
    t_xy = openmc.Tally(name="fission_xy")
    t_xy.filters = [openmc.MeshFilter(mesh_xy)]
    t_xy.scores = ["fission"]

    mesh_z = openmc.RegularMesh()
    mesh_z.dimension = (1, 1, t["fission_z_bins"])
    mesh_z.lower_left, mesh_z.upper_right = ll, ur
    t_z = openmc.Tally(name="fission_z")
    t_z.filters = [openmc.MeshFilter(mesh_z)]
    t_z.scores = ["fission"]

    sp = t["spectrum"]
    edges = np.logspace(math.log10(sp["energy_min_eV"]),
                        math.log10(sp["energy_max_eV"]), sp["groups"] + 1)
    t_spec = openmc.Tally(name="core_spectrum")
    t_spec.filters = [openmc.CellFilter(core), openmc.EnergyFilter(edges)]
    t_spec.scores = ["flux"]

    t_bal = openmc.Tally(name="material_balance")
    t_bal.filters = [openmc.MaterialFilter(list(lib.cache.values()))]
    t_bal.scores = ["absorption", "fission", "nu-fission"]

    model = openmc.Model(
        geometry=geometry,
        materials=openmc.Materials(list(lib.cache.values())),
        settings=st,
        tallies=openmc.Tallies([t_xy, t_z, t_spec, t_bal]))

    info = dict(lib=lib, fuel_mat=fuel_mat, n_fuel=n_fuel, n_hp=n_hp,
                z_fb=z_fb, z_ft=z_ft, z_rb=z_rb, z_rt=z_rt, z_bot=z_bot,
                **dims)
    return model, info


def build_kinf_model(cfg):
    """Core lattice only, reflective on every face: a no-leakage k estimate.
        If k-inf is less than 1.0, no core size can reach criticality. Must 
        adjust materials or fuel fraction instead.
    """
    validate(cfg)
    g = cfg["geometry"]
    lib = MaterialLibrary(cfg["materials"])
    L = g["core"]["active_length"]
    R_core = _radial_dims(cfg)["R_core"]

    fuel_u = fuel_pin_universe(g, lib)
    hp_u = heat_pipe_universe(g, lib)
    lattice, n_fuel, n_hp = build_lattice(g, lib, fuel_u, hp_u)

    s_cyl = openmc.ZCylinder(r=R_core, boundary_type="reflective")
    s_lo = openmc.ZPlane(z0=-L / 2, boundary_type="reflective")
    s_hi = openmc.ZPlane(z0=L / 2, boundary_type="reflective")
    core = openmc.Cell(name="core (reflective)", fill=lattice,
                       region=-s_cyl & +s_lo & -s_hi)

    s = cfg["settings"]
    st = openmc.Settings()
    st.run_mode = "eigenvalue"
    st.batches, st.inactive = s["batches"], s["inactive"]
    st.particles = s["particles"]
    st.temperature = {"method": s["temperature"]["method"],
                      "range": tuple(s["temperature"]["range"])}
    box = openmc.stats.Box((-R_core, -R_core, -L / 2),
                           (R_core, R_core, L / 2))
    st.source = openmc.IndependentSource(space=box,
                                         constraints={"fissionable": True})
    model = openmc.Model(
        geometry=openmc.Geometry(openmc.Universe(cells=[core])),
        materials=openmc.Materials(list(lib.cache.values())),
        settings=st)
    return model


# --------------------------------------------------------------------------
# 6. BOOKKEEPING
# --------------------------------------------------------------------------
def print_summary(cfg, info):
    g = cfg["geometry"]
    lib, pi = info["lib"], math.pi
    L = g["core"]["active_length"]
    d = g["control_drums"]
    R_barrel, R_ref = info["R_barrel"], info["R_ref"]

    def kg(vol_cm3, mat):
        return vol_cm3 * mat.density / 1000.0

    refl = lib.get(g["reflector"]["material"])
    vol_ax = 2 * pi * R_ref ** 2 * g["reflector"]["axial_thickness"]
    vol_drums = d["count"] * pi * d["radius"] ** 2 * L
    vol_abs = (d["count"] * pi * (d["radius"] ** 2 -
               (d["radius"] - d["absorber_thickness"]) ** 2) * L
               * (2 * d["absorber_half_angle_deg"] / 360.0))
    vol_rad = pi * (R_ref ** 2 - R_barrel ** 2) * L - vol_drums

    fuel = info["fuel_mat"]
    print("=" * 62)
    print(f"Fuel material : {g['fuel_pin']['fuel_material']}")
    print(f"Fuel pins     : {info['n_fuel']}    Heat pipes: {info['n_hp']}")
    print(f"Core radius   : {info['R_core']:.1f} cm   "
          f"Active length: {L:.1f} cm   Reflector OD: {2 * R_ref:.1f} cm")
    print(f"Fuel mass     : {fuel.get_mass() / 1000:.1f} kg")
    print(f"U-235 mass    : {fuel.get_mass('U235') / 1000:.1f} kg")
    print(f"Reflector     : {kg(vol_ax + vol_rad, refl):.0f} kg "
          f"(radial + axial, drums excluded)")
    print(f"Drums         : {kg(vol_abs, lib.get(d['absorber_material'])):.1f}"
          f" kg absorber + "
          f"{kg(vol_drums - vol_abs, lib.get(d['body_material'])):.0f}"
          f" kg body")
    if g["lattice"].get("interstitial_material"):
        lat = g["lattice"]
        n_pos = info["n_fuel"] + info["n_hp"]
        vol = (n_pos * lat["pitch"] ** 2 * math.sqrt(3) / 2
               - info["n_fuel"] * pi * g["fuel_pin"]["r_clad"] ** 2
               - info["n_hp"] * pi * g["heat_pipe"]["r_outer"] ** 2) * L
        print(f"Interstitial  : "
              f"{kg(vol, lib.get(lat['interstitial_material'])):.0f} kg "
              f"{lat['interstitial_material']} between pins")
    if g["shield"]["enabled"]:
        for layer in g["shield"]["layers"]:
            mass = kg(pi * R_ref ** 2 * layer["thickness"],
                      lib.get(layer["material"]))
            print(f"Shield {layer['name']:<20s}: {mass:.0f} kg")
    print("Excludes: pin cladding, heat pipes, structure, converter, radiator")
    print("=" * 62)


# --------------------------------------------------------------------------
# 7. PLOTTING
# --------------------------------------------------------------------------
def plot_geometry(cfg, model, info):
    import matplotlib
    matplotlib.use("Agg")                 # no display needed inside WSL
    import matplotlib.pyplot as plt

    lib = info["lib"]
    colors = lib.colors(cfg["materials"])
    root = model.geometry.root_universe
    outdir = cfg["plot"]["outdir"]
    os.makedirs(outdir, exist_ok=True)

    R, Rb = info["R_ref"], info["R_barrel"]
    z_span = info["z_rt"] - info["z_bot"]
    z_c = 0.5 * (info["z_rt"] + info["z_bot"])
    views = {
        "xy_midplane": dict(basis="xy", origin=(0, 0, 0),
                            width=(2.1 * R, 2.1 * R)),
        "xy_core_zoom": dict(basis="xy", origin=(0, 0, 0),
                             width=(2.1 * Rb, 2.1 * Rb)),
        "xz_section": dict(basis="xz", origin=(0, 0, z_c),
                           width=(2.1 * R, 1.05 * z_span)),
    }
    for name, view in views.items():
        ax = root.plot(pixels=cfg["plot"]["pixels"], color_by="material",
                       colors=colors, legend=True,
                       legend_kwargs={"loc": "upper left",
                                      "bbox_to_anchor": (1.02, 1),
                                      "fontsize": 9,
                                      "title": "Materials"},
                       **view)
        plt.savefig(os.path.join(outdir, f"{name}.png"), dpi=150,
                    bbox_inches="tight")
        plt.close("all")
    print(f"Saved geometry plots to ./{outdir}/")


# --------------------------------------------------------------------------
# 8. RUN MODES
# --------------------------------------------------------------------------
def run_case(cfg, cwd):
    model, info = build_model(cfg)
    sp_path = model.run(cwd=cwd)
    with openmc.StatePoint(sp_path) as sp:
        return sp.keff, info


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--config", default=DEFAULT_CONFIG)
    ap.add_argument("--mode", choices=["plot", "run", "sweep", "kinf"],
                    default="plot")
    ap.add_argument("--set", action="append", dest="overrides",
                    metavar="PATH=VALUE",
                    help="override a config value, e.g. geometry.lattice.pitch=2.0")
    args = ap.parse_args()

    cfg = apply_overrides(load_config(args.config), args.overrides)

    if args.mode == "plot":
        model, info = build_model(cfg)
        print_summary(cfg, info)
        plot_geometry(cfg, model, info)

    elif args.mode == "run":
        k, info = run_case(cfg, cwd=cfg["run"]["base_dir"])
        print_summary(cfg, info)
        print(f"k-eff = {k}")

    elif args.mode == "kinf":
        model = build_kinf_model(cfg)
        sp_path = model.run(cwd="run_kinf")
        with openmc.StatePoint(sp_path) as sp:
            print(f"k-inf (approx, fully reflected core) = {sp.keff}")

    elif args.mode == "sweep":
        results = {}
        for deg in cfg["run"]["drum_sweep_deg"]:
            c = deepcopy(cfg)
            c["geometry"]["control_drums"]["angle_deg"] = deg
            k, _ = run_case(c, cwd=f"run_drum_{int(deg)}")
            results[deg] = k
            print(f"drum angle = {deg:6.1f} deg  ->  k-eff = {k}")
        k_hi = max(results.values(), key=lambda k: k.n)
        k_lo = min(results.values(), key=lambda k: k.n)
        worth = (1 / k_lo - 1 / k_hi) * 1e5
        print(f"Drum worth (lowest-k -> highest-k case): "
              f"{worth.n:.0f} +/- {worth.s:.0f} pcm")


if __name__ == "__main__":
    main()
