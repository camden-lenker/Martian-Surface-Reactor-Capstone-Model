# Martian Surface Fission Reactor: OpenMC Model

A parametric Monte Carlo model of a 10-40 kWe, HALEU-fueled, heat-pipe-cooled
surface reactor for Mars, built for our nuclear engineering senior design
capstone. The model covers neutronics and shielding geometry; thermal
hydraulics is handled separately and fed in as fixed temperatures.

> **Status: v0 starting point.** The default dimensions are not an optimized or
> validated design, and the model has not yet been benchmarked against
> experiments. Early runs indicate the default fast-spectrum core is
> leakage-limited, which is why moderator options are included. Update this
> section as the team's results come in.

I am learning how github pulls work and branch stuff works -Logan
---

## Contents

- [How it works](#how-it-works)
- [Repository layout](#repository-layout)
- [Setup](#setup)
- [Running the model](#running-the-model)
- [Configuration reference](#configuration-reference)
- [Outputs](#outputs)
- [Assumptions and limitations](#assumptions-and-limitations)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)

---

## How it works

All materials, dimensions, and simulation settings live in
`reactor_config.json`. `reactor_model.py` reads that file and builds the
OpenMC model, so a design change is a data edit, not a code edit.

```
reactor_config.json  --->  reactor_model.py  --->  OpenMC  --->  k-eff, tallies, plots
 (materials, geometry,      (builds model,          (Monte Carlo
  settings)                  validates, runs)         transport)
```

### Default core concept

| Region | Description |
|---|---|
| Fuel | HALEU UN (19.75 wt% U-235, N-15 enriched), 95% theoretical density. UO2 is also defined. |
| Core lattice | Hexagonal lattice of fuel pins and sodium heat pipes. Every Nth slot is a heat pipe (placeholder pattern). |
| Pins | Fuel pellet, void gap, Nb-1Zr cladding. |
| Heat pipes | Void vapor core, liquid sodium wick, Nb-1Zr wall. |
| Between pins | Void by default. Can be filled with a material (BeO, YH2, ...) to model a moderating matrix. |
| Barrel | Nb-1Zr, surrounds the core. |
| Reflector | BeO radial reflector plus BeO axial end caps. |
| Control | Rotating control drums in the radial reflector: BeO bodies with a B4C absorber arc (90% B-10). |
| Shield | Optional stack of disks on the bottom end (default: tungsten for gammas, Li-7 hydride for neutrons). |

### Key modeling choices

- **Materials are built lazily.** Only materials referenced by the geometry are
  exported to OpenMC, so unused entries (for example the UO2 fuel when UN is
  selected) cost nothing.
- **Control drums use global planes.** Each absorber arc is defined with two
  vertical planes through the drum axis, so there are no rotated universes to
  get wrong. The `angle_deg` setting rotates the arc: `0` points the absorber at
  the core (shutdown), `180` points it away (full out).
- **Temperatures are fixed inputs** set per component in the JSON. Cross
  sections are interpolated between library temperatures.
- **The model validates inputs** before running: missing material references,
  overlapping drums, overlapping pins, and bad radii raise a clear error.
- **Tallies** include fission rate maps (radial and axial), a core flux
  spectrum, and a per-material absorption and fission balance.

---

## Repository layout

```
.
├── README.md
├── reactor_config.json     # all materials, geometry, and settings
├── reactor_model.py        # reads the config, builds and runs the OpenMC model
├── environment.yml         # conda environment (see Setup)
└── .gitignore
```

Run output (`run_*/`, `plots/`, `*.h5`) is generated locally and is not
committed.

---

## Setup

### 1. Operating system

OpenMC's conda-forge packages support **Linux and macOS only**.

- **Linux / macOS:** nothing extra needed.
- **Windows:** use WSL2 with Ubuntu. In an administrator PowerShell, run
  `wsl --install -d Ubuntu`, restart, and do everything below inside the Ubuntu
  terminal. Keep the project in the Linux filesystem (for example `~/capstone`),
  not under `/mnt/c/`, since file access across that boundary is slow.

### 2. Conda environment

Install [Miniforge](https://github.com/conda-forge/miniforge), then:

```bash
conda create -n openmc-env -c conda-forge --override-channels python=3.11 openmc
conda activate openmc-env
conda install -c conda-forge --override-channels \
    numpy scipy pandas matplotlib h5py uncertainties jupyterlab coolprop
```

Optional, for interactive geometry viewing:

```bash
conda install -c conda-forge --override-channels openmc-plotter
```

Verify the install:

```bash
python -c "import openmc; print(openmc.__version__)"
```

**OpenMC 0.14 or newer is required** (the model uses source constraints
introduced in 0.14).

| Package | Needed for | Required? |
|---|---|---|
| `openmc` | Transport solver and Python API | Yes |
| `numpy` | Energy grids and array handling | Yes |
| `matplotlib` | Geometry plots (`--mode plot`) | Yes |
| `uncertainties` | k-eff with statistical uncertainty (OpenMC depends on it) | Yes (installed with OpenMC) |
| `h5py` | Reading statepoint files directly | Recommended |
| `scipy`, `pandas` | Post-processing, transient and trade-study scripts | Recommended |
| `jupyterlab` | Exploratory analysis | Optional |
| `coolprop` | Gas and CO2 properties for radiator work (not liquid metals) | Optional |
| `openmc-plotter` | Interactive geometry viewer | Optional |

To share the exact environment with the team, export it and commit the result:

```bash
conda env export --from-history > environment.yml
```

Teammates then run `conda env create -f environment.yml`.

### 3. Nuclear data (not stored in this repo)

Download the pre-processed **ENDF/B-VIII.0** HDF5 library from the Official Data
Libraries page at <https://openmc.org/data>, extract it, and point OpenMC at it:

```bash
cd ~
wget -O endfb80.tar.xz "<link from the data page>"
tar -xf endfb80.tar.xz
echo 'export OPENMC_CROSS_SECTIONS=$HOME/endfb-viii.0-hdf5/cross_sections.xml' >> ~/.bashrc
source ~/.bashrc
ls $OPENMC_CROSS_SECTIONS      # should print the path, not an error
```

Adjust the folder name to match what you extracted. A depletion chain file from
the same page is needed later for burnup runs, but not for the steps below.

**Check the thermal scattering names** used by the config exist in your
library:

```bash
grep -o 'c_[A-Za-z0-9]*_in_BeO' $OPENMC_CROSS_SECTIONS | sort -u
grep -o 'c_[A-Za-z0-9]*_in_YH2' $OPENMC_CROSS_SECTIONS | sort -u
```

If the names differ, update the `s_alpha_beta` lists in the JSON.

---

## Running the model

Run from the repository folder with the environment active.

```bash
conda activate openmc-env
python reactor_model.py --mode plot     # geometry check (start here)
python reactor_model.py --mode kinf     # no-leakage diagnostic
python reactor_model.py --mode run      # k-eigenvalue run
python reactor_model.py --mode sweep    # control drum worth
```

| Mode | What it does |
|---|---|
| `plot` | Prints a mass summary and saves geometry images to `plots/`. Always check these first, especially that the drum absorber arcs face away from the core in the default (drums out) state. |
| `kinf` | Runs the core lattice with reflective boundaries on all faces. A k-infinity below 1.0 means no core size can reach criticality with these materials. Slightly approximate (includes the BeO corner filler). |
| `run` | Full model k-eigenvalue calculation. Prints the mass summary and k-eff with its uncertainty. |
| `sweep` | Runs the angles listed in `run.drum_sweep_deg` (default 0 and 180 degrees) and reports the total drum worth in pcm. |

### Quick test run

Use low statistics to confirm everything works before committing real compute
time:

```bash
python reactor_model.py --mode run \
    --set settings.particles=2000 --set settings.batches=60 --set settings.inactive=20
```

### Overriding any config value

Use `--set path=value`. The value is parsed as JSON, so strings need inner
quotes:

```bash
python reactor_model.py --mode run --set geometry.lattice.n_rings=11
python reactor_model.py --mode kinf --set 'geometry.lattice.interstitial_material="BeO"'
python reactor_model.py --mode run \
    --set 'geometry.lattice.interstitial_material="YH2"' \
    --set geometry.lattice.interstitial_temperature_K=800
python reactor_model.py --mode run --set geometry.shield.enabled=false
```

Use `--config other_design.json` to run a different config file. Give runs
separate output folders with `--set run.base_dir=run_name` so they don't
overwrite each other.

### Parameter scan example

```bash
for n in 9 11 13; do
  python reactor_model.py --mode run \
    --set geometry.lattice.n_rings=$n \
    --set settings.particles=5000 --set settings.batches=100 --set settings.inactive=30 \
    --set run.base_dir=run_n$n
done
```

### Reading results later

```bash
python -c "import openmc; sp=openmc.StatePoint('run_base/statepoint.300.h5'); print(sp.keff)"
```

Change `300` to match the number of batches you ran.

### Reading the output

- OpenMC prints a running k-eff each batch; early values swing, which is normal.
- Check that **Shannon entropy** has settled by the end of the inactive batches.
  If it is still drifting, increase `settings.inactive`.
- Compare k-eff to k-inf. The ratio is the fraction of neutrons that do not
  leak, and tells you whether a design is leakage-limited.

---

## Configuration reference

Lengths in cm, densities in g/cm3, temperatures in K. Keys beginning with `_`
are notes and are ignored.

| Section | Key fields |
|---|---|
| `materials` | Per material: `density_g_cm3`, optional `fraction_of_theoretical`, `components` (elements or nuclides with `percent`, `percent_type`, optional enrichment), optional `s_alpha_beta`, `plot_color`. |
| `geometry.fuel_pin` | `fuel_material`, `clad_material`, `r_fuel`, `gap`, `r_clad`, fuel and clad temperatures. |
| `geometry.heat_pipe` | `wall_material`, `wick_material`, `r_vapor`, `r_wick`, `r_outer`, temperature. |
| `geometry.lattice` | `pitch`, `n_rings`, `heat_pipe_every`, `orientation`, `filler_material`, `interstitial_material` (null = void), `interstitial_temperature_K`. |
| `geometry.core` | `active_length`, `margin_to_barrel`, `barrel_thickness`, `barrel_material`. |
| `geometry.reflector` | `material`, `temperature_K`, `axial_thickness`, `outboard_margin`. |
| `geometry.control_drums` | `count`, `radius`, `gap_to_barrel`, body and absorber materials, `absorber_thickness`, `absorber_half_angle_deg`, `angle_deg` (0 = absorber toward core, 180 = away). |
| `geometry.shield` | `enabled`, and a `layers` list (`name`, `material`, `thickness`, `temperature_K`) stacked downward from the bottom reflector. |
| `settings` | `batches`, `inactive`, `particles`, `seed`, temperature interpolation range, entropy mesh dimensions. |
| `tallies` | Fission map resolution, axial bin count, spectrum energy range and group count. |
| `run` | `base_dir`, `drum_sweep_deg`. |
| `plot` | `outdir`, `pixels`. |

Materials defined in the starter config: `fuel_UN`, `fuel_UO2`, `clad_NbZr`,
`sodium`, `BeO`, `B4C`, `tungsten`, `LiH`, `YH2`. Add a new material by adding
an entry under `materials` and referencing its key from the geometry.

---

## Outputs

Each run writes to its own folder (`run_base/`, `run_kinf/`,
`run_drum_<angle>/`, or whatever `run.base_dir` is set to) containing the
OpenMC input XML files, the statepoint `.h5` file, and tally results. Geometry
images go to `plots/` (`xy_midplane.png`, `xy_core_zoom.png`,
`xz_section.png`).

The mass summary printed by `plot` and `run` lists fuel and U-235 mass,
reflector, drum, interstitial, and shield masses. It excludes cladding, heat
pipes, structure, power conversion, and radiators.

---

## Assumptions and limitations

- Fuel-clad gaps, pin interstitials (unless a material is set), and heat pipe
  vapor cores are void.
- Heat pipes end at the core. There are no penetrations through the top
  reflector.
- Control drums are BeO cylinders with a B4C arc, with no sleeves, bearings, or
  actuators.
- Temperatures are fixed inputs, not coupled to a thermal-hydraulics solution.
- The shield is a flat disk stack, not a shadow cone sized to a dose limit.
- The fast, moderated, and hydrogenous cases have not been benchmarked. Run a
  published benchmark (for example KRUSTY) before trusting absolute k-eff.
- **YH2 density and H:Y ratio are assumptions** in the config. Real yttrium
  hydride is often sub-stoichiometric, and hydrogen retention limits usable
  temperature. Verify against a source before using results.
- The heat pipe pattern (every Nth slot) is a placeholder, not a designed
  layout.
- No Mars environment (regolith berm, dose to a habitat) is modeled yet.
  Adding a fixed-source dose run is a planned next step.

---

## Troubleshooting

| Problem | Likely cause and fix |
|---|---|
| `No match found for: openmc` during conda install | Native Windows. Use WSL2 (see Setup). |
| `cross_sections.xml not found` | `OPENMC_CROSS_SECTIONS` is not set in this terminal. Check with `echo $OPENMC_CROSS_SECTIONS`. |
| Error about a thermal scattering table | The name in `s_alpha_beta` does not match your library. Grep the library (see Setup) and edit the JSON. |
| `invalid choice` for `--mode` | You are running an older copy of `reactor_model.py`. Pull the latest version. |
| `no key ... in config` from `--set` | A typo in the path, or an older JSON missing that key. Compare against the current `reactor_config.json`. |
| Control drums overlap error | Reduce `control_drums.count` or `radius`, or increase the core size. |
| Very low k-eff | Run `--mode kinf` first. If k-inf is above 1 the core is leaking; consider a larger core, a moderating matrix, or a thicker reflector. Also confirm in the plots that the drums are in the intended position. |
| Slow runs under WSL | Keep the project in the Linux filesystem, not `/mnt/c/`. You can raise WSL memory and CPU limits in `.wslconfig`. |

---

## Contributing

- Work on a branch and open a pull request; changes to `main` should be
  reviewed. A mistyped density or enrichment in the JSON can silently change
  results.
- Describe what you changed and why, and include the k-eff (with uncertainty)
  before and after when a change affects physics.
- Do not commit the nuclear data library or run output. Suggested `.gitignore`:

  ```
  run_*/
  plots/
  *.h5
  *.xml
  __pycache__/
  ```

- Cite sources for any new material property or dimension, in the JSON `_notes`
  field or in the pull request.
- **Use only publicly available sources.** Do not add export-controlled,
  proprietary, or sponsor-restricted information to this repository. If you are
  unsure whether something can be shared, ask your advisor or your university's
  export control office first.

## Team and license

Add team members, advisor, course, and a license here.
