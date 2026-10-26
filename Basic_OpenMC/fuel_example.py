import openmc

#FUEL MATERIALS
fuel = openmc.Material(name="HALEU UO2")
fuel.add_element("U", 1.0, enrichment=19.75, percent_type="wo")  # or add_nuclide
fuel.add_element("O", 2.0)
fuel.set_density("g/cm3", 10.4)

mats = openmc.Materials([fuel])

#BUILD FUEL CELL
r = openmc.Sphere(r=20.0, boundary_type="vacuum")
cell = openmc.Cell(fill=fuel, region=-r)
geom = openmc.Geometry([cell])

#SETTINGS
sets = openmc.Settings()
sets.batches, sets.inactive, sets.particles = 150, 50, 20000
sets.source = openmc.IndependentSource(space=openmc.stats.Point())

#OPENMC RUN
openmc.Model(geom, mats, sets).run()