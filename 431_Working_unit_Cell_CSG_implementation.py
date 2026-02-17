# Test file for comparing DAGMC to OpenMC CGS 
# This file based on the OpenMC Pincell juptyr notebook, modified for a direct python input. 

# <<< import appropriate modules:
import openmc 
import datetime
import matplotlib.pyplot as py

now = datetime.datetime.now()
print('\n'*25) # space the terminal
print('Simulation run at:', now)

# <<< Define material properties:

# uo2
uo2 = openmc.Material(1,"uo2") # has ID = 1
uo2.add_nuclide('U235', 0.03)
uo2.add_nuclide('U238',0.97)
uo2.add_nuclide('O16', 2.0)
uo2.set_density('g/cm3', 10.0)

# zirconium
zirconium = openmc.Material(2,"zirconium") # has ID = 2
zirconium.add_element('Zr',1.0)
zirconium.set_density('g/cm3',6.6)

# water
water = openmc.Material(3,"h2o") # has ID = 3 
water.add_nuclide('H1', 2.0)
water.add_nuclide('O16', 1.0)
water.set_density('g/cm3', 1.0)
water.add_s_alpha_beta('c_H_in_H2O') # adds S_AlphaBeta table for bound CS at thermal energy

# Reminder: materials need at minimum nuclides to add from CS_xml and a density. 

# Collect and export to materals_xml: 
mats = openmc.Materials([uo2, zirconium,water])
mats.export_to_xml()

# <<< Define geometry properties: 

# Define cylindrical regions: 
fuel_or = openmc.ZCylinder(r=0.39)
clad_ir = openmc.ZCylinder(r=0.40)
clad_or = openmc.ZCylinder(r=0.46)

# make the regions: 
fuel_region = -fuel_or
gap_region = +fuel_or & -clad_ir
clad_region = +clad_ir & -clad_or

# make the appropriate CSG cells: 

# fuel:
fuel = openmc.Cell(1,'fuel') # has cell ID = 1
fuel.fill = uo2 
fuel.region = fuel_region

# gap: 
gap = openmc.Cell(2,'air gap') # has cell ID = 2
gap.region = gap_region

# clad: 
clad = openmc.Cell(3, 'clad') # has cell ID = 3
clad.fill = zirconium
clad.region = clad_region

# make bound geometry for coolant OUTSIDE of CSG: 
pitch = 1.26
left = openmc.XPlane(x0 = -pitch/2, boundary_type = 'reflective')
right = openmc.XPlane(x0 = pitch/2, boundary_type = 'reflective') 
bottom = openmc.YPlane(y0 = -pitch/2, boundary_type = 'reflective')
top = openmc.YPlane(y0 = pitch/2, boundary_type = 'reflective') 

water_region = +left & -right & +bottom &-top & +clad_or

moderator = openmc.Cell(4,'moderator')
moderator.fill = water
moderator.region = water_region

root = openmc.Universe(cells=(fuel,gap,clad,moderator))

geom = openmc.Geometry(root)
geom.export_to_xml()

point = openmc.stats.Point((0,0,0))
src = openmc.IndependentSource(space = point)

# <<< Define Run Settings: 

settings = openmc.Settings()
settings.source = src
settings.batches = 100
settings.inactive = 10
settings.particles = 1000
settings.export_to_xml()

# <<< define tallies: 
cell_filter = openmc.CellFilter(fuel)
t = openmc.Tally(1)
t.filters = [cell_filter]
t.nuclides = ['U235']
t.scores = ['total','fission','absorption','(n,gamma)']

tallies = openmc.Tallies([t])
tallies.export_to_xml()

# <<< Plot geometry to verify good import: 

plot = openmc.Plot(name = '431_CSG_plot') # creates the plot 
plot.basis = 'xy' # (viewed from above in Z direc)
plot.origin = (0,0,0) # Center point 
plot.width = (2, 2) # X by Y cm view area 
plot.pixels = (1080, 1080) # resolution
plot.color_by = 'material' # Auto generates based on random color selection
plot.filename = '431_CSG_Example' # writes the actual image file into the bin.
# Note: the images will persist, remember to delete them after each run, they won't update automatically, but plot_1 should.

plots = openmc.Plots([plot])
plots.export_to_xml()

# plot.colors = (# insert colors here)
# https://docs.openmc.org/en/stable/pythonapi/generated/openmc.Plot.html


# <<< Run OpenMC in plotting mode: 
# openmc.plot_geometry()

# <<< run OpenMC to initialize results
openmc.run()

# <<< print spacing to terminal: 
print('\n') 
print('End of simulation.')
print('\n')
