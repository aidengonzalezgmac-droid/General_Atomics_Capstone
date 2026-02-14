# An environment to perform semi-identical calculations of th OpenMC Pincell model.
# The workflow for this was provided by Coreform Inc. (Dr. Patrick Shriwise),
# as a part of their OpenMC Neutronics simulation workshop.

# Relevant commands for cubit/notes are also included in this document,

# This test file includes NEW material IDs to use in place of the pin cell to 
# verify that the materials are applied to the new model, though they are 
# identical in presentation to that of the Working CSG unit cell.


# Test file for comparing DAGMC to OpenMC CGS 
# This file based on the OpenMC Pincell juptyr notebook, modified for a direct python input. 

# Notes: 
# OpenMC interprets all dimensions from cubit as centimeters.
# Cubit models should be coincident if overlapping/interfacing
#   > this is to make simulations more efficient.
# imprint body all -> merge body all (might just need merge???)
# dagmc just needs to use the raw string of a material name, not "uo2" just uo2
# block id and name doesn't matter. remember to select when creating blocks.

# make volumes, then material names for openmc to reference,
# then make blocks on the volumes, then assign materials to those volumes,
# then assign boundary conditions to surfaces (sidesets), then mesh (SURFACE MESH NOT VOLUME) VOLUME FOR 3D FINITE ELEMENTS
# use coarse mesh for DAGMC for resolving (deviation angle) (selecting "all" in select surfaces means all surfaces, necessary for openmc raytracing)
# dagmc coarse setting is useful, reduces min. # of tris but with generally poor spatial distribution, using max surface graduation gets more familiar geometry.
# when mesh is complete, and all other steps complete: 
# free to export mesh to environment you call python file from. 
# with that you can then call it your environment. (here)
# the working filename following this (which should be uploaded to github) is:
# "431_DAGMC_Pincell.h5m"

# NOTE: ALL MODELS ARE EXPECTED TO BE WATERTIGHT ON CUBIT EXPORT
# CAN SEAL MODEL ON EXPORT FROM CUBIT IF NEEDED (No floating triangles)

# NOTE: use the create side-set to set boundary conditions, the SIDESET NAME needs to be = boundary:"INSERT NAME"
# i.e. (name = transmission, vacuum, reflective, so on, so forth)
# reflective = reflecting (the guide mentions reflecting, this isn't supported, use reflective.)
# NOTE: "reflecting" IS THE CORRECT CUBIT DAGMC RESULT. UNSURE IF VAC. CONDITIONS HAVE DIFFERENT NAMES, TRY STUFF.
# "Cubit> surface [id, specific ID or "all"] visibility on/off" to turn off / on visibility.
# to add sidesets to surfaces of interest:  sideset 1 add surface 9 10 11 12 where 9 10 11 and 12 are four independent surfaces.

# Begin actual command implementation: 

# All IDs have been arbitrarily promoted to the CSG model + 100, 
# i.e. (material ID/Cell ID in CSG + 100 = DAGMC material ID/CellID), 
# the model produced in Cubit has also been modeled appropriately for this purpose. 
# the Materials and geometry xml files should update appropriately to reflect the correct
# xml export was performed.



# <<< import appropriate modules:
import openmc 
import datetime
import matplotlib.pyplot as py

now = datetime.datetime.now()
print('\n'*25) # space the terminal
print('Simulation run at:', now)

# <<< Define material properties:

# uo2
uo2 = openmc.Material(101,"uo2") # has ID = 101
uo2.add_nuclide('U235', 0.03)
uo2.add_nuclide('U238',0.97)
uo2.add_nuclide('O16', 2.0)
uo2.set_density('g/cm3', 10.0)

# zirconium
zirconium = openmc.Material(102,"zirconium") # has ID = 102
zirconium.add_element('Zr',1.0)
zirconium.set_density('g/cm3',6.6)

# water
water = openmc.Material(103,"h2o") # has ID = 103 
water.add_nuclide('H1', 2.0)
water.add_nuclide('O16', 1.0)
water.set_density('g/cm3', 1.0)
water.add_s_alpha_beta('c_H_in_H2O') # adds S_AlphaBeta table for bound CS at thermal energy

# Reminder: materials need at minimum nuclides to add from CS_xml and a density. 

# Collect and export to materals_xml: 
mats = openmc.Materials([uo2, zirconium,water])
mats.export_to_xml()


# THE FOLLOWING GEOMETRY IS REPLACED WITH A DAGMC EQUIVALENT UNIVERSE.
# <<< Define geometry properties: 

dagmc_universe = openmc.DAGMCUniverse(filename = '/bin/431_DAGMC_Pincell1.h5m')
# This geometry can be inserted INSIDE of other CSG cells.

# on import, the only thing extra needed is to set the root universe as dagmc, or to add as appropriate.

root = dagmc_universe # sets the root to dagmc_universe.

geom = openmc.Geometry(root) # this is just added to the root universe, might be able to move it around.
geom.export_to_xml()

point = openmc.stats.Point((0,0,0))
src = openmc.IndependentSource(space = point)

# <<< Define Run Settings: 

settings = openmc.Settings()
settings.source = src
settings.batches = 100
settings.inactive = 10
settings.particles = 1000
settings.export_to_xml

# <<< define tallies:
# suppressed tallies since fuel is not present, working on method for how to extract dagmc geometry for tallying - JMD 

"""
cell_filter = openmc.CellFilter(fuel) # FIX THIS.
t = openmc.Tally(1)
t.filters = [cell_filter]
t.nuclides = ['U235']
t.scores = ['total','fission','absorption','(n,gamma)']

tallies = openmc.Tallies([t])
tallies.export_to_xml
"""

# <<< Plot geometry to verify good import: 

plot = openmc.Plot(name = '431_DAGMC_plot_Example') # creates the plot 
plot.basis = 'xy' # (viewed from above in Z direc)
plot.origin = (0,0,0) # Center point 
plot.width = (2, 2) # X by Y cm view area 
plot.pixels = (1080, 1080) # resolution
plot.color_by = 'material' # Auto generates based on random color selection
plot.filename = '431_DAGMC_Example' # writes the actual image file into the bin.
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

# final calibrated results for DAGMC: 

"""
k-effective (Collision)     = 1.40068 +/- 0.00132
 k-effective (Track-length)  = 1.40206 +/- 0.00161
 k-effective (Absorption)    = 1.40061 +/- 0.00095
 Combined k-effective        = 1.40061 +/- 0.00085
 Leakage Fraction            = 0.00000 +/- 0.00000
 """