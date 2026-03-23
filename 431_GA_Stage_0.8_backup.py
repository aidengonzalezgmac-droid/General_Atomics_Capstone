# import modules: 

import openmc 
import datetime
import matplotlib.pyplot as py # in case it is needed. 
import math
from openmc_plasma_source import fusion_ring_source # for point source 

R_major        = 187.5   # cm — plasma major radius
ring_radius_cm = R_major
z_plasma       = 119.4   # cm — tokamak Z offset in CAD coordinates

# space terminal by printing current system time to verify operability.
now = datetime.datetime.now()
print('\n'*25) # space the terminal
print('Simulation run at:', now)

# MATERIALS:
# Using the stage 0 cubit export: 

# Materials required: 
    # Wall_Concrete
    # Structural_Steel
    # inconel
    # Roof_Shielding
    # Wall_Shielding
    # Air
    # Void (For within the tokamak)

# Concrete: 
# Create ordinary concrete
concrete = openmc.Material(name='Concrete')
concrete.set_density('g/cm3', 2.3)
concrete.add_element('H', 0.010, percent_type='wo')
concrete.add_element('C', 0.001, percent_type='wo')
concrete.add_element('O', 0.529, percent_type='wo')
concrete.add_element('Na', 0.016, percent_type='wo')
concrete.add_element('Mg', 0.002, percent_type='wo')
concrete.add_element('Al', 0.034, percent_type='wo')
concrete.add_element('Si', 0.337, percent_type='wo')
concrete.add_element('K', 0.013, percent_type='wo')
concrete.add_element('Ca', 0.044, percent_type='wo')
concrete.add_element('Fe', 0.014, percent_type='wo')

# Barite: 
# Create barite (BaSO4) for heavy aggregate
barite = openmc.Material(name='Barite')
barite.set_density('g/cm3', 4.5)
barite.add_element('Ba', 0.589, percent_type='wo')
barite.add_element('S', 0.137, percent_type='wo')
barite.add_element('O', 0.274, percent_type='wo')

Wall_Concrete = openmc.Material(name='Wall_Concrete')
Wall_Concrete = openmc.Material.mix_materials([concrete, barite], [0.85, 0.15], 'wo')

#SA508 steel from steelprogroup.com
Structural_Steel = openmc.Material(name='Structural_Steel') 
Structural_Steel.set_density('g/cm3', 7.85)
Structural_Steel.add_element('C', 0.0025, percent_type='wo')   # 0.25% (taking max)
Structural_Steel.add_element('Mn', 0.0120, percent_type='wo')  # 1.20%
Structural_Steel.add_element('Si', 0.0060, percent_type='wo')  # 0.60%
Structural_Steel.add_element('P', 0.00035, percent_type='wo')  # 0.035% (taking max)
Structural_Steel.add_element('S', 0.00035, percent_type='wo')  # 0.035% (taking max)
Structural_Steel.add_element('Ni', 0.0025, percent_type='wo')  # 0.25%
Structural_Steel.add_element('Cr', 0.0025, percent_type='wo')  # 0.25%
Structural_Steel.add_element('Mo', 0.0008, percent_type='wo')  # 0.08%
Structural_Steel.add_element('Cu', 0.0035, percent_type='wo')  # 0.35%
Structural_Steel.add_element('V', 0.0005, percent_type='wo')   # 0.05%
Structural_Steel.add_element('B', 0.00003, percent_type='wo')  # 0.003% (taking max)
Structural_Steel.add_element('Fe', 0.97537, percent_type='wo')

# Inconel shell: 
Inconel = openmc.Material(name = "Inconel")
Inconel.set_density("g/cm3", 8.44) # from the public special materials document:
Inconel.add_element("Ni", 0.5922) # increased from 58.0 min to make up. 
Inconel.add_element("Cr", 0.21) # increased from 20.0 to 22
Inconel.add_element("Ni", 0.048) # decreased from 5 max
Inconel.add_element("Mo", 0.08)
Inconel.add_element("Nb", 0.0315)
Inconel.add_element("C", 0.01)
Inconel.add_element("Mn", 0.0050)
Inconel.add_element("Si", 0.0050)
Inconel.add_element("P", 0.00015)
Inconel.add_element("S", 0.00015)
Inconel.add_element("Al", 0.0040)
Inconel.add_element("Ti", 0.0040)
Inconel.add_element("Co", 0.01)

Roof_Shielding=openmc.Material(name='Roof_Shielding')
Roof_Shielding.set_density('g/cm3',1.12)
Roof_Shielding.add_element('H',0.133, percent_type='wo')
Roof_Shielding.add_element('C',0.817, percent_type='wo')
Roof_Shielding.add_element('B', 0.05, percent_type='wo')

# From Anderson et al. the density of the polyboron shield wall is 1.12 g/cc
Wall_Shielding=openmc.Material(name='Wall_Shielding')
Wall_Shielding.set_density('g/cm3',1.12)
Wall_Shielding.add_element('H',0.133, percent_type='wo')
Wall_Shielding.add_element('C',0.817, percent_type='wo')
Wall_Shielding.add_element('B', 0.05, percent_type='wo')

# air fill:
Air = openmc.Material(name="Air")
Air.set_density("g/cm3", 0.001225)
Air.add_element("N", 0.78)
Air.add_element("O", 0.21)
Air.add_element("Ar", 0.01)

# void fill (for within tokamak)
Void = openmc.Material(name="Void")
Void.set_density('g/cm3',0.1)
Void.add_element('H',1) # Arbitrary content, in reality there's a lot of plasma in there, it would be interesting to ask
# Igor and Chris about neutron-plasma interactions given the birth energies being so high.

mats = openmc.Materials([Wall_Concrete,Structural_Steel,Inconel,Roof_Shielding,Wall_Shielding,Air,Void])
# mats.cross_sections = "/bin/endfb-viii.0-hdf5/cross_sections.xml" #please rem to use your own file
# have you considered just *not?*
# references CS for docker version

mats.export_to_xml() # <<<<<<<<<<<>>>>>>>>>>>>

# bootstrap modification by re-naming Wall_Concrete in the file.

# GEOMETRY: 
# import the dagmc file after placed in here
Stage_0_Model = openmc.DAGMCUniverse(filename = '/bin/431_GA_Stage_0.8.h5m', auto_geom_ids = True)
room = openmc.Cell(name = 'Machine_Hall')
room.fill = Stage_0_Model

Boundcyl = openmc.Sphere(r=400,name ='bounds')
Boundcylreg = -Boundcyl # everything inside the cylinder.
Boundcell = openmc.Cell()
Boundcell.region = Boundcylreg
Boundcell.fill = Void

rootuni = openmc.Universe()
# rootuni.add_cell(Boundcell)
rootuni.add_cell(room)

root = openmc.Cell(name='root_uni')
root.fill = rootuni

geometry = openmc.Geometry([root]) # NOTE: FIGURE OUT HOW TO DO WHAT ROCCO SAID.
# h = geometry.get_all_universes()
# print(h)
# assumption: 
#  there should be an appended cell that is the implicit compliment for the geometry per 
# https://github.com/svalinn/DAGMC/issues/934, so if we assign the mesh to cell "9" we should be able to select
# the entire room to see what's bouncing.

# geometry = openmc.Geometry([room,Boundcell]) <- working originally
geometry.export_to_xml()
# I believe the Boundcell actually over-rides the void fill on the exterior. Neat.
# apply geometry

# add source 
# src_e = openmc.stats.Discrete(x=[12.0,], p=[1.0,])
# point = openmc.stats.Point((100, 100, 100))
# strength = 1, particle = ('neutron'),energy = src_e)

# CHANGE ME: 

my_source = fusion_ring_source(
    radius=R_major / 100,        # cm → m
    angles=(0.0, 2 * math.pi),
    z_placement=119.4 / 100,     # cm → m
    temperature=20000.0,
    fuel={"D": 1.0},
)
"""
# creates a simple isotropic neutron source in the center with 14MeV neutrons
my_source = openmc.Source()
# the distribution of radius is just a single value at the plasma major radius
radius = openmc.stats.Discrete([293.], [1])
# the distribution of source z values is just a single value
z_values = openmc.stats.Discrete([0], [1])
# the distribution of source azimuthal angles values is a uniform distribution between 0 and 0.5 Pi
# these angles must be the same as the reflective angles
angle = openmc.stats.Uniform(a=0., b=math.radians(90))
# this makes the ring source using the three distributions and a radius
my_source.space = openmc.stats.CylindricalIndependent(r=radius, phi=angle, z=z_values, origin=(0.0, 0.0, 0.0))
# sets the direction to isotropic
my_source.angle = openmc.stats.Isotropic()
# sets the energy distribution to a Muir distribution neutrons
my_source.energy = openmc.stats.Muir(e0=14080000.0, m_rat=5.0, kt=20000.0)
# sets energy to be super high at 14.08 MeV -> set to 2.345 for this test:

# https://github.com/fusion-energy/magnetic_fusion_openmc_dagmc_paramak_example/blob/main/2_run_openmc_dagmc_simulation.py

"""

# apply settings:
settings = openmc.Settings()
settings.source = my_source
settings.batches = 100
settings.inactive = 10
settings.particles = 1000
settings.run_mode = 'fixed source'
settings.export_to_xml()

squarecoord = 1250

mesh = openmc.RegularMesh(name = 'regmesh1',) # I want a mesh 
mesh.dimension = (squarecoord,squarecoord,1) # I want the mesh to have to have: 
# (N1 mesh cells in X direction, N2 mesh cells in Y direction, N3 mesh cells in Z direction)
# The total number of entries into the statepoint array is the product of N1,N2,N3
# A mesh of 0 is forbidden for any coordinate, they must be positive, a single value of 1 will give you one mesh in Z direc.
# SEE SHAPING FUNCTION BELOW:

# setting a single mesh.dimension to 1 zeros its contribution. THIS CHANGES THE ORIENTATION WE USE.

# limit seems to be about a 1250 side length oddly enough.

squarecoord2 = squarecoord # setting this to anything different zooms-in the plot. 
mesh.lower_left = (-squarecoord2,-squarecoord2,-1000)
# X, Y, Z inputs must be lower than the upper right coordinates.
mesh.upper_right = (squarecoord2,squarecoord2,1000) #square mesh with rel. arbitrary dimensions
# ref: https://docs.openmc.org/en/stable/pythonapi/generated/openmc.RegularMesh.html


# LL: -1250 -1250 400
# UR: 1250 1250 800 at CELL 9  square.
# lower left disappears when z is removed, no elevation present, this is an x-z
"""
mesh.lower_left = (-1500,-1500,-1000)
mesh.upper_right = (1500,1500,1000) #square mesh with rel. arbitrary dimensions - working dim.
working* with a mesh dim of 1000, 1,1000
"""

# adapted from openmc plotting of the tea-pot example: 
mesh_filter = openmc.MeshFilter(mesh)

inconel_filter = openmc.CellFilter([5]) # reference cell 5
inconel_tally = openmc.Tally()
inconel_tally.filters = [mesh_filter,inconel_filter]
inconel_tally.scores = ['flux']

implicit_comp_filter = openmc.CellFilter([26]) # Holdover from investigating.
implicit_comp_filter = openmc.MaterialFilter([Wall_Concrete]) # Structrual_Steel (variable) is good.

# highlights ALL the structrual_steel componets over one another, looking good and making progress.
# openmc can read Wall_Concrete appropriately. 
# You need ot handle each material filter separately it appears (otherwise your tally out will double)

implicit_comp_tally = openmc.Tally()
implicit_comp_tally.filters = [mesh_filter,implicit_comp_filter]
implicit_comp_tally.scores = ['flux']

# in XY config at -1250 to 1250 in XY direc, Z: -1000 to 1000
# cell 1: I guess the room, yields a single source at center 
# cell 2: seg fault no mapping
# cell 3: Identical to the original room, yields single source at center.
# cell 4: runs, don't know what its going to spit out.
# CELL 4: WRITTEN PERSPECTIVE AS IS SHOWS ROOM FROM SOUTH TO NORTH 
# X GOES RIGHT, Y NORTH, Z GOES UP. 
# CELL 4 CORRESPONDS TO BLOCK 7
# CELL 4 = N+1 FOR CELL ID
# cell 5: nothing appears, but it does run. 
# cell 6: TOP ROOF CORRESPONDS TO CELL 6.(I'm pretty sure)
# cell 7: INCONEL SHELL 
# cell 8: SHIELD BLOCKS SURROUNDING ROOM 
# Cell 9: Full shield extension? north shield block??
# cell 10: Unknown east-facing shield component?
# cell 11: unknown west-facing shield component? Volume 8 on DagMC
# cell 12: suspected horizontal beam at volume 9 on DagMC (free-floating beam)
# cell 12: doesn't even show up 
# cell 13: matches to an even THINNER spot?
# cell 14: corresponds to volume 12 in geometry xml -> I geuss expected???? it looks like a flat beam
# cell 21: skipping to cell 21, 21-2 = 19 = volume 19 in the coreform volume list, which is a shield
# block with known geometry that should appear. 
# cell 21: does NOT conform to expected, produces parts in top components.
# cell 22: CORRESPONDS TO THE SHIELD BLOCK FROM EARLIER CELL 19.
# cell 23: No expectatons: nothing shows up ()
# cell 24: expecting segfault for anything outside of cell boundaries. -> loops and gives a point source
# cell 25: seg fault, out of bounds
# cell 26: seg fault, out of bounds. 


#IMPORTANT: GIVE CELL IDs TO THE BLOCKS FOR PLOTTING.
# IMPORTANT: THE THINGS WE SEE IN THIS PLOT ARE 3D OVERLAYS OF THE PLOT, THEY ARE <<<NOT>>> SLICE PLOTS.



concrete_filter = openmc.CellFilter([3]) # reference cell 3
concrete_tally = openmc.Tally()
concrete_tally.filters = [mesh_filter,concrete_filter]
concrete_tally.scores = ['flux']

# tallies = openmc.Tallies([inconel_tally,concrete_tally])
tallies = openmc.Tallies([implicit_comp_tally])
tallies.export_to_xml()

# plot:  
plot = openmc.Plot(name = '431_GA_Stage_0YZ') # creates the plot 
plot.basis = 'yz' # (viewed from above in Z direc)
plot.origin = (0,0,0) # Center point 
plot.width = (3200, 3200) # X by Y cm view area  # A 31x by 31x cross sectional view of the model,
# centered at the origin # originally 3200x3200
plot.pixels = (1080, 1080) # resolution
plot.color_by = 'material'

plot.filename = '/bin/431_GA_Stage_0YZ' # writes the actual image file into the bin. 
plots = openmc.Plots([plot])
plots.export_to_xml()

# openmc.plot_geometry()
openmc.run()



# verify simulation ended: 
sp = openmc.StatePoint("statepoint.100.h5")

"""
inconel_tally = sp.get_tally(scores=['flux'],id = inconel_tally.id)
# inconel_tally = sp.get_tally(scores=['flux'],id = inconel_tally.id)
inconel_flux = inconel_tally.mean
inconel_flux.shape = (1000,1000)
inconel_flux = inconel_flux[::-1,:]
"""

air_tally = sp.get_tally(scores = ['flux'], id = implicit_comp_tally.id)
air_flux = air_tally.mean 
air_flux.shape = (squarecoord,squarecoord)
# The shape function on the air_flux value must match with the size of the array, it cannot be smaller or larger, it needs to be 1:1
# match this with the number of mesh cells in a given direction for best results
# if you have bounces in your geometry, you should be good.
air_flux = air_flux[::-1,:]

from matplotlib import pyplot as plt
fig = plt.figure(figsize=(20,20))
sub_plot1 = plt.subplot(121, title='inconel_flux')
sub_plot1.imshow(air_flux) # replace with the most recent "flux" variable name.

plt.savefig(fname='/bin/431_GA_Stage_0.8.1_proto_Map.png')
plt.show

print('\n') 
print('End of simulation.')
print('\n')

# report at about 3am on monday: 
    #: so all attempts to use this dagmc universe are reporting a segmentation fault error regardless of 
    # source type, setting type (k eigenvaule estimation versus fixed source)
    # from some of the reports it appears that there's a part of the model which is actively fighting and 
    # or overlapping some other part, likely this was caused by poor planning on my part. 
    # However, did learn how to effectively import solidworks models into dagmc and then to openmc,
    # so to fix this: 

    # We're going to need to import each model one by one into cubit and fix it like that, however hazarding
    # a guess, the problem child is either the inconel shell or the walls. 
    # Don't have a good enough background to find the error and try and brute force a fix, but also don't
    # want to re-compile on this machine out of fear of breaking my openmc install. 

    # To replicate in next test: 
        # Import the inconel shell into Cubit: 
            # Translate back to origin
            # Set the inconel shell up for 
            # hmmmmmmmm...
                # When I set up the air fill volume, I boolean'd the fuck out of the enviro, but I might 
                # have missed the inconel shell interior, which remained untouched
                # for the hell of it, lets see if we can fix this right here and now. 
                # can't too tired. this was a waste of time. 
        # import the inconel shell into cubit 
        # surround it with a box
        # take a volume subtraction out of it and apply appropriate conditions
        # (air, void, inconel)
        # set the appropriate surface boundaries
        # try and import this to openmc with a standard source behavior.

# so the geometry imports correctly, but there's no way of knowing what's causing the fatal fucking memory leak.
