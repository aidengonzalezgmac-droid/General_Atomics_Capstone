# Note for V-1.0+: 
# higher particle counts allow higher mesh qualities to be used 
# NOTE: Many working variables simply have an extra letter attached to them and don't have very descriptive names. 
# At one point the shorthand variable initials meant something to me, but now I don't recall - JMD 6/12/2026
 
# === IMPORT MODULES & PRINT TERMINAL SIM START ===

import openmc 
import datetime
import numpy as np
import matplotlib.pyplot as py # in case it is needed. 
import math
from openmc_plasma_source import fusion_ring_source # for point source 

# space terminal by printing current system time to verify operability.
now = datetime.datetime.now()
# print('\n'*25) # space the terminal
print('Simulation run at:', now)

# === MATERIAL DEFINITIONS ===

# Temperatures
AmbientT = 294 # in Kelvin, # ambient 70 deg Fahrenheit
SupercondT = 294 # in Kelvin # Nuclear data doesn't have CS for copper interactions at superconducting temps.
# for S_alpha_beta components: https://deepwiki.com/openmc-dev/openmc_mcnp_adapter/4.2.3-thermal-scattering-data

# PFCs (Poloidal Field Coil):
pfc = openmc.Material(name='pfc')
pfc.set_density("g/cm3", 1.31)
# pfc.add_element("C", 0.02)
pfc.add_element("H", 0.06) # changed hydrgoen content to 0.06 from 0.05 because C0 doesn't have 4 kelvin cross secitons.
pfc.add_element("O", 0.04)
pfc.add_element("Cu",0.9)
pfc.temperature = SupercondT

# TFCs (Torodial Field Coil):
tfc = openmc.Material(name="tfc")
tfc.set_density("g/cm3", 1.31)
# tfc.add_element("C", 0.02)
tfc.add_element("H", 0.06)
tfc.add_element("O", 0.04)
tfc.add_element("Cu",0.9)
tfc.temperature = SupercondT

# Neutral Beams: 
nbl = openmc.Material(name="nbl")
nbl.set_density("g/cm3", 1) # base results are 6, lowering for comparisons.
nbl.add_element("Al",0.5)
nbl.add_element("Cu",0.2)
nbl.add_element("Fe",0.1)
nbl.add_element("C",0.1)
nbl.add_element("H",0.1)
nbl.temperature = 310 # increased above ambient 294 k to 310 for average operation of equipment

# Concrete: 
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
Wall_Concrete.temperature = AmbientT

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
Structural_Steel.temperature = AmbientT

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
Inconel.temperature = 294 # Decreased due to assumed local air cooling, APPROXIMATION. 

Roof_Shielding=openmc.Material(name='Roof_Shielding')
Roof_Shielding.set_density('g/cm3',1.12)
Roof_Shielding.add_element('H',0.133, percent_type='wo')
Roof_Shielding.add_element('C',0.817, percent_type='wo')
Roof_Shielding.add_element('B', 0.05, percent_type='wo')
Roof_Shielding.temperature = AmbientT

# From Anderson et al. the density of the polyboron shield wall is 1.12 g/cc
Wall_Shielding=openmc.Material(name='Wall_Shielding')
Wall_Shielding.set_density('g/cm3',1.12)
Wall_Shielding.add_element('H',0.133, percent_type='wo')
Wall_Shielding.add_element('C',0.817, percent_type='wo')
Wall_Shielding.add_element('B', 0.05, percent_type='wo')
Wall_Shielding.temperature = AmbientT

# air fill:
Air = openmc.Material(name="Air")
Air.set_density("g/cm3", 0.001225)
Air.add_element("N", 0.78)
Air.add_element("O", 0.21)
Air.add_element("Ar", 0.01)
Air.temperature = AmbientT

# Air alternate: (water fill quick)
Air1 = openmc.Material(name="Air")
Air1.set_density("g/cm3",1)
Air1.add_element("H",0.66)
Air1.add_element("O",0.34)
Air1.add_s_alpha_beta('c_H_in_H2O')
Air1.temperature = AmbientT

# Void fill (for within tokamak)
Void = openmc.Material(name="Void")
Void.set_density('g/cm3',0.1)
Void.add_element('H',1) # Arbitrary content, in reality there's a lot of plasma in there, it would be interesting to ask
# Igor and Chris about neutron-plasma interactions given the birth energies being so high.

# CR39 compound for OSL chips (uses free atom cross sections)
# Disregard global definition of this material, it's not meant to be "global, but was convenient designation for my typing" - JMD
CR39 = openmc.Material(name="CR39")
CR39.set_density("g/cm3", 1.31)
# CR39 element proportions based on single polymer structure (repeated to form the actual chain)
CR39.add_element("C", 0.324)
CR39.add_element("H", 0.486)
CR39.add_element("O", 0.190) # 12 + 18 + 7 = 37, 12/37 = 0.324, 18/37 = 0.486, remainder =  0,190
CR39.temperature = AmbientT

# Graphite Armor Tiling:
GraphiteArmor = openmc.Material(name="GraphiteArmor")
GraphiteArmor.set_density("g/cm3",2.2)
GraphiteArmor.add_element("C",1)
GraphiteArmor.add_s_alpha_beta('c_Graphite')
GraphiteArmor.temperature = 340 # Upped, but consider expanding this wildly because these are-plasma facing after all.

# Anti-Torque Plate (lower)
LowATM = openmc.Material(name="LowATM")
LowATM.set_density('g/cm3', 7.85)
LowATM.add_element('C', 0.0025, percent_type='wo')   # 0.25% (taking max)
LowATM.add_element('Mn', 0.0120, percent_type='wo')  # 1.20%
LowATM.add_element('Si', 0.0060, percent_type='wo')  # 0.60%
LowATM.add_element('P', 0.00035, percent_type='wo')  # 0.035% (taking max)
LowATM.add_element('S', 0.00035, percent_type='wo')  # 0.035% (taking max)
LowATM.add_element('Ni', 0.0025, percent_type='wo')  # 0.25%
LowATM.add_element('Cr', 0.0025, percent_type='wo')  # 0.25%
LowATM.add_element('Mo', 0.0008, percent_type='wo')  # 0.08%
LowATM.add_element('Cu', 0.0035, percent_type='wo')  # 0.35%
LowATM.add_element('V', 0.0005, percent_type='wo')   # 0.05%
LowATM.add_element('B', 0.00003, percent_type='wo')  # 0.003% (taking max)
LowATM.add_element('Fe', 0.97537, percent_type='wo')
LowATM.temperature = AmbientT

# Anti-Torque Plate (Upper)
HighATM = openmc.Material(name="HighATM")
HighATM.set_density('g/cm3', 7.85)
HighATM.add_element('C', 0.0025, percent_type='wo')   # 0.25% (taking max)
HighATM.add_element('Mn', 0.0120, percent_type='wo')  # 1.20%
HighATM.add_element('Si', 0.0060, percent_type='wo')  # 0.60%
HighATM.add_element('P', 0.00035, percent_type='wo')  # 0.035% (taking max)
HighATM.add_element('S', 0.00035, percent_type='wo')  # 0.035% (taking max)
HighATM.add_element('Ni', 0.0025, percent_type='wo')  # 0.25%
HighATM.add_element('Cr', 0.0025, percent_type='wo')  # 0.25%
HighATM.add_element('Mo', 0.0008, percent_type='wo')  # 0.08%
HighATM.add_element('Cu', 0.0035, percent_type='wo')  # 0.35%
HighATM.add_element('V', 0.0005, percent_type='wo')   # 0.05%
HighATM.add_element('B', 0.00003, percent_type='wo')  # 0.003% (taking max)
HighATM.add_element('Fe', 0.97537, percent_type='wo')
HighATM.temperature = AmbientT

mats = openmc.Materials([Wall_Concrete,Structural_Steel,Inconel,Roof_Shielding,Wall_Shielding,Air,Void,CR39, pfc, tfc, nbl, GraphiteArmor, LowATM, HighATM])
# bootstrap modification by re-naming Wall_Concrete in the file, the aggregate doesn't take on the same name on export for some reason as requested, this should work though. 
mats.cross_sections = "/bin/endfb-viii.0-hdf5/cross_sections.xml" 
# SWITCH CROSS SECTION FILE TO LOCAL FILE


# <<<<<<<<<<<< MANUALLY SET TO CARBON 12 >>>>>>>>>>>> IF: Your file spits a "No Nuclide found with C0 on X database"

# mats.export_to_xml() # <<<<<<<<<<<>>>>>>>>>>>> <-- UN-COMMENT THIS IF YOU WANT TO RE-EXPORT AND CHANGE NEW MATERIAL FILE

# <<<<<<<<<<<< MANUALLY SET TO CARBON 12 >>>>>>>>>>>> IF: Your file spits a "No Nuclide found with C0 on X database"


# === GEOMETRY DEFINITIONS ===

# import the dagmc file after placed in here
Stage_0_Model = openmc.DAGMCUniverse(filename = '/bin/431_GA_Stage_2.7.h5m', auto_geom_ids = True)
room = openmc.Cell(name = 'Machine_Hall')
room.fill = Stage_0_Model

rootuni = openmc.Universe()
# rootuni.add_cell(Boundcell)
rootuni.add_cell(room)
root = openmc.Cell(name='root_uni')
root.fill = rootuni
geometry = openmc.Geometry([root]) 
geometry.export_to_xml()

# Note: This is probably a redundant organization for cell structures.


# === SOURCE SPECIFICATIONS ===

# Parameters for source setup:
R_major        = 187   # cm — plasma major radius
ring_radius_cm = R_major
z_plasma       = 0 # 119.4   # cm — tokamak Z offset in CAD coordinates
HSz = 35 # cm offset for high source
LSz = -35 # cm offset for low source

# ========== Manual Neutron Source addition ========== 
# based on math derived from the Fusion Neutronics Group - They reference a paper in their source code for setting up these ring sources.
# The Z offset has been benchmarked against this, the purpose of this source method is to format it in a
# method such that OpenMC can be given many of these sources and to manually change the relative strength and more
# realistically pick from a larger number of spatial points to generate neutrons (rather than a discrete source)

# Manual neutron source placement:

# NOTE: DESIGN KNOB: Choose Z_offset values to place sources within.
Zoff = [LSz,z_plasma,HSz] # list of Z offsets, this is a triple-ring source method, More of these can be added.
Roff = [160,187.5,160] # MUST MATCH Z_OFF ARRAY LENGTH, that is to say, 1:1 corresponding for each entry for clean loop implementation. 

Zoffa = [] # Empty array(s) to append with information
Roffa = []
for i in range(len(Zoff)): 
    Zoffa_1 = openmc.stats.Discrete([Zoff[i]],[1])
    Zoffa.append(Zoffa_1)
    Roffa_1 = openmc.stats.Discrete([Roff[i]],[1])
    Roffa.append(Roffa_1)
# ^ turns the input Zoff and R_offset values into discrete components. 
# I wonder if the neutral ion distribution for emission is really like a ring
# I also wonder if you could use this as a tracking system for D-D fusion reactions.

# construciton from: 
# https://fusion-energy.github.io/neutronics-workshop/tasks/task_04_make_sources/2_ring_source.html 
# Just ask for neutron generation in a ring-like pattern, not the source placement from the workshop: 

radius = openmc.stats.Discrete([187], [1])
z_values = openmc.stats.Discrete([Zoff[2]], [1])
angle = openmc.stats.Uniform(a=0., b=2 * math.pi)

# ========== TRS construction area =========== 
# This section contains many of the parameters fed to OpenMC for creating the ring source that the fusion neutronics group developed, most are 1:1 copies
# though some of the parameters have been changed to pure averages for this fluency estimate.

ion_temperature_kev = 20
a_1 = 4.69515
a_2 = -0.040729
a_3 = 0.47
a_4 = 0.81844

mean_delta = (
        a_1
        * ion_temperature_kev ** (2.0 / 3.0)
        / (1.0 + a_2 * ion_temperature_kev**a_3)
        + a_4 * ion_temperature_kev
    )
mean_delta *= 1e3
mean = 2.4495e6
neutron_energy_mean = mean + mean_delta

w_0 = 82.542
a_11 = 1.7013e-3
a_21 = 0.16888
a_31 = 0.49
a_41 = 7.9460e-4

delta = (
        a_11
        * ion_temperature_kev ** (2.0 / 3.0)
        / (1.0 + a_21 * ion_temperature_kev**a_31)
        + a_41 * ion_temperature_kev
    )

    # 2.3548200450309493 on the line below comes from equation 2* math.sqrt(math.log(2)*2) # from fusion source code.
variance = ((w_0 * (1 + delta)) ** 2 * ion_temperature_kev) / 2.3548200450309493**2
variance *= 1e6  # converting keV^2 back to eV^2
std_dev = np.sqrt(variance)
dd_source = openmc.stats.Normal(mean_value=neutron_energy_mean, std_dev=std_dev)

# ========== End TRS construction area =========== (MANUAL EXAMPLE)

HSzo = openmc.stats.Discrete([Zoff[2]], [1])
LSzo = openmc.stats.Discrete([Zoff[0]], [1]) # uses discrete inputs from Zoff collection

trsH = openmc.IndependentSource(
    space = openmc.stats.CylindricalIndependent(r=radius, phi=angle, z=HSzo, origin=(0.0, 0.0, 0.0)),
    angle = openmc.stats.Isotropic(),
    energy = dd_source,
    strength = 0.2
) # trs = test radial source

trsL = openmc.IndependentSource(
    space = openmc.stats.CylindricalIndependent(r=radius, phi=angle, z=LSzo, origin=(0.0, 0.0, 0.0)),
    angle = openmc.stats.Isotropic(),
    energy = dd_source,
    strength = 0.2
) # trs = test radial source

sources = [trsH, trsL]

universal_strength_modifier = 0.02 # should not exceed 1, sets relative particle generation points.
    # in short, a design knob for us to tune our model to. 

sources_main = []

for offset in range(len(Zoff)): 
    TRisrc = openmc.IndependentSource(
        space = openmc.stats.CylindricalIndependent(r=Roffa[offset],
                                                    phi=angle,
                                                    z=Zoffa[offset], 
                                                    origin=(0.0, 0.0, 0.0)),
        angle = openmc.stats.Isotropic(),
        energy = dd_source,
        strength = universal_strength_modifier
    )
    sources_main.append(TRisrc) 

# === SETTINGS ===
settings = openmc.Settings()
settings.source = sources_main # From main sources, uses input of Zoff and Roffa components.

No_Particles = 12500 # set for 100,000

settings.batches = 20 # Removed inactive batches.
settings.inactive = 10
settings.particles = No_Particles # changed to use the No_Particles variable so it can be called upon elsewhere.
settings.run_mode = 'fixed source'
settings.export_to_xml()
settings.output = {"tallies": False} # stops tally output for disk savings
# NOTE: THIS SETTINGS OUTPUT CHANGE DOESN'T SEEM TO WORK IN LOCAL DOCKER ENVIRONMENTS, probably doing something wrong - JMD


# === VISUALIZATION AND MESH CONSTRUCTION ===
# NOTE: A lot of redundant variables in here, be careful deleting some of them or consolidating the values.

VisMesh = 100 # Controls both the visualization tallies and the dose calculations.
squarecoord = VisMesh
VisZcoord2 = 592

mesh = openmc.RegularMesh(name = 'regmesh1',) # I want a mesh 
mesh.dimension = (squarecoord,squarecoord,1) # mesh cells in a given direction, position of 1 changes area.
squarecoord2 = 1500 # setting this to anything different zooms-in the plot.  
mesh.lower_left = (-squarecoord2,-squarecoord2,-VisZcoord2)
# X, Y, Z inputs must be lower than the upper right coordinates.
mesh.upper_right = (squarecoord2,squarecoord2,VisZcoord2) #square mesh with rel. arbitrary dimensions
# ref: https://docs.openmc.org/en/stable/pythonapi/generated/openmc.RegularMesh.html

# Standard settings: 
# Mesh SQ1: 1250
# Mesh dim = squarecoord,1,squarecoord
# mesh SQ2: 1500 
# Mesh LL = -squarecoord2,-squarecoord2,-1000
# Mesh UR = squarecoord2,squarecoord2,1000

mesh_filter = openmc.MeshFilter(mesh)
# Tally: - OLD CARRYOVER FROM 1.0 STAGE

ntally_fil = openmc.CellFilter([1])
# ntally_fil = openmc.MaterialFilter(Inconel)
ntally = openmc.Tally()
ntally.filters = [mesh_filter,ntally_fil]
ntally.scores = ['flux']



# ===== AG Mesh and Tallies: =====

mesh2xy = VisMesh # originally 40
mesh2z = VisMesh # 120 # sets the z component, originally 20 for a 40x40x20 dim.
mesh2 = openmc.RegularMesh(name='regmesh2')
mesh2.dimension   = (mesh2xy, mesh2xy, mesh2z) # 40 in X, 40 in Y, 20 in Z

sq3 = 1500  # manual overwrite for different dimensions if needed 
VZcoord3 = 1500 # manual overwrite for different dimensions if needed.
# NOTE: CHANGE VZCOORD 3 TO VISMESH FOR COMPLETELY CUBIC MESH GEOMETRY.  
mesh2.lower_left  = (-sq3, -sq3, -VZcoord3)
mesh2.upper_right = ( sq3,  sq3,  VZcoord3)
mesh_filter2 = openmc.MeshFilter(mesh2)

# calculate volume to correct source particle interaction
meshvol = mesh2.volumes[0][0][0]

# 3D flux tally over full mesh
flux_tallyAG = openmc.Tally(name="flux_3d")
flux_tallyAG.filters = [mesh_filter2]
flux_tallyAG.scores  = ["flux"]

energy_bins_n, dose_coeffs_n = openmc.data.dose_coefficients(
    particle='neutron',
    geometry='ISO'
)

energy_function_filter_n = openmc.EnergyFunctionFilter(
    energy=energy_bins_n,
    y=dose_coeffs_n,
    interpolation="cubic" ) # cubic interpolation is recommended by ICRP

neutron_particle_filter = openmc.ParticleFilter("neutron")
# adapted from fusion neutronics workshop ^

dose_tally = openmc.Tally(name="dose_3d")
dose_tally.filters = [mesh_filter2,
                      neutron_particle_filter, energy_function_filter_n] 
                      # energy_function_filter_n] # adding in ntally fill can work but kinda breaks dose?
dose_tally.scores  = ["flux"]   # flux x energy function filter = dose
# scored flux is particles-cm / source particles

# ===== end AG Tallies =====

# sq3 is(1500)
# VZcoord3 is (592) # for internal room dimensions!, for cubic geometry, changed to 1500.
# Vismesh is the geometry 
# mesh2xy is the dimension of the mesh in x-y coordinates (symmetric)
# mesh2z is the dimension of the mesh in z coordinates (can vary but we need to handle this case.)



# Get the relative coordinates for defining the upper and lower boundaries: 
def get_rel_coord(MeshX,MeshY,MeshZ,Xm,Ym,Zm,Coord): 
    # "N"m is the actual coordinate dimension normally passed 
    if len(Coord) > 3 or len(Coord) < 3:
        raise TypeError("Insufficient Dimension Passed to get_rel_coord command.")
        # ensure the coordinate is a 3D coordinate. 
    Xspan = 2*Xm 
    Yspan = 2*Ym
    Zspan = 2*Zm
    # ^ gets the total span of the dimensioned room 
    X_psp = Xspan/MeshX
    Y_psp = Yspan/MeshY
    Z_psp = Zspan/MeshZ
    # ^ gets the INDIVIDUAL SIDE LENGTHS of a single component
        # Coord has [X,Y,Z] arguments. 
        # if you want the respective bounding boxes, just multiple or add based on Psp value which will give you a coordinate to work with. 
    llx = (Coord[0]-X_psp)*0.5
    lly = (Coord[1]-Y_psp)*0.5
    llz = (Coord[2]-Z_psp)*0.5
    urx = (Coord[0]+X_psp)*0.5
    ury = (Coord[1]+Y_psp)*0.5
    urz = (Coord[2]+Z_psp)*0.5
    # ^ get individual coordinates
    # pass to a single component: 
    dimensions = [llx,lly,llz,urx,ury,urz]
    return dimensions 
# NOTE: I don't know what the rel_coord function does anymore, it was probably for an older version of the  code - JMD 6/13/2026


esbin = 10 # number of requested bins (last excluded) CANNOT exceed 10 because of logarihmic space. 
energy_bins   = np.logspace(3, 7, esbin+1)   # 1 keV to 10 MeV
print(energy_bins)
print(np.log10(energy_bins))
energy_filter = openmc.EnergyFilter(energy_bins)
spectrum_tally = openmc.Tally(name="energy_spectrum")
spectrum_tally.filters = [energy_filter,mesh_filter2] # throws mesh-filter into the mix, directly to spectrum. 
spectrum_tally.scores  = ["flux"]

# ===== Tally export: =====
tallies = openmc.Tallies([ntally, dose_tally,spectrum_tally]) # removed spectrum tally
# ntally for visualization
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


# <<<<<<<<<<<<<<<<<<<< RUN CONDITIONS >>>>>>>>>>>>>>>>>>>>>


# openmc.plot_geometry()
openmc.run()


# <<<<<<<<<<<<<<<<<<<< RUN CONDITIONS >>>>>>>>>>>>>>>>>>>>> 

# ===== AG Post Processing: =====
# NORMALIZATION — 14AUG2025 shot data

# Neutron totals from neutron cameras (shots 205034-205050)
#
# benchmark = 'all'   -> S_n = sum of all shots (~1.3e17)
#                        use this when comparing against Landauer badge totals
#                        (badges integrated dose over the full day)
# benchmark = <int>   -> S_n = single shot fluence
#                        use when comparing against a per-shot measurement
# ============================================================
shotnum_neutrons = {
    205034: 2.7e15,
    205035: 9.6e15,
    205036: 1.0e16,
    205037: 1.0e16,
    205038: 1.1e16,  # highest flux shot
    205039: 1.1e16,
    205042: 1.0e16,
    205043: 9.0e15,
    205044: 9.5e15,
    205045: 8.1e15,
    205046: 1.1e16,
    205047: 4.1e15,
    205048: 9.9e15,
    205049: 7.4e15,
    205050: 7.7e15,
}

benchmark = 'all'   # if you want a single shot, just put the number plain, not as a string or nada

if benchmark == 'all':
    S_n = sum(shotnum_neutrons.values())
    print(f"Benchmarking against full Aug 14 session: S_n = {S_n:.2e} n")
else:
    S_n = shotnum_neutrons[benchmark]
    print(f"Benchmarking against shot {benchmark}: S_n = {S_n:.2e} n")

pSv_to_mrem = 1*(10**(-7))   # 1 pSv = 1e-12 Sv = 1e-7 mrem # doubled to 1*10^-14

# POST-PROCESSING: 

# <<< STATE POINT DEF>>>
sp = openmc.StatePoint("statepoint.20.h5") 
# <<< STATE POINT DEF>>>

dose_data    = sp.get_tally(name="dose_3d")
#dose_3d_pSv  = dose_data.mean.reshape(nx, ny, nz)
dose_3d_pSv_slice = dose_data.get_slice(scores=['flux'])
dose_3d_pSv = dose_3d_pSv_slice.get_reshaped_data(expand_dims=True, value = 'mean')  # pSv/source-n
flub_factor = 1 # to be used for low particle simulation ONLY, set to 1 for lower results.
dose_mrem_3d = dose_3d_pSv*S_n*(pSv_to_mrem)/(meshvol)     # mrem integrated
# dose_mrem_3d is in psv-cm^3 / neutron source particle, multiplying by the 10^-7 factor and then dividing by mesh vol gives a result
# the product of pSv_to_mrem * mesh_vol needs to be on the order of 10^-14 for reasonable (all particle) results that we expect, so the mesh needs to be re-shaped to be smaller

# dose_mrem_3d: 
#   physical dose data in PicoSieverts times number of neutrons emitted that WE feed to it,
# times multiplier, divided by meshvolume. 

# dose_3d_psv intakes flux units and psv units: 
#   [pSv*cm^2]*[particle-cm/source_particle]
# multiply by # source_particles (SETTINGS PARTICLES)
    # [pSv*cm^2]*[particle-cm/source_particle]*[source_particles] = 
    # [pSv*cm^2]*[particle-cm]
    # [pSv*cm^3]*[particle]
# multiply by # of actual reference particles to get the total amount of neutrons simualted: 
    # [pSv*cm^3]*[(simulated)particle]*[Real_particles/1 simulated particle]
    # CORRECTION: particle -> source particle would be just "one simulated particle represents X real particles"
    # this would cancel with itself and would go away, so we just multiply by S_n (goodness I hope)
    # pSv*cm^3 / mesh volume = pSv
    # pSv to mrem by 1*10^-7, which gives us our results. 
# needs additional conversion to get rid of cm^2 which we get from flux units: 



fsd = sp.get_tally(name="energy_spectrum") # fsd = fluency spectrum data
fsd_slice = fsd.get_slice(scores=['flux'])
fsd_r = fsd_slice.get_reshaped_data(expand_dims=True, value ='mean') # returns slice plot flux values at each 
h = range(esbin)
"""
print('\n')
print('esbin range = ',h)
print(fsd_r[0])
print(fsd_r[9])
print('\n')
"""
def get_energy_data(xyz_cm,
                       ll=(-squarecoord2, -squarecoord2, -VisZcoord2),
                       ur=( squarecoord2,  squarecoord2,  VisZcoord2), 
                       dims=(mesh2xy, mesh2xy, mesh2z)): # consider changing these to sq3 and VisCoord3
    # ur, ll squarecoord2 and VisZcoord 2 used on vsize calculation.
    vsize = (np.array(ur) - np.array(ll)) / np.array(dims) 
    idx   = ((np.array(xyz_cm) - np.array(ll)) / vsize).astype(int) 
    idx   = np.clip(idx, 0, np.array(dims)-1) # clips dataset to get to closest mesh 
    testidx = fsd_r[0] 
    testidx2 = testidx[idx[0]-1,idx[1]-1,idx[2]-1] # diagnostic variable, not references in final version.
    edata = []
    for bin in range(esbin):
        workdata = fsd_r[bin] # this accesses the specific energy cell, then the index pulls the counts in 
        # each energy bin out. # but fsd_r[bin] is saying that, no, there are only 10 bins.
        # so something's not right with the fsd_data, is it the re-shaping? 
        target_data = [workdata[idx[0]-1,idx[1]-1,idx[2]-1]]
        # print(target_data) <- diagnostic to confirm against tallies.
        edata.append(target_data)
    hda = edata # I don't know why I chose hda, its a place holder name, it holds the energy data.
    final_data = [hda,idx]
    return final_data 



def get_mesh_dose_mrem(xyz_cm,
                       ll=(-squarecoord2, -squarecoord2, -VisZcoord2),
                       ur=( squarecoord2,  squarecoord2,  VisZcoord2), 
                       dims=(mesh2xy, mesh2xy, mesh2z)):
    # ur, ll squarecoord2 and VisZcoord 2 used on vsize calculation.
    """Return ICRP-116 H*(10) dose (mrem) at a point in CAD coords (cm)."""
    vsize = (np.array(ur) - np.array(ll)) / np.array(dims)
    # [dx,dy,dz]/[dimx,dimy,dimz] 
    idx   = ((np.array(xyz_cm) - np.array(ll)) / vsize).astype(int) 
    # subtracts the coordinates 
    idx   = np.clip(idx, 0, np.array(dims) - 1) # clips dataset to get to closest mesh
    print(idx) # DELETEME  
    testidx = dose_mrem_3d[idx[0], idx[1], idx[2]]
    # print('testidx = ', testidx) #prints in results the actual 4 dimension nested thingy.
    return dose_mrem_3d[idx[0], idx[1], idx[2]]

# ============================================================
# C/E BENCHMARKING vs Aug 14 OSL/TLD (Landauer) Data
# ============================================================
# Detector positions (x, y, z) in CAD cm — estimated from floor plan (slide 2)
# Z ~ 0 cm = lower level floor (tokamak midplane is at z=119.4 cm)
# *** REPLACE with surveyed CAD coordinates when available ***
#
# Landauer fast neutron (N-F) DDE from Q3 2025 report (mrem), period shown.
# Benchmarking FN only — photon component is activation gammas, not simulated.
# ============================================================
# test positions: 
# (0,0,0) xyz
rs = 400

# DETECTORS 1: Original diagnostic version, used to create multiple position reports about the origin of the room
detectors1 = {
    # Name                   xyz (cm) (tbd)            FN_mrem(from report)   Notes
    '==01== x,0,0'  :          ((rs,0,0),     1,      'Transrex area, far wall'),
    '==02== 0,y,0'  :          ((0,rs,0),     1,    'No neutron signal — below MDL'),
    '==03== 0,0,z'  :          ((0,0,rs),     1,      'Main entrance, well shielded'),
    '==04== 1,1,1'  :          ((rs,rs,rs),   1,     'Shielded'),
    '==05== x,y,0'  :          ((rs,rs,0),    1,    'LOS'),
    '==06== -y,-x,0':          ((-rs, -rs,0), 1,    'Outside Deyong shield box'),
    '==07== 0,y,z'  :          ((0,rs,rs),    1,    'Inside Deyong box — 2374 (2380) extrapolated)'),
    '==08== 0,-y,-z':          ((0,-rs,-rs),  1,    '150V+1 near vessel'),
    '==09== x,0,z'  :          ((rs,0,rs),    1,   '150V-1 near vessel — highest FN'),
    '==10== -x,0,-z':          ((-rs,0,-rs),  1,   'Second highest FN location'),
}

# DETECTORS: Diagnostic value 
detectors = {
    # Name                   xyz (cm) (tbd)            FN_mrem(from report)   Notes
    'Locator':          ((0,     0,     0),    10000,      'CenterEnviro'),
}

# DETECTORS 3: An initial diagnostic version
detectors3 = {
    # Name                   xyz (cm) (tbd)            FN_mrem(from report)   Notes
    'EXP_6_Transrex':          ((-0,     1140,     0),    30,      'Transrex area, far wall - XZ centered'),
    'EXP_7_Transrex':          ((-457.2, 1140,     0),    None,    'No neutron signal — below MDL - E Put aligned'),
    'EXP_8_MainEntrance':      ((-914.4,  1183,     0),    30,      'Main entrance, well shielded'),
    'EXP_9_NorthDoor_Non_LOS': ((-874.4,  -940,     00),    860,     'Lower N door WITH line-of-sight to vessel'),
    'EXP_10_NorthDoor_LOS':    ((-582.9,  -1040,     00),    8910,    'Lower N door WITHOUT vessel view — verify vs EXP_9'),
    'EXP_11_DeyongOut':        ((-557.2, -487.6,    80),    5060,    'Outside Deyong shield box'),
    'EXP_12_DeyongIn':         ((-557.2, -487.6,    40),    None,    'Inside Deyong box — FN exceeded detector cap'),
    'EXP_13_150V+1':           (( 50,      100,   200),    6620,    '150V+1 near vessel'),
    'EXP_14_150V-1':           (( 50,       50,  -200),    19540,   '150V-1 near vessel — highest FN'),
    'EXP_15':                  (( 200,   557.2,     0),    18210,   'Second highest FN location'),
}
# USE CUBIC 10 VISMESH DIMENSION FOR FILE VERIFICATION. + NP.LOG DATA ON FUNCTION. 

# DETECTORS 3: An inversion of Detectors 3, which trades X and Y coordinates between entries as we suspected 
# that we were checking different regions. 
detectors4 = {
    # Name                   Y-X-Z (cm) (tbd)            FN_mrem(from report)   Notes
    'EXP_6_Transrex':          ((1140,     0,     0),    30,      'Transrex area, far wall - XZ centered'),
    'EXP_7_Transrex':          ((1140, 457.2,     0),    None,    'No neutron signal — below MDL - E Put aligned'),
    'EXP_8_MainEntrance':      ((-1183,  -914.4,     0),    30,      'Main entrance, well shielded'),
    'EXP_9_NorthDoor_Non_LOS': ((-940,  -874.4,     00),    860,     'Lower N door WITH line-of-sight to vessel'),
    'EXP_10_NorthDoor_LOS':    ((-1040,  -582.9,     00),    8910,    'Lower N door WITHOUT vessel view — verify vs EXP_9'),
    'EXP_11_DeyongOut':        ((-487.6, -557.2,    80),    5060,    'Outside Deyong shield box'),
    'EXP_12_DeyongIn':         ((-487.6, -557.2,    40),    None,    'Inside Deyong box — FN exceeded detector cap'),
    'EXP_13_150V+1':           (( 100,      50,   200),    6620,    '150V+1 near vessel'),
    'EXP_14_150V-1':           (( 50,       50,  -200),    19540,   '150V-1 near vessel — highest FN'),
    'EXP_15':                  (( 557.2,   200,     0),    18210,   'Second highest FN location'),
} # switched X and Y positions to verify correct behavior 

# DETECTORS 5(Backup): A copied version of the original Detectors_5, containes the MOST accurate reference conditions for detectors based on our 
# best visual estimates.
detectors5_backup = {
    # Name                   X-Y-Z (cm) (tbd)            FN_mrem(from report)   Notes
    'EXP_6_Transrex':          ((0,     -1070,     0),    30,      'Transrex area, far wall - XZ centered'),
    'EXP_7_Transrex':          ((457.2, -1070,     0),    None,    'No neutron signal — below MDL - E Put aligned'),
    'EXP_8_MainEntrance':      ((-914.4,  940,     0),    30,      'Main entrance, well shielded'),
    'EXP_9_NorthDoor_Non_LOS': ((-874.4,  940,     0),    860,     'Lower N door WITH line-of-sight to vessel'),
    'EXP_10_NorthDoor_LOS':    ((-682.9,  975,     0),    8910,    'Lower N door WITHOUT vessel view — verify vs EXP_9'),
    'EXP_11_DeyongOut':        ((-557.2, 487.6,    80),    5060,    'Outside Deyong shield box'),
    'EXP_12_DeyongIn':         ((-557.2, 487.6,    40),    None,    'Inside Deyong box — FN exceeded detector cap'),
    'EXP_13_150V+1':           (( 50,      50,   300),    6620,    '150V+1 near vessel'),
    'EXP_14_150V-1':           (( 50,       50,  -300),    19540,   '150V-1 near vessel — highest FN'),
    'EXP_15':                  (( 200,   557.2,     0),    18210,   'Second highest FN location'),
} 
# Detectors 5_backup was a more-or-less exact position reference from center of tokamak in SLDWRKS

# DETECTORS 5: Final iteration, contains final reference positions for analysis, 
# NOTE: EXP_9a through 9i march through different X coordinates (though 9a and 9e are identical)
# to analyze where we have elevated dose reports, we didn't know for certain how far our geometry was off as compared to the real
# OSL/NeuTrak dosimeter positions, so we had to inspect a large amount of cells. 
detectors5 = {
    # Name                   X-Y-Z (cm) (tbd)            FN_mrem(from report)   Notes
    'EXP_6_Transrex':          ((120,     -970,     -60),    30,      'Transrex area, far wall - XZ centered'),
    'EXP_7_Transrex':          ((457.2, -1070,     0),    None,    'No neutron signal — below MDL - E Put aligned'),
    'EXP_8_MainEntrance':      ((-914,  940,     0),    30,      'Main entrance, well shielded'),
    'EXP_9a_NorthDoor_Non_LOS': ((-870,    930,     30),    860,     'Lower N door WITH line-of-sight to vessel'),
    'EXP_9b_NorthDoor_Non_LOS': ((-960,    930,     30),    860,     'Lower N door WITH line-of-sight to vessel'),
    'EXP_9c_NorthDoor_Non_LOS': ((-930,    930,     30),    860,     'Lower N door WITH line-of-sight to vessel'),
    'EXP_9d_NorthDoor_Non_LOS': ((-900,    930,     30),    860,     'Lower N door WITH line-of-sight to vessel'),
    'EXP_9e_NorthDoor_Non_LOS': ((-870,    930,     30),    860,     'Lower N door WITH line-of-sight to vessel'),
    'EXP_9f_NorthDoor_Non_LOS': ((-840,    930,     30),    860,     'Lower N door WITH line-of-sight to vessel'),
    'EXP_9g_NorthDoor_Non_LOS': ((-810,    930,     30),    860,     'Lower N door WITH line-of-sight to vessel'),
    'EXP_9h_NorthDoor_Non_LOS': ((-780,    930,     30),    860,     'Lower N door WITH line-of-sight to vessel'),
    'EXP_9i_NorthDoor_Non_LOS': ((-750,    930,     30),    860,     'Lower N door WITH line-of-sight to vessel'),
    'EXP_10_NorthDoor_LOS':    ((-650,  800,     0),    8910,    'Lower N door WITHOUT vessel view — verify vs EXP_9'),
    'EXP_11_DeyongOut':        ((-757.2, 200,    80),    5060,    'Outside Deyong shield box'),
    'EXP_12_DeyongIn':         ((-557.2, 487.6,    40),    None,    'Inside Deyong box — FN exceeded detector cap'),
    'EXP_13_150V+1':           (( 150,      150,   240),    6620,    '150V+1 near vessel'),
    'EXP_14_150V-1':           (( 70,       70,  -120),    19540,   '150V-1 near vessel — highest FN'),
    'EXP_15':                  (( 200,   600,     0),    18210,   'Second highest FN location'),
}
# Modified detector locations slightly, increased rel. height of Exp13 and 14 and moved LOS closer to wall
# -582.9 - LOS X position
# X+ POINTS SOUTH, Y+ POINTS EAST, 
# exp 9 -874.4

# note, index positions of about 150 are forbidden for some weird reason.

# Neutrak/OSL 12 was extrapolated from OSL/Neutrak 11 measurements (assumed gamma/neutron ration was sufficient.)
# if a cell exists
print("\n" + "=" * 95)
print(f"C/E BENCHMARKING — Aug 14 2025 | S_n={S_n:.2e} n | ICRP-116 ISO Fast Neutron Dose (mrem)")
print("=" * 95)
print(f"  {'Detector':<28} {'xyz (cm)':<26} {'C (mrem)':<12} {'E (mrem)':<12} {'C/E':<7} Flag")
print("-" * 95)

ce_vals = {}
for name, (xyz, E_mrem, note) in detectors5.items():
    if E_mrem is None:
        print(f"  {name:<28} {str(xyz):<26} {'---':<12} {'MDL/CAP':<12} {'---':<7} SKIP -- {note}")
        continue

        # delete this stuff below:
    C_mrem = get_mesh_dose_mrem(xyz)
    C_mrem =  C_mrem[0][0][0][0] # np.log(C_mrem[0][0][0][0])
    if C_mrem < 1e-6:
        flag = "WARNING: C~0 -- low stats or position outside active mesh"
        print(f"  {name:<28} {str(xyz):<26} {C_mrem:<12.2e} {E_mrem:<12} {'---':<7} {flag}")
    else:
        ratio = C_mrem / E_mrem
        ce_vals[name] = ratio
        if   ratio < 0.1:  flag = "FAIL  C << E  (>10x under)"
        elif ratio < 0.5:  flag = "WARN  Under by >2x"
        elif ratio <= 1.3: flag = "PASS  Within +/-30%"
        elif ratio <= 2.0: flag = "OK    Within factor of 2"
        else:              flag = "FAIL  C >> E  (>2x over)"
        print(f"  {name:<28} {str(xyz):<26} {C_mrem:<12.2f} {E_mrem:<12} {ratio:<7} {flag}")

print("-" * 95)
if ce_vals:
    r = list(ce_vals.values())
    in_band = sum(0.7 <= v <= 1.3 for v in r)
    print(f"\n  C/E summary:  mean={np.mean(r):.2f}  min={np.min(r):.2f}  "
          f"max={np.max(r):.2f}  |  Within +/-30%: {in_band}/{len(r)}")

#  Old raw flux data plotting:

air_tally = sp.get_tally(scores = ['flux'], id = ntally.id)
air_flux = air_tally.mean 
air_flux.shape = (squarecoord,squarecoord)
# declare shape of the data set by these coordinates 
air_flux = air_flux[::-1,:]

slice_value = int((VisMesh/2))
slice_value = 1 # <- manual override of slice value 
dose_mrem_3dlog = np.log10(dose_mrem_3d)  # dose_mrem_3d <- for just using raw linear input 
dose_plot_data = dose_mrem_3dlog[:,:,slice_value] # change to log or non-log for plotting purposes.
# THE SLICE VALUE YOU USE MUST BE WITHIN THE BOUNDS OF THE MESH SIZE - DEFAULTS TO VISMESH/2 FOR Z CENTER>
# changing the position of the slice index changes the 
# local X-Y-Z component to choose to slice from, props to the neutronics team for their slicing
# funciton, saves a lot of work for us. 
dose_plot_shaped = dose_plot_data.shape = (squarecoord,squarecoord)
dps_ndarray = np.array(dose_plot_data)
air_flux = dps_ndarray # <<<Delete this to get full flux spectrum plotting rather than mrem dose equivalent.

from matplotlib import pyplot as plt
from matplotlib import colormaps

# switch to log scale for plotting or not: 
logset = True 

# NEED TO REQEUST SPECIFIC OSL POSITONS BUT HERE WE GO: AT CENTER OF TOKAMAK 
Diagnostic_position = [0,  0,     0] # CONTROLS VIEWING COORDINATE 
Diagnostic_position1 = [0,  0,     0]
# Diagnostic_position = [56,70,50] # Diagonostic_position = [28, 76, 50]

esd = (get_energy_data(Diagnostic_position)) # the target cell is accessed by the input to this function.
# This also uses the same behavior of mesh access as the regular mesh. 
# the output of g_e_d function is all of the data, with the end being the array access.
cell_dim_array = esd[-1] # accesses the last element of the esd array (numpy array which has raw clipped values from g_e_d.)

slice_valuez = cell_dim_array[2]
slice_valuey = cell_dim_array[1]
slice_valuex = cell_dim_array[0]

# request which plotting data to use 
if logset == True: 
    fxyta = dose_mrem_3dlog[:,:,slice_valuez]
    fyzta = dose_mrem_3dlog[slice_valuex,:,:]
    fxzta = dose_mrem_3dlog[:,slice_valuey,:]
else: 
    fxyta = dose_mrem_3d[:,:,slice_valuez]
    fyzta = dose_mrem_3d[slice_valuex,:,:]
    fxzta = dose_mrem_3d[:,slice_valuey,:]

fxytb = fxyta.shape = (squarecoord,squarecoord)
fyztb = fyzta.shape = (squarecoord,squarecoord)
fxztb = fxzta.shape = (squarecoord,squarecoord)

fxyt = fxyta # np.array(fxytb)
fyzt = fyzta # np.array(fyztb)
fxzt = fxzta # np.array(fxztb)

# NOTE: See older file revisions (2.7, 3.0 for commented "XCONFIRMED/YCONFIRMED" Code)
# Depending on the air-flux of choice the data can provide an entire overlay image.  
# Future version: consider 6 different overlay mesh plots to go along with these positions based on the air_flux tally.

# setup plots: 
fig, (ax1,ax2,ax3,ax4) = plt.subplots(ncols=4,figsize=(30,10)) # creates 2x2 subplots for plotting multiple traces 

# ax5

# create xy at top left corner: 
pfxy = ax1.imshow(fxyt, cmap = 'plasma') 
ax1.set_title('Axis 0-0 fxyt')
ax1.set_xlabel("phY")
ax1.set_ylabel("phX") 
cbarxyt = plt.colorbar(pfxy, ax = ax1, cmap = 'plasma') 
cbarxyt.set_label('Neutron Dose along particular <X-Y> trace [mrem]') 

"""pfxy.set_clim(vmin = 0, vmax = 1e7) # does this work????""" # spoiler: it doesn't and breaks the color.

# create yz at top right corner: 
pfyz = ax2.imshow(fyzt, cmap = 'plasma') 
ax2.set_title('Axis 0-1 fyzt') 
ax2.set_xlabel("phZ")
ax2.set_ylabel("phY") 
cbaryzt = plt.colorbar(pfyz, ax = ax2, cmap = 'plasma') 
cbaryzt.set_label('Neutron dose along particular <Y-Z> trace [mrem]') 

# create xz at bottom left corner: 
pfxz = ax3.imshow(fxzt, cmap = 'plasma') 
ax3.set_title('Axis 1-0 fyzt') 
ax3.set_xlabel("phZ")
ax3.set_ylabel("phX")
cbarxzt = plt.colorbar(pfxz, ax = ax3, cmap = 'plasma')
cbarxzt.set_label('Neutron dose along a particular <X-Z> trace [mrem]') 

# make data for fluency estimates: 
echo = []
for i in range(len(esd[0])): 
    extracted_data = esd[0][i][0][0][0]
    echo.append(extracted_data)

# creates the step plot at bottom right corner: 
ax4.plot(range(esbin), echo) 
ax4.set_title('Neutron energy at particular mesh point of interest') 
ax4.set_xlabel("energy bin [eV]") 
ax4.set_ylabel("relative fluency [particles/cm]") 

# plot_af = ax5.imshow(air_flux) (build for total non-log scale traces)

# save figure data: 
plt.savefig(fname = '/bin/431_GA_Stage_3.0.png')

print('\n') 
print('End of simulation.')
print('\n')
