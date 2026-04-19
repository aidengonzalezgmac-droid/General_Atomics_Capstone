# Note for V-1.0+: 
# higher particle counts allow higher mesh qualities to be used 

# import modules: 

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

# S2 Materials: 
    # PFC
    # TFC # (D SHAPED)
    # NB 

# S3 Materials: 
    # Anti-Torque (Low and High) -> HighATM, LowATM
    # GraphiteArmor (For internal graphite paneling)

# PFCs:
pfc = openmc.Material(name='pfc')
pfc.set_density("g/cm3", 1.31)
pfc.add_element("C", 0.02)
pfc.add_element("H", 0.04)
pfc.add_element("O", 0.04)
pfc.add_element("Cu",0.9)

# TFCs:
tfc = openmc.Material(name="tfc")
tfc.set_density("g/cm3", 1.31)
tfc.add_element("C", 0.02)
tfc.add_element("H", 0.04)
tfc.add_element("O", 0.04)
tfc.add_element("Cu",0.9)

# NB: 
nbl = openmc.Material(name="nbl")
nbl.set_density("g/cm3", 6)
nbl.add_element("Al",0.5)
nbl.add_element("Cu",0.2)
nbl.add_element("Fe",0.1)
nbl.add_element("C",0.1)
nbl.add_element("H",0.1)


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

# Air alternate: (water fill quick)
Air1 = openmc.Material(name="Air")
Air1.set_density("g/cm3",1)
Air1.add_element("H",0.66)
Air1.add_element("O",0.34)

# void fill (for within tokamak)
Void = openmc.Material(name="Void")
Void.set_density('g/cm3',0.1)
Void.add_element('H',1) # Arbitrary content, in reality there's a lot of plasma in there, it would be interesting to ask
# Igor and Chris about neutron-plasma interactions given the birth energies being so high.

# CR39 compound for OSL chips (uses free atom cross sections)
CR39 = openmc.Material(name="CR39")
CR39.set_density("g/cm3", 1.31)
# CR39 element proportions based on single polymer structure (repeated to form the actual chain)
CR39.add_element("C", 0.324)
CR39.add_element("H", 0.486)
CR39.add_element("O", 0.190) # 12 + 18 + 7 = 37, 12/37 = 0.324, 18/37 = 0.486, remainder =  0,190

GraphiteArmor = openmc.Material(name="GraphiteArmor")
GraphiteArmor.set_density("g/cm3",2.2)
GraphiteArmor.add_element("C",1)

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

mats = openmc.Materials([Wall_Concrete,Structural_Steel,Inconel,Roof_Shielding,Wall_Shielding,Air,Void,CR39, pfc, tfc, nbl, GraphiteArmor, LowATM, HighATM])
# mats.cross_sections = "/bin/endfb-viii.0-hdf5/cross_sections.xml" #please rem to use your own file
# have you considered just *not?*
# references CS for docker version

mats.export_to_xml() # <<<<<<<<<<<>>>>>>>>>>>>

# bootstrap modification by re-naming Wall_Concrete in the file.

# GEOMETRY: 
# import the dagmc file after placed in here
Stage_0_Model = openmc.DAGMCUniverse(filename = '/bin/431_GA_Stage_2.7.h5m', auto_geom_ids = True)
room = openmc.Cell(name = 'Machine_Hall')
room.fill = Stage_0_Model


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

R_major        = 187   # cm — plasma major radius
ring_radius_cm = R_major
z_plasma       = 0 # 119.4   # cm — tokamak Z offset in CAD coordinates
HSz = 35 # cm offset for high source
LSz = -35 # cm offset for low source


# Manual neutron source placement:
Zoff = [LSz,z_plasma,HSz] # list of Z offsets 
Roff = [160,187.5,160] # MUST MATCH Z_OFF LENGTH. 

Zoffa = []
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


# ========== Manual Neutron Source addition based on math derived from the link above and from source code for 
# ring_source. Z offset has been benchmarked against this, the purpose of this source method is to format it in a
# method such that OpenMC can be given many of these sources and to manually change the relative strength.  

radius = openmc.stats.Discrete([187], [1])
z_values = openmc.stats.Discrete([Zoff[2]], [1])
angle = openmc.stats.Uniform(a=0., b=2 * math.pi)

#TRS math:
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

    # 2.3548200450309493 on the line below comes from equation 2* math.sqrt(math.log(2)*2)
variance = ((w_0 * (1 + delta)) ** 2 * ion_temperature_kev) / 2.3548200450309493**2
variance *= 1e6  # converting keV^2 back to eV^2
std_dev = np.sqrt(variance)

dd_source = openmc.stats.Normal(mean_value=neutron_energy_mean, std_dev=std_dev)

"""
trs = openmc.IndependentSource(
    space = openmc.stats.CylindricalIndependent(r=radius, phi=angle, z=z_values, origin=(0.0, 0.0, 0.0)),
    angle = openmc.stats.Isotropic(),
    energy = openmc.stats.muir(e0=2.4495e6, m_rat=4, kt=20000.0)
"""
trs = openmc.IndependentSource(
    space = openmc.stats.CylindricalIndependent(r=radius, phi=angle, z=z_values, origin=(0.0, 0.0, 0.0)),
    angle = openmc.stats.Isotropic(),
    energy = dd_source
) # trs = test radial source


# ========== End TRS construction area ===========

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

universal_strength_modifier = 0.015 # should not exceed 1, sets relative particle generation points.
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


my_source = fusion_ring_source(
    radius=R_major,        # cm → m
    angles=(0.0, 2 * math.pi),
    z_placement= z_plasma,     # cm → m
    temperature=20000.0,
    fuel={"D": 1.0},
)


HSource = fusion_ring_source(
    radius=R_major,        # cm → m
    angles=(0.0, 2 * math.pi),
    z_placement= HSz,     # cm → m
    temperature=20000.0,
    fuel={"D": 1.0}, 
)



LSource = fusion_ring_source(
    radius=R_major,        # cm → m
    angles=(0.0, 2 * math.pi),
    z_placement= LSz,     # cm → m
    temperature=20000.0,
    fuel={"D": 1.0}, 
)



# apply settings:
settings = openmc.Settings()
settings.source = sources_main # HSource # my_source

settings.batches = 20 # oriignally 100 but for testing we shouldnd't need this many.
settings.inactive = 10
settings.particles = 10000
settings.run_mode = 'fixed source'
settings.export_to_xml()
settings.output = {"tallies": False} # stops tally output for disk savings

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
# Tally:

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

# mesh2.lower_left  = (-1500, -1500, -592)
# mesh2.upper_right = ( 1500,  1500,  592)

mesh2.lower_left  = (-1500, -1500, -592)
mesh2.upper_right = ( 1500,  1500,  592)
mesh_filter2 = openmc.MeshFilter(mesh2)

# calculate volume to correct source particle interaction
meshvol = mesh2.volumes[0][0][0]

# 3D flux tally over full mesh
flux_tallyAG = openmc.Tally(name="flux_3d")
flux_tallyAG.filters = [mesh_filter2]
flux_tallyAG.scores  = ["flux"]

# Energy spectrum tally (global, no mesh filter — total flux vs energy)
energy_bins   = np.logspace(3, 7, 200)   # 1 keV to 10 MeV
energy_filter = openmc.EnergyFilter(energy_bins)
spectrum_tally = openmc.Tally(name="energy_spectrum")
spectrum_tally.filters = [energy_filter]
spectrum_tally.scores  = ["flux"]

# ICRP-116 energy-dependent dose tally
# dose_coefficients() returns (energies eV, coefficients pSv*cm^2)
# EnergyFunctionFilter applies H*(10) curve during scoring
# Result units: pSv per source neutron -> convert to mrem in post-processing
# geometry='ISO': isotropic irradiation — correct for area badges that see
#   neutrons from all directions (AP assumes a directed parallel beam)
# W finds 
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

# ===== Tally export: =====
tallies = openmc.Tallies([ntally, dose_tally]) # removed spectrum tally
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

# openmc.plot_geometry()
# openmc.run()


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

 # DELETE THIS WHEN FINISHING THE RIGHT VALUES.
# POST-PROCESSING

# <<< STATE POINT DEF>>>
sp = openmc.StatePoint("statepoint.20.h5") 
# <<< STATE POINT DEF>>>

nx, ny, nz = mesh2xy, mesh2xy, mesh2z   # must match mesh.dimension


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
# 

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
    print('testidx = ', testidx)
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

detectors = {
    # Name                   xyz (cm) (tbd)            FN_mrem(from report)   Notes
    'Locator':          ((0,     0,     0),    10000,      'CenterEnviro'),
}

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

detectors5 = {
    # Name                   X-Y-Z (cm) (tbd)            FN_mrem(from report)   Notes
    'EXP_6_Transrex':          ((0,     -1070,     0),    30,      'Transrex area, far wall - XZ centered'),
    'EXP_7_Transrex':          ((457.2, -1070,     0),    None,    'No neutron signal — below MDL - E Put aligned'),
    'EXP_8_MainEntrance':      ((-914.4,  940,     0),    30,      'Main entrance, well shielded'),
    'EXP_9_NorthDoor_Non_LOS': ((-874.4,  940,     00),    860,     'Lower N door WITH line-of-sight to vessel'),
    'EXP_10_NorthDoor_LOS':    ((-682.9,  975,     00),    8910,    'Lower N door WITHOUT vessel view — verify vs EXP_9'),
    'EXP_11_DeyongOut':        ((-557.2, 487.6,    80),    5060,    'Outside Deyong shield box'),
    'EXP_12_DeyongIn':         ((-557.2, 487.6,    40),    None,    'Inside Deyong box — FN exceeded detector cap'),
    'EXP_13_150V+1':           (( 50,      50,   300),    6620,    '150V+1 near vessel'),
    'EXP_14_150V-1':           (( 50,       50,  -300),    19540,   '150V-1 near vessel — highest FN'),
    'EXP_15':                  (( 200,   557.2,     0),    18210,   'Second highest FN location'),
} # Modified detector locations slightly, increased rel. height of Exp13 and 14 and moved LOS closer to wall
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
    # h = type(C_mrem)
    # q = C_mrem[0][0]
    # print(h)
    # print(q)
    # print(xyz)
    # print(C_mrem)
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

# Post Processing error: 
    # Passing dimensions which are OUTSIDE of the mesh bounding box gives values WITHIN the mesh bounding box,
    # and does NOT clip the values. 
    # The code used to scale and retrieve this data is working on this dataset relative to the dimensions of the plot, 
    # NOT The actual model dimensions. 
    # If we re-scale the mesh bounding box to the room we could cirucmvent this, but this is BAD for analysis.
        # - 4/19/2026
    # For now, just work on adding new source positions. 



#  Old raw flux data plotting:

air_tally = sp.get_tally(scores = ['flux'], id = ntally.id)
air_flux = air_tally.mean 
air_flux.shape = (squarecoord,squarecoord)
# declare shape of the data set by these coordinates 
air_flux = air_flux[::-1,:]


slice_value = int((VisMesh/2))
slice_value = 50 # <- manual override of slice value 
dose_mrem_3dlog = np.log(dose_mrem_3d)  # dose_mrem_3d <- for just using raw linear input 
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
import matplotlib as mpl

fig,ax = plt.subplots(figsize=(20,20))
sub_plot_1 = ax.imshow(air_flux, cmap = 'plasma')
ax.set_xlabel("YCONFIRMED Mesh Coordinates") 
# if using dps_ndarray then the this should be Y (if using the full plot this is )
ax.set_ylabel("XCONFIRMED Mesh Coordinates")
# if using the dps_ndaray then this should be X
cbar = plt.colorbar(sub_plot_1, ax=ax, cmap = 'plasma')
cbar.set_label('Total Neutron Dose [mrem/day]')


"""
sub_plot1 = plt.subplot(121, title='Testing_Plot')
sub_plot1.imshow(air_flux, cmap = 'plasma') # replace with the most recent "flux" variable name.
sub_plot1.set_xlabel("X Mesh Coordinate")
sub_plot1.set_ylabel("Y Mesh Coordinate")

"""

plt.savefig(fname='/bin/431_GA_Stage_0.9.png')
plt.show

print('\n') 
print('End of simulation.')
print('\n')

# measuring space decrease by some volumetric factor, particle counts must be increased accordingly to produce reasonable
# results, which is a fight between precision and measurement ability
# very fine mesh means very low counting in specific areas - you're a nuce engineer, come the fuck on you ought to know that
# collect everything in an area, get a good impact
# the mesh - particle balance is what is problematic. 