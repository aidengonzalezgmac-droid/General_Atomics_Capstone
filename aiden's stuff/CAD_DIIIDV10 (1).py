# DIII-D OpenMC + DagMC — V9
# Adds: relative error reporting per detector, MAGIC weight windows (2-pass)


import os
os.environ['OPENMC_CROSS_SECTIONS'] = '/storage/work/ajg7072/Capstone/endf/cross_sections.xml'

import openmc
import numpy as np
import math
from openmc_plasma_source import fusion_ring_source

# ============================================================
# MATERIALS — names must match H5M mat: tags exactly
# Stage 2.7 tags: 
"""
mat: Wall_Shielding
mat:Roof_Shielding
mat:Structural_Steel
mat:Concrete(0.85)-Barite(0.15)
mat:CR39
mat:CR39
mat:Graveyard
mat:Inconel
mat:nbl
mat:tfc
mat:pfc
mat:LowATM
mat:HighATM
mat:GraphiteArmor
No data for these special temps, only 250, 294, 600, 900, 1200, 2500 K
"""

AmbientT = 294 # in Kelvin, # ambient 70 deg Fahrenheit
SupercondT = 250 # in Kelvin # Nuclear data doesn't have CS for copper interactions at superconducting temps.

# PFCs:
pfc = openmc.Material(name='pfc')
pfc.set_density("g/cm3", 1.31)
pfc.add_element("C", 0.02)
pfc.add_element("H", 0.04) # changed hydrgoen content to 0.06 from 0.05 because C0 doesn't have 4 kelvin cross secitons.
pfc.add_element("O", 0.04)
pfc.add_element("Cu",0.9)
pfc.temperature = SupercondT

# TFCs:
tfc = openmc.Material(name="tfc")
tfc.set_density("g/cm3", 1.31)
tfc.add_element("C", 0.02)
tfc.add_element("H", 0.04)
tfc.add_element("O", 0.04)
tfc.add_element("Cu",0.9)
tfc.temperature = SupercondT

# NB: 
nbl = openmc.Material(name="nbl")
nbl.set_density("g/cm3", 1) # base results are 6, lowering for comparisons.
nbl.add_element("Al",0.5)
nbl.add_element("Cu",0.2)
nbl.add_element("Fe",0.1)
nbl.add_element("C",0.1)
nbl.add_element("H",0.1)
nbl.temperature = AmbientT # increased above ambient 294 k to 310 for average operation of equipment


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
Wall_Concrete = openmc.Material.mix_materials(
    [concrete, barite], [0.85, 0.15], 'wo',
    name='Concrete(0.85)-Barite(0.15)'
)
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
Inconel.add_element("Ni", 0.6402)  # 0.5922 + 0.048
Inconel.add_element("Cr", 0.21) # increased from 20.0 to 22
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
Air1 = openmc.Material(name="Air1")
Air1.set_density("g/cm3",1)
Air1.add_element("H",0.66)
Air1.add_element("O",0.34)
Air1.add_s_alpha_beta('c_H_in_H2O')
Air1.temperature = AmbientT

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
CR39.temperature = AmbientT

GraphiteArmor = openmc.Material(name="GraphiteArmor")
GraphiteArmor.set_density("g/cm3",2.2)
GraphiteArmor.add_element("C",1)
GraphiteArmor.add_s_alpha_beta('c_Graphite')
GraphiteArmor.temperature = AmbientT #No data for these special temps, only 250, 294, 600, 900, 1200, 2500 K

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

"""                       
Graveyard = openmc.Material(name="Graveyard")
Graveyard.set_density('g/cm3', 1e-10)
Graveyard.add_element('H', 1.0)
"""                     
mats = openmc.Materials([Wall_Concrete,Structural_Steel,Inconel,Roof_Shielding,Wall_Shielding,Air,Void,CR39, pfc, tfc, nbl, GraphiteArmor, LowATM, HighATM])
mats.export_to_xml()
# ============================================================
# GEOMETRY
# ============================================================
Stage_3_Model = openmc.DAGMCUniverse(
    filename='/storage/work/ajg7072/Capstone/Ring_source/cad_modeling/431_GA_Stage_2.7.h5m',
    auto_geom_ids=True,
)
geometry = openmc.Geometry(Stage_3_Model)
geometry.export_to_xml()

# ============================================================
# SOURCE
# Z centroid of Stage 2.7 H5M = 24.4 cm (verified via h5py query)
# fusion_ring_source takes METRES — divide cm by 100
# R_major = 187 cm = 1.87 m
# ============================================================
R_major  = 187    # cm
z_plasma = 24.4   # cm — updated for Stage 2.7 (was 119.4 for Stage 0.6)

my_source = fusion_ring_source(
    radius=R_major / 100,        # cm -> m
    angles=(0.0, 2 * math.pi),
    z_placement=z_plasma / 100,  # cm -> m
    temperature=20000.0,
    fuel={"D": 1.0},
)


# SETTINGS — Pass 1 (no weight windows, used to generate WW map)
# For Pass 2 (production), weight windows are loaded below
# Toggle PASS variable to switch between runs

PASS = 2   # 1 = generate weight windows,  2 = production run with WW

settings = openmc.Settings()
settings.source = my_source
settings.run_mode = 'fixed source'
print(my_source) 

if PASS == 1:
    settings.batches   = 20
    settings.particles = 50000
    print("PASS 1: Generating weight windows via MAGIC method")
else:
    settings.batches   = 20
    settings.particles = 100000
    wws = openmc.WeightWindowsList.from_hdf5('weight_windows.h5')
    settings.weight_windows = wws
    settings.weight_windows_on = True
    settings.weight_window_checkpoints = {'collision': True, 'surface': True}
    settings.survival_biasing = False
    print("PASS 2: Production run with MAGIC weight windows loaded")

#settings.max_splits = 10  # prevent runaway splitting
settings.max_lost_particles = 1000
settings.rel_max_lost_particles = 0.05
settings.max_write_lost_particles = 1
# MESH AND TALLIES
# Bounds set to Stage 2.7 H5M extents: X/Y ±1550 cm, Z ±950 cm
VisMesh    = 100
squarecoord2 = 1500
VisZcoord2   = 592

# Define weight window spatial mesh
ww_mesh = openmc.RegularMesh()
ww_mesh.dimension = (VisMesh, VisMesh, 1)
ww_mesh.lower_left = (-squarecoord2, -squarecoord2, -VisZcoord2)
ww_mesh.upper_right = (squarecoord2,  squarecoord2,  VisZcoord2)

# Create weight window object and adjust parameters
if PASS == 1:
    wwg = openmc.WeightWindowGenerator(
        method='magic',
        mesh=ww_mesh,
        max_realizations=settings.batches,
        on_the_fly=False
    )
    settings.weight_window_generators = wwg


# mesh1 — 2D slice for quick visualization (ntally)
mesh = openmc.RegularMesh(name='regmesh1')
mesh.dimension   = (VisMesh, VisMesh, 1)
mesh.lower_left  = (-squarecoord2, -squarecoord2, -VisZcoord2)
mesh.upper_right = ( squarecoord2,  squarecoord2,  VisZcoord2)
mesh_filter = openmc.MeshFilter(mesh)

ntally = openmc.Tally()
ntally.filters = [mesh_filter]
ntally.scores  = ['flux']

# mesh2 — 3D dose/flux mesh
mesh2xy = VisMesh
mesh2z  = VisMesh
mesh2 = openmc.RegularMesh(name='regmesh2')
mesh2.dimension   = (mesh2xy, mesh2xy, mesh2z)
mesh2.lower_left  = (-squarecoord2, -squarecoord2, -VisZcoord2)
mesh2.upper_right = ( squarecoord2,  squarecoord2,  VisZcoord2)
mesh_filter2 = openmc.MeshFilter(mesh2)

vox_x = (2 * squarecoord2) / mesh2xy   # 30 cm
vox_y = (2 * squarecoord2) / mesh2xy   # 30 cm  
vox_z = (2 * VisZcoord2)   / mesh2z    # 11.84 cm
meshvol = vox_x * vox_y * vox_z        # ~10,656 cm³

# 3D flux tally
flux_tally = openmc.Tally(name="flux_3d")
flux_tally.filters = [mesh_filter2]
flux_tally.scores  = ["flux"]

# ICRP-116 dose tally — ISO geometry, cubic interpolation
energy_bins_n, dose_coeffs_n = openmc.data.dose_coefficients(
    particle='neutron',
    geometry='ISO',
    data_source='icrp116',
)
dose_energy_filter = openmc.EnergyFunctionFilter(
    energy=energy_bins_n,
    y=dose_coeffs_n,
    interpolation="cubic",
)
neutron_particle_filter = openmc.ParticleFilter("neutron")

dose_tally = openmc.Tally(name="dose_3d")
dose_tally.filters = [mesh_filter2, neutron_particle_filter, dose_energy_filter]
dose_tally.scores  = ["flux"]

tallies = openmc.Tallies([ntally, flux_tally, dose_tally])
tallies.export_to_xml()
settings.export_to_xml()
############# RUN
# After Pass 1 completes, generate weight windows then re-run with PASS=2
# Only create the auto-generator in Pass 1
openmc.run()

# POST-PROCESSING (runs in both passes, meaningful in Pass 2)
shotnum_neutrons = {
    205034: 2.7e15,
    205035: 9.6e15,
    205036: 1.0e16,
    205037: 1.0e16,
    205038: 1.1e16,
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

benchmark = 'all'
if benchmark == 'all':
    S_n = sum(shotnum_neutrons.values())
    print(f"Benchmarking against full Aug 14 session: S_n = {S_n:.2e} n")
else:
    S_n = shotnum_neutrons[benchmark]
    print(f"Benchmarking against shot {benchmark}: S_n = {S_n:.2e} n")

pSv_to_mrem = 1e-7   # 1 pSv = 1e-12 Sv = 1e-7 mrem

sp = openmc.StatePoint(f"statepoint.{settings.batches}.h5")

nx, ny, nz = mesh2xy, mesh2xy, mesh2z

# --- Flux ---
flux_data  = sp.get_tally(name="flux_3d")
flux_3d    = flux_data.mean.reshape((nx, ny, nz))
flux_err   = flux_data.get_values(
                scores=['flux'], value='rel_err').reshape((nx, ny, nz))

# --- Dose ---
# Units: pSv per source neutron (EnergyFunctionFilter applied during scoring)
# Multiply by S_n -> pSv total; multiply by pSv_to_mrem -> mrem
# NO division by meshvol — OpenMC flux tallies are already track_length/volume
dose_data    = sp.get_tally(name="dose_3d")
dose_3d_pSv  = dose_data.mean.reshape((nx, ny, nz))      # pSv/source-n
dose_mrem_3d = dose_3d_pSv * S_n * pSv_to_mrem/meshvol           # mrem integrated, might need to divide by meshvol
dose_err_3d  = dose_data.get_values(
                scores=['flux'], value='rel_err').reshape((nx, ny, nz))

def get_mesh_idx(xyz_cm,
                 ll=(-squarecoord2, -squarecoord2, -VisZcoord2),
                 ur=( squarecoord2,  squarecoord2,  VisZcoord2),
                 dims=(mesh2xy, mesh2xy, mesh2z)):
    """Return voxel indices for a point in CAD coords (cm)."""
    vsize = (np.array(ur) - np.array(ll)) / np.array(dims)
    idx   = ((np.array(xyz_cm) - np.array(ll)) / vsize).astype(int)
    return np.clip(idx, 0, np.array(dims) - 1)

# ============================================================
# C/E BENCHMARKING
# *** xyz are floor plan estimates — replace with surveyed CAD coords ***
# ============================================================
detectors5 = {
    'EXP_6_Transrex':          ((0,       -1070,   0),  30,    'Transrex area, far wall'),
    'EXP_7_Transrex':          ((457.2,   -1070,   0),  None,  'Below MDL'),
    'EXP_8_MainEntrance':      ((-914.4,   940,    0),  30,    'Main entrance'),
    'EXP_9_NorthDoor_Non_LOS': ((-874.4,   940,    0),  860,   'N door no LOS'),
    'EXP_10_NorthDoor_LOS':    ((-682.9,   975,    0),  8910,  'N door with LOS'),
    'EXP_11_DeyongOut':        ((-557.2,   487.6,  80), 5060,  'Outside Deyong box'),
    'EXP_12_DeyongIn':         ((-557.2,   487.6,  40), None,  'Inside Deyong box — FN cap exceeded'),
    'EXP_13_150V+1':           ((50,       50,     300), 6620, '150V+1 near vessel'),
    'EXP_14_150V-1':           ((50,       50,    -300), 19540,'150V-1 near vessel'),
    'EXP_15':                  ((200,      557.2,   0),  18210,'Second highest FN'),
}

print("\n" + "=" * 105)
print(f"C/E BENCHMARKING — Aug 14 2025 | S_n={S_n:.2e} n | ICRP-116 ISO | PASS {PASS}")
print("=" * 105)
print(f"  {'Detector':<28} {'C (mrem)':<12} {'E (mrem)':<12} {'C/E':<7} {'Flux RelErr':<13} {'Dose RelErr':<13} Flag")
print("-" * 105)

ce_vals = {}
for name, (xyz, E_mrem, note) in detectors5.items():
    if E_mrem is None:
        print(f"  {name:<28} {'---':<12} {'MDL/CAP':<12} {'---':<7} {'---':<13} {'---':<13} SKIP -- {note}")
        continue

    idx = get_mesh_idx(xyz)
    ix, iy, iz = idx

    C_mrem      = float(dose_mrem_3d[ix, iy, iz])
    flux_re     = float(flux_err[ix, iy, iz])
    dose_re     = float(dose_err_3d[ix, iy, iz])

    flux_re_str = f"{flux_re:.1%}" if flux_re > 0 else "---"
    dose_re_str = f"{dose_re:.1%}" if dose_re > 0 else "---"

    if C_mrem < 1e-6:
        flag = "WARNING: C~0 -- low stats or outside mesh"
        print(f"  {name:<28} {C_mrem:<12.2e} {E_mrem:<12} {'---':<7} {flux_re_str:<13} {dose_re_str:<13} {flag}")
    else:
        ratio = C_mrem / E_mrem
        ce_vals[name] = ratio
        if   ratio < 0.1:  flag = "FAIL  C << E  (>10x under)"
        elif ratio < 0.5:  flag = "WARN  Under by >2x"
        elif ratio <= 1.3: flag = "PASS  Within +/-30%"
        elif ratio <= 2.0: flag = "OK    Within factor of 2"
        else:              flag = "FAIL  C >> E  (>2x over)"
        print(f"  {name:<28} {C_mrem:<12.1f} {E_mrem:<12} {ratio:<7.2f} {flux_re_str:<13} {dose_re_str:<13} {flag}")

print("-" * 105)
if ce_vals:
    r = list(ce_vals.values())
    in_band = sum(0.7 <= v <= 1.3 for v in r)
    print(f"\n  C/E summary:  mean={np.mean(r):.2f}  min={np.min(r):.2f}  "
          f"max={np.max(r):.2f}  |  Within +/-30%: {in_band}/{len(r)}")

print("""
  Notes:
    - Detector xyz are FLOOR PLAN ESTIMATES -- replace with surveyed CAD coords
    - Flux/Dose RelErr > 20% means that detector's C value is statistically unreliable
    - EXP 7: below MDL -- excluded; EXP 12: FN cap exceeded -- excluded
    - Photon dose not compared -- activation gammas not modelled
    - C/E targets: +/-30% near-field, factor of 2 in shielded regions
    - Run PASS=1 first to generate weight_windows.h5, then PASS=2 for production
""")
print("=" * 105)
print("ALL DONE")

