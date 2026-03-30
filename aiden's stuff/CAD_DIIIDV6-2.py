"""
DIII-D OpenMC + DagMC neutron transport
D-D ring source, 2.45 MeV
CAD geometry: 431_GA_Stage_0.9.h5m
"""

import os
os.environ['OPENMC_CROSS_SECTIONS'] = '/storage/work/ajg7072/Capstone/endf/cross_sections.xml'

import openmc
import numpy as np
import math
from openmc_plasma_source import fusion_ring_source

# DIII-D CAD coordinate system (from H5M centroid check)
R_major  = 187.5  # cm — plasma major radius
z_plasma = 119.4  # cm — tokamak Z offset in CAD coordinates

# MATERIALS
inconel = openmc.Material(name='Inconel')
inconel.set_density('g/cm3', 8.44)
inconel.add_element('Ni', 0.580, percent_type='wo')
inconel.add_element('Cr', 0.215, percent_type='wo')
inconel.add_element('Mo', 0.090, percent_type='wo')
inconel.add_element('Fe', 0.050, percent_type='wo')
inconel.add_element('Nb', 0.037, percent_type='wo')
inconel.add_element('Co', 0.010, percent_type='wo')
inconel.add_element('Mn', 0.005, percent_type='wo')
inconel.add_element('Si', 0.005, percent_type='wo')
inconel.add_element('Al', 0.004, percent_type='wo')
inconel.add_element('Ti', 0.004, percent_type='wo')

wall_shielding = openmc.Material(name='Wall_Shielding')
wall_shielding.set_density('g/cm3', 1.12)
wall_shielding.add_element('H', 0.133, percent_type='wo')
wall_shielding.add_element('C', 0.817, percent_type='wo')
wall_shielding.add_element('B', 0.05,  percent_type='wo')

structural_steel = openmc.Material(name='Structural_Steel')
structural_steel.set_density('g/cm3', 7.85)
structural_steel.add_element('Fe', 0.97537, percent_type='wo')
structural_steel.add_element('C',  0.0025,  percent_type='wo')
structural_steel.add_element('Mn', 0.0120,  percent_type='wo')
structural_steel.add_element('Si', 0.0060,  percent_type='wo')
structural_steel.add_element('Ni', 0.0025,  percent_type='wo')
structural_steel.add_element('Cr', 0.0025,  percent_type='wo')
structural_steel.add_element('Mo', 0.0008,  percent_type='wo')
structural_steel.add_element('Cu', 0.0035,  percent_type='wo')
structural_steel.add_element('P',  0.00035, percent_type='wo')
structural_steel.add_element('S',  0.00035, percent_type='wo')
structural_steel.add_element('V',  0.0005,  percent_type='wo')
structural_steel.add_element('B',  0.00003, percent_type='wo')

#Concrete to be mixed with Barite
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


# Placeholder until GA provides roof spec
roof_shielding = openmc.Material(name='Roof_Shielding')
roof_shielding.set_density('g/cm3', 2.3)
roof_shielding.add_element('H',  0.010, percent_type='wo')
roof_shielding.add_element('O',  0.529, percent_type='wo')
roof_shielding.add_element('Si', 0.337, percent_type='wo')
roof_shielding.add_element('Ca', 0.044, percent_type='wo')
roof_shielding.add_element('Al', 0.034, percent_type='wo')
roof_shielding.add_element('Fe', 0.014, percent_type='wo')
roof_shielding.add_element('Na', 0.016, percent_type='wo')
roof_shielding.add_element('K',  0.013, percent_type='wo')
roof_shielding.add_element('Mg', 0.002, percent_type='wo')
roof_shielding.add_element('C',  0.001, percent_type='wo')

air = openmc.Material(name='Air')
air.set_density('g/cm3', 0.001205)
air.add_element('N',  0.784, percent_type='wo')
air.add_element('O',  0.210, percent_type='wo')
air.add_element('Ar', 0.006, percent_type='wo')

void = openmc.Material(name='Void')
void.set_density('g/cm3', 1e-10)
void.add_element('He', 1.0, percent_type='wo')

# CR39 compound for OSL chips (uses free atom cross sections)
CR39 = openmc.Material(name="CR39")
CR39.set_density("g/cm3", 1.31)
# CR39 element proportions based on single polymer structure (repeated to form the actual chain)
CR39.add_element("C", 0.324)
CR39.add_element("H", 0.486)
CR39.add_element("O", 0.190) # 12 + 18 + 7 = 37, 12/37 = 0.324, 18/37 = 0.486, remainder =  0,190

mats = openmc.Materials([Wall_Concrete,
                         structural_steel,
                         inconel,
                         roof_shielding,
                         wall_shielding,
                         air,
                         void,
                         CR39])
mats.export_to_xml()


# GEOMETRY — DagMC, graveyard in H5M handles outer boundary

Stage_0_Model = openmc.DAGMCUniverse(
    filename='/storage/work/ajg7072/Capstone/Ring_source/cad_modeling/431_GA_Stage_0.9.h5m',
    auto_geom_ids=True,
)
geometry = openmc.Geometry(Stage_0_Model)
geometry.export_to_xml()


# SOURCE — D-D ring, 2.45 MeV
# fusion_ring_source uses METRES — always divide cm values by 100

my_source = fusion_ring_source(
    radius=R_major / 100,       # cm -> m  (187.5 cm = 1.875 m)
    angles=(0.0, 2 * math.pi),
    z_placement=z_plasma / 100, # cm -> m  (119.4 cm = 1.194 m)
    temperature=20000.0,
    fuel={"D": 1.0},            # pure D-D -> 2.45 MeV
)

 
# SETTINGS
settings = openmc.Settings()
settings.batches   = 100
settings.particles = 500000
settings.run_mode  = 'fixed source'
settings.sources   = [my_source]
settings.export_to_xml()

# MESH AND TALLIES
mesh = openmc.RegularMesh(name='regmesh1')
mesh.dimension   = (40, 40, 20)
mesh.lower_left  = (-1575, -1575, -1075)
mesh.upper_right = ( 1575,  1575,  1075)
mesh_filter = openmc.MeshFilter(mesh)

# 3D flux tally over full mesh
flux_tally = openmc.Tally(name="flux_3d")
flux_tally.filters = [mesh_filter]
flux_tally.scores  = ["flux"]

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
energies_icrp, dose_cf = openmc.data.dose_coefficients(
    particle='neutron',
    geometry='ISO',
    data_source='icrp116',
)
dose_energy_filter = openmc.EnergyFunctionFilter(energies_icrp, dose_cf)
dose_tally = openmc.Tally(name="dose_3d")
dose_tally.filters = [mesh_filter, dose_energy_filter]
dose_tally.scores  = ["flux"]   # flux x energy function filter = dose

tallies = openmc.Tallies([flux_tally, dose_tally, spectrum_tally])
tallies.export_to_xml()

# RUN
openmc.run()

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

pSv_to_mrem = 1e-7   # 1 pSv = 1e-12 Sv = 1e-7 mrem


# POST-PROCESSING
sp = openmc.StatePoint("statepoint.100.h5")

flux_data = sp.get_tally(name="flux_3d")
nx, ny, nz = 40, 40, 20   # must match mesh.dimension
flux_3d = flux_data.mean.reshape((nx, ny, nz))

dose_data    = sp.get_tally(name="dose_3d")
dose_3d_pSv  = dose_data.mean.reshape((nx, ny, nz))  # pSv/source-n
dose_mrem_3d = dose_3d_pSv * S_n * pSv_to_mrem       # mrem integrated

def get_mesh_dose_mrem(xyz_cm,
                       ll=(-1575, -1575, -1075),
                       ur=( 1575,  1575,  1075),
                       dims=(40, 40, 20)):
    """Return ICRP-116 H*(10) dose (mrem) at a point in CAD coords (cm)."""
    vsize = (np.array(ur) - np.array(ll)) / np.array(dims)
    idx   = ((np.array(xyz_cm) - np.array(ll)) / vsize).astype(int)
    idx   = np.clip(idx, 0, np.array(dims) - 1)
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
detectors = {
    # Name                   xyz (cm) (tbd)            FN_mrem(from report)   Notes
    'EXP_6_Transrex':          ((-900, -900,    0),    30,      'Transrex area, far wall'),
    'EXP_7_Transrex':          ((-950, -950,    0),    None,    'No neutron signal — below MDL'),
    'EXP_8_MainEntrance':      ((-800,  200,    0),    30,      'Main entrance, well shielded'),
    'EXP_9_NorthDoor_w_view':  ((-500,  700,    0),    860,     'Lower N door WITH line-of-sight to vessel'),
    'EXP_10_NorthDoor_no_view':((-600,  750,    0),    8910,    'Lower N door WITHOUT vessel view — verify vs EXP_9'),
    'EXP_11_DeyongOut':        ((-300, -200,    0),    5060,    'Outside Deyong shield box'),
    'EXP_12_DeyongIn':         ((-300, -200,    0),    None,    'Inside Deyong box — FN exceeded detector cap'),
    'EXP_13_150V+1':           (( 100,  100,    0),    6620,    '150V+1 near vessel'),
    'EXP_14_150V-1':           (( 100,   50,    0),    19540,   '150V-1 near vessel — highest FN'),
    'EXP_15':                  (( 200, -100,    0),    18210,   'Second highest FN location'),
}

print("\n" + "=" * 95)
print(f"C/E BENCHMARKING — Aug 14 2025 | S_n={S_n:.2e} n | ICRP-116 ISO Fast Neutron Dose (mrem)")
print("=" * 95)
print(f"  {'Detector':<28} {'xyz (cm)':<26} {'C (mrem)':<12} {'E (mrem)':<12} {'C/E':<7} Flag")
print("-" * 95)

ce_vals = {}
for name, (xyz, E_mrem, note) in detectors.items():
    if E_mrem is None:
        print(f"  {name:<28} {str(xyz):<26} {'---':<12} {'MDL/CAP':<12} {'---':<7} SKIP -- {note}")
        continue

    C_mrem = get_mesh_dose_mrem(xyz)

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
        print(f"  {name:<28} {str(xyz):<26} {C_mrem:<12.1f} {E_mrem:<12} {ratio:<7.2f} {flag}")

print("-" * 95)
if ce_vals:
    r = list(ce_vals.values())
    in_band = sum(0.7 <= v <= 1.3 for v in r)
    print(f"\n  C/E summary:  mean={np.mean(r):.2f}  min={np.min(r):.2f}  "
          f"max={np.max(r):.2f}  |  Within +/-30%: {in_band}/{len(r)}")

print("""
  Notes:
    - Detector xyz are FLOOR PLAN ESTIMATES -- replace with surveyed CAD coords
    - EXP 7: below minimum detectable level -- excluded from C/E
    - EXP 12: fast neutron component exceeded Landauer cap -- excluded (lower bound only)
    - EXP 10 has HIGHER dose than EXP 9 despite 'no vessel view' -- verify badge
      placement in field records (10 and 9 may be swapped)
    - Photon dose NOT compared -- activation gammas not in OpenMC model
    - C/E targets: +/-30% near-field, factor of 2 acceptable in deep shielded regions
""")
print("=" * 95)

print("=" * 50)
print("ALL DONE -- job completed successfully")
print("=" * 50)
