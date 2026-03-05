"""
This is the ringsource "downgrade". We dont need all the plasma parameters listed in the tokamak.py, GA said a ring source is adequate
This is for cad modeling
"""

import os
os.environ['OPENMC_CROSS_SECTIONS'] = '/storage/work/ajg7072/Capstone/endf/cross_sections.xml'

import openmc
import numpy as np
import matplotlib.pyplot as plt
import math
from openmc_plasma_source import fusion_ring_source

#
# DIII-D CAD coordinate system (from H5M centroid check)
R_major        = 170.0   # cm — plasma major radius
ring_radius_cm = R_major
z_plasma       = 119.4   # cm — tokamak Z offset in CAD coordinates

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

# Intermediate ingredients for the mix — not added to mats directly
_concrete = openmc.Material()
_concrete.set_density('g/cm3', 2.3)
_concrete.add_element('H',  0.010, percent_type='wo')
_concrete.add_element('O',  0.529, percent_type='wo')
_concrete.add_element('Si', 0.337, percent_type='wo')
_concrete.add_element('Ca', 0.044, percent_type='wo')
_concrete.add_element('Al', 0.034, percent_type='wo')
_concrete.add_element('Fe', 0.014, percent_type='wo')
_concrete.add_element('Na', 0.016, percent_type='wo')
_concrete.add_element('K',  0.013, percent_type='wo')
_concrete.add_element('Mg', 0.002, percent_type='wo')
_concrete.add_element('C',  0.001, percent_type='wo')

_barite = openmc.Material()
_barite.set_density('g/cm3', 4.5)
_barite.add_element('Ba', 0.589, percent_type='wo')
_barite.add_element('S',  0.137, percent_type='wo')
_barite.add_element('O',  0.274, percent_type='wo')

wall_concrete = openmc.Material.mix_materials([_concrete, _barite], [0.85, 0.15], 'wo', name='Wall_Concrete')

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

mats = openmc.Materials([inconel, wall_shielding, structural_steel, wall_concrete, roof_shielding])
mats.export_to_xml()

# GEOMETRY — surfaces first, then regions, then cells
Stage_0_Model = openmc.DAGMCUniverse(filename = '/storage/work/ajg7072/Capstone/Ring_source/cad_modeling/431_GA_Stage_0.6.h5m', auto_geom_ids = True)
geometry = openmc.Geometry(Stage_0_Model)
geometry.export_to_xml()

# SOURCE 
my_source = fusion_ring_source(
    radius=R_major / 100,        # cm → m
    angles=(0.0, 2 * math.pi),
    z_placement=119.4 / 100,     # cm → m
    temperature=20000.0,
    fuel={"D": 1.0},
)
# SETTINGS
settings = openmc.Settings()
settings.batches = 10
settings.particles = 50000
settings.run_mode = 'fixed source'
settings.sources = [my_source]
settings.export_to_xml()

# MESH AND TALLIES
mesh = openmc.RegularMesh(name = 'regmesh1',)
mesh.dimension = (50,50,25)
mesh.lower_left = (-1575,-1575,-1075)
mesh.upper_right = (1575,1575,1075) #square mesh with rel. arbitrary dimensions
mesh_filter = openmc.MeshFilter(mesh)

mesh_filter = openmc.MeshFilter(mesh)
flux_tally  = openmc.Tally(name="flux_3d")
flux_tally.filters = [mesh_filter]
flux_tally.scores  = ["flux"]

energy_bins  = np.logspace(3, 7, 200)
energy_filter  = openmc.EnergyFilter(energy_bins)
spectrum_tally = openmc.Tally(name="energy_spectrum")
spectrum_tally.filters = [energy_filter]
spectrum_tally.scores  = ["flux"]

tallies = openmc.Tallies([flux_tally, spectrum_tally])
tallies.export_to_xml()


# RUN
openmc.run()

#Put normalzing factors here to get accurate fluxes
##
###
####
#####
######
#####
####
##
# POST-PROCESSING: Energy Spectrum
sp = openmc.StatePoint("statepoint.10.h5")

spec_tally = sp.get_tally(name="energy_spectrum")
energy_bins_eV = energy_filter.values
energy_bins_MeV = energy_bins_eV / 1e6
flux_1d = spec_tally.mean[:, 0, 0]
dE = np.diff(energy_bins_eV)
spectrum = flux_1d / dE

fig_spec, ax = plt.subplots(figsize=(10, 6))
ax.step(energy_bins_MeV[:-1], spectrum * 1e6, where="post",
        color='#E8432A', linewidth=1.8, label='D-D Ring Source')
ax.axvline(x=2.45, color='gold', linestyle='--', linewidth=1.5,
           label='D-D peak (2.45 MeV)')
ax.set_xlabel("Energy (MeV)", fontsize=13)
ax.set_ylabel("Flux per unit energy [n/cm²/MeV/source]", fontsize=13)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_title("Neutron Energy Spectrum — DIII-D D-D Ring Source\n"
             "R = 1.7 m, a = 0.6 m, Inconel 625 Vacuum Vessel", fontsize=13)
ax.legend(fontsize=11)
ax.grid(True, which="both", alpha=0.3)
fig_spec.tight_layout()
fig_spec.savefig("DIIID_Energy_Spectrum.png", dpi=300)
plt.close()
print("Energy spectrum saved.")
###This needs to be revisted, the spectrum doesnt look accurate at all and its probably because of how to flux is measured, or rather where 


# POST-TRANSPORT: Export flux mesh to VTK for ParaView
from pyevtk.hl import gridToVTK

flux_data = sp.get_tally(name="flux_3d")
nx, ny, nz = 50, 50, 25
flux_3d = flux_data.mean.reshape((nx, ny, nz))
flux_err = flux_data.get_values(scores=['flux'], value='rel_err').reshape((nx, ny, nz))

x_vtk = np.linspace(-1575, 1575, nx+1, dtype=np.float64) / 100  # cm -> m
y_vtk = np.linspace(-1575, 1575, ny+1, dtype=np.float64) / 100
z_vtk = np.linspace(-1075, 1075, nz+1, dtype=np.float64) / 100

gridToVTK(
    './DIIID_flux',
    x_vtk, y_vtk, z_vtk,
    cellData={
        'neutron_flux': np.ascontiguousarray(flux_3d, dtype=np.float64),
        'rel_error':    np.ascontiguousarray(flux_err, dtype=np.float64),
    }
)
print("VTK saved as DIIID_flux.vts — open in ParaView")


print("="*50)
print("ALL DONE — job completed successfully")
print("="*50)