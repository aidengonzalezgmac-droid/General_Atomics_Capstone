"""
This is the ringsource "downgrade". We dont need all the plasma parameters listed in the tokamak.py, GA send a ring source is adequate
Since no cad models have been sent, this is just for proof of concept and understanding
"""

import os
os.environ['OPENMC_CROSS_SECTIONS'] = '/storage/work/ajg7072/Capstone/endf/cross_sections.xml'

import openmc
import numpy as np
import matplotlib.pyplot as plt
import math
from openmc_plasma_source import fusion_ring_source

# ============================================================
# DIII-D GEOMETRY PARAMETERS (cm)
# ============================================================
R_major        = 170.0
a_minor        = 60.0
t_wall         = 3.0
a_outer        = a_minor + t_wall  # 63 cm
ring_radius_cm = R_major

# ============================================================
# MATERIALS
# ============================================================
inconel = openmc.Material(name='inconel625')
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

air = openmc.Material(name='air')
air.set_density('g/cm3', 0.001225)
air.add_element('N',  0.78)
air.add_element('O',  0.21)
air.add_element('Ar', 0.01)

pbwall = openmc.Material(name='pbwall')
pbwall.set_density('g/cm3', 1.12)
pbwall.add_element('H', 0.133, percent_type='wo')
pbwall.add_element('C', 0.817, percent_type='wo')
pbwall.add_element('B', 0.05,  percent_type='wo')

steel = openmc.Material(name='steel')
steel.set_density('g/cm3', 7.85)
steel.add_element('Fe', 0.97537, percent_type='wo')
steel.add_element('C',  0.0025,  percent_type='wo')
steel.add_element('Mn', 0.0120,  percent_type='wo')
steel.add_element('Si', 0.0060,  percent_type='wo')
steel.add_element('Ni', 0.0025,  percent_type='wo')
steel.add_element('Cr', 0.0025,  percent_type='wo')
steel.add_element('Mo', 0.0008,  percent_type='wo')
steel.add_element('Cu', 0.0035,  percent_type='wo')
steel.add_element('P',  0.00035, percent_type='wo')
steel.add_element('S',  0.00035, percent_type='wo')
steel.add_element('V',  0.0005,  percent_type='wo')
steel.add_element('B',  0.00003, percent_type='wo')

concrete = openmc.Material(name='concrete')
concrete.set_density('g/cm3', 2.3)
concrete.add_element('H',  0.010, percent_type='wo')
concrete.add_element('O',  0.529, percent_type='wo')
concrete.add_element('Si', 0.337, percent_type='wo')
concrete.add_element('Ca', 0.044, percent_type='wo')
concrete.add_element('Al', 0.034, percent_type='wo')
concrete.add_element('Fe', 0.014, percent_type='wo')
concrete.add_element('Na', 0.016, percent_type='wo')
concrete.add_element('K',  0.013, percent_type='wo')
concrete.add_element('Mg', 0.002, percent_type='wo')
concrete.add_element('C',  0.001, percent_type='wo')

barite = openmc.Material(name='barite')
barite.set_density('g/cm3', 4.5)
barite.add_element('Ba', 0.589, percent_type='wo')
barite.add_element('S',  0.137, percent_type='wo')
barite.add_element('O',  0.274, percent_type='wo')

floor_mat = openmc.Material.mix_materials([concrete, barite], [0.85, 0.15], 'wo')

mats = openmc.Materials([inconel, air, pbwall, steel, floor_mat])
mats.export_to_xml()

# ============================================================
# GEOMETRY — surfaces first, then regions, then cells
# ============================================================
room_l  = 1830
room_h  = 856
room_w  = 1830
wall_t  = 45
roof_t  = 36
floor_t = 30

# --- Torus surfaces ---
inner_torus = openmc.ZTorus(x0=0, y0=0, z0=60, a=R_major, b=a_minor, c=a_minor)
outer_torus = openmc.ZTorus(x0=0, y0=0, z0=60, a=R_major, b=a_outer, c=a_outer)

# --- Room bounding planes ---
x_min_outer = openmc.XPlane(-room_l/2 - wall_t, boundary_type='vacuum')
x_max_outer = openmc.XPlane( room_l/2 + wall_t, boundary_type='vacuum')
y_min_outer = openmc.YPlane(-room_w/2 - wall_t, boundary_type='vacuum')
y_max_outer = openmc.YPlane( room_w/2 + wall_t, boundary_type='vacuum')
z_min_outer = openmc.ZPlane(-floor_t, boundary_type='vacuum')
z_max_outer = openmc.ZPlane( room_h + roof_t, boundary_type='vacuum')

x_min_inner = openmc.XPlane(-room_l/2)
x_max_inner = openmc.XPlane( room_l/2)
y_min_inner = openmc.YPlane(-room_w/2)
y_max_inner = openmc.YPlane( room_w/2)
z_min_inner = openmc.ZPlane(0)
z_max_inner = openmc.ZPlane(room_h)

# --- Room region (now defined after its surfaces) ---
room_region = (
    +x_min_inner & -x_max_inner &
    +y_min_inner & -y_max_inner &
    +z_min_inner & -z_max_inner
)

# 
plasma_cell = openmc.Cell(name='plasma')
plasma_cell.region = -inner_torus & room_region


vessel_cell = openmc.Cell(name='inconel_vessel')
vessel_cell.region = +inner_torus & -outer_torus & room_region
vessel_cell.fill = inconel

# --- Interior air: room minus the vessel ---
interior = openmc.Cell(name='interior')
interior.region = room_region & +outer_torus
interior.fill = air

# --- Walls ---
west_wall = openmc.Cell(name='west_wall')
west_wall.region = +x_min_outer & -x_min_inner & +y_min_outer & -y_max_outer & +z_min_inner & -z_max_inner
west_wall.fill = pbwall

east_wall = openmc.Cell(name='east_wall')
east_wall.region = +x_max_inner & -x_max_outer & +y_min_outer & -y_max_outer & +z_min_inner & -z_max_inner
east_wall.fill = pbwall

north_wall = openmc.Cell(name='north_wall')
north_wall.region = +x_min_inner & -x_max_inner & +y_max_inner & -y_max_outer & +z_min_inner & -z_max_inner
north_wall.fill = pbwall

south_wall = openmc.Cell(name='south_wall')
south_wall.region = +x_min_inner & -x_max_inner & +y_min_outer & -y_min_inner & +z_min_inner & -z_max_inner
south_wall.fill = pbwall

floor_cell = openmc.Cell(name='floor')
floor_cell.region = +x_min_outer & -x_max_outer & +y_min_outer & -y_max_outer & +z_min_outer & -z_min_inner
floor_cell.fill = floor_mat

roof_cell = openmc.Cell(name='roof')
roof_cell.region = +x_min_outer & -x_max_outer & +y_min_outer & -y_max_outer & +z_max_inner & -z_max_outer
roof_cell.fill = steel

universe = openmc.Universe(cells=[
    plasma_cell, vessel_cell,
    interior,
    west_wall, east_wall, north_wall, south_wall,
    floor_cell, roof_cell
])

geometry = openmc.Geometry(universe)
geometry.export_to_xml()

# SOURCE
my_source = fusion_ring_source(
    radius=ring_radius_cm,
    angles=(0.0, 2 * math.pi),
    temperature=20000.0,
    fuel={"D": 1.0},
)

# SETTINGS
settings = openmc.Settings()
settings.batches = 100
settings.particles = 5000
settings.run_mode = 'fixed source'
settings.sources = [my_source]
settings.export_to_xml()

# GEOMETRY PLOTS

plot_xy = openmc.Plot(name='xy_midplane')
plot_xy.basis = 'xy'
plot_xy.origin = (0, 0, 0)
plot_xy.width = (2 * (room_l/2 + wall_t), 2 * (room_w/2 + wall_t))
plot_xy.pixels = (1920, 1920)
plot_xy.color_by = 'material'
plot_xy.filename = 'DIIID_xy_midplane'

plot_xz = openmc.Plot(name='xz_midplane')
plot_xz.basis = 'xz'
plot_xz.origin = (0, 0, room_h/2)
plot_xz.width = (2 * (room_l/2 + wall_t), room_h + roof_t + floor_t)
plot_xz.pixels = (1920, 960)
plot_xz.color_by = 'material'
plot_xz.filename = 'DIIID_xz_midplane'

plots = openmc.Plots([plot_xy, plot_xz])
plots.export_to_xml()
openmc.plot_geometry()

# ============================================================
# MESH AND TALLIES
# ============================================================
mesh = openmc.RegularMesh()
mesh.dimension = (150, 150, 75)
mesh.lower_left  = (-room_l/2, -room_w/2, 0)
mesh.upper_right = ( room_l/2,  room_w/2, room_h)

mesh_filter = openmc.MeshFilter(mesh)
flux_tally  = openmc.Tally(name="flux_3d")
flux_tally.filters = [mesh_filter]
flux_tally.scores  = ["flux"]

energy_bins    = np.logspace(3, 7, 200)
energy_filter  = openmc.EnergyFilter(energy_bins)
spectrum_tally = openmc.Tally(name="energy_spectrum")
spectrum_tally.filters = [energy_filter]
spectrum_tally.scores  = ["flux"]

tallies = openmc.Tallies([flux_tally, spectrum_tally])
tallies.export_to_xml()

# ============================================================
# RUN
# ============================================================
openmc.run()

# ============================================================
# POST-PROCESSING: Energy Spectrum
# ============================================================
sp = openmc.StatePoint("statepoint.100.h5")

spec_tally      = sp.get_tally(name="energy_spectrum")
energy_bins_eV  = energy_filter.values
energy_bins_MeV = energy_bins_eV / 1e6
flux_1d  = spec_tally.mean[:, 0, 0]
dE       = np.diff(energy_bins_eV)
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

# ============================================================
# POST-PROCESSING: 2D Flux Slices
# ============================================================
flux_data = sp.get_tally(name="flux_3d")
nx, ny, nz = 150, 150, 75
flux_3d = flux_data.mean.reshape((nx, ny, nz))

x_arr = np.linspace(-room_l/2, room_l/2, nx) / 100
y_arr = np.linspace(-room_w/2, room_w/2, ny) / 100
z_arr = np.linspace(0, room_h, nz) / 100

# XY midplane
mid_z   = nz // 2
flux_xy = flux_3d[:, :, mid_z]
fig2d, ax2d = plt.subplots(figsize=(10, 10))
im = ax2d.pcolormesh(x_arr, y_arr, flux_xy.T,
                     norm=plt.matplotlib.colors.LogNorm(
                         vmin=np.percentile(flux_xy[flux_xy > 0], 5),
                         vmax=flux_xy.max()),
                     cmap='inferno', shading='auto')
fig2d.colorbar(im, ax=ax2d, label='Neutron Flux [n/cm²/source particle]')
ax2d.set_xlabel('X (m)', fontsize=13)
ax2d.set_ylabel('Y (m)', fontsize=13)
ax2d.set_title('Neutron Flux — XY Midplane\nDIII-D CSG Model, D-D Ring Source, Inconel 625 Vessel',
               fontsize=13)
ax2d.set_aspect('equal')
fig2d.tight_layout()
fig2d.savefig("DIIID_flux_xy_slice.png", dpi=300)
plt.close()
print("XY flux slice saved.")

# XZ slice
mid_y   = ny // 2
flux_xz = flux_3d[:, mid_y, :]
fig_xz, ax_xz = plt.subplots(figsize=(12, 6))
im_xz = ax_xz.pcolormesh(x_arr, z_arr, flux_xz.T,
                          norm=plt.matplotlib.colors.LogNorm(
                              vmin=np.percentile(flux_xz[flux_xz > 0], 5),
                              vmax=flux_xz.max()),
                          cmap='inferno', shading='auto')
fig_xz.colorbar(im_xz, ax=ax_xz, label='Neutron Flux [n/cm²/source particle]')
ax_xz.set_xlabel('X (m)', fontsize=13)
ax_xz.set_ylabel('Z (m)', fontsize=13)
ax_xz.set_title('Neutron Flux — XZ Slice (y = 0)\nDIII-D CSG Model, D-D Ring Source, Inconel 625 Vessel',
                fontsize=13)
fig_xz.tight_layout()
fig_xz.savefig("DIIID_flux_xz_slice.png", dpi=300)
plt.close()
print("XZ flux slice saved.")

# ============================================================
# VTK EXPORT FOR PARAVIEW
# ============================================================
try:
    from pyevtk.hl import gridToVTK

    x_vtk = np.linspace(-room_l/2, room_l/2, nx + 1, dtype=np.float64) / 100
    y_vtk = np.linspace(-room_w/2, room_w/2, ny + 1, dtype=np.float64) / 100
    z_vtk = np.linspace(0, room_h,           nz + 1, dtype=np.float64) / 100

    gridToVTK(
        "./DIIID_flux",
        x_vtk, y_vtk, z_vtk,
        cellData={"neutron_flux": np.ascontiguousarray(flux_3d, dtype=np.float64)}
    )
    print("VTK file saved as DIIID_flux.vts — open in ParaView")

except ImportError:
    print("pyevtk not found. Install with: conda install -c conda-forge pyevtk -y")
    np.save("flux_3d_data.npy", flux_3d)
    np.save("flux_coords_x.npy", x_arr)
    np.save("flux_coords_y.npy", y_arr)
    np.save("flux_coords_z.npy", z_arr)
    print("Flux arrays saved as .npy files instead.")

print('\nAll outputs complete.')