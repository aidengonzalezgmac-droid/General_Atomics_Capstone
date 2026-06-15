# 431_Gif Generation code 

# NOTE: Code requires a valid statepoint file generated to perform this work automatically.

# import modules: 
import openmc 
import datetime
import numpy as np
import matplotlib.pyplot as py # in case it is needed. 
import math

# print beginning of file: 
# space terminal by printing current system time to verify operability.
now = datetime.datetime.now()
print('\n') # space the terminal
print('image generation begins at:', now)
print('\n') # space the terminal


# set conditions/data used in plot generation files on 3.9 file:
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



VisMesh = 100
squarecoord = VisMesh
squarecoord2 = 1500
VisZcoord2 = 592
mesh2xy = VisMesh
mesh2z = VisMesh
esbin = 10
pSv_to_mrem = 1*(10**(-7))

sp = openmc.StatePoint("statepoint.20.h5")

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

meshvol = mesh2.volumes[0][0][0]

dose_data    = sp.get_tally(name="dose_3d")
#dose_3d_pSv  = dose_data.mean.reshape(nx, ny, nz)
dose_3d_pSv_slice = dose_data.get_slice(scores=['flux'])
dose_3d_pSv = dose_3d_pSv_slice.get_reshaped_data(expand_dims=True, value = 'mean')  # pSv/source-n
flub_factor = 1 # to be used for low particle simulation ONLY, set to 1 for lower results.
dose_mrem_3d = dose_3d_pSv*S_n*(pSv_to_mrem)/(meshvol)

fsd = sp.get_tally(name="energy_spectrum") # fsd = fluency spectrum data
fsd_slice = fsd.get_slice(scores=['flux'])
fsd_r = fsd_slice.get_reshaped_data(expand_dims=True, value ='mean') # returns slice plot flux values at each 
h = range(esbin)



def get_energy_data(xyz_cm,
                       ll=(-squarecoord2, -squarecoord2, -VisZcoord2),
                       ur=( squarecoord2,  squarecoord2,  VisZcoord2), 
                       dims=(mesh2xy, mesh2xy, mesh2z)): # consider changing these to sq3 and VisCoord3
    # ur, ll squarecoord2 and VisZcoord 2 used on vsize calculation.
    vsize = (np.array(ur) - np.array(ll)) / np.array(dims) 
    idx   = ((np.array(xyz_cm) - np.array(ll)) / vsize).astype(int) 
    idx   = np.clip(idx, 0, np.array(dims)-1) # clips dataset to get to closest mesh 
    testidx = fsd_r[0] 
    testidx2 = testidx[idx[0]-1,idx[1]-1,idx[2]-1] # diagnostic variable
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

# import statepoint file: 

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
logset = False


# setup loops: 

# I want a single image coordinate that spans each vismesh which is 100 total images,

bellerophon = np.linspace(0,100,num = 100).astype(int)
# bellerophon = [0,1,2,3]

for number in range(len(bellerophon)):
    # NEED TO REQEUST SPECIFIC OSL POSITONS BUT HERE WE GO: AT CENTER OF TOKAMAK 
    Diagnostic_position = [0,  0,     0] # CONTROLS VIEWING COORDINATE 
    Diagnostic_position1 = [bellerophon[number],  bellerophon[number],     bellerophon[number]]
    esd = (get_energy_data(Diagnostic_position1)) # the target cell is accessed by the input to this function.
    # This also uses the same behavior of mesh access as the regular mesh. 
    # the output of g_e_d function is all of the data, with the end being the array access.
    cell_dim_array = esd[-1] # accesses the last element of the esd array (numpy array which has raw clipped values from g_e_d.)

    slice_valuez = bellerophon[number] # cell_dim_array[2]
    slice_valuey = bellerophon[number] # cell_dim_array[1]
    slice_valuex = bellerophon[number] # cell_dim_array[0]

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
    plt.savefig(fname = f'/bin/431_combined_images/431_GA_Stage_3.0_{bellerophon[number]}.png')
    now2 = datetime.datetime.now()
    print(f'file 431_GA_stage_3.0_{bellerophon[number]} saved to folder at',now2)
    fig.clf()

# print end of file use.
now = datetime.datetime.now()
print('\n') # space the terminal
print('image generation ended at:', now)
print('\n') # space the terminal