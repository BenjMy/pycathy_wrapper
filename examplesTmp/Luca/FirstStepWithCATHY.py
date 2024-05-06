#!/usr/bin/env python
# coding: utf-8

import sys
from IPython import embed

import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from pyCATHY import cathy_tools
from pyCATHY.plotters import cathy_plots as cplt

# define project folder
path2prj = "."  # add your local path here
simu = cathy_tools.CATHY(
    dirName=path2prj,
    prj_name="my_first_project4"
)

simu.run_preprocessor(verbose=True)

# Change the preprocessing inputs (mesh size, ... vegetation map)

# simu.DEM
# simu.show_input('dem')
# plt.show()
# plt.close()

# dem number of elements
my_new_DEM = np.ones([40, 20])
# set the needed outlet point
my_new_DEM[-1, -1] = 1 - 1e-3




maxdepth = 4
zb = np.linspace(0, maxdepth, 10)   

nstr = len(zb)
zr = list((np.ones(len(zb))) / (nstr))
len(zr)


# zb = np.geomspace(1e-1, maxdepth, num=20)
# nstr = len(zb)
# # zr = [abs(zb[0] / maxdepth)]
# # zr.extend(list(abs(np.diff(zb) / maxdepth)))
# zr.extend(list(abs(np.diff(zb) / maxdepth)))

zr = np.geomspace(0.1, 1, 15)
print(zr)
zr /= np.sum(zr)
nstr = len(zr)
np.sum(zr)


simu.update_prepo_inputs(
    DEM=my_new_DEM,
    delta_x=0.5,
    delta_y=0.5,
    # xllcorner=0,
    # yllcorner=4e3,
    # nstr=nstr,
    zratio=list(zr),
    base=max(zb),
)



simu.create_mesh_vtk()

sys.exit()

DEM, DEM_header = simu.read_inputs('dem')

root_map, root_map_header = simu.read_inputs('root_map')
root_map_new = root_map
root_map_new[0: int(len(root_map_new) / 2), :] = 2
# simu.update_veg_map?
simu.update_veg_map(root_map_new)

simu.show_input('root_map')
plt.show()

# update soil with two vegetations Feddes
root_map_dict = {
    'PCANA': [0.0, 0.0],
    'PCREF': [-4.0, -4.0],
    'PCWLT': [-150.0, -150.0],
    'ZROOT': [-1.5, -1.5],
    'PZ': [1.0, 1.0],
    'OMGC': [1.0, 1.0],
}
simu.update_soil(FP_map=root_map_dict)

print(simu)

# Within the soil file we describe soil physical properties and root properties

# From github files:
# VGN,   - parameters for van Genuchten and extended van Genuchten 
# VGM,     moisture curves (other 'VG' parameters - specific storage,
# VGRMC,   porosity, and VGPNOT - are assigned nodally). VGM is 
# VGPSAT   derived from VGN. VGRMC is residual moisture content.
#          For IVGHU=0, VGPNOT is (porosity - VGRMC)/porosity,
#          or (1 - residual water saturation).
#          For IVGHU=1, VGPNOT is a continuity parameter, derived by
#          imposing a continuity requirement on the derivative of 
#          moisture content with respect to pressure head.
##
# retention curves parameters VGN, VGRMC, and VGPSAT
# - 'VGNCELL' (NSTR, NZONE): van Genuchten curve exponent  = n
# - 'VGRMCCELL' (NSTR, NZONE): residual moisture content = \thetaR
# - 'VGPSATCELL' (NSTR, NZONE): van Genuchten curve exponent -->
#    VGPSAT == -1/alpha (with alpha expressed in [L-1]);

SPP, FP = simu.read_inputs('soil', NVEG=2)

FP['ZROOT'].iloc[1] = 0.5
FP['ZROOT'].iloc[0] = 1

FP_map = FP.to_dict(orient='list')


# this would update the number of vegetation so that the new parameters can be added
# simu.update_cathyH(MAXVEG=2)

simu.update_soil(
    SPP=SPP,
    FP_map=FP_map,
    MAXVEG=2,
)

simu.update_veg_map(show=True)
plt.show()


# simu.show_input('soil')
# plt.show()


simu.run_preprocessor(verbose=True)


simu.create_mesh_vtk()

# pl = pv.Plotter(notebook=True)
# mymesh = pv.read('./my_first_project/vtk/my_first_project.vtk')
# pl.add_mesh(mymesh)
# pl.show_bounds()
# pl.show()



# simu.update_ic(INDP=0, pressure_head_ini=-5)  # pressure head in meter
# simu.update_ic(INDP=3, WTHEIGHT=1)  # pressure head in meter
# simu.update_ic(INDP=2, pressure_head_ini=-0.1)  # pressure head in meter
# Flag for pressure head initial conditions (all nodes). The default is 2.
#  - =0 for input of uniform initial conditions (one value read in)
#  - =1 for input of non-uniform IC's (one value read in for each node)
#  - =2 for calculation of fully saturated vertical hydrostatic equilibrium IC's
#    (calculated in subroutine ICVHE). In the case of IPOND>0, the fully saturated
#    hydrostatic IC is calculated (in subroutine ICVHEPOND) starting from the ponding head
#    values at the surface nodes, rather than surface pressure heads of 0.
#  - =3 for calculation of partially saturated vertical hydrostatic equilibrium IC's
#    (calculated in subroutine ICVHWT) with the water table height (relative to the base
#    of the 3‐d grid) given by parameter WTHEIGHT
#  - =4 for calculation of partially saturated vertical hydrostatic equilibrium IC's
#    (calculated in subroutine ICVDWT) with the water table depth (relative to the surface
#    of the 3‐d grid) given by parameter WTHEIGHT

df_atmbc = simu.read_inputs('atmbc')
print(df_atmbc.time)

simu.update_atmbc(
    HSPATM=1,
    IETO=0,
    time=list(df_atmbc.time),
    # VALUE=[Precipitation,ET]
    netValue=[1e-7, 5e-7]  # m/s
)
print(simu.atmbc)



# 86400.0/(60*60*24)

simu.parm['DTMIN'] = 1e-2  # in seconds
simu.parm['DELTAT'] = 1
vtk_times = np.linspace(0, 24, 25, endpoint=True) * 60 * 60
simu.update_parm(TIMPRTi=vtk_times.tolist())

# number of outputs (pressure, water content, some extra) to be included in the output vtk
# TIMPRTi selected output times
# simu.update_parm(VTKF=4)
# simu.update_parm(
#     TIMPRTi=list_of_outputs_times,
#      VTKF=2 # both saturation and pressure head
# )

simu.run_processor(
    IPRT1=2,
    DTMIN=1e-2,
    DTMAX=1e2,
    DELTAT=5,
    TRAFLAG=0,
    verbose=False,
    VTKF=4
)

pl = pv.Plotter()
cplt.show_vtk(
    unit="pressure", 
    timeStep=1, 
    path=simu.workdir + "/my_first_project/vtk/",
    ax=pl,
)
pl.show()


# pl = pv.Plotter()
# plt.show_vtk_TL()


# NOTE
## how to use modularity
# simu.show_input('name_of the file')
# simu.show_name_file
# out = simu.read_name_file()
# plot(out)
