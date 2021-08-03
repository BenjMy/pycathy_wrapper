#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 16:36:58 2021

@author: ben
"""

# cathy project rhizo lab
from cathy_tools import CATHY
# from rhizo_tools import create_infitration # create_DA

import plottools as cplt
#import cathy_plot as ctplot

import os
import numpy as np


# %% Run cathy processor
# this function is a wrapper. It builds a makefile and the executable.
# Then the exe is run and use core CATHY file to solve the problem

path2prj = os.getcwd()

simu = CATHY(dirName=path2prj,prjName='rhizo_prj',clear_outputs=True)


xyzb=np.loadtxt('/home/ben/Documents/CATHY/rhizoLab/CATHY_additional_files_rhizo/elecsXYZ_Rhizo_LAB.csv',skiprows=1,delimiter=',')
xyzb

xll=0
yll=0
zll=0



delta_x = 0.025
delta_y = 0.01
# # delta_z = 0.025
xb = np.arange(xll,0.45+delta_x,delta_x)
yb = np.arange(xll,0.03+delta_y,delta_y)
nstr=20
# zb = np.linspace(0,0.5,nstr+1)
zb = np.arange(zll,0.5+delta_x,delta_x)
zr=list(np.ones(len(zb)-1)/(len(zb)-1))
len(zr)
sum(zr)

N = len(xb)
M = len(yb)
dem_rhizo =  np.ones([M,N])
np.shape(dem_rhizo)

N_celle=N*M
simu.update_prepo_inputs(DEM=dem_rhizo,
                          xllcornery=xll,yllcorner=yll,
                          delta_x=delta_x,delta_y=delta_y,
                          N=N,M=M,N_celle=N_celle,
                          nstr=nstr,base=0.5,n1=20,
                          zratio=zr)

with open(os.path.join('/home/ben/Documents/GitHub/pycathy_wrapper/pyCATHY/rhizo_prj/prepro/src/dtm_13.val'), 'w+') as f:
    np.savetxt(f, dem_rhizo, fmt='%1.4e')   # use exponential notation

# dem_test =  np.ones([20,5])
# simu.update_prepo_inputs(DEM=dem_test,
#                       N=np.shape(dem_test)[0],
#                       M=np.shape(dem_test)[1],
#                       delta_x=10.,delta_y=10.,
#                       nstr=20,zratio=zr,base=5,m1=20,
#                       N_celle=np.shape(dem_test)[0]*np.shape(dem_test)[1],
#                       nzone=1)


# simu.update_prepo_inputs(DEM=dem_test,
#                       N=np.shape(dem_test)[0],
#                       M=np.shape(dem_test)[1],
#                       N_celle=np.shape(dem_test)[0]*np.shape(dem_test)[1])

# simu.hapin
# simu.dem_parameters

# # remove outlet as it is a lab experiment wirth no flux outside the rhizotron
# simu.run_preprocessor(verbose=True,KeepOutlet=False) 


simu.run_preprocessor(verbose=True,KeepOutlet=False) # remove outlet as it is a lab experiment with no flux outside the rhizotron

# simu.create_parm(NPRT=4,TIMPRTi=[1800.,3600.,7200.,80000.])
# simu.parm




# simu.set_elecs(xyz=xyzb[:,0:2])

# # xyz_drip=np.loadtxt('/home/ben/Documents/CATHY/rhizoLab/CATHY_additional_files_rhizo/drippers.txt', delimiter=',')
# # simu.set_drippers(xyz=xyz_drip)

# simu.create_infitration('inputs_rhizo')
# simu.drippers
# # create_infitration(drippersPos=[],RWU=False)
# # # rhizo.create_DA(drippersPos=[],RWU=False)


# generate xyz , grid2d.exp grid3d
# simu.run_processor(verbose=True,IPRT1=3)

# generate vtk files
simu.run_processor(verbose=True,IPRT1=2,TRAFLAG=0,TMAX=7200)



# %% Explore outputs

cplt.showvtk('./rhizo_prj/vtk/100.vtk')

# os.getcwd()
# import time
# import glob
# import pyvista as pv


# cplt.showvtkTL('./rhizo_prj/vtk/100.vtk')

# cplt.showvtk('./rhizo_prj/vtk/103.vtk')

# cplt.showvtk('./rhizo_prj/vtk/cele200.vtk')
# cplt.showvtk('./rhizo_prj/vtk/cnod100.vtk')
# import pyvista as pv
# filename='./rhizo_prj/vtk/100.vtk'
# mesh = pv.read(filename)
# plotter = pv.Plotter()
# _ = plotter.add_mesh(mesh,show_edges=True)

# legend_entries = []
# legend_entries.append(['Time='+ str(mesh['TIME']), 'w'])
# _ = plotter.add_legend(legend_entries)
# plotter.show_grid()
# #plotter.add_mesh_clip_box(mesh, color='white')
# cpos = plotter.show()

# mesh.array_names


# %% test
