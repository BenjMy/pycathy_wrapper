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

# xll=0
# yll=0
# delta_x = 0.025
# delta_y = 0.0025
# delta_z = 0.025

# xb = np.arange(xll,0.5,delta_x)
# yb = np.arange(xll,0.03,delta_y)
# nstr=20
# zb = np.linspace(0,0.5,nstr+1)
# zr = list(zb[1:]/sum(zb))
# N = len(xb)
# M = len(yb)
# dem_rhizo =  np.ones([N,M])
# N_celle=N*M

# simu.update_prepo_inputs(DEM=dem_rhizo,
#                          xllcornery=xll,yllcorner=yll,
#                          delta_x=delta_x,delta_y=delta_y,
#                          N=N,M=M,N_celle=N_celle,
#                          nstr=nstr,
#                          zratio=zr)


simu.update_prepo_inputs(DEM=dem_rhizo,
                         xllcornery=xll,yllcorner=yll,
                         delta_x=delta_x,delta_y=delta_y,
                         N=N,M=M,N_celle=N_celle,
                         nstr=nstr,
                         zratio=zr)



simu.hapin
simu.dem_parameters

simu.run_preprocessor(verbose=True,KeepOutlet=False) # remove outlet as it is a lab experiment wirth no flux outside the rhizotron


# simu.create_parm(NPRT=4,TIMPRTi=[1800.,3600.,7200.,80000.])
# simu.parm




# simu.set_elecs(xyz=xyzb[:,0:2])

# # xyz_drip=np.loadtxt('/home/ben/Documents/CATHY/rhizoLab/CATHY_additional_files_rhizo/drippers.txt', delimiter=',')
# # simu.set_drippers(xyz=xyz_drip)

# simu.create_infitration('inputs_rhizo')
# simu.drippers
# # create_infitration(drippersPos=[],RWU=False)
# # # rhizo.create_DA(drippersPos=[],RWU=False)



simu.run_processor(verbose=True)

# %% Explore outputs
# os.getcwd()
# import time
# import glob
# import pyvista as pv


# cplt.showvtkTL('./rhizo_prj/vtk/100.vtk')

cplt.showvtk('./rhizo_prj/vtk/100.vtk')
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
