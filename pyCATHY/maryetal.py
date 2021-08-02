#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 16:36:58 2021

@author: ben
"""

# cathy project rhizo lab
from cathy_tools import CATHY
from rhizo_tools import create_infitration # create_DA

import plottools as cplt
#import cathy_plot as ctplot

import os
import numpy as np


# %% Run cathy processor
# this function is a wrapper. It builds a makefile and the executable.
# Then the exe is run and use core CATHY file to solve the problem
            
path2prj = os.getcwd()

simu = CATHY(dirName=path2prj,prjName='rhizo_prj')

dem_test =  np.ones([15,3])
simu.create_DEM_mesh(DEM=dem_test,
                      N=np.shape(dem_test)[0],
                      M=np.shape(dem_test)[1],
                      N_celle=np.shape(dem_test)[0]*np.shape(dem_test)[1])


simu.run_preprocessor(verbose=True)


simu.create_parm(NPRT=4,TIMPRTi=[1800.,3600.,7200.,80000.])

# xyzb=np.loadtxt('/home/ben/Documents/CATHY/rhizoLab/CATHY_additional_files_rhizo/elecsXYZ_Rhizo_LAB.csv',skiprows=1,delimiter=',')
# simu.set_elecs(xyz=xyzb[:,0:2])

# xyz_drip=np.loadtxt('/home/ben/Documents/CATHY/rhizoLab/CATHY_additional_files_rhizo/drippers.txt', delimiter=',')
# simu.set_drippers(xyz=xyz_drip)

# create_infitration()
# create_infitration(drippersPos=[],RWU=False)
# # rhizo.create_DA(drippersPos=[],RWU=False)

    
    
simu.run_processor(verbose=True)


# parm_file = open('./rhizo_prj/input/parm', 'r')
# Lines = parm_file.readlines()

# x = Lines[0].strip().split(" ")
# print(Lines[0])
# tmp_id2rmv=[]
# for j in enumerate(list(x)):
#     print(j)
#     if len(j[1])<1:
#         tmp_id2rmv.append(j[0])
    
# line0=np.delete(x,tmp_id2rmv)
# testDIGIT = ([x for x in line0 if not str(x).isdigit()])



# print(line0)
# x2 = line0.strip().split(" ")

    
# import re
# m = re.search(r"\d", str(line0))
# print(m)

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


# %% test       



    
    