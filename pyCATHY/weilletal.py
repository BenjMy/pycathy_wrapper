#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 16:36:58 2021

@author: ben
"""

# cathy project rhizo lab
import matplotlib.pyplot as plt
from cathy_tools import CATHY
import plottools as cplt
#import cathy_plot as ctplot

import os
import numpy as np


# %% Run cathy processor
# this function is a wrapper. It builds a makefile and the executable.
# Then the exe is run and use core CATHY file to solve the problem

path2prj ='test_weilletal'

simu = CATHY(dirName=path2prj,clear_src=False, clear_outputs=True)
simu.project_name = 'my_cathy_prj'



# nstr=2
# zr=list(np.ones(nstr)/nstr)


# The values are those of a 200-min rainfall event at a uniform intensity of 
# 3.3·10-4 m/min, followed by 100 min of drainage.


# here we can change the hap.in and dem_13.val files
# simu.update_prepo_inputs(delta_x=1,delta_y=1)

simu.run_preprocessor(verbose=True)


simu_time_max = 30

t_atmbc=[0.,12e3,18e3]
v_atmbc=[5.5e-06, 0.0,  0.0]
# simu.update_atmbc()
# simu.atmbc
simu.update_atmbc(HSPATM=1,IETO=1,TIME=t_atmbc,VALUE=v_atmbc)   
fig, ax = plt.subplots()
ax.plot(t_atmbc,v_atmbc)


# generate xyz , grid2d.exp grid3d
simu.run_processor(verbose=True,IPRT1=3)


grid3d_file = open('./my_cathy_prj/output/grid3d', 'r')
nnod,nnod3,nel = np.loadtxt('./my_cathy_prj/output/grid3d',max_rows=1)
mesh3d_nodes = np.loadtxt('./my_cathy_prj/output/grid3d',skiprows=1,max_rows=int(nel))
mesh_tetra = np.loadtxt('./my_cathy_prj/output/grid3d',skiprows=1+int(nel))

xmesh = mesh_tetra[:,0]
ymesh = mesh_tetra[:,1]
zmesh = mesh_tetra[:,2]


# nodi nella mesh 2-D (NNOD)
# numero di nodi totale in tutta la mesh 3-D (N)
# numero totale di elementi tetraedrici (NT)


#    c     Read grid3d
#  read(28,*) nnod,nnod3,nel
#  do i=1,nel
#     read(28,*) 
#  enddo
#  a=0
#  do i=1,nnod3
#     read(28,*) x(i),y(i),z(i)
#     if ((x(i).eq.0).or.(x(i).eq.5).or.(y(i).eq.0).or.
# 1       (y(i).eq.5))then
#     a=a+1
#     endif
#  enddo

         
Lines = grid3d_file.readlines()
       
x = Lines[0].strip().split(" ")


simu.update_ic(INDP=2,IPOND=0,WTHEIGHT=0)
simu.ic


simu.update_soil(INDP=2,IPOND=0,WTHEIGHT=0)
simu.soil




ax.set(xlabel='time (s)', ylabel='Q (m/s)',
       title='atmbc')
ax.grid()


# here we can update the parm file directly passing new arguments
simu.run_processor(verbose=True,TMAX=simu_time_max,
                   TIMPRTi=[1,3600.,7200.])



# %% Explore outputs

cplt.showvtk('./my_cathy_prj/vtk/100.vtk',notebook=True)
cplt.showvtk(unit='saturation',timeStep=1,notebook=True, path='./my_cathy_prj/vtk/')
cplt.showvtk(unit='saturation',timeStep=2,notebook=True, path='./my_cathy_prj/vtk/')

cplt.showvtk(unit='pressure',timeStep=2,notebook=True, path='./my_cathy_prj/vtk/')

# cumflowvol
# dtcoupling
# grid2d.exp
# risul


# %% test
