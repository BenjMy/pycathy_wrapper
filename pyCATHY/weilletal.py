#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 16:36:58 2021

@author: ben
"""

# cathy project rhizo lab
from cathy_tools import CATHY
import plottools as cplt
#import cathy_plot as ctplot

import os
import numpy as np


# %% Run cathy processor
# this function is a wrapper. It builds a makefile and the executable.
# Then the exe is run and use core CATHY file to solve the problem
            
path2prj ='test_weilletal'

simu = CATHY(dirName=path2prj)
simu.project_name = 'my_cathy_prj'
# survey.create_DEM_mesh(delta_x=100)
    
simu.run_preprocessor(verbose=True)

simu.run_processor(verbose=True)



# %% Explore outputs
                    
cplt.showvtk('./my_cathy_prj/vtk/100.vtk')


# %% test       



    
    