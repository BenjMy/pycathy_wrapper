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


# %% Explore outputs



os.chdir('./test_weilletal/my_cathy_prj/vtk')

cplt.showvtk(unit='saturation',timeStep=1,notebook=True)

# cplt.showvtkTL(unit='saturation',notebook=True)

cplt.showvtk('100.vtk',notebook=True)
