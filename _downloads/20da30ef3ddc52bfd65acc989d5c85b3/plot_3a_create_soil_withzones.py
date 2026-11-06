#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Soil 3d from a Digital Elevation Model (DEM)
============================================

Weill, S., et al. « Coupling Water Flow and Solute Transport into a Physically-Based Surface–Subsurface Hydrological Model ». 
Advances in Water Resources, vol. 34, no 1, janvier 2011, p. 128‑36. DOI.org (Crossref), 
https://doi.org/10.1016/j.advwatres.2010.10.001.

This example shows how to use pyCATHY object to build a 3d soil properties from a DEM and run the hydrological model.

*Estimated time to run the notebook = 5min*

"""


# !! run preprocessor change the DEM shape !
# dtm_13 does not have the same shape anymore!

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import pyCATHY.meshtools as mt
from pyCATHY import cathy_tools
from pyCATHY.importers import cathy_inputs as in_CT
from pyCATHY.importers import cathy_outputs as out_CT
from pyCATHY.plotters import cathy_plots as cplt

#%% Init CATHY model
# ------------------------
path2prj = "../SSHydro/"  # add your local path here
simu = cathy_tools.CATHY(dirName=path2prj, 
                         prj_name="soil_withzones", 
                         clear_src=False
                         )

rootpath = os.path.join(simu.workdir + simu.project_name)
simu.run_preprocessor(verbose=False)
simu.run_processor(IPRT1=3,verbose=True)

#%% Define zones and create 

simu.DEM
zones = np.ones(np.shape(simu.DEM))
zones[:,0:2] = 2
zones[:,2:4] = 3
zones[:,4:6] = 4

simu.update_zone(zones)
simu.show_input('zone')

#%% Create an empty dataframe of SPP and set default SPP properties 

df_SPP_map = simu.init_soil_SPP_map_df(nzones=4,nstr=15)
SPP_map = simu.set_SOIL_defaults(SPP_map_default=True)

#%% Update soil file

simu.update_soil(SPP_map=SPP_map)


