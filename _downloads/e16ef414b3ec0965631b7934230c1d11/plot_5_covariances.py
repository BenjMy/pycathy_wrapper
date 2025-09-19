"""
Observations covariances, Localisation matrix
===========================================

The notebook illustrate how to read EM sensors dataset to be prepare for DA

*Estimated time to run the notebook = 2min*

"""
import numpy as np
from pyCATHY.DA.cathy_DA import DA
import pandas as pd
import matplotlib.pyplot as plt
from pyCATHY.DA.cathy_DA import DA, dictObs_2pd
from pyCATHY.DA.observations import read_observations, prepare_observations, make_data_cov
from pathlib import Path
import pyvista as pv

#%% Create a CATHY project
# -----------------------
