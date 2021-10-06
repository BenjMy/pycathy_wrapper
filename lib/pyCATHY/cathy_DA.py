#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 12:12:41 2021

@author: ben
"""

import os
import matplotlib.pyplot as plt
import numpy as np

# -------------------------------------------------------------------#
#%% DATA ASSIMILATION


def create_archie():
    pass


# NENS
# Number of realizations in the ensemble kalman filter (EnKF) and the particle filter (SIR). If a realization does not converge, NENS is decreased by one.


def perturbate_var(var, parameter, mean, sd, per_type, sampling_type = 'lognormal',
                        ensemble_size = 128):
    """
    Perturbate variable for the generation of the ensemble condition

    Parameters
    ----------
    var : dict
        DESCRIPTION.
    parameter : str
        specify the parameter type.
    mean : TYPE
        DESCRIPTION.
    sd : TYPE
        DESCRIPTION.
    per_type : TYPE
        DESCRIPTION.
    sampling_type : TYPE, optional
        DESCRIPTION. The default is 'lognormal'.
    ensemble_size : TYPE, optional
        DESCRIPTION. The default is 128.

    Returns
    -------
    var_per : np.array([])
        Perturbated variable

    """


    var_sampling = np.random.lognormal(mean, sigma=sd, size=ensemble_size)
    
    count, bins, ignored = plt.hist(var_sampling, ensemble_size, density=True, align='mid')
    
    x = np.linspace(min(bins), max(bins), ensemble_size)
    pdf = (np.exp(-(np.log(x) - mean)**2 / (2 * sd**2))
           / (x * sd * np.sqrt(2 * np.pi)))
    
    plt.plot(x, pdf, linewidth=2, color='r')
    plt.axis('tight')
    plt.show()

    var_mat = np.ones(ensemble_size)*var['kss']
    
    if per_type == 'multiplicative':
        var_per = var_mat*var_sampling
    elif per_type == 'additive':
        var_per = var_mat+var_sampling

    
    
    count, bins, ignored = plt.hist(var_per, ensemble_size, density=True, align='mid')
    


    return var_per 





def create_elec_nodes():

    pass
