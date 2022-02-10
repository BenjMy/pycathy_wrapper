#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 12:12:41 2021

@author: ben
"""

import pyvista as pv
import glob
import time 
import os
import matplotlib.pyplot as plt
import pandas as pd





def transform2_time_delta(t,x_units):
    '''
    Time to time delta

    Returns
    -------
    None.

    '''

    delta_t = pd.to_timedelta(t,unit=x_units) 

    return delta_t


# t = [0]
# dt = transform2_time_delta(t,'s')
# print(dt[0])


def convert_time_units(t, x_units):
    '''
    convert time units

    Parameters
    ----------
    t : TYPE
        DESCRIPTION.
    x_units : str
        DESCRIPTION.

    Returns
    -------
    xlabel : TYPE
        DESCRIPTION.
    t_new : TYPE
        DESCRIPTION.

    '''

    xlabel = " (s)"
    if x_units == "days":
        xlabel = " (days)"
        t_new = [x / (24 * 60 * 60) for x in t]
        t_new = round(t_new[0], 2)
    if x_units == "hours":
        xlabel = " (h)"
        t_new = [x / (60 * 60) for x in t]
        t_new = round(t_new[0], 1)

    return xlabel, t_new


def label_units(units,**kwargs):
    '''
    label units


    Parameters
    ----------
    units : TYPE
        DESCRIPTION.

    Returns
    -------
    label

    '''

    if units == "SW":
        label = "Soil Water Content \n ($m^{3}/m^{3}$)"
    elif units == "PH":
        label = "Pressure head (m)"
    elif units == "CKRW":
        label = "Relative hydraulic conductivity"
    elif units == "QTRANIE":
        label = "Root‐zone water uptake"
    else:
        label = ''                
        
    return label




def infer_VGP_literature(soil_type):
    
    SPP = []
    # SPP = {
    #     "PERMX": PERMX,
    #     "PERMY": PERMY,
    #     "PERMZ": PERMZ,
    #     "ELSTOR": ELSTOR,
    #     "POROS": POROS,
    #     "VGNCELL": VGNCELL,
    #     "VGRMCCELL": VGRMCCELL,
    #     "VGPSATCELL": VGPSATCELL,
    # }
                
    return SPP
    
    



