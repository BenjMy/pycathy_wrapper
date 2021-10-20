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




def convert_time_units(t, x_units):
    '''
    convert time units

    Parameters
    ----------
    t : TYPE
        DESCRIPTION.
    x_units : TYPE
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








# def RWU(**kwargs):
#     """
#     Short summary.

#     """

#     # import numpy as np

#     # PCANA=1.0
#     # PCREF=-4.0
#     # PCWLT=-150.0
#     # PSI=np.arange(-200,10.0,0.1)

        
#     # ALFA= np.zeros(len(PSI))       
#     # for i in range(len(PSI)):
#     #     S1 = PCANA
#     #     S2 = max([0.0, (PCANA+1.0E-03)])
#     #     SH2O = PSI[i]
#     #     GX1 = (SH2O-PCWLT)/(PCREF-PCWLT)
#     #     GX1 = min([1.0, max([0.0, GX1])])
#     #     GX2 = 1.0 - (SH2O-S1)/(S2-S1)
#     #     GX2 = min([1.0, max([0.0, GX2])])
#     #     GX = min([GX1, GX2])
#     #     ALFA[i] = GX



#     return
