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


# These are the van Genuchten (1980) equations
# The input is matric potential, psi and the hydraulic parameters.
# psi must be sent in as a numpy array.
# The pars variable is like a MATLAB structure.
# https://github.com/simpeg/simpeg/blob/main/SimPEG/flow/richards/empirical.py#L55-L69

def thetaFun(psi,pars):
  Se=(1+abs(psi*pars['alpha'])**pars['n'])**(-pars['m'])
  Se[psi>=0]=1.
  return pars['thetaR']+(pars['thetaS']-pars['thetaR'])*Se
  
def CFun(psi,pars):
  Se=(1+abs(psi*pars['alpha'])**pars['n'])**(-pars['m'])
  Se[psi>=0]=1.
  dSedh=pars['alpha']*pars['m']/(1-pars['m'])*Se**(1/pars['m'])*(1-Se**(1/pars['m']))**pars['m']
  # dSedh(psi>=0)=0;
  return Se*pars['Ss']+(pars['thetaS']-pars['thetaR'])*dSedh
  
def KFun(psi,pars):
  Se=(1+abs(psi*pars['alpha'])**pars['n'])**(-pars['m'])
  Se[psi>=0]=1.
  return pars['Ks']*Se**pars['neta']*(1-(1-Se**(1/pars['m']))**pars['m'])**2
  
def setpars():
  pars={}
  pars['thetaR']=float(raw_input("thetaR = "))
  pars['thetaS']=float(raw_input("thetaS = "))
  pars['alpha']=float(raw_input("alpha = "))
  pars['n']=float(raw_input("n = "))
  pars['m']=1-1/pars['n']
  pars['Ks']=float(raw_input("Ks = "))
  pars['neta']=float(raw_input("neta = "))
  pars['Ss']=float(raw_input("Ss = "))
  return pars



def PH2SW_VGN(SW):
    '''
    VGN
    '''
    
    SW = [] # saturation water


    return SW


def SW2PH_VGN(SW):
    '''
    VGN
    '''
    
    PH = [] # pressure heads

    return PH


def Leij_etal_1996():
    '''
    Description of soil physical parameters
    
    .. seealso::  Leij, F. J., W. J. Alves, and M. T. van Genuchten, 1996. The UNSODA unsaturated soil hydraulic
                    database: user’s manual, volume 96. National Risk Management Research Laboratory, Oﬃce of
                    Research and Development, US Environmental Protection Agency.

    Returns
    -------
    None.

    '''
    
def Carsel_and_Parrish_1988():
    '''
    Ref Rediel thesis
    Carsel_and_Parrish_1988
    '''

    medium = ['loam','loamy_sand','sandy_loam', 'silt_loam']
    names = {'medium', 'K0','alpha','n', 'tau', 'thetar', 'thetas'}

    K0= [2.9e-6,4.1e-5,1.2e-5,1.3e-6] # m/s
    alpha= [-3.6,-12.4,-7.5,-2.0] # 1/m
    n= [1.56,2.28,1.89,1.41] # []
    tau= [0.078,0.057,0.065,0.067] # []
    thetar= [0.5,0.5,0.5,0.5] # []
    thetas= [0.078,0.057,0.065,0.067] # []
    data = [medium, K0,alpha,n,tau,thetar,thetas]


    VG_CP_1988 = pd.DataFrame(data, columns=names)

    return VG_CP_1988


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
    
    



