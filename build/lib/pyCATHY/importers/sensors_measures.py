"""Reader for sensors measured data
""" 

import os
import numpy as np
import pandas as pd



def read_ERT(filename, **kwargs):
    '''
    Reader for csv exported resipy data container type

    Returns
    -------
    None.

    '''
    
    
    df_ERT = pd.read_csv(filename, sep=",", header='infer')    
    # PH: pressure head
    # SW: Soil Water
    # CKRW: Relative hydraulic conductivity output at all nodes
    # QTRANIE Root‐zone water uptake at current time level always positive, changes sign in BCPIC
    # dict_vp = {'PH': 'pressure head'}


    # transform a numpy array into panda df
    # ------------------------------------------------------------------------    
    
    
    dict_ERT = {}
    dict_ERT['electrodesXYZ'] = []
    dict_ERT['sequenceABMN'] = []

    return df_ERT   


def read_discharge(filename, **kwargs):
    '''

    Returns
    -------
    None.

    '''
    
    # discharge_file = open(filename, "r")
    # discharge = pd.read(discharge_file, skiprows=0, usecols=range(2))    
    df_discharge = pd.read_csv(filename, sep="\t", header='infer')
    # discharge_file.close()
    
    

    
    return df_discharge   



def read_tensiometers(filename, **kwargs):
    '''

    Returns
    -------
    None.

    '''
    
    tensio_file = open(filename, "r")

    
    tensio_file.close()
    
    

    
    return df_tensio   



def read_TDR_prob(filename, **kwargs):
    '''

    Returns
    -------
    None.

    '''
    
    TDR_file = open(filename, "r")

    
    TDR_file.close()
    
    

    
    return df_TDR  


