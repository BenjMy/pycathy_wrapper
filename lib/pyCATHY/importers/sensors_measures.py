"""Reader for sensors measured data
""" 

import os
import numpy as np
import pandas as pd



def read_ERT(filename, **kwargs):
    '''

    Returns
    -------
    None.

    '''
    
    
    df_ERT = pd.read_csv(filename, sep="\t", header='infer')

    
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


