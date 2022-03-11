"""Reader for sensors measured data
""" 

import os
import numpy as np
import pandas as pd
import pygimli as pg



def read_ERT(filename, data_format, **kwargs):
    '''
    Reader for csv exported resipy data container type

    Returns
    -------
    None.

    '''
    
    dict_ERT = {}
    
    if 'resipy' in data_format:
        df_ERT = pd.read_csv(filename, sep=",", header='infer')       
        # dict_ERT = {}
        # dict_ERT['electrodesXYZ'] = []
        # dict_ERT['sequenceABMN'] = []
    elif 'pygimli' in data_format:
        df_ERT = pg.load(filename)  
        # df_ERT['a']
        dict_ERT['elecs'] = np.array(df_ERT.sensorPositions())
        # df_ERT.sensors()
        # np.array(df_ERT.sensorPositions())
        
        
    else:
        raise ValueError('ERT data format not recognized')
        

    return df_ERT, dict_ERT 


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


