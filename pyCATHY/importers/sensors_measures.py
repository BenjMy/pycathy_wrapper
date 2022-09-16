"""Reader for sensors measured data
""" 

import os
import numpy as np
import pandas as pd

try: 
    import pygimli as pg
except ImportError: 
    pygimli = None



def read_ERT(filename, data_format, **kwargs):
    '''
    Reader for csv exported resipy data container type

    Returns
    -------
    None.

    '''
    
    dict_ERT = {}
    
    if 'resipy' in data_format:
        if resipy is None:
            raise ValueError(f'resipy module not imported. Please pip install.')

        df_ERT_new = pd.read_csv(filename, sep=",", header='infer')       
    elif 'pygimli' in data_format:
        if pygimli is None:
            raise ValueError(f'pygimli module not imported. Please pip install.')

        df_ERT = pg.load(filename)  
        
        df_ERT_new = pd.DataFrame([df_ERT['a'],df_ERT['b'],
                              df_ERT['k'],df_ERT['m'],
                              df_ERT['n'],df_ERT['r'],
                              df_ERT['rhoa'],df_ERT['valid']
                              ]
                              )
        df_ERT_new = df_ERT_new.T
        df_ERT_new.columns=['a', 'b', 'k', 'm', 'n', 'r', 'rhoa', 'valid']
        dict_ERT['elecs'] = np.array(df_ERT.sensorPositions())
        
    elif 'custum' in data_format:
        df_ERT = pg.load(filename)  
      
 
    else:
        raise ValueError('ERT data format not recognized')
        

    return df_ERT_new, dict_ERT 


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


