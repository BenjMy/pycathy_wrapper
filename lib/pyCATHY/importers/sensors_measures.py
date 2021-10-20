"""Reader for sensors measured data
""" 

import os
import numpy as np



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


