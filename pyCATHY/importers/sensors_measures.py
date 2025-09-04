"""Reader for sensors measured data
"""

import os

import numpy as np
import pandas as pd

try:
    import pygimli as pg
    from pygimli.physics import ert
except ImportError:
    pg = None

try:
    from resipy import Project
except ImportError:
    Project = None

try:
    import emagpy
    from emagpy import Problem
except ImportError:
    emagpy = None



def read_EM(filename, data_format, **kwargs):
    """
    Reader for emagpy data container type

    Returns
    -------
    None.

    """
    
    timeLapse = False
    if timeLapse in kwargs:
        timeLapse = kwargs.pop('timeLapse')
        
    k = Problem()
    k.createSurvey(str(filename))
    df_EM = k.surveys[0].df
    
    return df_EM,  k 


def read_ERT(filename, data_format, **kwargs):
    """
    Reader for csv exported resipy data container type

    Returns
    -------
    None.

    """

    dict_ERT = {}

    if "resipy" in data_format:
        if Project is None:
            raise ValueError(f"resipy module not imported. Please pip install.")

        df_ERT_new = pd.read_csv(filename, sep=",", header="infer")
        
    elif "pygimli" in data_format:
        if pg is None:
            raise ValueError(f"pygimli module not imported. Please pip install.")

        if '.csv' in filename:
            df_ERT_new = pd.read_csv(filename, sep=",", header="infer")
            
            if 'rhoa' not in df_ERT_new.columns:
                df_ERT_new['rhoa'] = df_ERT_new['K']*df_ERT_new['resist']
                
            # df_ERT_new['Rho']
            # df_ERT_new['resist']

        else:
            df_ERT = pg.load(filename)
            if np.sum(df_ERT['rhoa'])==0:
                if np.sum(df_ERT['k'])==0:
                    # df_ERT['k'] = ert.createGeometricFactors(df_ERT, 
                    #                                          numerical=True
                    #                                              )
                    df_ERT['k'] = ert.createGeometricFactors(df_ERT)
                    df_ERT['rhoa'] = df_ERT['k']*df_ERT['r']
                    # np.min(df_ERT['rhoa'])
                    # pg.show(df_ERT)
                    
            header = ["a", "b", "m", "n", "k", "r", "rhoa",
                      'err', 'rec_err','valid']
            dict_ERT_new = {}
            for hh in header:
                try:
                    dict_ERT_new[hh] = df_ERT[hh]
                except:
                    pass

            df_ERT_new = pd.DataFrame(dict_ERT_new)
            dict_ERT["elecs"] = np.array(df_ERT.sensorPositions())

    elif "custum" in data_format:
        df_ERT = pg.load(filename)

    else:
        raise ValueError("ERT data format not recognized")

    return df_ERT_new, dict_ERT


def read_discharge(filename, **kwargs):
    """

    Returns
    -------
    None.

    """

    # discharge_file = open(filename, "r")
    # discharge = pd.read(discharge_file, skiprows=0, usecols=range(2))
    df_discharge = pd.read_csv(filename, sep="\t", header="infer")
    # discharge_file.close()

    return df_discharge


def read_tensiometers(filename, **kwargs):
    """

    Returns
    -------
    None.

    """

    tensio_file = open(filename, "r")
    tensio_file.close()

    return df_tensio


def read_TDR_prob(filename, **kwargs):
    """

    Returns
    -------
    None.

    """

    TDR_file = open(filename, "r")

    TDR_file.close()

    return df_TDR
