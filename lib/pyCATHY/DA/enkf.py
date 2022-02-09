"""Class managing data Assimilation process analysis using enkf
"""
# Credit to:     https://github.com/hickmank/pyda/blob/7a2f04bd752e9c75bc8dcd2a45b21ff549736fa6/pyda/analysis_generator/kf/enkf1.py

import os
import matplotlib.pyplot as plt
from matplotlib import pyplot
import numpy as np
import shutil


def enkf_analysis(data,data_cov,param,ensemble,observation):
    '''

    parameters
    ----------
    data : TYPE
        DESCRIPTION.
    data_cov : TYPE
        DESCRIPTION.
    param : TYPE
        DESCRIPTION.
    ensemble : TYPE
        DESCRIPTION.
    observation : TYPE
        DESCRIPTION.

    Returns
    -------
    list
        DESCRIPTION.

    '''
    # Collect data sizes.
    ens_size = ensemble.shape[1]
    sim_size = ensemble.shape[0] 
    meas_size = np.array([data]).T.shape[0]

    # First combine the ensemble and param arrays.
   
    if len(param)>0:
        A = np.vstack([ensemble, param.transpose()])
    else:
        A = ensemble


    # Calculate mean of the ensemble
    Amean = (1./float(ens_size))*np.tile(A.sum(1), (ens_size,1)).transpose()
                    
    # Calculate ensemble perturbation from mean
    # Apert should be (sim_size+ParSize)x(ens_size)
    dA = A - Amean

    # data perturbation from ensemble measurements
    # dD should be (MeasSize)x(ens_size)
    dD = data - observation

    # ensemble measurement perturbation from ensemble measurement mean.
    # S is (MeasSize)x(ens_size)
    MeasAvg = (1./float(ens_size))*np.tile(observation.reshape(meas_size,ens_size).sum(1), (ens_size,1))
    S = observation - MeasAvg

    # Set up measurement covariance matrix
    # COV is (MeasSize)x(MeasSize)
    # print(data_cov)
    # data_cov = []
    if np.shape(data_cov)[1]>1:
        COV = data_cov
    else:
        # print(np.shape( (1./float(ens_size-1))*np.dot(S,S.transpose())))
        COV = (1./float(ens_size-1))*np.dot(S,S.transpose()) + data_cov
        # print(np.shape(COV))

        
    # Compute inv(COV)*dD
    # Should be (MeasSize)x(ens_size)
    B = np.linalg.solve(COV,dD.T)

    # Adjust ensemble perturbations
    # Should be (sim_size+ParSize)x(MeasSize)
    dAS = (1./float(ens_size-1))*np.dot(dA,S)

    # Compute analysis
    # Analysis is (sim_size+ParSize)x(ens_size)
    Analysis = A + np.dot(dAS,B)

    # Separate and return Analyzed ensemble and Analyzed parameters.
    Analysisparam = Analysis[sim_size:,:].transpose()
    Analysis = Analysis[0:sim_size,:]
    
            
    # return [Analysis,Analysisparam]
    return [A, Amean, dA, dD, MeasAvg, S, COV, B, dAS, Analysis, Analysisparam]

    

def enkf_analysis_inflation(data,data_cov,param,ensemble,observation):

    alpha = 1  #<alpha> = (scalar) covariance inflation parameter. Usually alpha >= 1.
    
    # Collect data sizes.
    ens_size = ensemble.shape[1]
    sim_size = ensemble.shape[0] 
    meas_size = data.shape[0]

    # First combine the ensemble and param arrays.
    # A is (sim_size+ParSize)x(ens_size)
    A = np.vstack([ensemble, param.transpose()])

    # Calculate ensemble mean
    Amean = (1./float(ens_size))*np.tile(A.sum(1), (ens_size,1)).transpose()

    # Calculate ensemble observation mean
    MeasAvg = (1./float(ens_size))*np.tile(observation.reshape(meas_size,ens_size).sum(1), (ens_size,1)).transpose()
    
    # Inflate only the simulation ensemble, not the parameter ensemble 
    A[:(sim_size+1),:] = np.sqrt(alpha)*(A[:(sim_size+1),:] - Amean[:(sim_size+1),:]) + Amean[:(sim_size+1),:]

    # Inflate the ensemble observations
    observation = np.sqrt(alpha)*(observation - MeasAvg) + MeasAvg

    # Calculate ensemble perturbation from mean
    # Apert should be (sim_size+ParSize)x(ens_size)
    dA = A - Amean

    # data perturbation from ensemble measurements
    # dD should be (MeasSize)x(ens_size)
    dD = data - observation

    # ensemble measurement perturbation from ensemble measurement mean.
    # S is (MeasSize)x(ens_size)
    S = observation - MeasAvg

    # Set up measurement covariance matrix
    # COV is (MeasSize)x(MeasSize)
    try:
        np.shape(data_cov)[1]
        COV = data_cov
    except:
        COV = (1./float(ens_size-1))*np.dot(S,S.transpose()) + data_cov

    # Compute inv(COV)*dD
    # Should be (MeasSize)x(ens_size)
    B = np.linalg.solve(COV,dD)

    # Adjust ensemble perturbations
    # Should be (sim_size+ParSize)x(MeasSize)
    dAS = (1./float(ens_size-1))*np.dot(dA,S.transpose())

    # Compute analysis
    # Analysis is (sim_size+ParSize)x(ens_size)
    Analysis = A + np.dot(dAS,B)

    # Separate and return Analyzed ensemble and Analyzed parameters.
    Analysisparam = Analysis[sim_size:,:].transpose()
    Analysis = Analysis[0:sim_size,:]
            
    return  [A, Amean, dA, dD, MeasAvg, S, COV, B, dAS, Analysis, Analysisparam]


