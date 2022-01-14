"""Class managing Data Assimilation process analysis using enkf
"""
# Credit to:     https://github.com/hickmank/pyda/blob/7a2f04bd752e9c75bc8dcd2a45b21ff549736fa6/pyda/analysis_generator/kf/enkf1.py

import os
import matplotlib.pyplot as plt
from matplotlib import pyplot
import numpy as np
import shutil


def enkf_analysis(Data,DataCov,Param,Ensemble,Observation):
    '''

    Parameters
    ----------
    Data : TYPE
        DESCRIPTION.
    DataCov : TYPE
        DESCRIPTION.
    Param : TYPE
        DESCRIPTION.
    Ensemble : TYPE
        DESCRIPTION.
    Observation : TYPE
        DESCRIPTION.

    Returns
    -------
    list
        DESCRIPTION.

    '''
    # Collect data sizes.
    EnSize = Ensemble.shape[1]
    SimSize = Ensemble.shape[0] 
    MeaSize = Data.shape[0]

    # First combine the Ensemble and Param arrays.
   
    if len(Param)>0:
        A = np.vstack([Ensemble, Param.transpose()])
    else:
        A = Ensemble


    # Calculate mean of the ensemble
    Amean = (1./float(EnSize))*np.tile(A.sum(1), (EnSize,1)).transpose()
                    
    # Calculate ensemble perturbation from mean
    # Apert should be (SimSize+ParSize)x(EnSize)
    dA = A - Amean

    # Data perturbation from ensemble measurements
    # dD should be (MeasSize)x(EnSize)
    dD = Data - Observation.T
    
    print('dD')
    print(dD)

    # Ensemble measurement perturbation from ensemble measurement mean.
    # S is (MeasSize)x(EnSize)
    MeasAvg = (1./float(EnSize))*np.tile(Observation.reshape(MeaSize,EnSize).sum(1), (EnSize,1)).transpose()
    S = Observation - MeasAvg

    # Set up measurement covariance matrix
    # COV is (MeasSize)x(MeasSize)
    # print(DataCov)
    # DataCov = []
    if np.shape(DataCov)[1]>1:
        COV = DataCov
    else:
        print(np.shape( (1./float(EnSize-1))*np.dot(S,S.transpose())))
        COV = (1./float(EnSize-1))*np.dot(S,S.transpose()) + DataCov
        print(np.shape(COV))

        
    # Compute inv(COV)*dD
    # Should be (MeasSize)x(EnSize)
    B = np.linalg.solve(COV,dD.T)

    # Adjust ensemble perturbations
    # Should be (SimSize+ParSize)x(MeasSize)
    dAS = (1./float(EnSize-1))*np.dot(dA,S.transpose())

    # Compute analysis
    # Analysis is (SimSize+ParSize)x(EnSize)
    Analysis = A + np.dot(dAS,B)

    # Separate and return Analyzed Ensemble and Analyzed Parameters.
    AnalysisParam = Analysis[SimSize:,:].transpose()
    Analysis = Analysis[0:SimSize,:]
    
            
    # return [Analysis,AnalysisParam]
    return [A, Amean, dA, dD, MeasAvg, S, COV, B, dAS, Analysis, AnalysisParam]

    

def enkf_analysis_inflation(Data,DataCov,Param,Ensemble,Observation):

    alpha = 1  #<alpha> = (scalar) covariance inflation parameter. Usually alpha >= 1.
    
    # Collect data sizes.
    EnSize = Ensemble.shape[1]
    SimSize = Ensemble.shape[0] 
    MeaSize = Data.shape[0]

    # First combine the Ensemble and Param arrays.
    # A is (SimSize+ParSize)x(EnSize)
    A = np.vstack([Ensemble, Param.transpose()])

    # Calculate ensemble mean
    Amean = (1./float(EnSize))*np.tile(A.sum(1), (EnSize,1)).transpose()

    # Calculate ensemble observation mean
    MeasAvg = (1./float(EnSize))*np.tile(Observation.reshape(MeaSize,EnSize).sum(1), (EnSize,1)).transpose()
    
    # Inflate only the simulation ensemble, not the parameter ensemble 
    A[:(SimSize+1),:] = np.sqrt(alpha)*(A[:(SimSize+1),:] - Amean[:(SimSize+1),:]) + Amean[:(SimSize+1),:]

    # Inflate the ensemble observations
    Observation = np.sqrt(alpha)*(Observation - MeasAvg) + MeasAvg

    # Calculate ensemble perturbation from mean
    # Apert should be (SimSize+ParSize)x(EnSize)
    dA = A - Amean

    # Data perturbation from ensemble measurements
    # dD should be (MeasSize)x(EnSize)
    dD = Data - Observation

    # Ensemble measurement perturbation from ensemble measurement mean.
    # S is (MeasSize)x(EnSize)
    S = Observation - MeasAvg

    # Set up measurement covariance matrix
    # COV is (MeasSize)x(MeasSize)
    try:
        np.shape(DataCov)[1]
        COV = DataCov
    except:
        COV = (1./float(EnSize-1))*np.dot(S,S.transpose()) + DataCov

    # Compute inv(COV)*dD
    # Should be (MeasSize)x(EnSize)
    B = np.linalg.solve(COV,dD)

    # Adjust ensemble perturbations
    # Should be (SimSize+ParSize)x(MeasSize)
    dAS = (1./float(EnSize-1))*np.dot(dA,S.transpose())

    # Compute analysis
    # Analysis is (SimSize+ParSize)x(EnSize)
    Analysis = A + np.dot(dAS,B)

    # Separate and return Analyzed Ensemble and Analyzed Parameters.
    AnalysisParam = Analysis[SimSize:,:].transpose()
    Analysis = Analysis[0:SimSize,:]
            
    return  [A, Amean, dA, dD, MeasAvg, S, COV, B, dAS, Analysis, AnalysisParam]


