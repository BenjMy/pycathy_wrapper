"""Class managing Data Assimilation process analysis using enkf
"""

import os
import matplotlib.pyplot as plt
from matplotlib import pyplot
import numpy as np
import shutil


def enkf_analysis(Data,DataCov,Param,Ensemble,Observation):
    '''
    https://github.com/hickmank/pyda/blob/7a2f04bd752e9c75bc8dcd2a45b21ff549736fa6/pyda/analysis_generator/kf/enkf1.py

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
    # A is (SimSize+ParSize)x(EnSize)
    A = np.vstack([Ensemble, Param.transpose()])

    # Calculate mean
    Amean = (1./float(EnSize))*np.tile(A.sum(1), (EnSize,1)).transpose()

    # Calculate ensemble perturbation from mean
    # Apert should be (SimSize+ParSize)x(EnSize)
    dA = A - Amean

    # Data perturbation from ensemble measurements
    # dD should be (MeasSize)x(EnSize)
    dD = Data - Observation.T

    # Ensemble measurement perturbation from ensemble measurement mean.
    # S is (MeasSize)x(EnSize)
    MeasAvg = (1./float(EnSize))*np.tile(Observation.reshape(MeaSize,EnSize).sum(1), (EnSize,1)).transpose()
    S = Observation - MeasAvg

    # Set up measurement covariance matrix
    # COV is (MeasSize)x(MeasSize)
    COV = (1./float(EnSize-1))*np.dot(S,S.transpose()) + DataCov

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


def ENKF(self,A,dAS,B):
    '''
    THIS SHOULD BE MOVED TO DA CLASS 

    ENKF from Sakov and Eversen (2009)
    
    .. math::

        X^{u} = x^{u} + A^{u}

    x: ensemble average
    A: matrix of ensemble anomalies
    superscript u means updated
    
    Returns
    -------
    None.

    '''
    
    self.console.print('ENKF')
    
    # ensemble mean of the simulated observations
    
    

    # Compute analysis
    # Analysis is (SimSize+ParSize)x(EnSize)
    X_augm_analysis = A + np.dot(dAS,B)

    # Separate and return Analyzed Ensemble and Analyzed Parameters.
    # AnalysisParam = Analysis[SimSize:,:].transpose()
    # Analysis = Analysis[0:SimSize,:]
    
    
    
    
    # %*% is matrix multiplication
    # C = A.dot(B) in python
    
    # X be an ensemble matrix of M rows and N columns, where
    # N is the number of realizations and M is the state dimen-
    # sion, i.e., the number of nodes in the finite element grid, aug-
    # mented by the number of parameters that are subject to up-
    # date.
    
    # x: ensemble average
    # A: matrix of ensemble anomalies
    # suscript u means updated
    
    

    # R: measurement error covariance matrix
    # D: difference between the measurements,
    # Hx: ensemble mean of the simulated observations
    # S is the matrix of scaled ensemble innovation anomalies
    # HA being the simulated measurement anomalies  
    
    # When updating the states only, the elements of X are the 
    # pressure heads at each node of the finite element grid,

    #STEP 1 Estimation matrix [A] which represents the difference 
    # between matrix [X] and the mean values estimated over all the 
    # scenarios
    
    # update_ic
    # X_u = []
    
    
    return X_a

