"""Class managing data Assimilation process analysis using enkf
"""
# Credit to:     https://github.com/hickmank/pyda/blob/7a2f04bd752e9c75bc8dcd2a45b21ff549736fa6/pyda/analysis_generator/kf/enkf1.py

import os
import matplotlib.pyplot as plt
from matplotlib import pyplot
import numpy as np
import shutil
import rich.console
from rich.progress import track
from rich import print


def enkf_analysis(data,data_cov,param,ensemble,observation,**kwargs):
    '''
    Case with non linear observation operator:
        
    ..note::
        $ K_{t} = P^{f}_{t}H^{T}(HP^{f}_{t}H^{T} + R^{T})^{-1} $
        
        H is the linear observation operator matrix $[N_{obs}\times N_{u}]$, 
        $P^f$ is the forecast error covariance matrix $[N_{u}\times N_{u}]$, 
        and R the measurement error covariance matrix $[N_{obs}\times N_{obs}]$.

    parameters
    ----------
    data : np.array([])
        stacked measured data.
    data_cov : TYPE
        measured data covariance matrice.
    param : TYPE
        model parameters (perturbated and to update).kwrags
    ensemble : np.array([])
        state values (can be either pressure heads or saturation water).
    observation : TYPE
        stacked predicted observation (after mapping) in the same physical quantity than data.

    Returns
    -------
    list
        DESCRIPTION.

    '''
    
    # get kwargs optionnal arguments
    # --------------------------------
    localize = False
    if 'localize' in kwargs:
        localize = True
        
        
    # Collect data sizes.
    #- -------------------------------
    ens_size = ensemble.shape[1]
    sim_size = ensemble.shape[0] 
    meas_size = data.shape[0]

    # observation = observation # (MeasSize)x(ens_size)
    # data = data # (MeasSize)

    # First combine the ensemble and param arrays.
    # ------------------------------------------------------------------------
    ensemble_mean = (1./float(ens_size))*np.tile(ensemble.sum(1), (ens_size,1)).transpose()
    
    if len(param)>0:
        
        # Take the mean line by line i.e. each line is an independant parameter to update
        
        param_mean = []
        for l in range(len(param)):
            param_mean.append((1./float(ens_size))*np.tile(param[l,:].sum(0), (ens_size,1)).transpose())
            # param_mean.append((1./float(ens_size))*np.tile(param[:,:].sum(0), (ens_size,1)).transpose())
    
        param_mean = np.vstack(param_mean)
        
        # print('mean pRam')
        # print(param_mean)
        
    if len(param)>0:
        augm_state_mean = np.vstack([ensemble_mean, param_mean])
        augm_state = np.vstack([ensemble, param])
        # augm_data =  np.vstack([data, param.T])
        # np.shape(augm_state)
        # np.shape(data)
        # np.shape(observation)
    else:
        augm_state_mean = ensemble_mean
        augm_state = ensemble


    # Calculate ensemble perturbation from mean
    # augm_state_pert should be (sim_size+ParSize)x(ens_size)
    augm_state_pert = augm_state - augm_state_mean
    
    # data perturbation from ensemble measurements
    # data_pert should be (MeasSize)x(ens_size)
    data_pert = (data - observation.T).T
    
    
    # print('data_pert')
    # print(data_pert)
        
        
        
    if np.max(abs(data_pert))>1e3:
        raise ValueError('predictions are too far from observations')
    # S = ensemble measurement perturbation from ensemble measurement mean.
    # S is (MeasSize)x(ens_size)        
    meas_avg = (1./float(ens_size))*np.tile(observation.reshape(meas_size,ens_size).sum(1), (ens_size,1)).transpose()
    obs_pert = observation - meas_avg
    
    if obs_pert.mean() == 0:
        print('Ensemble measurement perturbation from ensemble'
                'measurement mean is too small (=0):'
                '- Increase perturbation!'
                '- Check if the system is not in steady state'
                )
        
    # Set up measurement covariance matrix
    # COV is (MeasSize)x(MeasSize)
    # if np.shape(data_cov)[1]>1:
    #     COV = data_cov
    # else:
    # print(np.shape( (1./float(ens_size-1))*np.dot(S,S.transpose())))
    
    Sakov = False
    if Sakov:
        COV = data_cov.transpose()
    else:
        COV = (1./float(ens_size-1))*np.dot(obs_pert,obs_pert.transpose()) + data_cov.transpose()
        # COV = data_cov.transpose()

    # print(np.shape(COV))

    # np.shape(data_pert)
    # Compute inv(COV)*dD
    # Should be (MeasSize)x(ens_size)
    inv_data_pert = np.linalg.solve(COV,data_pert)

    # Adjust ensemble perturbations
    # Should be (sim_size+ParSize)x(MeasSize)
    pert = (1./float(ens_size-1))*np.dot(augm_state_pert,obs_pert.T)

    # Compute analysis
    # Analysis is (sim_size+ParSize)x(ens_size)
    analysis = augm_state + np.dot(pert,inv_data_pert)

    # Separate and return Analyzed ensemble and Analyzed parameters.
    analysis_param = analysis[sim_size:,:].transpose()
    analysis = analysis[0:sim_size,:]
    

    # print('analysis_param')
    # print(analysis_param)
    
    
    # return [Analysis,Analysisparam]
    return [augm_state, augm_state_mean, augm_state_pert, data_pert, meas_avg, 
            obs_pert, COV, inv_data_pert, pert, analysis, analysis_param]

    

def enkf_analysis_inflation(data,data_cov,param,ensemble,observation,**kwargs):

    alpha = 1  #<alpha> = (scalar) covariance inflation parameter. Usually alpha >= 1.
    
    if 'alpha' in kwargs:
        alpha = kwargs['alpha']
        
    
    
    # Collect data sizes.
    ens_size = ensemble.shape[1]
    sim_size = ensemble.shape[0] 
    
    
    observation = observation.T
    # data = np.array([data]).T
    meas_size = data.shape[0]

    # First combine the ensemble and param arrays.
    # A is (sim_size+ParSize)x(ens_size)
    A = np.vstack([ensemble, param])

    # np.shape(ensemble)
    # np.shape(param)
    
    # Calculate ensemble mean
    Amean = (1./float(ens_size))*np.tile(A.sum(1), (ens_size,1)).transpose()

    # Calculate ensemble observation mean
    MeasAvg = (1./float(ens_size))*np.tile(observation.reshape(meas_size,ens_size).sum(1), (ens_size,1)).transpose()
    
    # Inflate only the simulation ensemble, not the parameter ensemble 
    # A[:(sim_size+1),:] = np.sqrt(alpha)*(A[:(sim_size+1),:] - Amean[:(sim_size+1),:]) + Amean[:(sim_size+1),:]
    A[:,:] = np.sqrt(alpha)*(A[:,:] - Amean[:,:]) + Amean[:,:]

    # Inflate the ensemble observations
    observation = np.sqrt(alpha)*(observation - MeasAvg.transpose()) + MeasAvg.transpose()

    # Calculate ensemble perturbation from mean
    # Apert should be (sim_size+ParSize)x(ens_size)
    dA = A - Amean
    #np.shape(dA)
    # data perturbation from ensemble measurements
    # dD should be (MeasSize)x(ens_size)
    dD = (data - observation).T
    # np.shape(dD)

    # ensemble measurement perturbation from ensemble measurement mean.
    # S is (MeasSize)x(ens_size)
    S = observation.T - MeasAvg

    # Set up measurement covariance matrix
    # COV is (MeasSize)x(MeasSize)
    COV = (1./float(ens_size-1))*np.dot(S,S.transpose()) + data_cov

    # Compute inv(COV)*dD
    # Should be (MeasSize)x(ens_size)
    B = np.linalg.solve(COV,dD)
    # np.shape(B)
    
    # Adjust ensemble perturbations
    # Should be (sim_size+ParSize)x(MeasSize)
    dAS = (1./float(ens_size-1))*np.dot(dA,S.transpose())

    # Compute analysis
    # Analysis is (sim_size+ParSize)x(ens_size)
    Analysis = A + np.dot(dAS,B)

    # Separate and return Analyzed ensemble and Analyzed parameters.
    Analysisparam = Analysis[sim_size:,:].transpose()
    Analysis = Analysis[0:sim_size,:]
    
    print(np.shape(Analysisparam))
            
    return  [A, Amean, dA, dD, MeasAvg, S, COV, B, dAS, Analysis, Analysisparam]


# [augm_state, augm_state_mean, augm_state_pert, data_pert, meas_avg, 
#         obs_pert, COV, inv_data_pert, pert, analysis, analysis_param]