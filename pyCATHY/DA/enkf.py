"""Class managing data Assimilation process analysis using enkf
"""
# Credit to:     https://github.com/hickmank/pyda/blob/7a2f04bd752e9c75bc8dcd2a45b21ff549736fa6/pyda/analysis_generator/kf/enkf1.py

import os
import shutil

import matplotlib.pyplot as plt
import numpy as np
import rich.console
from matplotlib import pyplot
from rich import print
from rich.progress import track


def enkf_analysis(data, data_cov, param, ensemble, observation, **kwargs):
    """
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
    """

    # get kwargs optionnal arguments
    # --------------------------------
    # localize = False
    # if "localize" in kwargs:
    #     localize = True

    # print(kwargs)

    Sakov = False
    if 'Sakov' in kwargs:
        Sakov = kwargs.pop('Sakov')
    print('Sakov' + str(Sakov))
    
    # Collect data sizes.
    # - -------------------------------
    ens_size = ensemble.shape[1]
    sim_size = ensemble.shape[0]
    meas_size = data.shape[0]

    # First combine the ensemble and param arrays (augmented state matrice)
    # -------------------------------------------------------------------------
    ensemble_mean = np.mean(ensemble, axis=1, keepdims=True)
    ensemble_mean = np.tile(ensemble_mean, (1, ens_size))

    
    if len(param) > 0:
        # Take the mean line by line i.e. each line is an independent parameter to update
        param_mean = np.mean(param, axis=1, keepdims=True)
        param_mean = np.tile(param_mean, (1, ens_size))

    # Combine the state and parameter means to form the augmented state mean
    if len(param) > 0:
        augm_state_mean = np.vstack([ensemble_mean, param_mean])
        augm_state = np.vstack([ensemble, param])
    else:
        augm_state_mean = ensemble_mean
        augm_state = ensemble
        
    # Calculate ensemble perturbation from mean
    # -------------------------------------------------------------------------
    # augm_state_pert should be (sim_size+ParSize)x(ens_size)
    augm_state_pert = augm_state - augm_state_mean
    # augm_state_pert[sim_size:, :][0]

    
    # Calculate data perturbation from ensemble measurements
    # -------------------------------------------------------------------------
    # data_pert should be (MeasSize)x(ens_size)
    if isinstance(observation, list):
        if len(list) == 1:
            observation = np.array(observation[0])
        else:
            print(observation)
            print("observation is a list? should be a numpy array")
            
    if  data.ndim>0:
        # print(f"Number of dimensions >0:
        # Case of dual analysis where data are within an ensemble
        # OR case where data are perturbated (add noise and create an ensemble)
        data_pert = (data.T - observation.T).T
    else:
        data_pert = (data - observation.T).T


    # print(np.shape(data))
    # print(np.shape(observation.T))
    # data_pert = (data - observation).T

    if np.max(abs(data_pert)) > 1e3:
        print(f'data mean: {np.mean(data)}, min: {np.min(data)}, max: {np.max(data)}')
        print(f'data obs: {np.mean(observation)}, min: {np.min(observation)}, max: {np.max(observation)}')
        print('!predictions are too far from observations!')

    # Calculate S = ensemble observations perturbation from ensemble observation mean.
    # -------------------------------------------------------------------------
    # S is (MeasSize)x(ens_size)
    obs_avg = (1.0 / float(ens_size)) * np.tile(
        observation.reshape(meas_size, ens_size).sum(1), (ens_size, 1)
    ).transpose()
    
    # np.max(obs_avg)
    # np.min(obs_avg)
    
    obs_pert = observation - obs_avg
    # np.max(obs_pert)
    # np.min(obs_pert)
    
    if obs_pert.mean() == 0:
        print(
            "Ensemble observations perturbation from ensemble"
            " measurement mean is too small (=0):"
            "- Increase perturbation!"
            "- Check if the system is not in steady state"
        )

    # Set up observations covariance matrix
    # -------------------------------------------------------------------------
    if Sakov:
        COV = data_cov.transpose()
    else:
        COV = ( (1.0 / float(ens_size - 1)) * np.dot(obs_pert, obs_pert.transpose()) + data_cov.transpose())
    # Compute inv(COV)*dD
    # -------------------------------------------------------------------------
    # Should be (MeasSize)x(ens_size)
    inv_data_pert = np.linalg.solve(COV, data_pert)
    print("Compute inv(COV)*dD")
    print(np.shape(COV))


    # Adjust ensemble perturbations
    # -------------------------------------------------------------------------
    # Should be (sim_size+ParSize)x(MeasSize)
    pert = (1.0 / float(ens_size - 1)) * np.dot(augm_state_pert, 
                                                obs_pert.T
                                                )

    # Compute analysis
    # -------------------------------------------------------------------------
    # Analysis is (sim_size+ParSize)x(ens_size)
    analysis = augm_state + np.dot(pert, inv_data_pert)
    print("Analysis ...")

    # Separate and return Analyzed ensemble and Analyzed parameters.
    # -------------------------------------------------------------------------
    analysis_param = analysis[sim_size:, :].transpose()
    analysis = analysis[0:sim_size, :]

    return [
        augm_state,
        augm_state_mean,
        augm_state_pert,
        data_pert,
        obs_avg,
        obs_pert,
        COV,
        inv_data_pert,
        pert,
        analysis,
        analysis_param,
    ]


def enkf_dual_analysis(data, data_cov, param, ensemble, observation, **kwargs):
    
    # 1st step update state but not parameters
    # -----------------------------------------
    first_step_analysis = enkf_analysis(data, data_cov, param, ensemble, observation, **kwargs)
    cov_data_first_step = first_step_analysis[6]   
    
    # 2nd step parameters, are updated by assimilating new states as observations
    # -----------------------------------------
    cov_data_second_step = np.diag([np.mean(np.diag(cov_data_first_step))] * len(first_step_analysis[9]))   
    observation_2nd_step = ensemble
    return  enkf_analysis(first_step_analysis[9], 
                          cov_data_second_step, 
                          param, ensemble, 
                          observation_2nd_step, 
                          **kwargs
                          )
    
# def enkf_analysis_Sakov(data, data_cov, param, ensemble, observation,Sakov=True):
    
#     [
#        augm_state,
#        augm_state_mean,
#        augm_state_pert,
#        data_pert,
#        obs_avg,
#        obs_pert,
#        COV,
#        inv_data_pert,
#        pert,
#        analysis,
#        analysis_param,
#    ] = enkf_analysis(data, data_cov, param, ensemble, observation, **kwargs)
    
#     return [
#         augm_state,
#         augm_state_mean,
#         augm_state_pert,
#         data_pert,
#         obs_avg,
#         obs_pert,
#         COV,
#         inv_data_pert,
#         pert,
#         analysis,
#         analysis_param,
#     ]


def enkf_analysis_inflation(data, data_cov, param, ensemble, observation, **kwargs):

    alpha = 1  # <alpha> = (scalar) covariance inflation parameter. Usually alpha >= 1.

    if "alpha" in kwargs:
        alpha = kwargs["alpha"]

    # Collect data sizes.
    ens_size = ensemble.shape[1]
    sim_size = ensemble.shape[0]

    observation = observation.T
    # data = np.array([data]).T
    meas_size = data.shape[0]

    # First combine the ensemble and param arrays.
    # A is (sim_size+ParSize)x(ens_size)
    A = np.vstack([ensemble, param])

    # Calculate ensemble mean
    Amean = (1.0 / float(ens_size)) * np.tile(A.sum(1), (ens_size, 1)).transpose()

    # Calculate ensemble observation mean
    MeasAvg = (1.0 / float(ens_size)) * np.tile(
        observation.reshape(meas_size, ens_size).sum(1), (ens_size, 1)
    ).transpose()

    # Inflate only the simulation ensemble, not the parameter ensemble
    # A[:(sim_size+1),:] = np.sqrt(alpha)*(A[:(sim_size+1),:] - Amean[:(sim_size+1),:]) + Amean[:(sim_size+1),:]
    A[:, :] = np.sqrt(alpha) * (A[:, :] - Amean[:, :]) + Amean[:, :]

    # Inflate the ensemble observations
    observation = (
        np.sqrt(alpha) * (observation - MeasAvg.transpose()) + MeasAvg.transpose()
    )

    # Calculate ensemble perturbation from mean
    # Apert should be (sim_size+ParSize)x(ens_size)
    dA = A - Amean
    # np.shape(dA)
    # data perturbation from ensemble measurements
    # dD should be (MeasSize)x(ens_size)
    dD = (data - observation).T
    # np.shape(dD)

    # ensemble measurement perturbation from ensemble measurement mean.
    # S is (MeasSize)x(ens_size)
    S = observation.T - MeasAvg

    # Set up measurement covariance matrix
    # COV is (MeasSize)x(MeasSize)
    COV = (1.0 / float(ens_size - 1)) * np.dot(S, S.transpose()) + data_cov

    # Compute inv(COV)*dD
    # Should be (MeasSize)x(ens_size)
    B = np.linalg.solve(COV, dD)
    # np.shape(B)

    # Adjust ensemble perturbations
    # Should be (sim_size+ParSize)x(MeasSize)
    dAS = (1.0 / float(ens_size - 1)) * np.dot(dA, S.transpose())

    # Compute analysis
    # Analysis is (sim_size+ParSize)x(ens_size)
    Analysis = A + np.dot(dAS, B)

    # Separate and return Analyzed ensemble and Analyzed parameters.
    Analysisparam = Analysis[sim_size:, :].transpose()
    Analysis = Analysis[0:sim_size, :]

    # print(np.shape(Analysisparam))

    return [A, Amean, dA, dD, MeasAvg, S, COV, B, dAS, Analysis, Analysisparam]



def enkf_analysis_inflation_multiparm(data, data_cov, param, ensemble, observation, **kwargs):

    alpha = 1  # <alpha> = (scalar) covariance inflation parameter. Usually alpha >= 1.

    if "alpha" in kwargs:
        alpha = kwargs["alpha"]

    # Collect data sizes.
    ens_size = ensemble.shape[1]
    sim_size = ensemble.shape[0]

    observation = observation.T
    # data = np.array([data]).T
    meas_size = data.shape[0]

    # calculate parameter mean
    if len(param) > 0:
        # Take the mean line by line i.e. each line is an independant parameter to update
        param_mean = []
        for l in range(len(param)):
            param_mean.append(
                (1.0 / float(ens_size))
                * np.tile(param[l, :].sum(0), (ens_size, 1)).transpose()
            )
        param_mean = np.vstack(param_mean)

    else:
        param_mean = param        
        
    # First combine the ensemble and param arrays.
    # A is (sim_size+ParSize)x(ens_size)
    
    A = np.vstack([ensemble, param])

    # Calculate ensemble mean
    ensemble_mean = (1.0 / float(ens_size)) * np.tile(ensemble.sum(1), 
                                                      (ens_size, 1)).transpose()
    
    # Logic ensemble mean
    Amean =  np.vstack([ensemble_mean, param_mean])

    # Calculate ensemble observation mean
    MeasAvg = (1.0 / float(ens_size)) * np.tile(
        observation.reshape(meas_size, ens_size).sum(1), (ens_size, 1)
    ).transpose()

    # Inflate only the simulation ensemble, not the parameter ensemble
    # A[:(sim_size+1),:] = np.sqrt(alpha)*(A[:(sim_size+1),:] - Amean[:(sim_size+1),:]) + Amean[:(sim_size+1),:]
    A[:, :] = np.sqrt(alpha) * (A[:, :] - Amean[:, :]) + Amean[:, :]

    # Inflate the ensemble observations
    observation = (
        np.sqrt(alpha) * (observation - MeasAvg.transpose()) + MeasAvg.transpose()
    )

    # Calculate ensemble perturbation from mean
    # Apert should be (sim_size+ParSize)x(ens_size)
    dA = A - Amean
    # np.shape(dA)
    # data perturbation from ensemble measurements
    # dD should be (MeasSize)x(ens_size)
    dD = (data - observation).T
    # np.shape(dD)

    # ensemble measurement perturbation from ensemble measurement mean.
    # S is (MeasSize)x(ens_size)
    S = observation.T - MeasAvg

    # Set up measurement covariance matrix
    # COV is (MeasSize)x(MeasSize)
    COV = (1.0 / float(ens_size - 1)) * np.dot(S, S.transpose()) + data_cov

    # Compute inv(COV)*dD
    # Should be (MeasSize)x(ens_size)
    B = np.linalg.solve(COV, dD)
    # np.shape(B)

    # Adjust ensemble perturbations
    # Should be (sim_size+ParSize)x(MeasSize)
    dAS = (1.0 / float(ens_size - 1)) * np.dot(dA, S.transpose())

    # Compute analysis
    # Analysis is (sim_size+ParSize)x(ens_size)
    Analysis = A + np.dot(dAS, B)

    # Separate and return Analyzed ensemble and Analyzed parameters.
    Analysisparam = Analysis[sim_size:, :].transpose()
    Analysis = Analysis[0:sim_size, :]


    return [A, Amean, dA, dD, MeasAvg, S, COV, B, dAS, Analysis, Analysisparam]


