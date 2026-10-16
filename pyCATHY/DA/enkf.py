"""Class managing data Assimilation process analysis using enkf
"""
# Credit to:     https://github.com/hickmank/pyda/blob/7a2f04bd752e9c75bc8dcd2a45b21ff549736fa6/pyda/analysis_generator/kf/enkf1.py

import logging
import os
import shutil

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import pyplot
from rich.progress import track

logger = logging.getLogger(__name__)


def enkf_analysis(data, data_cov, param, ensemble, predict_obs, **kwargs):
    """
    Ensemble Kalman Filter (EnKF) analysis step with detailed math,
    physical meaning, and debug prints.

    -----
    Theory
    -----
    EnKF approximates the Kalman analysis step by operating directly
    on an ensemble of states and parameters.

    **Kalman analysis equation**:

        x_a = x_f + K (y - H x_f)

    where:
        - x_f : forecast state (prior)
        - x_a : analysis state (posterior)
        - y   : observations
        - H   : observation operator
        - K   : Kalman gain

    **Kalman gain**:

        K = P^f H^T (H P^f H^T + R)^(-1)

    where:
        - P^f : forecast error covariance
        - R   : measurement error covariance

    In EnKF:
        - P^f is approximated from ensemble statistics.
        - H can be nonlinear → we use ensemble-predicted observations.

    -----
    Physical meaning of key terms
    -----
    - Innovation:   (y - H x_f) = mismatch between observations and model predictions.
    - P^f:          spread of ensemble (uncertainty in forecast).
    - H P^f H^T:    uncertainty of predictions in observation space.
    - Cross-cov:    correlation between errors in states and errors in observations.
                    This term tells the filter *how changes in states affect obs*.
    - Gain K:       balance between trusting model vs. trusting observations.
    - Analysis:     updated states/params nudged toward observations,
                    with correction proportional to cross-correlation.
    """

    # --- Options ---
    Sakov = kwargs.pop('Sakov', False)

    # --- Sizes ---
    ens_size = ensemble.shape[1]   # number of ensemble members
    sim_size = ensemble.shape[0]   # number of state variables
    meas_size = data.shape[0]      # number of observations

    logger.info(
        "EnKF analysis (Sakov=%s): N_ens=%d, N_state=%d, N_obs=%d",
        Sakov, ens_size, sim_size, meas_size,
    )

    # ============================================================
    # Step 1: Ensemble mean
    # ------------------------------------------------------------
    # x̄_f = (1/N) Σ x_f^(i)
    # Physical meaning: the "best guess" of the system before assimilation.
    # ============================================================
    ensemble_mean = np.mean(ensemble, axis=1, keepdims=True)
    ensemble_mean = np.tile(ensemble_mean, (1, ens_size))

    if len(param) > 0:
        param_mean = np.mean(param, axis=1, keepdims=True)
        param_mean = np.tile(param_mean, (1, ens_size))

    if len(param) > 0:
        augm_state_mean = np.vstack([ensemble_mean, param_mean])
        augm_state      = np.vstack([ensemble, param])
    else:
        augm_state_mean = ensemble_mean
        augm_state      = ensemble

    # ============================================================
    # Step 2: Perturbations
    # ------------------------------------------------------------
    # X' = X - x̄
    # Physical meaning: deviations around the mean, encode ensemble spread.
    # Spread = measure of uncertainty.
    # ============================================================
    augm_state_pert = augm_state - augm_state_mean
    logger.debug("state+param perturbations: var=%.4g", np.var(augm_state_pert))

    # ============================================================
    # Step 3: Innovation
    # ------------------------------------------------------------
    # d = y - Hx_f
    # Physical meaning: mismatch between obs and model predictions.
    # Large d → model diverges from reality.
    # ============================================================
    if isinstance(predict_obs, list):
        if len(predict_obs) == 1:
            predict_obs = np.array(predict_obs[0])
        else:
            raise ValueError("predict_obs should be numpy array")

    if data.ndim > 0:
        data_pert = (data.T - predict_obs.T).T
    else:
        data_pert = (data - predict_obs.T).T
    logger.debug(
        "innovation (obs - pred): mean=%.4g, std=%.4g",
        np.mean(data_pert), np.std(data_pert),
    )

    max_residual = 1e2
    if np.max(abs(data_pert)) > max_residual:
        logger.warning(
            "Very large innovations (max |d| = %.3g > %g): check model bias or obs errors",
            np.max(abs(data_pert)), max_residual,
        )

    # ============================================================
    # Step 4: Predicted observation mean and perturbations
    # ------------------------------------------------------------
    # ŷ̄ = (1/N) Σ Hx_f^(i)
    # S = Y - ŷ̄
    # Physical meaning:
    # - obs_avg: best guess in obs space
    # - obs_pert: ensemble spread in obs space (model’s uncertainty in obs)
    # ============================================================
    obs_avg = (1.0 / ens_size) * np.tile(
        predict_obs.reshape(meas_size, ens_size).sum(1), (ens_size, 1)
    ).T
    obs_pert = predict_obs - obs_avg
    logger.debug("predicted obs perturbations: var=%.4g", np.var(obs_pert))

    # (the mean of perturbations around the ensemble mean is 0 by construction:
    # the spread is what tells about a collapse)
    if np.std(obs_pert) <= 1e-12 * max(1.0, float(np.max(np.abs(obs_avg)))):
        logger.warning("Obs perturbations ~0 (low ensemble spread: risk of collapse)")

    # ============================================================
    # Step 5: Observation covariance
    # ------------------------------------------------------------
    # COV = (1/(N-1)) S S^T + R
    # Physical meaning: total uncertainty in obs space =
    #                   model-prediction uncertainty + measurement noise.
    # ============================================================
    if Sakov:
        COV = data_cov.T
        inv_data_pert = data_pert / np.diag(COV)[:, None]
        logger.debug("Sakov simplified covariance form")
    else:
        COV = (1.0 / (ens_size - 1)) * (obs_pert @ obs_pert.T) + data_cov.T
        inv_data_pert = np.linalg.solve(COV, data_pert)
        logger.debug("full covariance matrix COV = SS^T + R")

    if logger.isEnabledFor(logging.DEBUG):  # the condition number is an SVD of COV
        logger.debug("COV shape=%s, cond#=%.2e", COV.shape, np.linalg.cond(COV))

    # ============================================================
    # Step 6: Cross-covariance (state vs obs)
    # ------------------------------------------------------------
    # P_xo = (1/(N-1)) X' S^T
    # Physical meaning: how errors in states are correlated
    # with errors in observations.
    # This drives how much each state is corrected by each observation.
    # ============================================================
    ensemble_pert = (1.0 / (ens_size - 1)) * (augm_state_pert @ obs_pert.T)
    logger.debug("state-observation cross covariance: norm=%.4g", np.linalg.norm(ensemble_pert))
#     ensemble_pert = (1.0 / float(ens_size - 1)) * np.dot(augm_state_pert, 
#                                                 obs_pert.T
#                                                 )

    # ============================================================
    # Step 7: Analysis update
    # ------------------------------------------------------------
    # X_a = X_f + P_xo COV^{-1} d
    # Physical meaning: ensemble is nudged towards observations
    # proportionally to correlations.
    # If strongly correlated → strong update.
    # If weakly correlated → weak update.
    # ============================================================
    analysis = augm_state + (ensemble_pert @ inv_data_pert)
    logger.debug("analysis update applied")

    # ============================================================
    # Step 8: Split results
    # ------------------------------------------------------------
    # Separate updated states and parameters.
    # ============================================================
    analysis_param = analysis[sim_size:, :].T
    analysis       = analysis[0:sim_size, :]


    # covariance diagnostics (result not used here): only when debugging
    if logger.isEnabledFor(logging.DEBUG):
        analyze_covariances(data_cov, obs_pert, augm_state_pert)

    return [
        augm_state,
        augm_state_mean,
        augm_state_pert,
        data_pert,
        obs_avg,
        obs_pert,
        COV,
        inv_data_pert,
        ensemble_pert,
        analysis,
        analysis_param,
    ]
def enkf_analysis_localized_with_inflation(
    data,
    data_cov,
    ensemble,
    param,
    predict_obs,
    L=None,  # localization matrix (state x obs)
    **kwargs
):
    """
    EnKF analysis with an augmented state (states + parameters) and
    localization applied only to states (state x obs).

    Parameters
    ----------
    data : array
        Observations (y)
    data_cov : array
        Observation error covariance (R)
    ensemble : array
        Ensemble of states (X_f)
    param : array
        Ensemble of parameters
    predict_obs : array
        Predicted observations from ensemble
    L : array or None
        Localization matrix (state x obs)
    kwargs : dict
        Options:
            - Sakov : bool (default False)
            - inflate_states : float
            - inflate_params : float
            - jitter_params : float

    Returns
    -------
    list of arrays
        [augm_state,
         augm_state_mean,
         augm_state_pert,
         data_pert,
         obs_avg,
         obs_pert,
         COV,
         inv_data_pert,
         ensemble_pert,
         analysis,
         analysis_param]
    """
    Sakov = kwargs.pop('Sakov')
    inflate_states = kwargs.pop('inflate_states')
    inflate_params = kwargs.pop('inflate_params')
    jitter_params = kwargs.pop('jitter_params')

    ens_size = ensemble.shape[1]
    sim_size = ensemble.shape[0]

    # Step 1: Ensemble mean
    ensemble_mean = np.mean(ensemble, axis=1, keepdims=True)
    param_mean = np.mean(param, axis=1, keepdims=True)
    augm_mean = np.vstack([ensemble_mean, param_mean])
    augm_state = np.vstack([ensemble, param])
    augm_state_pert = augm_state - np.tile(augm_mean, (1, ens_size))
    ensemble_pert = ensemble - np.tile(ensemble_mean, (1, ens_size))

    # Step 2: Innovation
    obs_avg = np.mean(predict_obs, axis=1, keepdims=True)
    obs_pert = predict_obs - obs_avg
    data_pert = data.reshape(-1,1) - predict_obs

    # Step 3: Observation covariance
    if Sakov:
        COV = data_cov
        inv_data_pert = data_pert / np.diag(COV)[:, None]
    else:
        COV = (obs_pert @ obs_pert.T) / (ens_size - 1) + data_cov
        inv_data_pert = np.linalg.solve(COV, data_pert)

    # Step 4: Cross-covariance (augmented)
    P_xo = (augm_state_pert @ obs_pert.T) / (ens_size - 1)

    # Apply localization only to states (top rows)
    if L is not None:
        if L.shape != P_xo[:sim_size, :].shape:
            raise ValueError(f"Localization matrix shape {L.shape} does not match P_xo[:states, :] {P_xo[:sim_size, :].shape}")
        P_xo[:sim_size, :] = P_xo[:sim_size, :] * L

    # Step 5: Analysis update
    analysis_augm = augm_state + P_xo @ inv_data_pert
    analysis = analysis_augm[:sim_size, :]
    analysis_param = analysis_augm[sim_size:, :]

    # Step 6: Inflation
    if inflate_states != 1.0:
        mean_s = np.mean(analysis, axis=1, keepdims=True)
        analysis = mean_s + inflate_states * (analysis - mean_s)

    if inflate_params != 1.0:
        mean_p = np.mean(analysis_param, axis=1, keepdims=True)
        analysis_param = mean_p + inflate_params * (analysis_param - mean_p)

    if jitter_params > 0.0:
        noise = np.random.normal(0, jitter_params, analysis_param.shape)
        analysis_param += noise

    return [
        augm_state,
        augm_mean,
        augm_state_pert,
        data_pert,
        obs_avg,
        obs_pert,
        COV,
        inv_data_pert,
        ensemble_pert,
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

    Sakov = False
    if 'Sakov' in kwargs:
        Sakov = kwargs.pop('Sakov')
    logger.debug('Sakov=%s', Sakov)
    
    
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
    # COV = (1.0 / float(ens_size - 1)) * np.dot(S, S.transpose()) + data_cov

    # Set up observations covariance matrix
    # -------------------------------------------------------------------------
    if Sakov:
        COV = data_cov.transpose()
        B = dD /  np.diag(COV)[:, None]
    else:
        COV = ( (1.0 / float(ens_size - 1)) * np.dot(S, S.transpose()) + data_cov.transpose())
        # Compute inv(COV)*dD
        # -------------------------------------------------------------------------
        # Should be (MeasSize)x(ens_size)
        B = np.linalg.solve(COV, dD)
        
        
        
    # Compute inv(COV)*dD
    # Should be (MeasSize)x(ens_size)
    # B = np.linalg.solve(COV, dD)
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




def enkf_assimilation_with_parm_cov(ensemble, param, data, observation, data_cov, **kwargs):
    Sakov = False
    if 'Sakov' in kwargs:
        Sakov = kwargs.pop('Sakov')
    logger.debug('Sakov=%s', Sakov)

    # Collect data sizes.
    # -------------------------------
    ens_size = ensemble.shape[1]
    sim_size = ensemble.shape[0]
    meas_size = data.shape[0]
    param_size = param.shape[0] if len(param) > 0 else 0

    # First combine the ensemble and param arrays (augmented state matrix)
    # -------------------------------------------------------------------------
    ensemble_mean = np.mean(ensemble, axis=1, keepdims=True)
    ensemble_mean = np.tile(ensemble_mean, (1, ens_size))

    if len(param) > 0:
        # Take the mean line by line for the parameters
        param_mean = np.mean(param, axis=1, keepdims=True)
        param_mean = np.tile(param_mean, (1, ens_size))

    # Combine the state and parameter means to form the augmented state mean
    if len(param) > 0:
        augm_state_mean = np.vstack([ensemble_mean, param_mean])
        augm_state = np.vstack([ensemble, param])
    else:
        augm_state_mean = ensemble_mean
        augm_state = ensemble

    # Calculate ensemble perturbation from the mean
    # -------------------------------------------------------------------------
    augm_state_pert = augm_state - augm_state_mean  # (sim_size+param_size)x(ens_size)

    # Calculate data perturbation from ensemble measurements
    # -------------------------------------------------------------------------
    if isinstance(observation, list):
        if len(observation) == 1:
            observation = np.array(observation[0])
        else:
            logger.warning(
                "observation is a list of %d arrays: should be a numpy array", len(observation)
            )
            
    if data.ndim > 0:
        # Case of dual analysis where data is within an ensemble
        data_pert = (data.T - observation.T).T
    else:
        data_pert = (data - observation.T).T

    if np.max(abs(data_pert)) > 1e3:
        logger.warning(
            "Predictions are too far from observations: data mean/min/max = %.4g/%.4g/%.4g, "
            "predictions mean/min/max = %.4g/%.4g/%.4g",
            np.mean(data), np.min(data), np.max(data),
            np.mean(observation), np.min(observation), np.max(observation),
        )

    # Calculate S = ensemble observation perturbation from ensemble observation mean
    # -------------------------------------------------------------------------
    obs_avg = (1.0 / float(ens_size)) * np.tile(
        observation.reshape(meas_size, ens_size).sum(1), (ens_size, 1)
    ).transpose()
    obs_pert = observation - obs_avg

    if np.std(obs_pert) == 0:
        logger.warning(
            "Ensemble observation perturbations from the ensemble mean are 0: "
            "increase the perturbation / check that the system is not in steady state"
        )

    # Define the parameter covariance matrix P_param (example values)
    # -------------------------------------------------------------------------
    # These values should be based on physical understanding or empirical data.
    sigma_Ks = 1.0e-5  # Std. dev. for hydraulic conductivity (Ks)
    sigma_ZROOT = 0.1  # Std. dev. for root depth (ZROOT)
    cov_Ks_ZROOT = 0.0  # Covariance between Ks and ZROOT (change if needed)

    if param_size > 0:
        P_param = np.array([[sigma_Ks**2, cov_Ks_ZROOT],
                            [cov_Ks_ZROOT, sigma_ZROOT**2]])
    else:
        P_param = np.zeros((0, 0))

    # Set up state covariance matrix P_xx (example values)
    # -------------------------------------------------------------------------
    # This should be replaced with the actual state covariance matrix.
    P_xx = (1.0 / float(ens_size - 1)) * np.dot(augm_state_pert[:sim_size, :], augm_state_pert[:sim_size, :].T)

    # Combine state and parameter covariances into an augmented covariance matrix
    # -------------------------------------------------------------------------
    if param_size > 0:
        P_augm = np.block([
            [P_xx, np.zeros((sim_size, param_size))],     # State to state, no cov. with parameters
            [np.zeros((param_size, sim_size)), P_param]   # Param to param covariances
        ])
    else:
        P_augm = P_xx

    # Set up observations covariance matrix
    # -------------------------------------------------------------------------
    if Sakov:
        COV = data_cov.transpose()
    else:
        COV = ( (1.0 / float(ens_size - 1)) *
                np.dot(obs_pert, obs_pert.transpose()) 
               + data_cov.transpose()
              )

    # Compute inv(COV) * dD
    # -------------------------------------------------------------------------
    inv_data_pert = np.linalg.solve(COV, data_pert)

    # Adjust ensemble perturbations
    # -------------------------------------------------------------------------
    pert = (1.0 / float(ens_size - 1)) * np.dot(augm_state_pert, obs_pert.T)

    # Compute analysis
    # -------------------------------------------------------------------------
    analysis = augm_state + np.dot(pert, inv_data_pert)

    # Separate and return analyzed ensemble and analyzed parameters
    # -------------------------------------------------------------------------
    if param_size > 0:
        analysis_param = analysis[sim_size:, :].transpose()
        analysis_state = analysis[0:sim_size, :]
    else:
        analysis_param = None
        analysis_state = analysis

    # return analysis_state, analysis_param
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


def enkf_analysis_inflation_with_parm_cov(data, data_cov, param, ensemble, observation, param_cov=None, **kwargs):
    alpha = 1  # Covariance inflation parameter. Usually alpha >= 1.

    if "alpha" in kwargs:
        alpha = kwargs["alpha"]

    # Collect data sizes.
    ens_size = ensemble.shape[1]
    sim_size = ensemble.shape[0]
    meas_size = data.shape[0]

    observation = observation.T

    # Combine the ensemble and parameter arrays.
    # A is (sim_size + param_size) x ens_size
    A = np.vstack([ensemble, param])

    # Calculate ensemble mean
    Amean = (1.0 / float(ens_size)) * np.tile(A.sum(1), (ens_size, 1)).transpose()

    # Calculate ensemble observation mean
    MeasAvg = (1.0 / float(ens_size)) * np.tile(
        observation.reshape(meas_size, ens_size).sum(1), (ens_size, 1)
    ).transpose()

    # Inflate the entire ensemble (state + parameters)
    A = np.sqrt(alpha) * (A - Amean) + Amean

    # Inflate the ensemble observations
    observation = (
        np.sqrt(alpha) * (observation - MeasAvg.transpose()) + MeasAvg.transpose()
    )

    # Calculate ensemble perturbation from mean
    dA = A - Amean

    # Data perturbation from ensemble measurements
    dD = (data - observation).T

    # Ensemble measurement perturbation from ensemble measurement mean.
    S = observation.T - MeasAvg

    # Set up measurement covariance matrix
    # COV is (MeasSize)x(MeasSize)
    COV = (1.0 / float(ens_size - 1)) * np.dot(S, S.transpose()) + data_cov

    # If parameter covariance is provided, create a combined covariance matrix for the parameters
    if param_cov is not None:
        # Check if param_cov is a square matrix and matches the size of parameters
        if param_cov.shape != (param.shape[0], param.shape[0]):
            raise ValueError("param_cov must be a square matrix of size (param_size, param_size).")

        # Combine parameter covariance with measurement covariance
        # COV_param is the parameter covariance matrix
        COV_param = param_cov
        
        # Combine the measurement and parameter covariance matrices
        # Assuming you want to create a block diagonal covariance matrix
        COV_combined = np.block([[COV, np.zeros((COV.shape[0], COV_param.shape[1]))],
                                  [np.zeros((COV_param.shape[0], COV.shape[1])), COV_param]])
    else:
        COV_combined = COV  # If no parameter covariance is provided, use only measurement covariance

    # Compute inv(COV_combined)*dD
    B = np.linalg.solve(COV_combined, dD)

    # Adjust ensemble perturbations
    # dAS should be (sim_size + param_size)x(MeasSize)
    dAS = (1.0 / float(ens_size - 1)) * np.dot(dA, S.transpose())

    # Compute analysis
    # Analysis is (sim_size + param_size)x(ens_size)
    Analysis = A + np.dot(dAS, B)

    # Separate and return Analyzed ensemble and Analyzed parameters.
    Analysisparam = Analysis[sim_size:, :].transpose()
    Analysis = Analysis[0:sim_size, :]

    logger.debug('Analysisparam shape: %s', np.shape(Analysisparam))

    return [A, Amean, dA, dD, MeasAvg, S, COV_combined, B, dAS, Analysis, Analysisparam]

# Example usage
# data = np.random.rand(meas_size)  # Replace with actual data
# data_cov = np.random.rand(meas_size, meas_size)  # Replace with actual data covariance
# param = np.random.rand(param_size, ens_size)  # Replace with actual parameter ensemble
# ensemble = np.random.rand(sim_size, ens_size)  # Replace with actual state ensemble
# observation = np.random.rand(meas_size, ens_size)  # Replace with actual observations
# param_cov = np.random.rand(param_size, param_size)  # Define your own parameter covariance
# results = enkf_analysis_inflation(data, data_cov, param, ensemble, observation, param_cov=param_cov)


def _cov_diag_offdiag_means(pert):
    """Mean of the diagonal and of (C - diag(C)) for C = P P^T / (N-1), P = ``pert``.

    Equal to ``np.mean(np.diag(C))`` and ``np.mean(C - np.diag(np.diag(C)))`` but
    computed from P without building C (O(n N) memory instead of O(n^2)):
    sum(C) = ||P^T 1||^2 / (N-1) and diag(C)_i = ||P_i||^2 / (N-1).
    """
    pert = np.asarray(pert, dtype=float)
    n, n_ens = pert.shape
    sum_diag = np.sum(pert**2) / (n_ens - 1)
    sum_all = np.sum(pert.sum(axis=0) ** 2) / (n_ens - 1)
    return sum_diag / n, (sum_all - sum_diag) / n**2


def analyze_covariances(data_cov, obs_pert, augm_state_pert):
    """
    Analyze the structure and meaning of covariances used in EnKF.

    Parameters
    ----------
    data_cov : np.ndarray, shape (N_obs, N_obs)
        Measurement error covariance matrix R.
    obs_pert : np.ndarray, shape (N_obs, N_ens)
        Perturbations of predicted observations.
    augm_state_pert : np.ndarray, shape (N_state+N_param, N_ens)
        Perturbations of augmented state.

    Returns
    -------
    summary : dict
        Dictionary summarizing covariance types and diagnostics.
    """

    n_ens = obs_pert.shape[1]
    diag_R = np.allclose(data_cov, np.diag(np.diag(data_cov)))
    logger.debug(
        "R: shape=%s, %s", data_cov.shape,
        "diagonal (independent measurement errors)" if diag_R
        else "full (correlated measurement errors)",
    )

    # means of the diagonal / off-diagonal of P P^T/(N-1), without forming the
    # (n x n) matrix (n_state can be ~4e4, i.e. ~14 GB)
    obs_diag, obs_off = _cov_diag_offdiag_means(obs_pert)
    logger.debug(
        "HPH^T: shape=(%d, %d), diagonal mean=%.4g, off-diagonal mean=%.4g",
        obs_pert.shape[0], obs_pert.shape[0], obs_diag, obs_off,
    )
    st_diag, st_off = _cov_diag_offdiag_means(augm_state_pert)
    logger.debug(
        "P^f: shape=(%d, %d), diagonal mean=%.4g, off-diagonal mean=%.4g",
        augm_state_pert.shape[0], augm_state_pert.shape[0], st_diag, st_off,
    )
    cross_cov_shape = (augm_state_pert.shape[0], obs_pert.shape[0])
    logger.debug("P_xo: shape=%s", cross_cov_shape)

    # Return a structured summary for later inspection
    return {
        "data_cov_diag": diag_R,
        "obs_cov_diag_mean": obs_diag,
        "obs_cov_offdiag_mean": obs_off,
        "state_cov_diag_mean": st_diag,
        "state_cov_offdiag_mean": st_off,
        "cross_cov_shape": cross_cov_shape,
    }
