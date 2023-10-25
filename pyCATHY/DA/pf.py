#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 10:19:40 2022

@author: ben
"""
# Credit to:     https://github.com/hickmank/pyda/blob/7a2f04bd752e9c75bc8dcd2a45b21ff549736fa6/pyda/analysis_generator/kf/enkf1.py
import numpy as np
import numpy.random as rn


def weight(Data, DataCov, Ensemble, Observation):

    # Collect data sizes.
    EnSize = Ensemble.shape[1]

    # First create weight array.
    # W is (1)x(EnSize)
    W = np.zeros(EnSize)

    # Calculate data perturbations from ensemble measurements
    # Dpert = (MeasSize)x(EnSize)
    Dpert = (Data - Observation.T).T
    # np.shape(Dpert)

    # import matplotlib.pyplot as plt

    # fig, ax = plt.subplots()
    # cax = ax.matshow(Dpert, aspect="auto")
    # cbar = fig.colorbar(cax, location="bottom")

    # Compute inv(DataCov)*Dpert
    # Should be (MeasSize)x(EnSize)
    B = np.linalg.solve(1 / DataCov, Dpert)

    # import matplotlib.pyplot as plt

    # fig, ax = plt.subplots()
    # cax = ax.matshow(B, aspect="auto")
    # cbar = fig.colorbar(cax, location="bottom")

    # Calculate un-normalized weight for each particle using observations
    NormArg = np.diag(np.dot(Dpert.T, B))
    W = np.exp(-(0.5) * (NormArg))

    # Now normalize weights
    W = W / np.sum(W)

    # Weight the ensemble
    # self.W = W

    return W


# Resampling functions use the Ensemble and Parameter arrays,
# along with weights calulated with one of the Particle filters,
# to generate an analysis Ensemble with equal weights. Ensemble members
# are resampled according to their weights and a random perturbation from a
# normal distribution with mean zero and standard deviation 'sigma' is added.
# This will reduce filter collapse.


def resample(Ensemble, Param, W, sigma):
    # Get ensemble size
    EnSize = Ensemble.shape[1]

    # Generate resampled indices
    index = range(EnSize)

    if np.all(np.isnan(W)) == False:
        resamp_index = np.random.choice(index, size=EnSize, replace=True, p=W)
    else:
        resamp_index = index

    # Create analysis ensembles
    AnalysisEnsemble = Ensemble[:, resamp_index] + sigma * rn.randn(
        Ensemble.shape[0], EnSize
    )
    AnalysisParams = Param[resamp_index] + sigma * rn.randn(EnSize, 1).T

    return [AnalysisEnsemble, AnalysisParams]


# Returns the analysis ensemble array and the analysis parameter array.
# Analysis = (Ntimestep*SimulationDimension)x(EnSize) numpy array
# AnalysisParam = (Parameter Size + Initialization Size)x(EnSize) numpy array


def pf_analysis(Data, DataCov, Param, Ensemble, Observation, sigma=1):
    # Weight the ensemble by data likelihood

    W = weight(Data, DataCov, Ensemble, Observation)

    # Resample ensemble by weights
    [Analysis, AnalysisParam] = resample(Ensemble, Param, W, sigma)

    return [Analysis, AnalysisParam]
