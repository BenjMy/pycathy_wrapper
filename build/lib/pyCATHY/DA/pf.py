#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 10:19:40 2022

@author: ben
"""
# Credit to:     https://github.com/hickmank/pyda/blob/7a2f04bd752e9c75bc8dcd2a45b21ff549736fa6/pyda/analysis_generator/kf/enkf1.py
import numpy as np
import numpy.random as rn


def weight(self,Data,DataCov,Ensemble,Observation):
    
    # Collect data sizes.
    EnSize = Ensemble.shape[1]

    # First create weight array.
    # W is (1)x(EnSize)
    W = np.zeros(EnSize)

    # Calculate data perturbations from ensemble measurements
    # Dpert = (MeasSize)x(EnSize)
    Dpert = Data - Observation

    # Compute inv(DataCov)*Dpert
    # Should be (MeasSize)x(EnSize)
    B = np.linalg.solve(DataCov,Dpert)

    # Calculate un-normalized weight for each particle using observations
    NormArg = np.diag(np.dot(Dpert.transpose(),B))
    W = np.exp(-(0.5)*(NormArg))

    # Now normalize weights
    W = W/np.sum(W)

    # Weight the ensemble
    self.W = W

# Resampling functions use the Ensemble and Parameter arrays,
# along with weights calulated with one of the Particle filters,
# to generate an analysis Ensemble with equal weights. Ensemble members 
# are resampled according to their weights and a random perturbation from a 
# normal distribution with mean zero and standard deviation 'sigma' is added. 
# This will reduce filter collapse.


def resample(self,Ensemble,Param):
    # Get ensemble size
    EnSize = Ensemble.shape[1]
    
    # Generate resampled indices
    index = range(EnSize)
    resamp_index = np.random.choice(index,size=EnSize,replace=True,p=self.W)

    # Create analysis ensembles
    AnalysisEnsemble = Ensemble[:,resamp_index] + self.sigma*rn.randn(Ensemble.shape[0],EnSize)
    AnalysisParams = Param[resamp_index,:] + self.sigma*rn.randn(EnSize,Param.shape[1])

    return [AnalysisEnsemble,AnalysisParams]

# Returns the analysis ensemble array and the analysis parameter array.
# Analysis = (Ntimestep*SimulationDimension)x(EnSize) numpy array
# AnalysisParam = (Parameter Size + Initialization Size)x(EnSize) numpy array


def pf_analysis(self,Data,DataCov,Param,Ensemble,Observation):
    # Weight the ensemble by data likelihood
    self.weight(Data,DataCov,Ensemble,Observation)
    
    # Resample ensemble by weights
    [Analysis,AnalysisParam] = self.resample(Ensemble,Param)
            
    return [Analysis,AnalysisParam]