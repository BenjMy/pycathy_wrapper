#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 12:12:41 2021

@author: ben
"""

import pyvista as pv
import glob
import time 
import os
import matplotlib.pyplot as plt


def RWU(**kwargs):
    """Short summary.

    """

    # import numpy as np

    # PCANA=1.0
    # PCREF=-4.0
    # PCWLT=-150.0
    # PSI=np.arange(-200,10.0,0.1)

        
    # ALFA= np.zeros(len(PSI))       
    # for i in range(len(PSI)):
    #     S1 = PCANA
    #     S2 = max([0.0, (PCANA+1.0E-03)])
    #     SH2O = PSI[i]
    #     GX1 = (SH2O-PCWLT)/(PCREF-PCWLT)
    #     GX1 = min([1.0, max([0.0, GX1])])
    #     GX2 = 1.0 - (SH2O-S1)/(S2-S1)
    #     GX2 = min([1.0, max([0.0, GX2])])
    #     GX = min([GX1, GX2])
    #     ALFA[i] = GX



    return

def dtcoupling(self):
    """processes the file dtcoupling and compares potential and actuale.

    Returns
    -------
    type
        Description of returned object.

    """
    dtcoupling_file = open(os.path.join(self.workdir ,'output' ,'dtcoupling'), 'r')
    Lines = dtcoupling_file.readlines()
    count = len(Lines)
    dtcoupling_file.close()
    
    # Using readline()
    dtcoupling_file = open(os.path.join(workdir ,'output' ,'dtcoupling'), 'r')
    nstep = count-31 # Number of timesteps
    # DT = np.(file,'%g',[22,nstep]); % Read data
    DT = np.loadtxt(dtcoupling_file,skiprows=22,max_rows=22+nstep)
    dtcoupling_file.close()
    
    DT[-1]
    
    fig, axs = plt.subplots(3, 2)
    
    axs[0,0].plot(DT[:,2],DT[:,8],'k:')
    axs[0,0].set(xlabel='Time (h)', ylabel='Pot. & act. atm. fluxes (m^3/h)')
    # axs[0,0].set_xlim([0,max(time)])
    # axs[0,0].set_ylim([1000*(min(min(DT(:,11),DT(:,15)))),0])
    axs[0,0].plot(DT[:,2],DT[:,13],'k-')
    axs[0,0].legend(['Potential','Actual'])
    
    axs[1,0].plot(DT[:,2],DT[:,16],'k:')
    axs[1,0].plot(DT[:,2],DT[:,17],'k--')
    axs[1,0].plot(DT[:,2],DT[:,18],'k-.')
    axs[1,0].plot(DT[:,2],DT[:,19],'k-')
    axs[1,0].set(xlabel='Time (h)', ylabel='Surface saturation fractions')
    axs[1,0].legend(['Horton','Dunne','Ponded','Saturated'])
    
    axs[2,0].plot(DT[:,2],DT[:,6],'k-')
    axs[2,0].plot(DT[:,2],DT[:,7],'k-')
    axs[2,0].set(xlabel='Time (h)', ylabel='Surface routing time steps')
    axs[2,0].legend(['No backstep','Backstep'])
    
    axs[0,1].plot(DT[:,2],DT[:,21]/DT[:,20],'k-')
    # axs[1,0].plot(DT[:,2],DT[:,21],'k-')
    axs[0,1].set(xlabel='Time (h)', ylabel='Surface/subsurface CPU')
    
    axs[1,1].plot(DT[:,2],DT[:,1],'k-')
    # axs[1,0].plot(DT[:,2],DT[:,21],'k-')
    axs[1,1].set(xlabel='Time (h)', ylabel='Subsurface step size (h)')
    
    axs[2,1].plot(DT[:,2],DT[:,4],'k-')
    axs[2,1].plot(DT[:,2],DT[:,5],'k:')
    axs[2,1].legend(['No backstep','Backstep'])
    axs[2,1].set(xlabel='Time (h)', ylabel='Nonlinear iterations')



    return


def COCumflowvol(self):
    """Processes Cumflowvol.

    Returns
    -------
    type
        Description of returned object.

    """
    cumflowvol_file = open(os.path.join(self.workdir ,'output' ,'cumflowvol'), 'r')
    Lines = cumflowvol_file.readlines()
    count = len(Lines)
    cumflowvol_file.close()
    
    nstep = count-3 # Number of timesteps
    
    cumflowvol_file = open(os.path.join(self.workdir ,'output' ,'cumflowvol'), 'r')
    CUMFLOWVOL = np.loadtxt(cumflowvol_file,skiprows=8,max_rows=8+nstep)
    cumflowvol_file.close()
    
    fig, ax = plt.subplots()
    ax.plot(CUMFLOWVOL[:,2],-CUMFLOWVOL[:,7])
    ax.plot(CUMFLOWVOL[:,2]/3600,CUMFLOWVOL[:,7])
    ax.set_title('Cumulative flow volume')
    ax.set(xlabel='Time (h)', ylabel='Net Flow Volume (m^3)')
    ax.legend(['Total flow volume','nansfdir flow volume'])


    return

