#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 29 12:12:41 2021

@author: ben
"""

import os
import matplotlib.pyplot as plt
import numpy as np
import shutil

from pyCATHY.cathy_tools import CATHY


# -------------------------------------------------------------------#
#%% DATA ASSIMILATION

class DA(CATHY): #         NO TESTED YET THE INHERITANCE with CATHY MAIN class
    def __init__(self, *args, **kwargs):

        self.vars_per = {} # dict of dict of perturbated variables parameters

        pass
    
    def perturbate_var(self, var, parameter, mean, sd, per_type, sampling_type = 'lognormal',
                            ensemble_size = 128):
        """
        Perturbate variable for the generation of the ensemble
    
        Parameters
        ----------
        var : dict
            DESCRIPTION.
        parameter : str
            specify the parameter type.
        mean : TYPE
            DESCRIPTION.
        sd : TYPE
            DESCRIPTION.
        per_type : TYPE
            DESCRIPTION.
        sampling_type : TYPE, optional
            DESCRIPTION. The default is 'lognormal'.
        ensemble_size : TYPE, optional
            DESCRIPTION. The default is 128.
    
        Returns
        -------
        var_per : dict
            Perturbated variable
    
        """
    
    
        var_sampling = np.random.lognormal(mean, sigma=sd, size=ensemble_size)
        
        count, bins, ignored = plt.hist(var_sampling, ensemble_size, density=True, align='mid')
        
        x = np.linspace(min(bins), max(bins), ensemble_size)
        pdf = (np.exp(-(np.log(x) - mean)**2 / (2 * sd**2))
               / (x * sd * np.sqrt(2 * np.pi)))
        
        plt.plot(x, pdf, linewidth=2, color='r')
        plt.axis('tight')
        plt.show()
    
        var_mat = np.ones(ensemble_size)*var[parameter+'_nominal']
        
        if per_type == 'multiplicative':
            var_per_array = var_mat*var_sampling
        elif per_type == 'additive':
            var_per_array = var_mat+var_sampling
    
        
        
        count, bins, ignored = plt.hist(var_per_array, ensemble_size, density=True, align='mid')
        
    
        # copy initiail variable dict and add 'sampling' and 'perturbated' attributes
        # -------------------------------------------------------------------------
        var_per = {}
        var_per[parameter] = var 
        
        key = 'sampling_type'
        var_per[parameter][key] = sampling_type
    
        key = 'sampling_mean'
        var_per[parameter][key] = mean
    
        key = 'sampling_sd'
        var_per[parameter][key] = sd
    
        key = 'sampling'
        var_per[parameter][key] = var_sampling
    
        key = 'per_type'
        var_per[parameter][key] = per_type
    
        key = 'perturbated'
        var_per[parameter][key] = var_per_array
    
        self._add_to_perturbated_dict(var_per)
        
        return self.vars_per 
    
    
    def _add_to_perturbated_dict(self,var_per_2add):
               
        self.vars_per = self.vars_per | var_per_2add
        
        return self.vars_per


    def _create_subfolders_ensemble(self,NENS):
        '''
        NO TESTED HERE YET THE INHERITANCE

        Parameters
        ----------
        NENS : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''

        if not os.path.exists(os.path.join(self.workdir, self.project_name, 
                                           "DA_Ensemble/cathy_origin")):
            os.makedirs(os.path.join(self.workdir, self.project_name, 
                                     "DA_Ensemble/cathy_origin"))
            
            # copy input, output and vtk dir
            for dir2copy in enumerate(['input','output','vtk']):                
                shutil.copytree(os.path.join(self.workdir, self.project_name, dir2copy[1]), 
                                os.path.join(self.workdir, self.project_name, 
                                             "DA_Ensemble/cathy_origin", dir2copy[1])
                )
        
        try:
            # copy exe into cathy_origin folder
            shutil.move(
                os.path.join(self.workdir, self.project_name, self.processor_name),
                os.path.join(self.workdir, self.project_name, 
                             "DA_Ensemble/cathy_origin", self.processor_name)
            )
            # copy cathy.fnames into cathy_origin folder
            shutil.copy(
                os.path.join(self.workdir, self.project_name, 
                             'cathy.fnames'),
                os.path.join(self.workdir, self.project_name, 
                             "DA_Ensemble/cathy_origin/cathy.fnames")
            ) 
            # copy prepro into cathy_origin folder            
            shutil.copytree(os.path.join(self.workdir,  self.project_name, 
                                         'prepro'),
                        os.path.join(self.workdir, self.project_name, 
                                     "DA_Ensemble/cathy_origin/prepro"))
                        
                        
        except:
            self.console.print(":worried_face: [b]processor exe not found[/b]")

            
        path_origin = os.path.join(self.workdir, self.project_name,
                                   "DA_Ensemble/cathy_origin")
        
        # copy origin folder to each ensemble subfolders
        for i in range(NENS):
            path_nudn_i =  os.path.join(self.workdir, self.project_name, 
                                        "DA_Ensemble/cathy_" + str(i+1))  
            
            if not os.path.exists(path_nudn_i):
                shutil.copytree(
                    path_origin, path_nudn_i
                )
        
    
