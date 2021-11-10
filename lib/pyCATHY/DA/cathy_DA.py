"""Class managing Data Assimilation process
"""

import os
import matplotlib.pyplot as plt
# from matplotlib import pyplot
import numpy as np
import shutil
# import matplotlib.pyplot as plt 
from pyCATHY.cathy_tools import CATHY



class DA(): #         NO TESTED YET THE INHERITANCE with CATHY MAIN class
    def __init__(self, *args, **kwargs):

        self.var_per_list = {} # dict of dict of perturbated variables parameters

        pass
    
    def perturbate_parm(self, parm, type_parm, mean, sd, per_type, sampling_type = 'lognormal',
                            ensemble_size = 128, show=False, **kwargs):
        """
        Perturbate parameter for the generation of the ensemble
        
        
        Possible variable to perturbate:
            - initial conditions
            - hyetograph
            - van Genuchten retention curves parameters
            - Feddes parameters
    
        Parameters
        ----------
        parm : dict
            DESCRIPTION.
        type_parm : str
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
    

        def Evensen2003():
            print('not yet implemented - see Botto 2018')
            
            pass
        # add Johnson1970 transformation in kwargs 
        
        # check if parameters in part of van Genuchten retention curves
        if type_parm in ['Alpha', 'n', 'thethaR']: #van Genuchten retention curves
            
            print('The parameters of the van Genuchten retention curves α,' + 
                  'n, and θ r are perturbed taking into account their mutual cor-' + 
                  'relation according to Carsel and Parrish (1988)')
            
            # transformed them into normally
            # distributed variables via the Johnson system (Johnson, 1970)
            def Johnson1970():
                print('not yet implemented - see Botto 2018')

                        
        # all parameters except Van Genuchten
        else:
    
            if sampling_type == 'lognormal':
                parm_sampling = np.random.lognormal(mean, sigma=sd, size=ensemble_size)
            elif sampling_type == 'normal':
                parm_sampling = np.random.normal(mean, sd, size=ensemble_size)
    
        
            parm_mat = np.ones(ensemble_size)*parm[type_parm+'_nominal']
            
            if per_type == 'multiplicative':
                parm_per_array = parm_mat*parm_sampling
            elif per_type == 'additive':
                parm_per_array = parm_mat+parm_sampling
        
        
            
            if show == True:
                
                fig = plt.figure(figsize=(6, 3), dpi=150)
                plt.hist(parm_sampling, ensemble_size, alpha=0.5, label='sampling')
                plt.hist(parm_per_array, ensemble_size, alpha=0.5, label='ini_perturbation')
                plt.legend(loc='upper right')
                plt.xlabel(parm[type_parm + '_units'])
                plt.ylabel('Probability')
                plt.title('Histogram of ' + type_parm)
                plt.axvline(x=parm[type_parm+'_nominal'],linestyle='--', color='red')
                plt.show()
                
                if 'savefig' in kwargs:
                    fig.savefig(kwargs['savefig'],dpi=350)
                

    
        # copy initiail variable dict and add 'sampling' and 'ini_perturbation' attributes
        # -------------------------------------------------------------------------
        var_per = {}
        var_per[type_parm] = parm 
        
        key = 'sampling_type'
        var_per[type_parm][key] = sampling_type
    
        key = 'sampling_mean'
        var_per[type_parm][key] = mean
    
        key = 'sampling_sd'
        var_per[type_parm][key] = sd
    
        key = 'sampling'
        var_per[type_parm][key] = parm_sampling
    
        key = 'per_type'
        var_per[type_parm][key] = per_type
    
        key = 'ini_perturbation'
        var_per[type_parm][key] = parm_per_array
    
        self._add_to_perturbated_dict(var_per)
        
        
        
        return self.var_per_list 
    
    
    def _add_to_perturbated_dict(self,var_per_2add):
        '''
        Dict of dict of perturbated variables

        Parameters
        ----------
        var_per_2add : dict
            dict to add to the root dict of perturbated variable.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
               
        self.var_per_list = self.var_per_list | var_per_2add
        
        return self.var_per_list


    def _create_subfolders_ensemble(self,NENS):
        '''
        NO TESTED HERE YET IN INHERITANCE
        SEE CATHY TOOLS for _create_subfolders_ensemble

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
        
    
    # def _prepare_measures(self,measure2prep):
    #     '''
    #     prepare measure before DA

    #     Parameters
    #     ----------
    #     measure2prep : dict
    #         dict containing measure data + metadata.

    #     Returns
    #     -------
    #     TYPE
    #         DESCRIPTION.

    #     '''
               
    #     # define measurement error covariance matrix, R
    #     #---------------------------------------------------------------------
        
    #     return self.vars_per