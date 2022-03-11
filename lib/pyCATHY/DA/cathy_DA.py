"""Class managing Data Assimilation process
"""

import os
import matplotlib.pyplot as plt
# from matplotlib import pyplot
import numpy as np
import shutil
# import matplotlib.pyplot as plt 
from pyCATHY.cathy_tools import CATHY

from pyCATHY.DA import enkf, pf



def run_analysis(DA_type,data,data_cov,param,ensembleX,prediction, default_state='psi'):
    '''
    Perform the DA analysis step 

    Parameters
    ----------
    DA_type : str
        type of analysis i.e ENKF or Particul filters.
    data : np.array([])
        measured data.
    data_cov : np.array([]) or int
        measured data covariance matrice.
    param : np.array([])
        model parameters to update.
    ensemble_state : np.array([])
        [pressure head, saturation water] at each nodes.
    prediction : np.array([])
        predicted observation (after mapping).
    rejected_ens : list
        list of ensemble to reject
    Returns
    -------
    A : TYPE
        DESCRIPTION.
    Amean : TYPE
        DESCRIPTION.
    dA : TYPE
        DESCRIPTION.
    dD : TYPE
        DESCRIPTION.
    MeasAvg : TYPE
        DESCRIPTION.
    S : TYPE
        DESCRIPTION.
    COV : TYPE
        DESCRIPTION.
    B : TYPE
        DESCRIPTION.
    dAS : TYPE
        DESCRIPTION.
    analysis : TYPE
        DESCRIPTION.
    analysis_param : TYPE
        DESCRIPTION.

    '''
    
    id_state = 0
    if default_state == 'sw':
        id_state = 1
        
    
    if DA_type=='enkf_Evensen2009_Sakov':
        
        [A, Amean, dA, 
         dD, MeasAvg, S, 
         COV, B, dAS, 
         analysis, 
         analysis_param] = enkf.enkf_analysis(data, 
                                              data_cov, 
                                              param, 
                                              ensembleX[id_state], 
                                              prediction
                                              )
                                              
        return [A, Amean, dA, 
         dD, MeasAvg, S, 
         COV, B, dAS, 
         analysis, 
         analysis_param] 
    
    
    if DA_type=='enkf_analysis_inflation':

        [A, Amean, dA, 
         dD, MeasAvg, S, 
         COV, B, dAS, 
         analysis, 
         analysis_param] = enkf.enkf_analysis_inflation(data,
                                                        data_cov,
                                                        param,
                                                        ensembleX[id_state],
                                                        prediction)
                                                        
        return [A, Amean, dA, 
         dD, MeasAvg, S, 
         COV, B, dAS, 
         analysis, 
         analysis_param] 
        
        
    elif DA_type=='pf':
        print('not yet implemented')
        
        [Analysis,AnalysisParam] = pf.pf_analysis(data,
                                                  data_cov,
                                                  param,
                                                  ensembleX[id_state],
                                                  prediction)
        
        
        return Analysis,AnalysisParam

        
    
    
class DA(): #         NO TESTED YET THE INHERITANCE with CATHY MAIN class
    def __init__(self, *args, **kwargs):

        self.var_per_list = {} # dict of dict of perturbated variables parameters

        pass
               
            
    def perturbate_parm(self, parm, type_parm, mean=[], sd=[], per_type=None, 
                        sampling_type = 'lognormal',
                        ensemble_size = 128, 
                        show=False, 
                        **kwargs):
        '''
        Perturbate parameter for the generation of the ensemble
        Possible variable to perturbate:
            - initial conditions
            - hyetograph
            - van Genuchten retention curves parameters
            - Feddes parameters
    

        Parameters
        ----------
        parm : TYPE
            DESCRIPTION.
        type_parm : TYPE
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
        show : TYPE, optional
            DESCRIPTION. The default is False.
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        var_per : dict
            Perturbated variabl

        '''

        var_per = {}


        # copy initiail variable dict and add 'sampling' and 'ini_perturbation' attributes
        # -------------------------------------------------------------------------
        var_per[type_parm] = parm 
        
        key = 'sampling_type'
        var_per[type_parm][key] = sampling_type
    
        key = 'sampling_mean'
        var_per[type_parm][key] = mean
    
        key = 'sampling_sd'
        var_per[type_parm][key] = sd
    
        key = 'per_type'
        var_per[type_parm][key] = per_type
    

        # Parameter perturbation rules
        #----------------------------------------------------------------------

        def Evensen2003(qk_0, wk,deltaT,Tau, time):
            '''
            Ensemble Generation of time-variable atmospheric forcing rates.

            Parameters
            ----------
            wk : np.array([])
                is a sequence of
                white noise drawn from the standard normal distribution.
            deltaT : float
                assimilation interval in sec
            Tau : np.array([])
                the specified time decorrelation length
            time : int
                assimilation time index
                
            Returns
            -------
            Ensemble of time-variable atmospheric forcing rates.

            '''
            
            gamma = 1 - deltaT/Tau
            print('not yet implemented - see Botto 2018')
            qk1 = gamma * qk_0 + np.sqrt(1-gamma*gamma) * wk
            
            return qk1
        
        
        
        # add Johnson1970 transformation in kwargs         
        # transformed them into normally
        # distributed variables via the Johnson system (Johnson, 1970)
        def Johnson1970():
            print('not yet implemented - see Botto 2018')


        # check if parameters in part of van Genuchten retention curves
        #----------------------------------------------------------------------
        if type_parm in ['Alpha', 'n', 'thethaR']: #van Genuchten retention curves
            
            print('The parameters of the van Genuchten retention curves α,' + 
                  'n, and θ r are perturbed taking into account their mutual cor-' + 
                  'relation according to Carsel and Parrish (1988)')
            

                        
        # all parameters except Van Genuchten
        #----------------------------------------------------------------------
        else:
    
            np.random.seed(1)
            if sampling_type == 'lognormal':
                parm_sampling = np.random.lognormal(mean, sigma=sd, size=ensemble_size)
            elif sampling_type == 'normal':
                # parm_sampling = np.random.normal(mean, sd, size=ensemble_size)
                parm_sampling = np.random.normal(mean,scale=sd, size=ensemble_size)
            elif sampling_type == 'uniform':
                minmax_uni = kwargs['minmax_uni']
                parm_sampling = np.random.uniform(minmax_uni[0],minmax_uni[1],ensemble_size)
                
                
            parm_mat = np.ones(ensemble_size)*parm[type_parm+'_nominal']
            
            if per_type == None:
                parm_per_array = parm_sampling
            if per_type == 'multiplicative':
                parm_per_array = parm_mat*parm_sampling
            elif per_type == 'additive':
                parm_per_array = parm_mat+parm_sampling

        key = 'ini_perturbation'
        var_per[type_parm][key] = parm_per_array
        
        
        key = 'sampling'
        var_per[type_parm][key] = parm_sampling
    
        

        # transf_type = '' 
        # Parameter tranformation
        # --------------------------------------------------------------------
        
        if 'transf_type' in kwargs:
            var_per[type_parm]['transf_type'] = kwargs['transf_type']
            if 'transf_bounds' in kwargs:
                var_per[type_parm]['transf_bounds'] = kwargs['transf_bounds']
        else:
            var_per[type_parm]['transf_type'] = None

        # Parameter spatial extension
        # --------------------------------------------------------------------
        if 'surf_param' in kwargs:
            nb_surf_nodes = kwargs['surf_param']
            parm_per_array = np.tile(parm_per_array,nb_surf_nodes)
            var_per[type_parm]['surf_param'] = kwargs['surf_param']
    
        self._add_to_perturbated_dict(var_per)
        

        if show == True:
            
            fig = plt.figure(figsize=(6, 3), dpi=150)
            
            if var_per[type_parm]['transf_type'] is not None:
                if 'log'.casefold() in var_per[type_parm]['transf_type'].casefold():
                    plt.hist(np.log10(parm_sampling), ensemble_size, alpha=0.5, label='sampling')
                    plt.hist(np.log10(parm_per_array), ensemble_size, alpha=0.5, label='ini_perturbation')
                    plt.axvline(x=np.log10(parm[type_parm+'_nominal']),linestyle='--', color='red')
            else:
                # plt.hist(parm_sampling, ensemble_size, alpha=0.5, label='sampling')
                plt.hist(parm_per_array, ensemble_size, alpha=0.5, label='ini_perturbation')
                plt.axvline(x=parm[type_parm+'_nominal'],linestyle='--', color='red')

            plt.legend(loc='upper right')
            plt.xlabel(parm[type_parm + '_units'])
            plt.ylabel('Probability')
            plt.title('Histogram of ' + type_parm)
            

            plt.title('Histogram of ' + type_parm)
            plt.tight_layout()
            
            
            if 'savefig' in kwargs:
                fig.savefig(os.path.join(os.getcwd(),
                                         kwargs['savefig']),
                            dpi=350
                            )
                
        
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
               
        self.var_per_list = self.var_per_list | var_per_2add # only for python 3.9
        #self.var_per_list.update(var_per_2add)
        
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