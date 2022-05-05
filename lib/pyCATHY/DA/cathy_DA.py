"""Class managing Data Assimilation process
"""

import os
import math
import matplotlib.pyplot as plt
# from matplotlib import pyplot
import numpy as np
import scipy.stats as stats

import shutil
# import matplotlib.pyplot as plt 
from pyCATHY.cathy_tools import CATHY

from pyCATHY.DA import enkf, pf



def run_analysis(DA_type,
                 data,data_cov,
                 param,list_update_parm,
                 ensembleX,prediction, 
                 default_state='psi',
                 **kwargs,
                 ):
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
                                                        prediction,
                                                        **kwargs)
                                                        
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
               
    
    
    
    # for dp in list_pert:
        
    #     # need to call perturbate_var as many times as variable to perturbate
    #     # return a dict merging all variable perturbate to parse into prepare_DA
    #     cathyDA = cathy_DA.DA()
    #     parm_per = cathyDA.perturbate_parm(parm=dp, 
    #                                         type_parm = dp['type_parm'], # can also be VAN GENUCHTEN PARAMETERS
    #                                         mean =  dp['mean'],
    #                                         sd =  dp['sd'],
    #                                         sampling_type =  dp['sampling_type'],
    #                                         ensemble_size =  dp['ensemble_size'], # size of the ensemble
    #                                         per_type= dp['per_type'],
    #                                         show= dp['show'],
    #                                         savefig= os.path.join(prj_name,
    #                                                               prj_name + dp['savefig'])
    #                                         )
        
            
    def perturbate_parm(self, parm, type_parm, mean=[], sd=[], per_type=None, 
                        sampling_type = 'lognormal',
                        ensemble_size = 128, 
                        show=False, seed=True,
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
    


        # def sampling_dist_trunc(sampling_type,mean,sd,ensemble_size):
        #     # sampling
        #     np.random.seed(1)
        #     if sampling_type == 'lognormal':
        #         parm_sampling = np.random.lognormal(mean, sigma=sd, size=ensemble_size)
        #     elif sampling_type == 'normal':
        #         # parm_sampling = np.random.normal(mean, sd, size=ensemble_size)
        #         parm_sampling = np.random.normal(mean,scale=sd, size=ensemble_size)
        #     elif sampling_type == 'uniform':
        #         minmax_uni = kwargs['minmax_uni']
        #         parm_sampling = np.random.uniform(minmax_uni[0],minmax_uni[1],ensemble_size)
        #     return parm_sampling
        

        def sampling_dist(sampling_type,mean,sd,ensemble_size):
            # sampling
            np.random.seed(1)
            if sampling_type == 'lognormal':
                parm_sampling = np.random.lognormal(mean, sigma=sd, size=ensemble_size)
            elif sampling_type == 'normal':
                # parm_sampling = np.random.normal(mean, sd, size=ensemble_size)
                parm_sampling = np.random.normal(mean,scale=sd, size=ensemble_size)
            elif sampling_type == 'uniform':
                minmax_uni = kwargs['minmax_uni']
                parm_sampling = np.random.uniform(minmax_uni[0],minmax_uni[1],ensemble_size)
            return parm_sampling
        
        def perturbate_dist(parm,per_type,parm_sampling,ensemble_size):
            # pertubate
            parm_mat = np.ones(ensemble_size)*parm['nominal']
            if per_type == None:
                parm_per_array = parm_sampling
            if per_type == 'multiplicative':
                parm_per_array = parm_mat*parm_sampling
            elif per_type == 'additive':
                parm_per_array = parm_mat+parm_sampling
            return parm_per_array
            
        # Parameter perturbation rules
        #----------------------------------------------------------------------

        def Evensen2003(qk_0, wk,deltaT,Tau):
            '''
            Ensemble Generation of time-variable atmospheric forcing rates.
            Parameters
            ----------
            wk : np.array([])
                is a sequence of white noise drawn from the standard normal distribution.
            deltaT : float
                assimilation interval in sec
            Tau : np.array([])
                the specified time decorrelation length
            time : int
                assimilation time index
                
            Returns
            -------
            qki : Ensemble of time-variable atmospheric forcing rates at time i 

            '''
            if Tau<deltaT:
                raise ValueError('Time decorrelation length is too small; should be at least>=' + str(deltaT))
            gamma = 1 - deltaT/Tau
            qki = gamma * qk_0 + np.sqrt(1-gamma*gamma) * wk
            
            return qki
        
        # check if parameters in part of van Genuchten retention curves
        #----------------------------------------------------------------------
        # if type_parm in ['Alpha', 'nVG', 'thethaR']: #van Genuchten retention curves
        
        def Carsel_Parrish_VGN_pert():

            cholesky_diag_mat = np.diag(3)
            # VGN_means =
            # VGN_parm_per = VGN_means + cholesky_diag_mat.T*z
        
        # add Johnson1970 transformation in kwargs         
        # transformed them into normally
        # distributed variables via the Johnson system (Johnson, 1970)
        def Johnson1970():
            print('not yet implemented - see Botto 2018')

        def sampling_dist_trunc(myclip_a,myclip_b,ensemble_size, **kwargs):
            # https://stackoverflow.com/questions/18441779/how-to-specify-upper-and-lower-limits-when-using-numpy-random-normal
            
            X = stats.truncnorm((myclip_a - kwargs['loc']) / kwargs['scale'], (myclip_b - kwargs['loc']) / kwargs['scale'],
                                loc=kwargs['loc'], scale=kwargs['scale'])
            
            
            # if len(np.where(X.rvs(ensemble_size)<=myclip_a)):
            #     raise ValueError('need to decrease sd')
            # elif len(np.where(X.rvs(ensemble_size)>=myclip_b)):
            #     raise ValueError('need to decrease sd')
                        
            # fig, ax = plt.subplots(2, sharex=True)
            # ax[0].hist(X.rvs(ensemble_size), normed=True)
            
                
            return X.rvs(ensemble_size)




        # perturbation
        # - parameters with constrainsts (truncated distribution)
        # - parameters time dependant
        # - other types of parameters 
        # --------------------------------------------------------------------
        if 'Archie' in type_parm:
            # a : TYPE, optional
            #     Tortuosity factor. The default is [1.0].
            # m : TYPE, optional
            #     Cementation exponent. The default is [2.0]. (usually in the range 1.3 -- 2.5 for sandstones)
            # n : TYPE, optional
            #     Saturation exponent. The default is [2.0].
           
            if 'rFluid' in type_parm:
                parm_sampling = sampling_dist_trunc(myclip_a=0,
                                                    myclip_b=np.inf,
                                                    ensemble_size=ensemble_size,
                                                    loc=mean,
                                                    scale=sd)
            elif 'a' in type_parm:
                parm_sampling = sampling_dist_trunc(myclip_a=0,
                                                    myclip_b=2.5,
                                                    ensemble_size=ensemble_size,
                                                    loc=mean,
                                                    scale=sd)
            elif 'm' in type_parm:
                parm_sampling = sampling_dist_trunc(myclip_a=1,
                                                    myclip_b=2.5,
                                                    ensemble_size=ensemble_size,
                                                    loc=mean,
                                                    scale=sd)
            elif 'n' in type_parm:
                parm_sampling = sampling_dist_trunc(myclip_a=2.5,
                                                    myclip_b=3,
                                                    ensemble_size=ensemble_size,
                                                    loc=mean,
                                                    scale=sd)
            else:
                parm_sampling = sampling_dist(sampling_type,mean,sd,ensemble_size)
                
                
            parm_per_array = perturbate_dist(parm,per_type,parm_sampling,ensemble_size)

            # import matplotlib.pyplot as plt
            # fig, ax = plt.subplots(1, 1)
            # ax.plot(x, truncnorm.pdf(x, a, b),
            #    'r-', lw=5, alpha=0.6, label='truncnorm pdf')
        
        elif 'porosity' in type_parm: 
            
            parm_sampling = sampling_dist_trunc(myclip_a=0,
                                                myclip_b=1,
                                                ensemble_size=ensemble_size,
                                                loc=mean,
                                                scale=sd)
            parm_per_array = perturbate_dist(parm,per_type,parm_sampling,ensemble_size)

            
            
        elif 'VGN' in type_parm: #van Genuchten retention curves
            
            print('The parameters of the van Genuchten retention curves α,' + 
                  'n, and θ r are perturbed taking into account their mutual cor-' + 
                  'relation according to Carsel and Parrish (1988)')
            
            Carsel_Parrish_VGN_pert()

        elif 'atmbc' in type_parm: 
            
            parm_sampling = sampling_dist(sampling_type,mean,sd,ensemble_size)
            parm_per_array = perturbate_dist(parm,per_type,parm_sampling,ensemble_size)

            var_per[type_parm] = parm                
            Tau = parm['time_decorrelation_len']
            wk0 = parm_per_array
            atmbc_times = parm['data2assimilate']['TIME']
            atmbc_values = parm['data2assimilate']['VALUE']
            
            parm_per_array_time_variable = []
            for i, t in enumerate(atmbc_times):
                if i==0:
                    qk_0 = wk0
                    parm_per_array_time_variable.append(qk_0)
                else:
                    qk_0 = parm_per_array_time_variable[i-1]
                    parm['nominal'] = atmbc_values[i]
                    wk = perturbate_dist(parm,per_type,parm_sampling,ensemble_size)
                    
                    deltaT = abs(atmbc_times[i] - atmbc_times[i-1])
                    qk_i = Evensen2003(qk_0, wk, deltaT, Tau)
                    parm_per_array_time_variable.append(qk_i)
            
            key = 'time_variable_perturbation'
            var_per[type_parm][key] = parm_per_array_time_variable
            
        else:
            # other types of parameters 
            #----------------------------------------------------------------------
                
            parm_sampling = sampling_dist(sampling_type,mean,sd,ensemble_size)
            parm_per_array = perturbate_dist(parm,per_type,parm_sampling,ensemble_size)



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
        
        if 'surf_zones_param' in kwargs:
            nb_surf_zones = kwargs['surf_zones_param']
            # parm_per_array = np.tile(parm_per_array,nb_surf_zones)
            var_per[type_parm]['surf_zones_param'] = kwargs['surf_zones_param']
    
        self._add_to_perturbated_dict(var_per)
        

        if show == True:
            
            fig = plt.figure(figsize=(6, 3), dpi=150)
            
            w = 0.2
            nbins = math.ceil((parm_per_array.max() - parm_per_array.min())/w)

            if var_per[type_parm]['transf_type'] is not None:
                if 'log'.casefold() in var_per[type_parm]['transf_type'].casefold():
                    # plt.hist(np.log10(parm_sampling), ensemble_size, alpha=0.5, label='sampling')
                    plt.hist(np.log10(parm_per_array), nbins, alpha=0.5, label='ini_perturbation')
                    plt.axvline(x=np.log10(parm[type_parm+'_nominal']),linestyle='--', color='red')
            else:
                # plt.hist(parm_sampling, ensemble_size/2, alpha=0.5, label='sampling')
                plt.hist(parm_per_array, nbins, alpha=0.5, label='ini_perturbation')
                plt.axvline(x=parm['nominal'],linestyle='--', color='red')

            plt.legend(loc='upper right')
            plt.xlabel(parm['units'])
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