"""Class managing Data Assimilation process
"""

import os
import math
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import pandas as pd
from collections import OrderedDict

import shutil
from pyCATHY.cathy_tools import CATHY
from pyCATHY.DA import enkf, pf
from pyCATHY.plotters import cathy_plots as plt_CT

import re
from pyCATHY.importers import cathy_outputs as out_CT
from pyCATHY.importers import cathy_inputs as in_CT
from pyCATHY.importers import sensors_measures as in_meas
from pyCATHY.ERT import petro_Archie as Archie
from pyCATHY import cathy_utils as utils_CT


    
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





# utils observations
# ------------------

def dictObs_2pd(dict_obs):
    '''dict of observation to dataframe of observation'''
    df_obs = pd.DataFrame.from_dict(dict_obs).stack().to_frame()
    df_obs = pd.DataFrame(df_obs[0].values.T.tolist(),
                            index=df_obs.index)
    df_obs.index.names= ['sensorNameidx','assimilation time']
    return df_obs

def resynchronise_times(data_measure,atmbc_df):
    ''' old key is elapsed time in second from the first observation,
    while new key is from the first atmbc time
    '''
    data_measure_sync = dict(data_measure)
    try:
        for d in range(len(data_measure_sync.keys())):
            # print(d)
            items_dict = list(data_measure_sync.items())
            # print(items_dict)
            # list(items_dict[d][1].keys())
            for sensor in list(items_dict[d][1].keys()):
                new_key = (atmbc_df[atmbc_df['sensor_name']==sensor]['diff'].dt.total_seconds().to_numpy()[d])
                old_key = list(data_measure_sync.keys())[d]
                data_measure_sync[old_key][sensor]['assimilation_times'] = new_key
                data_measure_sync[new_key] = data_measure_sync.pop(old_key)
    except:
        print('datetime first atmbc point')
        print(atmbc_df['datetime'][0])
        print('datetime first measurement')
        print(data_measure[0]['tensiometer']['datetime'])
        print('cant synchronise times - continue without') 
    return data_measure_sync



# SAMPLING distribution
# ----------------------

def sampling_dist_trunc(myclip_a,myclip_b,ensemble_size, **kwargs):
    # https://stackoverflow.com/questions/18441779/how-to-specify-upper-and-lower-limits-when-using-numpy-random-normal
    X = stats.truncnorm((myclip_a - kwargs['loc']) / kwargs['scale'], (myclip_b - kwargs['loc']) / kwargs['scale'],
                        loc=kwargs['loc'], scale=kwargs['scale'])
    return X.rvs(ensemble_size)

def sampling_dist(sampling_type,mean,sd,ensemble_size,**kwargs):
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


# PERTUBATE distribution
# ----------------------

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


def Carsel_Parrish_VGN_pert():
    cholesky_diag_mat = np.diag(3)
    pass

def Archie_pert_rules(parm,type_parm,ensemble_size,mean,sd,per_type,sampling_type):
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
                                            scale=sd
                                            )
    elif 'm' in type_parm:
        parm_sampling = sampling_dist_trunc(myclip_a=1.3,
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
    return parm_per_array


def VG_pert_rules(var_per_2add,parm,type_parm,ensemble_size,mean,sd,per_type,sampling_type,**kwargs):
    print('The parameters of the van Genuchten retention curves α,' +
          'n, and θ r are perturbed taking into account their mutual cor-' +
          'relation according to Carsel and Parrish (1988)')

    if 'Carsel_Parrish_VGN_pert' in kwargs:
        utils.Carsel_Parrish_1988(soilTexture=None)
        Carsel_Parrish_VGN_pert()
    else:
        if 'clip_min' in var_per_2add[type_parm].keys():
            parm_sampling = sampling_dist_trunc(myclip_a=var_per_2add[type_parm]['clip_min'],
                                                     myclip_b=var_per_2add[type_parm]['clip_max'],
                                                     ensemble_size=ensemble_size,
                                                     loc=mean,
                                                     scale=sd
                                                     )
        else:
            parm_sampling = sampling_dist(sampling_type,mean,sd,ensemble_size)
        parm_per_array = perturbate_dist(parm,per_type,parm_sampling,ensemble_size)
    return parm_per_array


def atmbc_pert_rules(var_per_2add,parm,type_parm,ensemble_size,mean,sd,per_type,sampling_type,**kwargs):
    
    parm_sampling = sampling_dist(sampling_type,mean,sd,ensemble_size)
    parm_per_array = perturbate_dist(parm,per_type,parm_sampling,ensemble_size)
    var_per_2add[type_parm] = parm
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
    var_per_2add[type_parm][key] = parm_per_array_time_variable
    return var_per_2add



def Johnson1970(self):
    print('not yet implemented - see Botto 2018')
        
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


def build_dict_attributes_pert(var_per_2add,type_parm,parm_per_array,parm_sampling,**kwargs):
    
    key = 'ini_perturbation'
    var_per_2add[type_parm][key] = parm_per_array
    key = 'sampling'
    var_per_2add[type_parm][key] = parm_sampling
    # Parameter tranformation
    # --------------------------------------------------------------------
    if 'transf_type' in kwargs:
        var_per_2add[type_parm]['transf_type'] = kwargs['transf_type']
        if 'transf_bounds' in kwargs:
            var_per_2add[type_parm]['transf_bounds'] = kwargs['transf_bounds']
    else:
        var_per_2add[type_parm]['transf_type'] = None
    # Parameter spatial extension
    # --------------------------------------------------------------------
    if 'surf_zones_param' in kwargs:
        nb_surf_zones = kwargs['surf_zones_param']
        # parm_per_array = np.tile(parm_per_array,nb_surf_zones)
        var_per_2add[type_parm]['surf_zones_param'] = kwargs['surf_zones_param']
    return var_per_2add



#%%  
# DA class
# ---------

class DA(): #         NO TESTED YET THE INHERITANCE with CATHY MAIN class
    # def __init__(self):
        #     self.var_per_dict = {} # dict of dict of perturbated variables parameters
        #     print(CATHY)
        #     C = CATHY()
        #     C.__init__()
        #     C.workdir
        #     print(self.workdir)
        #     pass
    
        #self.workdir
        #self.project_name
        #self.processor_name
        #self.console
        # self.stacked_data_cov = [] # merged data covariance matrice of all the observation data and all times
        # self.dict_obs = OrderedDict() # dictionnary containing all the observation data
        # self.var_per_dict = {} # list of variable perturbated
        # pass

    # -------------------------------------------------------------------#
    # %% DATA ASSIMILATION FCTS

    def _DA_init(self, NENS=[], ENS_times=[], parm_pert=[],update_parm_list='all'):
        """
        Initial preparation for DA


        .. note::

            1. update cathyH file given NENS, ENS_times + flag self.parm['DAFLAG']==True

            2. create the processor cathy_origin.exe

            3. create directory tree for the ensemble

            4. create panda dataframe for each perturbated variable

            5. update input files

        Parameters
        ----------
        NENS : int
            # Number of realizations in the ensemble
        ENS_times : int
            # Observation times

        Returns
        -------
        None.

        Note: for now only implemented for EnkF DA
        """

        self.NENS = NENS  # THIS IS TEMPORARY
        # updated_parm_list - update_parm_list


        if hasattr(self,'ens_valid') is False:
            self.ens_valid = list(np.arange(0,self.NENS))

        # create sub directories for each ensemble
        # ---------------------------------------------------------------------
        self._create_subfolders_ensemble(NENS)


        # update list of updated parameters based on problem heterogeneity
        # ---------------------------------------------------------------------

        updated_parm_list = ['St. var.']
        for d in self.dict_parm_pert.keys():
            for l in update_parm_list:
                if l in d:
                    updated_parm_list.append(d)

        return updated_parm_list


    def _DA_analysis(self, prediction,
                     DA_type='enkf_Evensen2009_Sakov',
                     list_update_parm=['St. var.'],
                     list_assimilated_obs='all',
                     ens_valid=[]):
        """
        THIS SHOULD BE MOVED TO DA CLASS

        Analysis ensemble using DA
        print(list_assimilated_obs)


        1. map state variable 2 Observations

        2. apply filter

        3. checkings

        Parameters
        ----------
        typ : str
            type of Data Assimilation
        NUDN : int
            #  NUDN  = number of observation points for nudging or EnKF (NUDN=0 for no nudging)
        ENS_times : np.array([])
            # Observation time for ensemble kalman filter (EnKF).
        rejected_ens : list
            list of ensemble to reject
        Returns
        -------
        None.

        Note: for now only implemented for EnkF DA
        """

        update_key = 'ini_perturbation'
        if self.count_DA_cycle > 0:
            update_key = 'update_nb' + str(self.count_DA_cycle)

        # find the method to map for this time step
        # --------------------------------------------------------
        obskey2map, obs2map = self._obs_key_select(list_assimilated_obs)

        # prepare states
        # ---------------------------------------------------------------------
        ensemble_psi, ensemble_sw, ens_size, sim_size = self._read_state_ensemble()

        # ---------------------------------------------------------------------
        # prepare data
        # ---------------------------------------------------------------------
        # select ONLY data to assimilate
        # ---------------------------------------------------------------------
        data, _ = self._get_data2assimilate(list_assimilated_obs)
        print('data size: ' + str(len(data)))

        # check size of data_cov
        # ---------------------------------------------------------------------
        # if len(self.stacked_data_cov)/len(self.dict_obs.items()) != len(data):
        # if 'ERT' in sensors:
        #     if np.shape(self.stacked_data_cov[self.count_DA_cycle])[0] != len(data[0]):
        #         raise ValueError('need to compute data covariance')
        # else:
        if np.shape(self.stacked_data_cov[self.count_DA_cycle])[0] != len(data):
            raise ValueError('need to compute data covariance')

        data_cov =self.stacked_data_cov[self.count_DA_cycle]



        # # filter non-valid ensemble before Analysis
        # # ---------------------------------------------------------------------
        prediction_valid = prediction
        ensemble_psi_valid = ensemble_psi[:,ens_valid]
        ensemble_sw_valid = ensemble_sw[:,ens_valid]
        # print('prediction valid size: ' + str(len(prediction_valid[:,0])))



        # prepare parameters
        # ---------------------------------------------------------------------
        # When updating the states only, the elements of X are the
        # pressure heads at each node of the finite element grid, while
        # the state augmentation technique is used when also updat-
        # ing the parameters
        param = self.transform_parameters(list_update_parm,update_key)
        print('parm size: ' + str(len(param)))

        # if update_key == 'ini_perturbation':
        #     param = self.spatialize_parameters(list_update_parm,param)

        param_valid = []
        if len(param)>0:
            param_valid = param[:, ens_valid]


        # run Analysis
        # ---------------------------------------------------------------------

        result_analysis = run_analysis(DA_type,
                                        data,
                                        data_cov,
                                        param_valid,
                                        list_update_parm,
                                        [ensemble_psi_valid,ensemble_sw_valid],
                                        prediction_valid,
                                        alpha=self.damping,
                                        )


        # plot ensemble covariance matrices and changes (only for ENKF)
        # ---------------------------------------------------------------------
        if len(result_analysis)>2:
            [A, Amean, dA,
             dD, MeasAvg, S,
             COV, B, dAS,
             analysis,
             analysis_param]  = result_analysis

            try:
                plt_CT.show_DA_process_ens(ensemble_psi_valid,
                                            data,COV,
                                            dD,dAS,B,analysis,
                                            savefig=True,
                                            savename= os.path.join(self.workdir,
                                                                   self.project_name,
                                                                   'DA_Ensemble',
                                                                   'DA_Matrices_t'
                                                                   +str(self.count_DA_cycle)
                                                                  ),
                                            # label_sensor=str([*self.dict_obs[0]])
                                            label_sensor=str([*obskey2map])
                                            )
            except:
                self.console.print("[b]Impossible to plot matrices[/b]")


        # particule filter
        # ----------------
        else:
            [analysis,
             analysis_param]  = result_analysis


        # Back-transformation of the parameters
        # ---------------------------------------------------------------------
        # for pp in enumerate(list_update_parm[:]):
        #     if 'St. var.' in pp[1]:
        #         pass
        #     else:

        # analysis_param_valid = []
        # if len(param)>0:
        #     analysis_param_valid = analysis_param[ens_valid]

        self.transform_parameters(list_update_parm,
                                  param=np.hstack(analysis_param),
                                  back=True)

                # if self.dict_parm_pert[pp[1]]['transf_type'] == 'log':
                #     print('back log transformation')
                #     analysis_param = np.exp(analysis_param)



        # return ensemble_psi, Analysis, AnalysisParam, Data, data_cov
        # ---------------------------------------------------------------------
        return(ensemble_psi, ensemble_sw, data,  analysis,
                 analysis_param)


    def _check_before_analysis(self, ens_valid=[],
                               threshold_rejected=10):
        '''
        Filter is applied only on selected ensemble

        Parameters
        ----------
        update_key : TYPE
            DESCRIPTION.

        Returns
        -------
        rejected_ens

        '''

        self.console.print(":white_check_mark: [b]check scenarii before analysis[/b]")

        ###CHECK which scenarios are OK and which are to be rejected and then create and applied the filter
        ###only on the selected scenarios which are to be considered.

        rejected_ens_new = []
        for n in range(self.NENS):

            if n in ens_valid:
                # print('ensemble_nb ' + str(n) + ' before update analysis')
                try:
                    df_mbeconv = out_CT.read_mbeconv(os.path.join(self.workdir, self.project_name,
                                                       "DA_Ensemble/cathy_" + str(n+1),
                                                       'output/mbeconv'))
                    if df_mbeconv['CUM.'].isnull().values.any():
                        rejected_ens_new.append(True)
                        self.console.print(":x: [b]df_mbeconv['CUM.'][/b]" + str(df_mbeconv['CUM.'].isnull().values.any())+ ', ens_nb:' + str(n+1))

                    elif (np.round(df_mbeconv['TIME'].iloc[-1]))<np.round(self.parm['(TIMPRT(I),I=1,NPRT)'][-1]) - self.parm['DELTAT']:
                        rejected_ens_new.append(True)
                        self.console.print(":x: [b]df_mbeconv['time'][/b]" + str(df_mbeconv['TIME'].iloc[-1]) + ', ens_nb:' + str(n+1))

                    elif len(df_mbeconv) == 0:
                        rejected_ens_new.append(True)
                        self.console.print(":x: [b]len(df_mbeconv)['TIME'][/b]" + str(len(df_mbeconv))+ ', ens_nb:' + str(n+1))

                    else:
                        rejected_ens_new.append(False)

                except:
                    rejected_ens_new.append(True)
                    print('cannot read mbeconv')
                    pass
                    # UnboundLocalError: local variable 'df_mbeconv' referenced before assignment

            else:
                rejected_ens_new.append(True)

        #check whether the number of discarded scenarios is major than the 10% of the total number [N];
        #if the answer is yes than the execution stops
        # --------------------------------------------------------------------
        if sum(rejected_ens_new)>=(threshold_rejected/100)*self.NENS:
            print('new parameters (if updated)')
            print(self.dict_parm_pert)

            raise ValueError('% number of rejected ensemble is too high:' + str((sum(rejected_ens_new)*100)/(self.NENS)))

        return  rejected_ens_new




    def _check_after_analysis(self,update_key,list_update_parm):
        '''
        CHECK which scenarios parameters are OK and which one to discard

        Returns
        -------
        None.

        '''
        self.console.print(":white_check_mark: [b]check scenarii post update[/b]")

        id_valid = [list(np.arange(0,self.NENS))]
        id_nonvalid = []


        def test_negative_values(update_key,pp):
            id_nonvalid.append(list(np.where(self.dict_parm_pert[pp[1]][update_key]<0)[0]))
            id_valid.append(list(np.where(self.dict_parm_pert[pp[1]][update_key]>0)[0]))

            for i in id_nonvalid:
                self.console.print(":x: [b]negative" + str(pp) + ":[/b]" +
                                   str(self.dict_parm_pert[pp[1]][update_key][i])+
                                   ', ens_nb:' + str(i))
            return id_nonvalid, id_valid

        def test_range_values(update_key,pp,min_r,max_r):
            id_nonvalid.append(list(np.where(self.dict_parm_pert[pp[1]][update_key]<min_r)[0]))
            id_nonvalid.append(list(np.where(self.dict_parm_pert[pp[1]][update_key]>max_r)[0]))
            id_valid.append(list(np.where(self.dict_parm_pert[pp[1]][update_key]>min_r)[0]))
            id_valid.append(list(np.where(self.dict_parm_pert[pp[1]][update_key]<max_r)[0]))
            for i in id_nonvalid:
                self.console.print(":x: [b]out of range" + str(pp) + ":[/b]" +
                                   str(self.dict_parm_pert[pp[1]][update_key][i])+
                                   ', ens_nb:' + str(i))
            return id_nonvalid, id_valid



        for pp in enumerate(list_update_parm[:]):
            if 'St. var.' in pp[1]:
                pass
            else:
                if 'rFluid'.casefold() in pp[1].casefold():
                    id_nonvalid, id_valid = test_negative_values(update_key,pp)
                elif 'porosity'.casefold() in pp[1].casefold():
                    id_nonvalid, id_valid = test_negative_values(update_key,pp)
                elif 'a_Archie'.casefold() in pp[1].casefold():
                    id_nonvalid, id_valid = test_range_values(update_key,pp, 0, 3)
                elif 'Ks'.casefold() in pp[1].casefold():

                    id_nonvalid.append(list(np.where(self.dict_parm_pert[pp[1]][update_key]<0)[0]))
                    id_valid.append(list(np.where(self.dict_parm_pert[pp[1]][update_key]>0)[0]))
                    for i in id_nonvalid:
                        self.console.print(":x: [b]negative Ks:[/b]" +
                                           str(self.dict_parm_pert[pp[1]][update_key][i])+
                                           ', ens_nb:' + str(i))
                elif 'PCREF'.casefold() in pp[1].casefold():
                    id_nonvalid.append(list(np.where(self.dict_parm_pert[pp[1]][update_key]>0)[0]))
                    id_valid.append(list(np.where(self.dict_parm_pert[pp[1]][update_key]<0)[0]))
                    for i in id_nonvalid:
                        self.console.print('Positive values encountered in PCREF:' +
                              str(self.dict_parm_pert[pp[1]][update_key][i])+
                                           ', ens_nb:' + str(i))
                # ckeck if new value of Zroot is feasible
                # ------------------------------------------------------------
                elif 'Zroot'.casefold() in pp[1].casefold():
                    id_nonvalid.append(list(np.where(self.dict_parm_pert[pp[1]][update_key]>
                                                abs(min(self.grid3d['nodes_idxyz'][:, -1])))[0]))
                    id_valid.append(list(np.where(self.dict_parm_pert[pp[1]][update_key]<
                                             abs(min(self.grid3d['nodes_idxyz'][:, -1])))[0])
                                    )
                    for i in id_nonvalid:
                        self.console.print(":x: [b]unfeasible root depth:[/b]" +
                                           str(self.dict_parm_pert[pp[1]][update_key][i])+ ', ens_nb:' + str(i))

        id_nonvalid_flat = [item for sublist in id_nonvalid for item in sublist]
        id_valid = set.difference(set(list(np.arange(0,self.NENS))), set(list(np.unique(id_nonvalid_flat))))
        return id_valid


        pass


    def spatialize_parameters(self,list_update_parm,parm):
        '''
        Extend the number of parameters according to the number of zones
        (useful for the heteregeneous case)
        '''
        i=0
        parm_extended = []
        for l in list_update_parm:
            print(l)
            if 'St. var.' in l:
                pass
            else:
                nb_of_zones = 1
                if 'surf_zones_param' in self.dict_parm_pert[l]:
                    nb_of_zones =  self.dict_parm_pert[l]['surf_zones_param']
                    print(nb_of_zones)
                    parm_extended.append(np.tile(parm[i,:],(nb_of_zones,1)))
                i = i + 1
        parm_extended = np.vstack(parm_extended)
        return parm_extended


    def transform_parameters(self,list_update_parm,
                             update_key=None,
                             back=False,
                             param=None):
        '''
        Parameter (back) transformation before analysis (log, bounded log, ...)

        ..note:: parameters are log transform so that they become Gaussian distributed
                 to be updated by the Ensemble Kalman Filter
        '''
        param_new = []
        for pp in enumerate(list_update_parm[:]):
            param_trans = []

            if 'St. var.' in pp[1]:
                pass
            else:
                if update_key:
                    param = self.dict_parm_pert[pp[1]][update_key]
                if self.dict_parm_pert[pp[1]]['transf_type'] == 'log':
                    print('log transformation')
                    if back:
                        param_trans = np.exp(param)
                    else:
                        param_trans = np.log(param)
                elif self.dict_parm_pert[pp[1]]['transf_type'] == 'log-ratio':
                   print('bounded log transformation')
                   range_tr = self.dict_parm_pert[pp[1]]['transf_bounds']
                   A = range_tr['A'] # min range
                   B = range_tr['B'] # max range
                   if back:
                       param_trans = (
                                       np.exp(  (param - A)/
                                                (B - param)
                                              )
                                    )
                   else:
                       param_trans = (
                                       np.log(  (param - A)/
                                                (B - param)
                                              )
                                    )
                elif self.dict_parm_pert[pp[1]]['transf_type'] == 'hyperbolic':
                   print('hyperbolic transformation')
                   # Y = sinh−1(U )
                   param_trans = np.log(self.dict_parm_pert[pp[1]][update_key])

                else:
                    param_trans = param
                param_new.append(param_trans)
            if len(param_new)>0:
                np.vstack(param_new)

        return np.array(param_new)


    def perturbate_parm(var_per, parm, type_parm, mean=[], sd=[], per_type=None,
                        sampling_type = 'lognormal',
                        ensemble_size = 128, seed=True,
                        **kwargs):
        '''
        Perturbate parameter for ensemble generation

        Possible variable to perturbate:
            - initial conditions
            - hyetograph atmbc
            - van Genuchten retention curves parameters
            - Feddes parameters

        Perturbation:
            - parameters with constrainsts (truncated distribution)
            - parameters time dependant
            - other types of parameters

        Returns
        -------
        var_per : dict
            Perturbated variable dictionnary

        '''

        var_per_2add = {}


        # copy initiail variable dict and add 'sampling' and 'ini_perturbation' attributes
        # -------------------------------------------------------------------------
        var_per_2add[type_parm] = parm

        key = 'sampling_type'
        var_per_2add[type_parm][key] = sampling_type
        key = 'sampling_mean'
        var_per_2add[type_parm][key] = mean
        key = 'sampling_sd'
        var_per_2add[type_parm][key] = sd
        key = 'per_type'
        var_per_2add[type_parm][key] = per_type

        # Contrainsted perturbation (bounded)
        # --------------------------------------------------------------------
        if 'Archie' in type_parm:
            parm_per_array = Archie_pert_rules(parm,
                                               type_parm,
                                               ensemble_size,
                                               mean,sd,
                                               per_type,
                                               sampling_type
                                               )
        elif 'porosity' in type_parm:
            parm_sampling = sampling_dist_trunc(myclip_a=0,
                                                myclip_b=1,
                                                ensemble_size=ensemble_size,
                                                loc=mean,
                                                scale=sd)
            parm_per_array = perturbate_dist(parm,per_type,parm_sampling,ensemble_size)

        elif 'zroot'.casefold() in type_parm.casefold():
            parm_sampling = sampling_dist_trunc(myclip_a=var_per_2add[type_parm]['clip_min'],
                                                     myclip_b=var_per_2add[type_parm]['clip_max'],
                                                     ensemble_size=ensemble_size,
                                                     loc=mean,
                                                     scale=sd
                                                     )
            parm_per_array = perturbate_dist(parm,per_type,parm_sampling,ensemble_size)



        # check if parameters in part of van Genuchten retention curves
        #----------------------------------------------------------------------
        # if type_parm in ['Alpha', 'nVG', 'thethaR']: #van Genuchten retention curves
        elif 'VG' in type_parm: #van Genuchten retention curves
            parm_per_array = VG_pert_rules(var_per_2add,parm,type_parm,ensemble_size,mean,sd,per_type,sampling_type,**kwargs)

        # Time dependant perturbation 
        # --------------------------------------------------------------------
        elif 'atmbc' in type_parm:
            var_per_2add = atmbc_pert_rules()
            
        # For all other types of perturbation 
        # --------------------------------------------------------------------
        else:
            if 'clip_min' in var_per_2add[type_parm].keys():
                parm_sampling = sampling_dist_trunc(myclip_a=var_per_2add[type_parm]['clip_min'],
                                                         myclip_b=var_per_2add[type_parm]['clip_max'],
                                                         ensemble_size=ensemble_size,
                                                         loc=mean,
                                                         scale=sd
                                                         )
            else:
                parm_sampling = sampling_dist(sampling_type,mean,sd,ensemble_size)
            parm_per_array = perturbate_dist(parm,per_type,parm_sampling,ensemble_size)


        # build dictionnary of perturbated variable
        # --------------------------------------------------------------------
        var_per_2add = build_dict_attributes_pert(var_per_2add,
                                                  type_parm,
                                                  parm_per_array,
                                                  parm_sampling,
                                                  **kwargs
                                                  )
                                                  
        # Add to var perturbated stacked dict
        # ----------------------------------
        var_per = var_per | var_per_2add


        if kwargs['savefig']:
            plt_CT.plot_hist_perturbated_parm(parm,var_per,type_parm,parm_per_array,
                                             **kwargs
                                             )

        return var_per


    def update_ENS_files(self, parm_pert, update_parm_list, **kwargs):
        '''
        THIS SHOULD BE MOVED TO DA CLASS

        Update by overwriting ensemble files (usually after analysis step or initially to build the ensemble)
        - update state
        - update model parameters:
            - initial conditions
            - Archie
            - Hydraulic conductivity ks
            - Feddes parameters
            - Atmbc boundary conditions
        Returns
        -------
        Overwritten files.

        '''
        update_key = 'ini_perturbation'
        if 'cycle_nb' in kwargs:
            if kwargs['cycle_nb'] > 0:
                update_key = 'update_nb' + str(kwargs['cycle_nb'])

        if 'analysis' in kwargs:
            analysis = kwargs['analysis']


        if update_parm_list=='all':
            update_parm_list = ['St. var.']
            for pp in parm_pert:
                update_parm_list.append(pp)

        # loop over ensemble files to update ic -->
        # this is always done since psi are state variable to update
        # ------------------------------------------------------------------
        for ens_nb in range(self.NENS):
            # change directory according to ensmble file nb
            os.chdir(os.path.join(self.workdir, self.project_name,
                                  './DA_Ensemble/cathy_' + str(ens_nb+1)))

            # state variable update use ANALYSIS update or not
            # --------------------------------------------------------------
            if 'St. var.' in update_parm_list:
                if kwargs['cycle_nb'] > 0:
                        df_psi = self.read_outputs(filename='psi',
                                                  path=os.path.join(os.getcwd(),
                                                                    self.output_dirname))
                        self.update_ic(INDP=1, IPOND=0,
                                       pressure_head_ini=analysis[:,ens_nb],
                                       filename=os.path.join(os.getcwd(), 'input/ic'),
                                       backup=True)
            else:
                if kwargs['cycle_nb'] > 0:
                    raise ValueError('no state variable update - use last iteration as initial conditions ?')
                    df_psi = self.read_outputs(filename='psi',
                                              path=os.path.join(os.getcwd(),
                                                                self.output_dirname))
                    if kwargs['cycle_nb'] > 0:
                            self.update_ic(INDP=1, IPOND=0,
                                           pressure_head_ini=df_psi[-1, :],
                                           filename=os.path.join(os.getcwd(), 'input/ic'),
                                           backup=False)


        # loop over dict of perturbated variable
        # ----------------------------------------------------------------------
        FeddesParam_mat_ens = [] # matrice of ensemble for Feddes parameters
        Archie_parms_mat_ens = [] # matrice of ensemble for Archie parameters
        VG_parms_mat_ens = [] # matrice of ensemble for VG parameters
        PERMX_het_ens = [] # matrice of ensemble for hydraulic conductivity parameters
        Archie_p_names = ['porosity', 'rFluid_Archie','a_Archie','m_Archie','n_Archie','pert_sigma_Archie']
        VG_p_possible_names = ['n_VG', 'thetar_VG', 'alpha_VG','VGPSATCELL_VG']
        VG_p_possible_names_positions_in_soil_table = [5,6,7,7]

        for parm_i, key in enumerate(update_parm_list):  # loop over perturbated variables dict
            key_root = re.split('(\d+)', key)
            if len(key_root)==1:
                key_root.append('0')

            # ------------------------------------------------------------------
            for ens_nb in range(self.NENS):
                # change directory according to ensmble file nb
                os.chdir(os.path.join(self.workdir, self.project_name,
                                      './DA_Ensemble/cathy_' + str(ens_nb+1)))

                # atmbc update
                # --------------------------------------------------------------
                if key == 'atmbc':
                    print('hietograph perturbated not yet implemented')
                    self.update_atmbc(verbose=True)

                # ic update (i.e. water table position update)
                # --------------------------------------------------------------
                elif key.casefold() in 'ic'.casefold():
                    if kwargs['cycle_nb']==0:
                        self.update_ic(INDP=0, IPOND=0,
                                       pressure_head_ini=parm_pert[key
                                           ]['ini_perturbation'][ens_nb],
                                       filename=os.path.join(os.getcwd(), 'input/ic'),
                                       backup=True)

                # kss update
                # --------------------------------------------------------------
                elif key_root[0].casefold() in 'ks'.casefold():
                    if len(PERMX_het_ens)==0:
                        SPP =  self.soil_SPP['SPP_map']
                        PERMX_het_ens = np.ones([len(SPP['PERMX']),self.NENS])
                        PERMX_het_ens[:] = np.nan

                    PERMX_het_ens[int(key_root[1]),ens_nb]=parm_pert[key][update_key][ens_nb]

                    SPP_ensi = []
                    for es in range(self.NENS):
                        spd = dict()
                        spd['PERMX'] = list(PERMX_het_ens[:,ens_nb])
                        spd['PERMY'] = list(PERMX_het_ens[:,ens_nb])
                        spd['PERMZ'] = list(PERMX_het_ens[:,ens_nb])
                        SPP_ensi.append(spd)

                    SPP.update(SPP_ensi[ens_nb])

                    # print(PERMX_het_ens[:,ens_nb])
                    if np.isnan(PERMX_het_ens[:,ens_nb]).any()==False:
                        self.update_soil(SPP=SPP,
                                          FP=self.soil_FP['FP_map'],
                                          verbose=True,
                                          filename=os.path.join(os.getcwd(), 'input/soil')
                                          )  # specify filename path
                # VG parameters update
                # --------------------------------------------------------------
                elif key_root[0] in VG_p_possible_names:

                    # equivalence_CATHY = {
                    #                      'thetar_VG':'VGRMCCELL',
                    #                      'alpha_VG':'-1/VGPSATCELL',
                    #                      'n_VG':'VGNCELL',
                    #                      }

                    # equivalence_CATHY
                    # retention curves parameters VGN, VGRMC, and VGPSAT
                    # - 'VGNCELL' (NSTR, NZONE): van Genuchten curve exponent  = n
                    # - 'VGRMCCELL' (NSTR, NZONE): residual moisture content = \thetaR
                    # - 'VGPSATCELL' (NSTR, NZONE): van Genuchten curve exponent -->
                    #                               VGPSAT == -1/alpha (with alpha expressed in [L-1]);

                    SPP =  self.soil_SPP['SPP_map']
                    if len(VG_parms_mat_ens)==0:
                        VG_parms_mat = self.soil_SPP['SPP']
                        np.shape(VG_parms_mat)
                        if len(SPP['PERMX'])<2:
                            VG_parms_mat_ens = np.matlib.repmat(VG_parms_mat[0], self.NENS, 1)
                        else:
                            # VG_parms_mat_ens = np.repeat(VG_parms_mat, self.NENS, axis=0)
                            print('Non homogeneous soil VG properties not implemented')


                    if len(SPP['PERMX'])<2:
                        # for k in parm_pert.keys():
                        for k in VG_p_possible_names:
                            match = [parm_incr for parm_incr in parm_pert.keys() if k in parm_incr]
                            if len(match)>0:
                                # print(match)
                                # print(k)
                                idVG = VG_p_possible_names_positions_in_soil_table[VG_p_possible_names.index(k)]
                                # print(idVG)
                                VG_parms_mat_ens[:,idVG] = parm_pert[match[0]][update_key]
                            else:
                                pass
                    else:
                        # VG_parms_mat_ens[:,:,idVG] = np.matlib.repmat(parm_pert[key][update_key],len(VG_parms_mat),1).T
                        print('Non homogeneous soil VG properties not implemented')

                    VG_parms_ensi = []
                    for es in range(self.NENS):
                        fed = dict()
                        for i, f in  enumerate(list(SPP.keys())):
                            if len(SPP['PERMX'])<2:
                                if f == 'alpha_VG': # convert to CATHY VGPSATCELL == -1/alpha
                                    fed[f] = -1/VG_parms_mat_ens[es,i]
                                else:
                                    fed[f] = [VG_parms_mat_ens[es,i]]
                            else:
                                print('Non homogeneous soil VG properties not implemented')
                                # fed[f] = list(VG_parms_mat_ens[es,:,i])
                        VG_parms_ensi.append(fed)

                    self.update_soil(SPP=VG_parms_ensi[ens_nb],
                                     FP=self.soil_FP['FP_map'],
                                     verbose=False,
                                     filename=os.path.join(os.getcwd(), 'input/soil'))  # specify filename path

                # FeddesParam update
                # --------------------------------------------------------------
                elif key_root[0] in ['PCANA', 'PCREF', 'PCWLT', 'ZROOT', 'PZ', 'OMGC']:
                    SPP =  self.soil_SPP['SPP_map']

                    if len(FeddesParam_mat_ens)==0:
                        FeddesParam_mat = self.soil_FP['FP']
                        FeddesParam_mat_ens = np.repeat(FeddesParam_mat, self.NENS, axis=0)
                        FeddesParam_mat_ens = np.reshape(FeddesParam_mat_ens,(self.NENS,self.cathyH["MAXVEG"],6))

                    idFeddes = ['PCANA', 'PCREF', 'PCWLT', 'ZROOT', 'PZ', 'OMGC'].index(key_root[0])


                    # if self.cathyH["MAXVEG"]<2: # Case where only 1 vegetation type
                    FeddesParam_mat_ens[:,int(key_root[1]),idFeddes] = parm_pert[key][update_key]
                    # FeddesParam_mat_ens[:,:,idFeddes] = parm_pert[key][update_key]
                    # else: # Case where only 1 vegetation type


                    FeddesParam_ensi = []
                    for es in range(self.NENS):
                        fed = dict()
                        for i, f in  enumerate(['PCANA', 'PCREF', 'PCWLT', 'ZROOT', 'PZ', 'OMGC']):
                            fed[f] = list(FeddesParam_mat_ens[es,:,i])
                        FeddesParam_ensi.append(fed)
                    self.update_soil(FP=FeddesParam_ensi[ens_nb],
                                     SPP=SPP,
                                     verbose=False,
                                     filename=os.path.join(os.getcwd(), 'input/soil'))  # specify filename path

                # Archie_p update
                # --------------------------------------------------------------
                elif key_root[0] in Archie_p_names:
                    # self.Archie_parms = {'rFluid':rFluid, 'a':a, 'm':m, 'n':n, 'pert_sigma':pert_sigma}
                    idArchie = Archie_p_names.index(key_root[0])
                    if len(Archie_parms_mat_ens)==0:
                        for p in Archie_p_names:
                            if len(self.Archie_parms[p]) != self.NENS:
                                self.Archie_parms[p]=list(self.Archie_parms[p]*np.ones(self.NENS))
                        Archie_parms_mat_ens = np.zeros([len(Archie_p_names),self.NENS])
                        Archie_parms_mat_ens[:] = np.nan
                    Archie_parms_mat_ens[idArchie,ens_nb] = parm_pert[key][update_key][ens_nb]

                    if np.isnan(Archie_parms_mat_ens[idArchie,:]).any()==False:
                        self.Archie_parms[key_root[0]]=list(Archie_parms_mat_ens[idArchie,:])

                    # print(self.Archie_parms)

                else:
                    # variable perturbated to ommit: state, Archie param
                    var_pert_to_ommit = ['St. var.']

                    if not key in var_pert_to_ommit:
                        raise  ValueError('key:' + str(key) + ' to update is not existing!')


        pass
    

    def _performance_assessement(self,list_assimilated_obs,
                                 data,
                                 prediction,
                                 t_obs,
                                 **kwargs):
        '''
        THIS SHOULD BE MOVED TO DA CLASS
        (Normalized) root mean square errors (NRMSEs)
        RMSE is compute separately for each observation assimilated


        ----------
        list_assimilated_obs : TYPE
            DESCRIPTION.
        Data : TYPE
            Refers to the measured data.
        Observation : TYPE
            Refers to the simulated data.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''


        OL_bool = False
        if 'openLoop' in kwargs:
            OL_bool = True


        if hasattr(self, 'df_performance') is False:
            # initiate if simulation just started
            # -------------------------------------------------------------
            self.df_performance = pd.DataFrame()
            self.RMSE_avg_stacked = []
            self.RMSE_sensor_stacked = []


            # NAs default
            # -------------------------------------------------------------
            # RMSE_avg = np.nan # root mean square error (RMSE)
            # NMRMSE = np.nan # time-averaged normalized root mror: unexpected indent



        # average differences over the number of ensemble
        # ALL OBSERVATIONS
        # ------------------------------


        all_Obs_diff_mat = np.zeros(np.shape(prediction))
        for i in range(len(data)):
            for j in range(np.shape(prediction)[1]): # Loop over ensemble collumns
                all_Obs_diff_mat[i,j] = abs(data[i]-prediction[i,j])
                # all_Obs_diff_mat[:,j] = abs(data[i]-prediction[i,j])

        all_Obs_diff_avg = np.nansum(all_Obs_diff_mat,axis=1)*(1/len(prediction))

        # average differences over all the observations
        all_Obs_RMSE_avg_ti = np.sum(all_Obs_diff_avg,axis=0)*(1/len(data))



        # # compute metrics for each observation variable
        # # ------------------------------------------
        # Loop trought observation dictionnary for a given assimilation time (count_DA_cycle)
        # -----------------------------------------------------------------
        obs2eval_key = [] # data type 2 eval
        obs2eval = []    # data dict 2 eval

        _ , obs2eval_key = self._get_data2assimilate(list_assimilated_obs)


        start_line_obs = 0
        for s, name_sensor in enumerate(obs2eval_key): # case where there is other observations than ERT
            obs2eval, _ = self._get_data2assimilate([name_sensor], match=True)

            # number of observation at a given time
            # for ERT number of observation = is the number of grid cells
            n_obs = len(obs2eval)
            prediction2eval = prediction[start_line_obs:start_line_obs+n_obs]
            obs2eval_diff_mat = np.zeros(np.shape(prediction2eval))

            for i in range(len(obs2eval)):
                for j in range(np.shape(prediction2eval)[1]): # Loop over ensemble collumns
                    obs2eval_diff_mat[i,j] = abs(obs2eval[i]-prediction2eval[i,j])

            obs2eval_diff_avg = np.nansum(obs2eval_diff_mat,axis=1)*(1/len(prediction2eval))

            start_line_obs = start_line_obs + n_obs

            if 'ERT' in name_sensor:
                RMSE_sensor_ti = np.sum(obs2eval_diff_avg,axis=0)*(1/len(data))
                # Plot here scatter Transfer resistance scatter plot between the observed and simulated data
                # plt.scatter(obs2eval[i],prediction2eval[i,j])
            else:
                RMSE_sensor_ti = obs2eval_diff_avg


            # compute normalised RMSE
            #---------------------------------------------------------

            self.RMSE_avg_stacked.append(all_Obs_RMSE_avg_ti)
            self.RMSE_sensor_stacked.append(RMSE_sensor_ti)

            # if hasattr(self, 'RMSE') is False:
            if t_obs == 0:
                NMRMSE_sensor_ti = RMSE_sensor_ti
                NMRMSE_avg_ti = all_Obs_RMSE_avg_ti
            else:
                NMRMSE_sensor_ti = (1/(t_obs+1))*np.sum(self.RMSE_sensor_stacked)
                NMRMSE_avg_ti = (1/(t_obs+1))*np.sum(self.RMSE_avg_stacked)




            # root names for the collumns name
            # -------------------------------------------------------------
            cols_root = ['time',
                         'ObsType',
                         'RMSE'+ name_sensor,
                         'RMSE_avg',
                         'NMRMSE'+ name_sensor,
                         'NMRMSE_avg',
                         'OL']

            # root data
            # -------------------------------------------------------------
            data_df_root = [[t_obs],
                            [name_sensor],
                            [RMSE_sensor_ti],
                            [all_Obs_RMSE_avg_ti],
                            [NMRMSE_sensor_ti],
                            [NMRMSE_avg_ti],
                            [OL_bool]]


            df_performance_ti = pd.DataFrame(np.array(data_df_root,dtype=object).T,
                                  columns=cols_root)


            # concatenate with main RMSE dataframe
            # -------------------------------------------------------------
            self.df_performance= pd.concat([self.df_performance, df_performance_ti],
                                       axis=0,
                                       ignore_index=True)




        return self.df_performance




    def _get_data2assimilate(self,list_assimilated_obs,time_ass=None,match=False):
        '''
        Loop over observation dictionnary and select corresponding data for a given assimilation time

        Parameters
        ----------
        list_assimilated_obs : list
            list of observation to assimilate.

        Returns
        -------
        data : list
            list of data values.
        '''

        # find the method to map for this time step
        # --------------------------------------------------------
        obskey2map, obs2map = self._obs_key_select(list_assimilated_obs)

        def extract_data(sensor, time_ass, data):
            data2add = []
            if 'ERT' in sensor:
                if 'pygimli' in items_dict[time_ass][1][sensor]['data_format']:
                # if 'pygimli' in obs2map[time_ass]['data_format']:
                # if 'pygimli' in self.dict_obs[time_ass]['ERT']['data_format']:
                    data2add.append(np.array(items_dict[time_ass][1][sensor]['data']['rhoa']))
                else:
                    data2add.append(items_dict[time_ass][1][sensor]['data']['resist'].to_numpy())
                data2add= np.hstack(data2add)
            else:
                data2add.append(items_dict[time_ass][1][sensor]['data'])

            return data2add


        if time_ass is None:
            time_ass = self.count_DA_cycle

        data = []
        # Loop trought observation dictionnary for a given assimilation time (count_DA_cycle)
        # -----------------------------------------------------------------
        items_dict = list(self.dict_obs.items())
        for sensor in items_dict[time_ass][1].keys():
            if 'all' in list_assimilated_obs:
                # obskey2map.append(sensor)
                data2add = extract_data(sensor, time_ass, data)
                data.append(data2add)
            else:
                for l in list_assimilated_obs:
                    if match:
                        if l == sensor:
                            # obskey2map.append(sensor)
                            data2add = extract_data(sensor, time_ass, data)
                            data.append(data2add)
                    else:
                        if l in sensor:
                            # obskey2map.append(sensor)
                            data2add = extract_data(sensor, time_ass, data)
                            data.append(data2add)
            # len(data2add)
        return np.hstack(data), obskey2map


    def _obs_key_select(self,list_assimilated_obs):
        ''' find data to map from dictionnary of observations '''
        # Loop trought observation dictionnary for a given assimilation time (count_DA_cycle)
        # -----------------------------------------------------------------
        obskey2map = [] # data type 2 map
        obs2map = []    # data dict 2 map
        items_dict = list(self.dict_obs.items())
        for sensor in items_dict[self.count_DA_cycle][1].keys():
            if 'all' in list_assimilated_obs:
                obskey2map.append(sensor)
                obs2map.append(items_dict[self.count_DA_cycle][1][sensor])
            else:
                for l in list_assimilated_obs:
                    if l in sensor:
                        obskey2map.append(sensor)
                        obs2map.append(items_dict[self.count_DA_cycle][1][sensor])
        return obskey2map,obs2map




    def map_states2Observations(self,
                                 list_assimilated_obs='all',
                                 parallel=False,
                                 default_state = 'psi',
                                 verbose = False,
                                 **kwargs):
        '''
        Translate (map) the state values (pressure head or saturation) to observations (or predicted) values


        In other words: apply H mapping operator to convert model to predicted value Hx

        In practice, for each assimilation time, loop over the observation
        dictionnary and infer observation type. Stack Hx.
        Hx = [Hx0,Hx1,Hx_nb_of_observation]

        .. note::

            - For ERT data, H is Archie law, SWC --> ER0 --> ERapp
            - For SWC data, H is a multiplicator (porosity)
            - For tensiometer data, no mapping needed or VGP model if model state used is saturation
            - For discharge, no mapping needed
            - For scale data: no mapping needed


        .. note::

            - For ERT data each ERT mesh node is an observation

        THIS SHOULD BE MOVED TO DA CLASS

        Parameters
        ----------
        state : np.array([])
            [pressure heads, saturation] for each mesh nodes. The default is [None, None].
        list_assimilated_obs : TYPE, optional
            DESCRIPTION. The default is 'all'.
        parallel : Bool, optional
            DESCRIPTION. The default is False.
        verbose : Bool, optional
            DESCRIPTION. The default is False.
        **kwargts : TYPE
            DESCRIPTION.

        Returns
        -------
        Hx : np.array
            Ensemble of the simulated (predicted) observations.
        '''


        if parallel:
            path_fwd_CATHY_list = [] # list of ensemble path of cathy output
            for ens_nb in self.ens_valid: # loop over ensemble files
                path_fwd_CATHY_list.append(os.path.join(self.workdir,
                                        self.project_name,
                                        'DA_Ensemble/cathy_' + str(ens_nb+1)))




        Hx_ens = [] # matrice of predicted observation for each ensemble realisation
        for ens_nb in self.ens_valid: # loop over ensemble files
            # print('ens nb:' + str(ens_nb))

            path_fwd_CATHY = os.path.join(self.workdir, self.project_name,
                                  './DA_Ensemble/cathy_' + str(ens_nb+1))

            df_psi = self.read_outputs(filename='psi',
                                       path=os.path.join(path_fwd_CATHY,
                                                         self.output_dirname)
                                       )
            df_sw = self.read_outputs(filename='sw',
                                      path=os.path.join(path_fwd_CATHY,
                                                        self.output_dirname)
                                      )

            # infer soil parameters properties
            # ---------------------------------
            porosity = self.soil_SPP['SPP'][:, 4][0]

            # find data to map with dictionnary of observations
            # --------------------------------------------
            obskey2map, obs2map = self._obs_key_select(list_assimilated_obs)
            state = [df_psi[-1],df_sw[-1]]


            Hx_stacked= [] # stacked predicted observation
            # Loop over observations to map
            # ---------------------------------------------------------------------
            for i, obs_key in enumerate(obskey2map):

                # print('obs_key nb:' + obs_key)

                if 'tensio' in obs_key:
                    # case 1: pressure head assimilation (Hx_PH)
                    # -------------------------------------------------------------
                    Hx_PH = state[0][obs2map[i]['mesh_nodes']]
                    Hx_stacked.append(Hx_PH)

                if 'swc' in obs_key:
                    # case 2: sw assimilation (Hx_SW)
                    # --------------------------------------------------------------------
                    Hx_SW = state[1][obs2map[i]['mesh_nodes']] * porosity
                    Hx_stacked.append(Hx_SW)
                    # note: the value of the porosity can be unique or not depending on the soil physical properties defined

                if 'scale' in obs_key:
                    #Atmpot-vf (9) : Potential atmospheric forcing (rain +ve / evap -ve) as a volumetric flux [L^3/T]
                    #Atmpot-v (10) : Potential atmospheric forcing volume [L^3] (See parm input file for units)
                    #Atmpot-r (11) : Potential atmospheric forcing rate [L/T]
                    #Atmpot-d (12) : Potential atmospheric forcing depth [L]
                    #Atmact-vf(13) : Actual infiltration (+ve) or exfiltration (-ve) at atmospheric BC nodes as a volumetric flux [L^3/T]
                    #Atmact-v (14) : Actual infiltration (+ve) or exfiltration (-ve) volume [L^3]
                    #Atmact-r (15) : Actual infiltration (+ve) or exfiltration (-ve) rate [L/T]
                    #Atmact-d (16) : Actual infiltration (+ve) or exfiltration (-ve) depth [L]

                    print('Not yet implemented')

                    df_dtcoupling = self.read_outputs(filename='dtcoupling',
                                              path=os.path.join(path_fwd_CATHY,
                                                                self.output_dirname)
                                              )

                    Hx_scale = state[1][obs2map[i]['mesh_nodes']] * porosity
                    Hx_stacked.append(Hx_scale)



                if 'discharge' in obs_key:
                    # case 3: discharge
                    # need to read the hgsfdet file (Hx_Q)
                    # --------------------------------------------------------------------
                    # if key[0] in 'discharge':
                    # derivation of the dircharge, Q from file 'hgsfdet'
                    # Hx_Q = []
                    # Hx.vstack(Hx_Q)
                    print('Not yet implemented')

                if 'ERT' in obs_key:
                    if parallel == False:
                        Hx_ERT = self._map_ERT(state,
                                                path_fwd_CATHY,
                                                ens_nb,
                                                **kwargs)
                        if 'pygimli' in obs2map[i]['data_format']:
                            Hx_stacked.append(Hx_ERT['rhoa'])
                        else:
                            Hx_stacked.append(Hx_ERT['resist'])
                        Hx_stacked = np.hstack(Hx_stacked)

            if len(Hx_stacked)>2:
                np.hstack(Hx_stacked)
            Hx_ens.append(Hx_stacked)

        if len(Hx_ens)>2:
            Hx_ens= np.hstack(Hx_ens)


        # special case of ERT // during sequential assimilation
        # ---------------------------------------------------------------------
        for i, obs_key in enumerate(obskey2map):
            if 'ERT' in obs_key:
                if parallel:
                    Hx_ens_ERT = self._map_ERT_parallel(path_fwd_CATHY_list,
                                                        savefig = True,
                                                        DA_cnb = self.count_DA_cycle,
                                                        )
                    if len(Hx_ens)>0:
                        Hx_ens = np.vstack([Hx_ens,Hx_ens_ERT])
                        # np.shape(Hx_ens)
                    else:
                        Hx_ens = Hx_ens_ERT

        return Hx_ens #meas_size * ens_size


    def _add_2_ensemble_Archie(self, df_Archie_2add):
        '''
        Store in a dataframe Archie relationship for all ensembles and all assimilation times

        Parameters
        ----------
        df_Archie_2add : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''

        if hasattr(self,'Archie') == False:
            self.Archie =  pd.DataFrame(columns=['time',
                                                 'ens_nb',
                                                 'sw',
                                                 'ER_converted',
                                                 'OL'])

        # print(self.Archie['time'].max())
        # print(self.Archie['ens_nb'].max())
        self.Archie = pd.concat([self.Archie,df_Archie_2add])
        # print(self.Archie['time'].max())
        # print(self.Archie['ens_nb'].max())




    def _add_2_ensemble_Hx(self, Hx, Hx_2add):
        '''
        Store in an array predicted value for all ensembles and all assimilation times
        '''
        try:
            Hx.append(Hx_2add)
        except:
            Hx = list(Hx)
            Hx.append(Hx_2add)


        return Hx


    def _read_state_ensemble(self):
        # read grid to infer dimension of the ensemble X
        # --------------------------------------------------------------------
        # self.grid3d = {}
        if len(self.grid3d) == 0:
            self.grid3d = in_CT.read_grid3d(
                os.path.join(self.workdir, self.project_name))

        M_rows = np.shape(self.grid3d['nodes_idxyz'])[0]
        N_col = self.NENS
        # N_col = len(ens_valid)

        # Ensemble matrix X of M rows and N  + fill with psi values
        # --------------------------------------------------------------------
        psi = np.zeros([M_rows, N_col])*1e99
        sw = np.zeros([M_rows, N_col])*1e99
        for j in range(self.NENS):
        # for j in range(len(ens_valid)):

            try:
                df_psi = out_CT.read_psi(os.path.join(self.workdir, self.project_name,
                                                   "DA_Ensemble/cathy_" + str(j+1),
                                                   'output/psi'))
                psi[:, j] = df_psi[-1, :]
            except:
                pass

            try:
                df_sw = out_CT.read_sw(os.path.join(self.workdir, self.project_name,
                                                   "DA_Ensemble/cathy_" + str(j+1),
                                                   'output/sw'))
                sw[:, j] = df_sw[-1, :]
            except:
                pass
        # check if there is still zeros
        if np.count_nonzero(psi==1e99) != 0:
            print('!!!unconsistent filled X, missing values!!!')
            # sys.exit()

        return psi, sw, N_col, M_rows










    def _update_input_ensemble(self, NENS, ENS_times,
                               parm_pert=[],
                               update_parm_list=['St. var.'],
                               analysis=[],
                               **kwargs):
        '''
        THIS SHOULD BE MOVED TO DA CLASS


        Select the time window of the hietograph.
        Write new variable perturbated/updated value into corresponding input file

        Parameters
        ----------
        NENS : int
            size of the ensemble.
        parm_pert : dict
            dict of all perturbated variables holding values and metadata.

        Returns
        -------
        New file written/overwritten into the input dir.

        '''

        self.console.print(":sponge: [b]update input ensemble[/b]")



        # resample atmbc file for the given DA window
        # ---------------------------------------------------------------------
        self.selec_atmbc_window(NENS, ENS_times)

        # ---------------------------------------------------------------------

        # analysis = []
        # if 'analysis' in kwargs:
        #     analysis = kwargs['analysis']

        self.update_ENS_files(parm_pert,
                              update_parm_list=update_parm_list,
                              cycle_nb=self.count_atmbc_cycle,
                              analysis=analysis)

        pass

    def selec_atmbc_window(self, NENS, ENS_times):
        '''
        THIS SHOULD BE MOVED TO DA CLASS

        Select the time window of the hietograph
        == time between two assimilation observation

        Parameters
        ----------
        NENS : TYPE
            DESCRIPTION.
        ENS_times : TYPE
            DESCRIPTION.
        '''

        if len(self.grid3d) == 0:
            self.run_processor(IPRT1=3, verbose=False, DAFLAG=0)
            self.grid3d = in_CT.read_grid3d(
                os.path.join(self.workdir, self.project_name))

        # read full simulation atmbc and filter time window
        # ----------------------------------------------------------------------
        df_atmbc, HSPATM, IETO = in_CT.read_atmbc(os.path.join(
                                                                self.workdir,
                                                                self.project_name,
                                                                'input', 'atmbc'
                                                               ),
                                                  grid=self.grid3d
                                                  )

        if self.count_atmbc_cycle is not None:
            try:
                time_window_atmbc = [
                                        df_atmbc.iloc[self.count_atmbc_cycle]['time'],
                                        df_atmbc.iloc[self.count_atmbc_cycle+1]['time']
                                    ]
            except:
                pass
        else:
                time_window_atmbc = [ENS_times[0], ENS_times[1]]

        df_atmbc_window = df_atmbc[(df_atmbc['time'] >= time_window_atmbc[0]) &
                      (df_atmbc['time'] <= time_window_atmbc[1])]


        for ens_nb in range(NENS):
            os.chdir(os.path.join(self.workdir, self.project_name,
                                  './DA_Ensemble/cathy_' + str(ens_nb+1)))
            diff_time = time_window_atmbc[1] - time_window_atmbc[0]
            self.update_parm(TIMPRTi=[0,diff_time],TMAX=diff_time,
                              filename=os.path.join(os.getcwd(), 'input/parm'),
                              backup=True)
            self.console.print(":warning: [b]Making the assumption that atmbc are homogeneous![/b]")
            VALUE = []
            for t in df_atmbc_window['time'].unique():
                # VALUE.append(df_atmbc_window[df_atmbc_window['time']==t]['value'].mean())
                VALUE.append(df_atmbc_window[df_atmbc_window['time']==t]['value'])
            if len(VALUE)>0:
                self.update_atmbc(HSPATM=1,IETO=0,
                                  time=[0,diff_time],
                                  VALUE=[VALUE[0]],
                                  filename=os.path.join(os.getcwd(), 'input/atmbc'))
            else:
                pass

        pass










    def _DA_df(self, state=[None,None], state_analysis=None, rejected_ens=[], **kwargs):
        '''
        THIS SHOULD BE MOVED TO DA CLASS

        Build a dataframe container to keep track of the changes before and after update for each ensemble.
        This dataframe is the support for DA results visualisation.


        Parameters
        ----------
        sw : np.array([])
            State variable mesh nodes values (before update).
        sw_analysis : np.array([])
            State variable mesh nodes values after DA analysis.

        Returns
        -------
        None.

        '''

        df_DA_ti_ni = pd.DataFrame() # nested df for a given ensemble/given time
        df_DA_ti = pd.DataFrame() # nested df for a given time

        cols_root = ['time',
                     'Ensemble_nb',
                     'psi_bef_update',
                     'sw_bef_update_',
                     'analysis',
                     'OL',# Open loop boolean
                     'rejected']


        # Open loop kwargs
        # --------------------------
        t_ass = []
        if 't_ass' in kwargs:
            t_ass = kwargs['t_ass']

        NENS = []
        if 'ens_nb' in kwargs:
            ens_nb = kwargs['ens_nb']
        # --------------------------


        if 'openLoop' in kwargs:

            data_df_root = [t_ass*np.ones(len(state[0])),
                            ens_nb*np.ones(len(state[0])),
                            state[0],
                            state[1],
                            state[0],
                            True*np.ones(len(state[0])),
                            False*np.ones(len(state[0]))]
            df_DA_ti = pd.DataFrame(np.transpose(data_df_root),
                                  columns=cols_root)

        else:

            for n in range(self.NENS):

                data_df_root = [self.count_atmbc_cycle*np.ones(len(state[0])),
                                n*np.ones(len(state[0])),
                                state[0][:,n],
                                state[1][:,n],
                                state_analysis[:,n],
                                False*np.ones(len(state[0])),
                                int(rejected_ens[n])**np.ones(len(state[0]))
                                ]
                df_DA_ti_ni = pd.DataFrame(np.transpose(data_df_root),
                                      columns=cols_root)

                df_DA_ti= pd.concat([df_DA_ti, df_DA_ti_ni], axis=0)


        self.df_DA= pd.concat([self.df_DA, df_DA_ti], axis=0, ignore_index=True)


        pass
    def _add_to_perturbated_dict(self,var_per_2add):
        '''
        Update dict of perturbated variable.
        '''
        self.var_per_dict = (self.var_per_dict | var_per_2add)
        return self.var_per_dict


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
