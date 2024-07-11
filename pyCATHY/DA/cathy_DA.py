"""
Class managing Data Assimilation process i.e.: 
    
    1. Parameters and data perturbation
    
        This step consit in the generation of the ensemble.
    
    DA class
    ---
    Steps 2 and 3 are controlled via a specific class
    
    2.  Open Loop

       Without data run the hydrological model for an ensemble of realisations
       using the `_run_hydro_DA_openLoop` module.
       
    3. Loop of DA:
        
        - Running hydro model
        
        Run iteratively (between two observations) the hydrological model
        using the `_run_ensemble_hydrological_model` module.
        
        
        - Map states 2 Observations 
         
        
        The module `map_states2Observations()` takes care of the mapping.
        
        For complex mapping such as ERT, it is done via a separate class.
        Example for the ERT we have a `class mappingERT()` which uses the 
        Archie's petrophysical relationship.
        
        
        - Running the analysis
        
        The analysis is controlled by `run_analysis()` which takes as argument
        the type of analysis to run i.e. EnkF or Pf.
        
        - Update ensemble
        
        The ensemble file update is controlled by `update_ENS_files()`
        
        
        - Evaluate performance


"""
import glob

import multiprocessing
import os
import re
import shutil
import subprocess
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rich.progress import track

from pyCATHY.cathy_tools import CATHY, subprocess_run_multi
from pyCATHY.DA import enkf, pf, mapper
from pyCATHY.importers import cathy_inputs as in_CT
from pyCATHY.importers import cathy_outputs as out_CT
# from pyCATHY.importers import sensors_measures as in_meas
# import pyCATHY.meshtools as mt
from pyCATHY.plotters import cathy_plots as plt_CT

warnings.filterwarnings("ignore", 
                        category=DeprecationWarning
                        )

# from functools import partial

REMOVE_CPU = 1

def run_analysis(
    DA_type,
    data,
    data_cov,
    param,
    list_update_parm,
    ensembleX,
    prediction,
    default_state="psi",
    **kwargs,
):
    """
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
    """

    id_state = 0
    if default_state == "sw":
        id_state = 1
        
        
    if DA_type == "enkf_Evensen2009":
        Analysis = enkf.enkf_analysis(data, 
                               data_cov, 
                               param, 
                               ensembleX[id_state], 
                               prediction,
                               )
        return Analysis
    
    if DA_type == "enkf_Evensen2009_Sakov":
        Analysis = enkf.enkf_analysis(data, 
                               data_cov, 
                               param, 
                               ensembleX[id_state], 
                               prediction,
                               Sakov=True,
                               )
        return Analysis
    
    if DA_type == "enkf_analysis_inflation":
        Analysis = enkf.enkf_analysis_inflation(
                                                data, 
                                                data_cov, 
                                                param, 
                                                ensembleX[id_state], 
                                                prediction, 
                                                **kwargs
                                                )
        return Analysis
    
    
    if DA_type == "enkf_analysis_inflation_multiparm":
        Analysis = enkf.enkf_analysis_inflation_multiparm(
                                                        data, 
                                                        data_cov, 
                                                        param, 
                                                        ensembleX[id_state], 
                                                        prediction, 
                                                        **kwargs
                                                        )
        return Analysis


    if DA_type == "enkf_dual_analysis":
        Analysis = enkf.enkf_dual_analysis(data, 
                                           data_cov, 
                                           param, 
                                           ensembleX[id_state], 
                                           prediction,
                                           Sakov=False,
                                           )
        return Analysis
    
    if DA_type == "enkf_dual_analysis_Sakov":
        Analysis = enkf.enkf_dual_analysis(data, 
                                           data_cov, 
                                           param, 
                                           ensembleX[id_state], 
                                           prediction,
                                           Sakov=True,
                                           )
        return Analysis
    
    
    if DA_type == "pf":
        print("not yet implemented")

        [Analysis, AnalysisParam] = pf.pf_analysis(
                                                    data, 
                                                    data_cov,
                                                    param, 
                                                    ensembleX[id_state], 
                                                    prediction
                                                    )
        return Analysis, AnalysisParam


# utils observations
# ------------------


def dictObs_2pd(dict_obs):
    """dict of observation to dataframe of observation"""
    df_obs = pd.DataFrame.from_dict(dict_obs).stack().to_frame()
    df_obs = pd.DataFrame(df_obs[0].values.T.tolist(), index=df_obs.index)
    df_obs.index.names = ["sensorNameidx", "assimilation time"]
    return df_obs


def resynchronise_times(data_measure, atmbc_df, time_offset):
    """resynchronise dict as old key is elapsed time in second from the first observation,
    while new key is from the first atmbc time
    """
    data_measure_sync = dict(data_measure)
    try:
        new_keys = []
        for d in range(len(data_measure_sync.keys())):
            items_dict = list(data_measure_sync.items())
            # new_key = (
            #     atmbc_df["diff"]
            #     .dt.total_seconds()
            #     .unique()[d]
            # )

            old_key = list(data_measure_sync.keys())[d]
            new_key = list(data_measure_sync.keys())[d] - time_offset  

            for sensor in list(items_dict[d][1].keys()):
                data_measure_sync[old_key][sensor]["assimilation_times"] = new_key
            new_keys.append(new_key)
        data_measure_sync = dict(zip(new_keys, data_measure_sync.values()))

    except:
        print("datetime first atmbc point")
        print(atmbc_df["datetime"][0])
        print("datetime first measurement")
        print(data_measure[0][sensor]["datetime"])
        print("cant synchronise times - continue without")
    return data_measure_sync


#%%
# DA class
# ---------

class DA(CATHY):

    # def __init__(self,**kwargs):
    #     self.damping = 1
    # -------------------------------------------------------------------#
    # %% DATA ASSIMILATION FCTS

    def _DA_init(self, dict_obs, dict_parm_pert, list_parm2update="all"):
        """
        Initial preparation for DA


        .. note::

            0. check dictionnaries validity

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

        # Initiate
        # -------------------------------------------------------------------
        update_key = "ini_perturbation"

        # check if dict_obs is ordered in time
        # ------------------------------------
        if (
            all(i < j for i, j in zip(list(dict_obs.keys()), list(dict_obs.keys())[1:]))
        ) is False:
            raise ValueError("Observation List is not sorted.")

        # if hasattr(self,'dict_obs') is False:
        self.dict_obs = dict_obs  # self.dict_obs is already assigned (in read observatio! change it to self.obs
        self.dict_parm_pert = dict_parm_pert
        self.df_DA = pd.DataFrame()
        self.df_Archie = pd.DataFrame()
        
        # backup all spatial ET iterations (self.backupfort777)
        # ET file are large!
        self.backupfort777 = False 
        
        # Infer ensemble size NENS from perturbated parameter dictionnary
        # -------------------------------------------------------------------
        for name in self.dict_parm_pert:
            self.NENS = len(self.dict_parm_pert[name]["ini_perturbation"])

        # Infer ensemble update times ENS_times from observation dictionnary
        # -------------------------------------------------------------------
        ENS_times = []
        for ti in self.dict_obs:
            ENS_times.append(float(ti))

        # start DA cycle counter
        # -------------------------------------------------------------------
        self.count_DA_cycle = 0
        self.count_atmbc_cycle = 0
        # (the counter is incremented during the update analysis)

        if hasattr(self, "ens_valid") is False:
            self.ens_valid = list(np.arange(0, self.NENS))

        # create sub directories for each ensemble
        # ---------------------------------------------------------------------
        self._create_subfolders_ensemble()

        # update list of updated parameters based on problem heterogeneity
        # ---------------------------------------------------------------------
        newlist_parm2update = ["St. var."]
        for d in self.dict_parm_pert.keys():
            for l in list_parm2update:
                if l in d:
                    newlist_parm2update.append(d)

        return ENS_times, newlist_parm2update, update_key

    def run_DA_sequential(
                            self,
                            DA_type,
                            dict_obs,
                            list_parm2update,
                            dict_parm_pert,
                            list_assimilated_obs,
                            open_loop_run=False,
                            parallel=True,
                            **kwargs,
                        ):

        """

        Run Data Assimilation

        .. note::

            Steps are:
            1. DA init (create subfolders)
            2a. run CATHY hydrological model (open loop)
            2b. run CATHY hydrological model recursively using DA times
                <- Loop-->
                    3. check before analysis
                    4. analysis
                    5. check after analysis
                    6. update input files
                <- Loop-->

        Parameters
        ----------
        callexe : str
            processor exe filename
        parallel : bool
            True for multiple ensemble folder simulations.
        DA_type : str
            type of Data Assimilation.
        list_assimilated_obs : TYPE
            list of names of observations to assimilate.
        dict_obs : dict
            observations dictionnary.
        list_update_parm : list
            list of names of parameters to update.
        dict_parm_pert : dict
            parameters perturbated dictionnary.


        """
        self.run_processor(DAFLAG=1,
                           runProcess=False,
                           **kwargs) # to create the processor exe
        callexe = "./" + self.processor_name

        self.damping = 1
        if "damping" in kwargs:
            self.damping = kwargs["damping"]

        threshold_rejected = 10
        if "threshold_rejected" in kwargs:
            threshold_rejected = kwargs["threshold_rejected"]

        self.verbose = False
        if "verbose" in kwargs:
            self.verbose = kwargs["verbose"]

        # initiate DA
        # -------------------------------------------------------------------
        ENS_times, list_parm2update, update_key = self._DA_init(
            dict_obs,
            dict_parm_pert,
            list_parm2update,
        )

        # initiate mapping petro
        # -------------------------------------------------------------------
        self._mapping_petro_init()

        # update the ensemble files with pertubarted parameters
        # -------------------------------------------------------------------
        self.update_ENS_files(
            dict_parm_pert,
            list_parm2update="all",  
            cycle_nb=self.count_DA_cycle,
        )

        all_atmbc_times = np.copy(self.atmbc["time"])

        # Open loop
        # -------------------------------------------------------------------
        if open_loop_run:
            self._DA_openLoop(ENS_times, list_assimilated_obs, parallel)

        # end of Open loop
        # -------------------------------------------------------------------

        # -------------------------------------------------------------------
        # update input files ensemble again (time-windowed)
        # -------------------------------------------------------------------
        self._update_input_ensemble(
            list(self.atmbc["time"]),
            list_parm2update=list_parm2update,         
        )

        # -----------------------------------
        # Run hydrological model sequentially
        # = Loop over atmbc times (including assimilation observation times)
        # -----------------------------------
        for t_atmbc in all_atmbc_times:  # atmbc times MUST include assimilation observation times

            self._run_ensemble_hydrological_model(parallel, callexe)
            os.chdir(os.path.join(self.workdir))

            # self.plot_ini_state_cov()
            # def plot_ini_state_cov():
            #     ensemble_psi, ensemble_sw, ens_size, sim_size = self._read_state_ensemble()

            # check scenario (result of the hydro simulation)
            # ----------------------------------------------------------------
            rejected_ens = self._check_before_analysis(
                self.ens_valid,
                threshold_rejected,
            )
            id_valid = ~np.array(rejected_ens)
            self.ens_valid = list(np.arange(0, self.NENS)[id_valid])
            ensemble_psi_valid_bef_update = []
            #!! define subloop here to optimize the inflation!!
            # if 'optimize_inflation' in DA_type:
            # threshold_convergence = 10
            # while self.df_performance < threshold_convergence:

            # check if there is an observation at the given atmbc time
            # --------------------------------------------------------
            if t_atmbc in ENS_times:

                #%%
                # map states to observation = apply H operator to state variable
                # ----------------------------------------------------------------
                prediction = self.map_states2Observations(
                                                            list_assimilated_obs,
                                                            default_state="psi",
                                                            parallel=parallel,
                                                        )

                # print(['Im done with mapping']*10)
                #%%
                # DA analysis
                # ----------------------------------------------------------------
                self.console.print(
                    r""":chart_with_upwards_trend: [b]Analysis[/b]:
                                       - DA type: {}
                                       - Damping: {}
                                   """.format(
                        DA_type, self.damping
                    )
                )
    
    
                (
                    ensemble_psi_valid_bef_update,
                    ensemble_sw_valid_bef_update,
                    ens_size,
                    sim_size,
                ) = self._read_state_ensemble()
                
                
                (
                    ensemble_psi,
                    ensemble_sw,
                    data,
                    analysis,
                    analysis_param,
                ) = self._DA_analysis(
                                    prediction,
                                    DA_type,
                                    list_parm2update,
                                    list_assimilated_obs,
                                    ens_valid=self.ens_valid,
                )
                # print(['Im done with _DA_analysis']*10)

                # DA mark_invalid_ensemble
                # ----------------------------------------------------------------
                (
                    prediction_valid,
                    ensemble_psi_valid,
                    ensemble_sw_valid,
                    analysis_valid,
                    analysis_param_valid,
                ) = self._mark_invalid_ensemble(
                                                self.ens_valid,
                                                prediction,
                                                ensemble_psi,
                                                ensemble_sw,
                                                analysis,
                                                analysis_param,
                )
                # print(['Im done with _invalid_ensemble']*10)

                #%%
                # check analysis quality
                # ----------------------------------------------------------------

                self.console.print(
                    ":face_with_monocle: [b]check analysis performance[/b]"
                )
                
                self._performance_assessement(
                                                list_assimilated_obs,
                                                data,
                                                prediction_valid,
                                                t_obs=self.count_DA_cycle,
                                            )
                # print(['Im done with _analysis performance']*10)

                #%%
                # the counter is incremented here
                # ----------------------------------------------------------------
                self.count_DA_cycle = self.count_DA_cycle + 1

                # # update parameter dictionnary
                # ----------------------------------------------------------------
                def check_ensemble(ensemble_psi_valid, ensemble_sw_valid):
                    if np.any(ensemble_psi_valid > 0):
                        print("!positive pressure heads observed!")

                check_ensemble(ensemble_psi_valid, ensemble_sw_valid)

                self.update_pert_parm_dict(
                                                update_key,
                                                list_parm2update, 
                                                analysis_param_valid
                                                )

            else:
                self.console.print(
                    ":confused: No observation for this time - run hydrological model only"
                )
  
                (
                    ensemble_psi_valid,
                    ensemble_sw_valid,
                    ens_size,
                    sim_size,
                ) = self._read_state_ensemble()
                
                print("!shortcut here ensemble are not validated!")

                analysis_valid = ensemble_psi_valid
                
                # DA mark_invalid_ensemble
                # ----------------------------------------------------------------
                # (
                #     prediction_valid,
                #     ensemble_psi_valid,
                #     ensemble_sw_valid,
                #     analysis_valid,
                #     analysis_param_valid,
                # ) = self._mark_invalid_ensemble(
                #                                 self.ens_valid,
                #                                 prediction,
                #                                 ensemble_psi,
                #                                 ensemble_sw,
                #                                 analysis,
                #                                 analysis_param,
                # )
                    
                print("!shortcut here ensemble are not validated!")
                # print(f'''max, min, mean psi: 
                #       {np.max(analysis_valid)}, 
                #       {np.min(analysis_valid)}, 
                #       {np.mean(analysis_valid)}
                #       ''')

            self.count_atmbc_cycle = self.count_atmbc_cycle + 1

            self.console.rule(
                ":round_pushpin: end of time step (s) "
                + str(int(t_atmbc))
                + "/"
                + str(int(all_atmbc_times[-1]))
                + " :round_pushpin:",
                style="yellow",
            )
            self.console.rule(
                ":round_pushpin: end of atmbc update # "
                + str(self.count_atmbc_cycle)
                + "/"
                + str(len(all_atmbc_times))
                + " :round_pushpin:",
                style="yellow",
            )

            # create dataframe _DA_var_pert_df holding the results of the DA update
            # ---------------------------------------------------------------------
            # self.DA_df
            
            if len(ensemble_psi_valid_bef_update)==0:
                (
                    ensemble_psi_valid_bef_update,
                    ensemble_sw_valid_bef_update,
                    ens_size,
                    sim_size,
                ) = self._read_state_ensemble()
            else:
                self._DA_df(
                    state=[ensemble_psi_valid_bef_update, 
                           ensemble_sw_valid_bef_update
                           ],
                    state_analysis=analysis_valid,
                    rejected_ens=rejected_ens,
                )
            # self.df_DA
            
            # print(['Im done with dataframe _DA_var_pert_df']*10)


            # export summary results of DA
            # ----------------------------------------------------------------
            meta_DA = {
                "listAssimilatedObs": list_assimilated_obs,
                "listUpdatedparm": list_parm2update,
                # '':,
                # '':,
                # '':,
            }
             
            # backup ET map for each iteration of ensemble numbers
            if self.backupfort777: 
                self.backup_fort777()
            self.backup_results_DA(meta_DA)
            self.backup_simu()
            

            # overwrite input files ensemble (perturbated variables)
            # ---------------------------------------------------------------------
            if (
                self.count_atmbc_cycle < len(all_atmbc_times) - 1
            ):  # -1 cause all_atmbc_times include TMAX
                # len(all_atmbc_times)
                self._update_input_ensemble(
                                                all_atmbc_times, 
                                                list_parm2update, 
                                                analysis=analysis_valid
                )
            else:
                self.console.rule(
                    ":red_circle: end of DA" ":red_circle:", style="yellow"
                )
                pass
            self.console.rule(
                                ":red_circle: end of DA update "
                                + str(self.count_DA_cycle)
                                + "/"
                                + str(len(ENS_times))
                                + " :red_circle:",
                                style="yellow",
                            )
            self.console.rule(
                                ":dart: % of valid ensemble: "
                                + str((len(self.ens_valid) * 100) / (self.NENS))
                                + ":dart:",
                                style="yellow",
                            )
            plt.close("all")

        # clean all vtk file produced
        # --------------------------
        directory = "./"
        pathname = directory + "/**/*converted*.vtk"
        files = glob.glob(pathname, recursive=True)


        remov_convert = False
        if remov_convert:
            [os.remove(f) for f in files]

        hard_remove = False
        if hard_remove:
            directory = "./"
            pathname = directory + "/**/*.vtk"
            files = glob.glob(pathname, recursive=True)
            [os.remove(f) for f in files]
        pass


    # backup ET map file at each DA iteration
    def backup_fort777(self):
       
        for ensi in range(self.NENS):
            dst_dir = os.path.join(
                                    self.workdir,
                                    self.project_name,
                                    'DA_Ensemble/cathy_'+ str(ensi+1),
                                    'DAnb' + str(self.count_DA_cycle) + 'fort.777' 
                                    )
    
            shutil.copy(os.path.join(self.workdir,
                                     self.project_name,
                                     'DA_Ensemble/cathy_'+ str(ensi+1),
                                     'fort.777'),
                        dst_dir)
        
    def _run_hydro_DA_openLoop(
        self, time_of_interest, nodes_of_interest, simu_time_max, ens_nb
    ):
        """
        Run openLoop and save result into DA dataframe (ensemble nb = 999)

        Parameters
        ----------
        time_of_interest : TYPE
            DESCRIPTION.
        nodes_of_interest : TYPE
            DESCRIPTION.
        simu_time_max : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """

        cwd = os.getcwd()

        self.console.print(":unlock: [b]open loop call[/b]")
        self.run_processor(
            recompile=True,
            TIMPRTi=time_of_interest,
            NODVP=nodes_of_interest,
            TMAX=simu_time_max,
            DAFLAG=0,
            path_CATHY_folder=cwd,
        )

    def _DA_openLoop(self, ENS_times, list_assimilated_obs, parallel, verbose=False):
        """
        Run open Loop (no update/assimilation) hydro simulation for an ensemble of realisation
        Evaluate the performance by comparison with measured data after mapping

        Parameters
        ----------
        ENS_times : list
            List of the time steps (in sec) where observations are assimilated.
        list_assimilated_obs : list
            List of assimilation observations.
        parallel : Bool, optional
            Parallelize. The default is True.
        verbose : Bool, optional
            Display prints. The default is False.

        Returns
        -------
        Prediction_OL.
        Performance dataframe.

        """

        # multi run CATHY hydrological model from the independant folders
        # composing the ensemble
        # ----------------------------------------------------------------
        if parallel == True:
            self.console.print(
                ":athletic_shoe: [b]Run hydrological model // :unlock: [/b]"
            )

            pathexe_list = []
            for ens_i in range(self.NENS):
                path_exe = os.path.join(
                    self.workdir,
                    self.project_name,
                    "DA_Ensemble/cathy_" + str(ens_i + 1),
                )
                pathexe_list.append(path_exe)
            with multiprocessing.Pool(processes=multiprocessing.cpu_count()-REMOVE_CPU) as pool:
                result = pool.map(subprocess_run_multi, pathexe_list)
                if verbose == True:
                    self.console.print(result)

        # Loop over ensemble realisations one by one
        # ------------------------------------------------------------
        else:
            for ens_i in range(self.NENS):
                os.chdir(
                    os.path.join(
                        self.workdir,
                        self.project_name,
                        "DA_Ensemble/cathy_" + str(ens_i + 1),
                    )
                )
                # len(list(self.atmbc['time']))
                self._run_hydro_DA_openLoop(
                    time_of_interest=list(self.atmbc["time"]),
                    nodes_of_interest=[],
                    simu_time_max=max(list(self.atmbc["time"])),
                    ens_nb=ens_i + 1,
                )

        # save into the DA_df dataframe
        # -----------------------------------------------------------------
        path_fwd_CATHY_list = []
        for ens_i in range(self.NENS):
            path_fwd_CATHY_list.append(
                os.path.join(
                    self.workdir,
                    self.project_name,
                    "DA_Ensemble/cathy_" + str(ens_i + 1),
                )
            )

            os.chdir(
                os.path.join(
                    self.workdir,
                    self.project_name,
                    "DA_Ensemble/cathy_" + str(ens_i + 1),
                )
            )
            df_psi = self.read_outputs(
                filename="psi", path=os.path.join(os.getcwd(), "output")
            )
            df_sw, _ = self.read_outputs(
                filename="sw", path=os.path.join(os.getcwd(), "output")
            )

            shift = len(df_psi) - self.parm["NPRT"]
            if shift < 0 | shift > 2:
                print("Error for the ensemble nb:" + str(ens_i))
                raise ValueError(
                    "Error on the simulation:"
                    "nb of times contained in the outputs files is too small;"
                    "Check "
                )
            for t in range(np.shape(df_psi)[0] - 2):
                self._DA_df(
                    state=[df_psi[t + shift, :], df_sw.iloc[t + shift, :]],
                    t_ass=t,
                    openLoop=True,
                    ens_nb=ens_i + 1,
                )

        os.chdir(
            os.path.join(
                self.workdir,
            )
        )
        prediction_OL = self._evaluate_perf_OL(
            parallel, list_assimilated_obs, path_fwd_CATHY_list, ENS_times
        )

        # ------------------------------------------------------
        # END of Open Loop simulation and performance evaluation
        # ------------------------------------------------------
        pass

    def _DA_analysis(
        self,
        prediction,
        DA_type="enkf_Evensen2009_Sakov",
        list_update_parm=["St. var."],
        list_assimilated_obs="all",
        ens_valid=[],
    ):
        """
        Analysis ensemble using DA

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

        update_key = "ini_perturbation"
        if self.count_DA_cycle > 0:
            update_key = "update_nb" + str(self.count_DA_cycle)

        # find the method to map for this time step
        # --------------------------------------------------------
        #%%
        obskey2map, obs2map = self._obs_key_select(list_assimilated_obs)
        
        #%%

        # prepare states
        # ---------------------------------------------------------------------
        ensemble_psi, ensemble_sw, ens_size, sim_size = self._read_state_ensemble()

        # ---------------------------------------------------------------------
        # prepare data
        # ---------------------------------------------------------------------
        # select ONLY data to assimilate
        # ---------------------------------------------------------------------
        data, _ = self._get_data2assimilate(list_assimilated_obs)
        self.console.print(
            r""":
                               - Data size: {}
                                 --> Observations --> {}
                           """.format(
                str(np.shape(data)), list_assimilated_obs
            )
        )

        # check size of data_cov
        # ---------------------------------------------------------------------
        if np.shape(self.stacked_data_cov[self.count_DA_cycle])[0] != len(data):
            raise ValueError("Wrong stacked_data_cov shape -- need to (re)compute data covariance")

        data_cov = self.stacked_data_cov[self.count_DA_cycle]

        # # filter non-valid ensemble before Analysis
        # # ---------------------------------------------------------------------
        prediction_valid = prediction
        ensemble_psi_valid = ensemble_psi[:, ens_valid]
        ensemble_sw_valid = ensemble_sw[:, ens_valid]

        # prepare parameters
        # ---------------------------------------------------------------------
        # When updating the states only, the elements of X are the
        # pressure heads at each node of the finite element grid, while
        # the state augmentation technique is used when also updat-
        # ing the parameters
        param = self.transform_parameters(list_update_parm, update_key)
        print("parm size: " + str(len(param)))

        param_valid = []
        if len(param) > 0:
            param_valid = param[:, ens_valid]

        # run Analysis
        # ---------------------------------------------------------------------
        result_analysis = run_analysis(
                                        DA_type,
                                        data,
                                        data_cov,
                                        param_valid,
                                        list_update_parm,
                                        [ensemble_psi_valid, ensemble_sw_valid],
                                        prediction_valid,
                                        alpha=self.damping,
                                    )
        
        _ , obs2eval_key = self._get_data2assimilate(list_assimilated_obs)

        # fig, ax = plt.subplots(figsize=(4,4))
        # # for i in range(len(prediction_valid.T)):
        # ax.scatter(data,np.mean(prediction_valid,axis=1),label='mean')
        # ax.scatter(data,np.median(prediction_valid,axis=1),label='median')
        # ax.scatter(data,np.min(prediction_valid,axis=1),color='k',alpha=0.1)
        # ax.scatter(data,np.max(prediction_valid,axis=1),color='k',alpha=0.1)
        # plt.legend()
        # # Add 1:1 line
        # lims = [min(min(data), min(np.mean(prediction_valid, axis=1))),
        #         max(max(data), max(np.mean(prediction_valid, axis=1)))]
        # ax.plot(lims, lims, 'k--', alpha=0.75, zorder=0)
        
        # ax.set_xlabel(f' {obs2eval_key[0]} Observations')
        # ax.set_ylabel(f' {obs2eval_key[0]} Predictions')
        # ax.set_xlim([min(data), max(data)])
        # ax.set_ylim([min(data), max(data)])
        # ax.set_aspect('equal')

        # fig.savefig(os.path.join(self.workdir,
        #                          self.project_name,
        #                          'T' + str(self.count_DA_cycle) +
        #                          '_11line.png'
        #                          )
        #             )
        
        # plot ensemble covariance matrices and changes (only for ENKF)
        # ---------------------------------------------------------------------
        if len(result_analysis) > 2:
            [
                A,
                Amean,
                dA,
                dD,
                MeasAvg,
                S,
                COV,
                B,
                dAS,
                analysis,
                analysis_param,
            ] = result_analysis

            try:
                plt_CT.show_DA_process_ens(
                    ensemble_psi_valid,
                    data,
                    COV,
                    dD,
                    dAS,
                    B,
                    analysis,
                    savefig=True,
                    savename=os.path.join(
                        self.workdir,
                        self.project_name,
                        "DA_Ensemble",
                        "DA_Matrices_t" + str(self.count_DA_cycle),
                    ),
                    # label_sensor=str([*self.dict_obs[0]])
                    label_sensor=str([*obskey2map]),
                )
            except:
                self.console.print("[b]Impossible to plot matrices[/b]")

        # particule filter
        # ----------------
        else:
            [analysis, analysis_param] = result_analysis

        # Back-transformation of the parameters
        # ---------------------------------------------------------------------
        # for pp in enumerate(list_update_parm[:]):
        #     if 'St. var.' in pp[1]:
        #         pass
        #     else:

        # analysis_param_valid = []
        # if len(param)>0:
        #     analysis_param_valid = analysis_param[ens_valid]

        self.transform_parameters(
            list_update_parm, param=np.hstack(analysis_param), back=True
        )
        # if self.dict_parm_pert[pp[1]]['transf_type'] == 'log':
        #     print('back log transformation')
        #     analysis_param = np.exp(analysis_param)

        # return ensemble_psi, Analysis, AnalysisParam, Data, data_cov
        # ---------------------------------------------------------------------
        return (ensemble_psi, ensemble_sw, data, analysis, analysis_param)

    def _mark_invalid_ensemble(
        self, ens_valid, prediction, ensemble_psi, ensemble_sw, analysis, analysis_param
    ):
        """mark invalid ensemble - invalid ensemble are filled with NaN values"""

        ensemble_psi_valid = np.empty(ensemble_psi.shape)
        ensemble_psi_valid[:] = np.NaN
        ensemble_psi_valid[:, ens_valid] = ensemble_psi[:, ens_valid]

        ensemble_sw_valid = np.empty(ensemble_psi.shape)
        ensemble_sw_valid[:] = np.NaN
        ensemble_sw_valid[:, ens_valid] = ensemble_sw[:, ens_valid]

        analysis_valid = np.empty(ensemble_psi.shape)
        analysis_valid[:] = np.NaN
        analysis_valid[:, ens_valid] = analysis

        prediction_valid = np.empty([prediction.shape[0], ensemble_psi.shape[1]])
        prediction_valid[:] = np.NaN
        prediction_valid[:, ens_valid] = prediction

        analysis_param_valid = []
        if len(analysis_param[0]) > 0:
            analysis_param_valid = np.empty(
                [analysis_param.shape[1], ensemble_psi.shape[1]]
            )
            analysis_param_valid[:] = np.NaN
            analysis_param_valid[:, ens_valid] = analysis_param.T

        return (
            prediction_valid,
            ensemble_psi_valid,
            ensemble_sw_valid,
            analysis_valid,
            analysis_param_valid,
        )

    def _run_ensemble_hydrological_model(self, parallel, callexe):
        """multi run CATHY hydrological model from the independant folders composing the ensemble"""

        self.console.print(":athletic_shoe: [b]Run hydrological model[/b]")
        # ----------------------------------------------------------------
        if parallel == True:
            pathexe_list = []
            # for ens_i in range(self.NENS):

            for ens_i in self.ens_valid:
                path_exe = os.path.join(
                    self.workdir,
                    self.project_name,
                    "DA_Ensemble/cathy_" + str(ens_i + 1),
                )
                pathexe_list.append(path_exe)
            with multiprocessing.Pool(processes=multiprocessing.cpu_count()-REMOVE_CPU) as pool:
                result = pool.map(subprocess_run_multi, pathexe_list)

                if self.verbose:
                    self.console.print(result)
        # process each ensemble folder one by one
        # ----------------------------------------------------------------
        else:
            # Loop over ensemble realisations
            # ------------------------------------------------------------
            for ens_i in track(
                range(self.NENS), description="Run hydrological fwd model..."
            ):

                self.console.print(
                    ":keycap_number_sign: [b]ensemble nb:[/b]"
                    + str(ens_i + 1)
                    + "/"
                    + str(self.NENS)
                )
                os.chdir(
                    os.path.join(
                        self.workdir,
                        self.project_name,
                        "DA_Ensemble/cathy_" + str(ens_i + 1),
                    )
                )
                p = subprocess.run([callexe], text=True, capture_output=True)

                if self.verbose:
                    self.console.print(p.stdout)
                    self.console.print(p.stderr)

        pass

    def _check_before_analysis(self, ens_valid=[], threshold_rejected=10):
        """
        Filter is applied only on selected ensemble

        Parameters
        ----------
        update_key : TYPE
            DESCRIPTION.

        Returns
        -------
        rejected_ens

        """

        self.console.print(":white_check_mark: [b]check scenarii before analysis[/b]")

        ###CHECK which scenarios are OK and which are to be rejected and then create and applied the filter
        ###only on the selected scenarios which are to be considered.

        rejected_ens_new = []
        for n in range(self.NENS):

            if n in ens_valid:
                # print('ensemble_nb ' + str(n) + ' before update analysis')
                try:
                    df_mbeconv = out_CT.read_mbeconv(
                        os.path.join(
                            self.workdir,
                            self.project_name,
                            "DA_Ensemble/cathy_" + str(n + 1),
                            "output/mbeconv",
                        )
                    )
                    if df_mbeconv["CUM."].isnull().values.any():
                        rejected_ens_new.append(True)
                        self.console.print(
                            ":x: [b]df_mbeconv['CUM.'][/b]"
                            + str(df_mbeconv["CUM."].isnull().values.any())
                            + ", ens_nb:"
                            + str(n + 1)
                        )

                    elif (np.round(df_mbeconv["TIME"].iloc[-1])) < np.round(
                        self.parm["(TIMPRT(I),I=1,NPRT)"][-1]
                    ) - self.parm["DELTAT"]:
                        rejected_ens_new.append(True)
                        self.console.print(
                            ":x: [b]df_mbeconv['time'][/b]"
                            + str(df_mbeconv["TIME"].iloc[-1])
                            + ", ens_nb:"
                            + str(n + 1)
                        )

                    elif len(df_mbeconv) == 0:
                        rejected_ens_new.append(True)
                        self.console.print(
                            ":x: [b]len(df_mbeconv)['TIME'][/b]"
                            + str(len(df_mbeconv))
                            + ", ens_nb:"
                            + str(n + 1)
                        )

                    else:
                        rejected_ens_new.append(False)

                except:
                    rejected_ens_new.append(True)
                    print("cannot read mbeconv")
                    pass
                    # UnboundLocalError: local variable 'df_mbeconv' referenced before assignment

            else:
                rejected_ens_new.append(True)

        # check whether the number of discarded scenarios is major than the 10% of the total number [N];
        # if the answer is yes than the execution stops
        # --------------------------------------------------------------------
        if sum(rejected_ens_new) >= (threshold_rejected / 100) * self.NENS:
            print("new parameters (if updated)")
            print(self.dict_parm_pert)

            raise ValueError(
                "% number of rejected ensemble is too high:"
                + str((sum(rejected_ens_new) * 100) / (self.NENS))
            )

        return rejected_ens_new

    def _check_after_analysis(self, update_key, list_update_parm):
        """
        Check which scenarios parameters are OK and which one to discard after analysis.
        Test the range of each physical properties.

        Parameters
        ----------
        update_key : str
            Update cycle key identifier.
        list_update_parm : list
            list of updated model parameters (after analysis).
        """

        self.console.print(":white_check_mark: [b]check scenarii post update[/b]")

        id_valid = [list(np.arange(0, self.NENS))]
        id_nonvalid = []

        def test_negative_values(update_key, pp):
            id_nonvalid.append(
                list(np.where(self.dict_parm_pert[pp[1]][update_key] < 0)[0])
            )
            id_valid.append(
                list(np.where(self.dict_parm_pert[pp[1]][update_key] > 0)[0])
            )

            for i in id_nonvalid:
                self.console.print(
                    ":x: [b]negative"
                    + str(pp)
                    + ":[/b]"
                    + str(self.dict_parm_pert[pp[1]][update_key][i])
                    + ", ens_nb:"
                    + str(i)
                )
            return id_nonvalid, id_valid

        def test_range_values(update_key, pp, min_r, max_r):
            id_nonvalid.append(
                list(np.where(self.dict_parm_pert[pp[1]][update_key] < min_r)[0])
            )
            id_nonvalid.append(
                list(np.where(self.dict_parm_pert[pp[1]][update_key] > max_r)[0])
            )
            id_valid.append(
                list(np.where(self.dict_parm_pert[pp[1]][update_key] > min_r)[0])
            )
            id_valid.append(
                list(np.where(self.dict_parm_pert[pp[1]][update_key] < max_r)[0])
            )
            for i in id_nonvalid:
                self.console.print(
                    ":x: [b]out of range"
                    + str(pp)
                    + ":[/b]"
                    + str(self.dict_parm_pert[pp[1]][update_key][i])
                    + ", ens_nb:"
                    + str(i)
                )
            return id_nonvalid, id_valid

        for pp in enumerate(list_update_parm[:]):
            if "St. var." in pp[1]:
                pass
            else:
                if "rFluid".casefold() in pp[1].casefold():
                    id_nonvalid, id_valid = test_negative_values(update_key, pp)
                elif "porosity".casefold() in pp[1].casefold():
                    id_nonvalid, id_valid = test_negative_values(update_key, pp)
                elif "a_Archie".casefold() in pp[1].casefold():
                    id_nonvalid, id_valid = test_range_values(update_key, pp, 0, 3)
                elif "Ks".casefold() in pp[1].casefold():

                    id_nonvalid.append(
                        list(np.where(self.dict_parm_pert[pp[1]][update_key] < 0)[0])
                    )
                    id_valid.append(
                        list(np.where(self.dict_parm_pert[pp[1]][update_key] > 0)[0])
                    )
                    for i in id_nonvalid:
                        self.console.print(
                            ":x: [b]negative Ks:[/b]"
                            + str(self.dict_parm_pert[pp[1]][update_key][i])
                            + ", ens_nb:"
                            + str(i)
                        )
                elif "PCREF".casefold() in pp[1].casefold():
                    id_nonvalid.append(
                        list(np.where(self.dict_parm_pert[pp[1]][update_key] > 0)[0])
                    )
                    id_valid.append(
                        list(np.where(self.dict_parm_pert[pp[1]][update_key] < 0)[0])
                    )
                    for i in id_nonvalid:
                        self.console.print(
                            "Positive values encountered in PCREF:"
                            + str(self.dict_parm_pert[pp[1]][update_key][i])
                            + ", ens_nb:"
                            + str(i)
                        )
                # ckeck if new value of Zroot is feasible
                # ------------------------------------------------------------
                elif "Zroot".casefold() in pp[1].casefold():
                    print('*'*26)
                    print(self.dict_parm_pert[pp[1]][update_key])
                    print(abs(min(self.grid3d["mesh3d_nodes"][:, -1])))

                    id_nonvalid.append(
                        list(
                            np.where(
                                self.dict_parm_pert[pp[1]][update_key]
                                > abs(min(self.grid3d["mesh3d_nodes"][:, -1]))
                            )[0]
                        )
                    )
                    id_valid.append(
                        list(
                            np.where(
                                self.dict_parm_pert[pp[1]][update_key]
                                < abs(min(self.grid3d["mesh3d_nodes"][:, -1]))
                            )[0]
                        )
                    )
                    for i in id_nonvalid:
                        self.console.print(
                            ":x: [b]unfeasible root depth:[/b]"
                            + str(self.dict_parm_pert[pp[1]][update_key][i])
                            + ", ens_nb:"
                            + str(i)
                        )

        id_nonvalid_flat = [item for sublist in id_nonvalid for item in sublist]
        id_valid = set.difference(
            set(list(np.arange(0, self.NENS))), set(list(np.unique(id_nonvalid_flat)))
        )
        return id_valid

    def spatialize_parameters(self, list_update_parm, parm):
        """
        Extend the number of parameters according to the number of zones
        (useful for the heteregeneous case)
        """
        i = 0
        parm_extended = []
        for l in list_update_parm:
            print(l)
            if "St. var." in l:
                pass
            else:
                nb_of_zones = 1
                if "surf_zones_param" in self.dict_parm_pert[l]:
                    nb_of_zones = self.dict_parm_pert[l]["surf_zones_param"]
                    print(nb_of_zones)
                    parm_extended.append(np.tile(parm[i, :], (nb_of_zones, 1)))
                i = i + 1
        parm_extended = np.vstack(parm_extended)
        return parm_extended

    def transform_parameters(
        self, list_update_parm, update_key=None, back=False, param=None
    ):
        """
        Parameter (back) transformation before analysis (log, bounded log, ...)

        ..note:: parameters are log transform so that they become Gaussian distributed
                 to be updated by the Ensemble Kalman Filter
        """
        param_new = []
        for pp in enumerate(list_update_parm[:]):
            param_trans = []

            if "St. var." in pp[1]:
                pass
            else:
                if update_key:
                    param = self.dict_parm_pert[pp[1]][update_key]
                if self.dict_parm_pert[pp[1]]["transf_type"] == "log":
                    print("log transformation")
                    if back:
                        param_trans = np.exp(param)
                    else:
                        param_trans = np.log(param)
                elif self.dict_parm_pert[pp[1]]["transf_type"] == "log-ratio":
                    print("bounded log transformation")
                    range_tr = self.dict_parm_pert[pp[1]]["transf_bounds"]
                    A = range_tr["A"]  # min range
                    B = range_tr["B"]  # max range
                    if back:
                        param_trans = np.exp((param - A) / (B - param))
                    else:
                        param_trans = np.log((param - A) / (B - param))
                elif self.dict_parm_pert[pp[1]]["transf_type"] == "hyperbolic":
                    print("hyperbolic transformation")
                    # Y = sinh−1(U )
                    param_trans = np.log(self.dict_parm_pert[pp[1]][update_key])

                else:
                    param_trans = param
                param_new.append(param_trans)
            if len(param_new) > 0:
                np.vstack(param_new)

        return np.array(param_new)

    def _evaluate_perf_OL(
        self, parallel, list_assimilated_obs, path_fwd_CATHY_list, ENS_times
    ):

        prediction_OL = []
        if parallel:
            if "ERT" in list_assimilated_obs:
                prediction_OL_ERT = mapper._map_ERT_parallel(
                    
                                                            self.dict_obs,
                                                            self.count_DA_cycle,
                                                            self.project_name,
                                                            self.Archie_parms,
                                                            path_fwd_CATHY_list,
                                                            list_assimilated_obs="all",
                                                            default_state="psi",
                                                            verbose=False,
                                                            savefig=True,
                                                            DA_cnb=self.count_DA_cycle,
                                                            ENS_times=self.ENS_times,
                                                        )
            else:
                write2shell_map = True
                for t in range(len(ENS_times)):
                    prediction_stacked = self.map_states2Observations(
                        list_assimilated_obs,
                        ENS_times=ENS_times,
                        savefig=False,
                        parallel=parallel,
                        write2shell_map=write2shell_map,
                    )
                    prediction_OL.append(prediction_stacked)
                    write2shell_map = False

        else:
            for t in range(len(ENS_times)):
                prediction_stacked = self.map_states2Observations(
                    list_assimilated_obs,
                    ENS_times=ENS_times,
                    savefig=False,
                    parallel=parallel,
                )
                # prediction_stacked is meas_size * ens_size
                prediction_OL.append(prediction_stacked)
                # prediction_OL is meas_size * ens_size * ENS_times size

        if parallel:
            if "ERT" in list_assimilated_obs:
                prediction_OL.append(prediction_OL_ERT)

                if np.shape(prediction_OL)[0] > 1:
                    prediction_OL = np.hstack(prediction_OL)

        prediction_OL = np.reshape(
            prediction_OL,
            [
                np.shape(prediction_OL)[1],
                self.NENS,
                len(ENS_times),
            ],
        )
        for t in range(len(ENS_times)):
            # print(str(t) + 't perf OL')
            data_t, _ = self._get_data2assimilate(
                list_assimilated_obs,
                time_ass=t,
            )
            self._performance_assessement(
                list_assimilated_obs,
                data_t,
                prediction_OL[:, :, t],
                t_obs=t,
                openLoop=True,
            )
        return prediction_OL

    def _mapping_petro_init(self):
        """
        Initiate Archie and VGP petro/pedophysical parameters

        If Archie and VGP dictionnaries are not existing fill with default values
        """

        # print('Initiate Archie and VGP')
        warnings_petro = []

        # check that porosity is a list
        # -------------------------------------
        
        # if hasattr(self, 'soil_SPP') is False:
            # self.read_inputs('soil')
            
        porosity = self.soil_SPP["SPP"][:, 4][0]
        if not isinstance(porosity, list):
            porosity = [porosity]

        # check if Archie parameters exists
        # ---------------------------------
        if hasattr(self, "Archie_parms") == False:
            warnings_petro.append(
                "Archie parameters not defined - Falling back to defaults"
            )
            print("Archie parameters not defined set defaults")
            self.Archie_parms = {
                "porosity": porosity,
                "rFluid_Archie": [1.0],
                "a_Archie": [1.0],
                "m_Archie": [2.0],
                "n_Archie": [2.0],
                "pert_sigma_Archie": [0],
            }

        # check if VGN parameters exists
        # ---------------------------------
        if hasattr(self, "VGN_parms") == False:
            warnings_petro.append(
                "VGN parameters not defined - Falling back to defaults"
            )
            self.VGN_parms = {}

        self.console.rule(":warning: warning messages above :warning:", style="yellow")
        for ww in warnings_petro:
            self.console.print(ww, style="yellow")
        self.console.rule("", style="yellow")

        pass


    def _read_dict_pert_update_ens_SPP(self,
                                       dict_parm_pert,
                                       df_SPP_2fill,
                                       key_root,
                                       key,
                                       update_key,
                                       ens_nb,
                                       shellprint_update
                                       ):
        '''
        Use dictionnary of perturbated parameters and parse it progressively to 
        df_SPP_2fill to build a matrice and update soil
        '''        
        # initiate SPP_2_fill if not existing
        # -----------------------------------------------------------
        if len(df_SPP_2fill) == 0:
            
            for i in range(self.NENS):

                df_SPP, _ = in_CT.read_soil(os.path.join(
                                                        self.workdir,
                                                        self.project_name,
                                                        'DA_Ensemble/cathy_'+ str(i+1),
                                                        'input/soil'
                                                        ),
                                            self.dem_parameters,
                                            self.cathyH["MAXVEG"]
                                            )           
                # df_SPP_2fill = [df_SPP.copy() for _ in range(self.NENS)]
                df_SPP_2fill.append(df_SPP.copy())
                
                if key_root[0].casefold() == 'ks':
                    df_SPP_2fill[i]['PERMX'] = np.nan
                    df_SPP_2fill[i]['PERMY'] = np.nan
                    df_SPP_2fill[i]['PERMZ'] = np.nan
                if key_root[0].casefold() == 'porosity':
                    df_SPP_2fill[i]['POROS'] = np.nan

            print('SPP ensemble update per zones is not yet implemented')

            print('''
                  Assuming that the dictionnary of perturbated parameters 
                  is build per layers i.e. ks0= layer 0, ks1=layer 1, etc...
                  '''
                  )        
        if key_root[0].casefold() == 'ks':
            df_SPP_2fill[ens_nb].loc[(slice(None), int(key_root[1])), 'PERMX'] = dict_parm_pert[key][update_key][ens_nb]
            df_SPP_2fill[ens_nb].loc[(slice(None), int(key_root[1])), 'PERMY'] = dict_parm_pert[key][update_key][ens_nb]
            df_SPP_2fill[ens_nb].loc[(slice(None), int(key_root[1])), 'PERMZ'] = dict_parm_pert[key][update_key][ens_nb]
        elif key_root[0].casefold() in 'porosity':
            df_SPP_2fill[ens_nb].loc[(slice(None), int(key_root[1])), 'POROS'] = dict_parm_pert[key][update_key][ens_nb]
            
        #%% check if contains nan (if not update soil file)
        contains_nan = df_SPP_2fill[ens_nb].isna().any().any()
        
        if hasattr(self, 'zone3d') is False:
            nstr = self.dem_parameters["nstr"]
            self.zone3d= [np.ones(np.shape(self.DEM))]*nstr
        
        if contains_nan==False:
            self.update_soil(
                SPP_map=df_SPP_2fill[ens_nb],
                FP_map=self.soil_FP["FP_map"],
                zone3d=self.zone3d,
                verbose=self.verbose,
                filename=os.path.join(os.getcwd(), "input/soil"),
                shellprint_update=shellprint_update,
                backup=True,
            )

        return df_SPP_2fill
                            
    def update_ENS_files(self, dict_parm_pert, list_parm2update, **kwargs):
        """
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

        """
        self.console.print(":arrows_counterclockwise: [b]Update ensemble[/b]")

        # replace dictionnary key with cycle nb
        # -------------------------------------
        update_key = "ini_perturbation"
        if "cycle_nb" in kwargs:
            if kwargs["cycle_nb"] > 0:
                update_key = "update_nb" + str(kwargs["cycle_nb"])

        if "analysis" in kwargs:
            analysis = kwargs["analysis"]

        if list_parm2update == "all":
            list_parm2update = ["St. var."]
            for pp in dict_parm_pert:
                list_parm2update.append(pp)

        # loop over ensemble files to update ic -->
        # this is always done since psi are state variable to update
        # ------------------------------------------------------------------
        shellprint_update = True
        for ens_nb in range(self.NENS):
            # change directory according to ensmble file nb
            os.chdir(
                os.path.join(
                    self.workdir,
                    self.project_name,
                    "./DA_Ensemble/cathy_" + str(ens_nb + 1),
                )
            )

            # state variable update use ANALYSIS update or not
            # --------------------------------------------------------------
            if "St. var." in list_parm2update:
                if kwargs["cycle_nb"] > 0:
                    df_psi = self.read_outputs(
                        filename="psi",
                        path=os.path.join(os.getcwd(), self.output_dirname),
                    )
                    self.update_ic(
                        INDP=1,
                        IPOND=0,
                        pressure_head_ini=analysis[:, ens_nb],
                        filename=os.path.join(os.getcwd(), "input/ic"),
                        backup=True,
                        shellprint_update=shellprint_update,
                    )
            else:
                if kwargs["cycle_nb"] > 0:
                    raise ValueError(
                        "no state variable update - use last iteration as initial conditions ?"
                    )
                    df_psi = self.read_outputs(
                        filename="psi",
                        path=os.path.join(os.getcwd(), self.output_dirname),
                    )
                    if kwargs["cycle_nb"] > 0:
                        self.update_ic(
                            INDP=1,
                            IPOND=0,
                            pressure_head_ini=df_psi[-1, :],
                            filename=os.path.join(os.getcwd(), "input/ic"),
                            backup=False,
                            shellprint_update=shellprint_update,
                        )
            shellprint_update = False
        
        # Initialise matrice of ensemble
        # ------------------------------
        FeddesParam_mat_ens = []  # matrice of ensemble for Feddes parameters
        VG_parms_mat_ens = []  # matrice of ensemble for VG parameters
        VG_p_possible_names = ["n_VG", "thetar_VG", "alpha_VG", "VGPSATCELL_VG"]
        VG_p_possible_names_positions_in_soil_table = [5, 6, 7, 7]
       
        # SPP_possible_names = ["porosity", "ks"]
        ks_het_ens = []  # matrice of ensemble for heterogeneous SPP parameters
        POROS_het_ens = []  # matrice of ensemble for heterogeneous SPP parameters
        
        # Archie
        # -------
        Archie_parms_mat_ens = []  # matrice of ensemble for Archie parameters
        Archie_p_names = [
            "porosity",
            "rFluid_Archie",
            "a_Archie",
            "m_Archie",
            "n_Archie",
            "pert_sigma_Archie",
        ]
        # Soil 3d if existing
        # -------------------
        zone3d = []
        if hasattr(self, 'zone3d'):
            zone3d = self.zone3d

        # loop over dict of perturbated variable
        # ----------------------------------------------------------------------
        for parm_i, key in enumerate(
            list_parm2update
        ):  # loop over perturbated variables dict
            # print('update parm' + key)


            key_root = re.split("(\d+)", key)
            if len(key_root) == 1:
                key_root.append("0")

            shellprint_update = True
            # ------------------------------------------------------------------
            for ens_nb in range(self.NENS):
                # print('update ensemble' + str(ens_nb))
                # change directory according to ensmble file nb
                os.chdir(
                    os.path.join(
                        self.workdir,
                        self.project_name,
                        "./DA_Ensemble/cathy_" + str(ens_nb + 1),
                    )
                )

                # atmbc update
                # --------------------------------------------------------------
                # if key == "atmbc":
                # if non-homogeneous atmbc, take only the first node perturbated parameters
                # tp avoid looping over all the perturbated atmbc nodes as the update is done
                # only once over the stacked perturbated atmbc parameters 
                if key_root[0].casefold() in "atmbc".casefold():
                    
                    if key.casefold() in "atmbc0".casefold():
                        tper_stacked = []
                        for kk in dict_parm_pert:
                            if 'atmbc' in kk:
                                tper_stacked.append(dict_parm_pert[kk]['time_variable_perturbation'])
                        tper_stacked = np.array(tper_stacked)
                        VALUE = tper_stacked[:,:,ens_nb]
                        times = dict_parm_pert[key]['data2perturbate']['time']                           
                        self.update_atmbc(
                                            HSPATM=self.atmbc['HSPATM'],
                                            IETO=self.atmbc['IETO'],
                                            time=times,
                                            netValue=VALUE.T,
                                            filename=os.path.join(os.getcwd(), "input/atmbc"),
                                            verbose=self.verbose,
                                        )
                        
                    else:
                        pass


                # ic update (i.e. water table position update)
                # --------------------------------------------------------------
                elif key_root[0].casefold() in "WTPOSITION".casefold():
                    if kwargs["cycle_nb"] == 0:
                        self.update_ic(
                            INDP=4,
                            IPOND=0,
                            WTPOSITION= dict_parm_pert[key]["ini_perturbation"][
                                ens_nb
                            ],
                            filename=os.path.join(os.getcwd(), "input/ic"),
                            backup=False,
                            shellprint_update=shellprint_update,
                        )
                        
                # ic update (i.e. initial pressure head update)
                # --------------------------------------------------------------
                elif key_root[0].casefold() in "ic".casefold():
                    match_ic_withLayers = re.search(r'ic\d+', key)
                        
                    if match_ic_withLayers:
                        DApath = os.getcwd()
                        saveMeshPath = os.path.join(DApath, "vtk", 
                                                     self.project_name + '.vtk'
                                                     ) 
                        ic_values_by_layers = dict_parm_pert[key]['ini_ic_withLayers'][:,ens_nb]
                        ic_nodes = self.map_prop_2mesh_markers('ic',
                                                                ic_values_by_layers, 
                                                                to_nodes=True,
                                                                saveMeshPath=saveMeshPath,
                                                                )
                        self.update_ic(
                            INDP=0,
                            IPOND=0,
                            pressure_head_ini=ic_nodes,
                            filename=os.path.join(DApath, "input/ic"),
                            backup=False,
                            shellprint_update=shellprint_update,
                        )
                    else:
                        if kwargs["cycle_nb"] == 0:
                            self.update_ic(
                                INDP=0,
                                IPOND=0,
                                pressure_head_ini=dict_parm_pert[key]["ini_perturbation"][ens_nb],
                                filename=os.path.join(os.getcwd(), "input/ic"),
                                backup=False,
                                shellprint_update=shellprint_update,
                            )

                # kss update
                # --------------------------------------------------------------
                elif key_root[0].casefold() in 'ks':                   
                    ks_het_ens = self._read_dict_pert_update_ens_SPP(
                                                       dict_parm_pert,
                                                       ks_het_ens,
                                                       key_root,
                                                       key,
                                                       update_key,
                                                       ens_nb,
                                                       shellprint_update
                                                       )
                    
                    
                elif key_root[0].casefold() in 'porosity':
                    
                    if hasattr(self, 'zone3d'):
                        POROS_het_ens = self._read_dict_pert_update_ens_SPP(
                                                           dict_parm_pert,
                                                           POROS_het_ens,
                                                           key_root,
                                                           key,
                                                           update_key,
                                                           ens_nb,
                                                           shellprint_update
                                                           )

                    else:
                        df_SPP, _ = in_CT.read_soil(os.path.join(
                                                                self.workdir,
                                                                self.project_name,
                                                                'DA_Ensemble/cathy_'+ str(ens_nb+1),
                                                                'input/soil'
                                                                ),
                                                    self.dem_parameters,
                                                    self.cathyH["MAXVEG"]
                                                    )    
                        
                        df_SPP_2fill = df_SPP.copy()
                        df_SPP_2fill['POROS'] = dict_parm_pert[key][update_key][ens_nb]
                        df_SPP_2fill['POROS']
                        
                        self.update_soil(
                            SPP_map=df_SPP_2fill,
                            FP_map=self.soil_FP["FP_map"],
                            verbose=self.verbose,
                            filename=os.path.join(os.getcwd(), "input/soil"),
                            shellprint_update=shellprint_update,
                            backup=True,
                        )
                        
                    
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

                    SPP_map = self.soil_SPP["SPP_map"]
                    if len(VG_parms_mat_ens) == 0:
                        VG_parms_mat = self.soil_SPP["SPP"]
                        np.shape(VG_parms_mat)
                        if len(SPP_map["PERMX"]) < 2:
                            VG_parms_mat_ens = np.matlib.repmat(
                                VG_parms_mat[0], self.NENS, 1
                            )
                        else:
                            # VG_parms_mat_ens = np.repeat(VG_parms_mat, self.NENS, axis=0)
                            print("Non homogeneous soil VG properties not implemented")

                    if len(SPP_map["PERMX"]) < 2:
                        # for k in parm_pert.keys():
                        for k in VG_p_possible_names:
                            match = [
                                parm_incr
                                for parm_incr in dict_parm_pert.keys()
                                if k in parm_incr
                            ]
                            if len(match) > 0:
                                idVG = VG_p_possible_names_positions_in_soil_table[
                                    VG_p_possible_names.index(k)
                                ]
                                VG_parms_mat_ens[:, idVG] = dict_parm_pert[match[0]][
                                    update_key
                                ]
                            else:
                                pass
                    else:
                        # VG_parms_mat_ens[:,:,idVG] = np.matlib.repmat(parm_pert[key][update_key],len(VG_parms_mat),1).T
                        print("Non homogeneous soil VG properties not implemented")

                    VG_parms_ensi = []
                    for es in range(self.NENS):
                        fed = dict()
                        for i, f in enumerate(list(SPP_map.keys())):
                            if len(SPP_map["PERMX"]) < 2:
                                if (
                                    f == "alpha_VG"
                                ):  # convert to CATHY VGPSATCELL == -1/alpha
                                    fed[f] = -1 / VG_parms_mat_ens[es, i]
                                else:
                                    fed[f] = [VG_parms_mat_ens[es, i]]
                            else:
                                print(
                                    "Non homogeneous soil VG properties not implemented"
                                )
                                # fed[f] = list(VG_parms_mat_ens[es,:,i])
                        VG_parms_ensi.append(fed)

                    self.update_soil(
                        SPP_map=VG_parms_ensi[ens_nb],
                        FP_map=self.soil_FP["FP_map"],
                        # zone3d=self.zone3d,
                        verbose=self.verbose,
                        filename=os.path.join(os.getcwd(), "input/soil"),
                        shellprint_update=shellprint_update,
                        backup=True,
                    )

                # FeddesParam update
                # --------------------------------------------------------------
                elif key_root[0] in ["PCANA", "PCREF", "PCWLT", "ZROOT", "PZ", "OMGC"]:
                    SPP_map = self.soil_SPP["SPP_map"]

                    if len(FeddesParam_mat_ens) == 0:
                        FeddesParam_mat = self.soil_FP["FP"]
                        FeddesParam_mat_ens = np.repeat(
                            FeddesParam_mat, self.NENS, axis=0
                        )
                        FeddesParam_mat_ens = np.reshape(
                            FeddesParam_mat_ens, (self.NENS, self.cathyH["MAXVEG"], 6)
                        )

                    idFeddes = ["PCANA", "PCREF", "PCWLT", "ZROOT", "PZ", "OMGC"].index(
                        key_root[0]
                    )

                    # if self.cathyH["MAXVEG"]<2: # Case where only 1 vegetation type
                    FeddesParam_mat_ens[:, int(key_root[1]), idFeddes] = dict_parm_pert[
                        key
                    ][update_key]
                    # FeddesParam_mat_ens[:,:,idFeddes] = parm_pert[key][update_key]
                    # else: # Case where only 1 vegetation type

                    FeddesParam_map_ensi = []
                    for es in range(self.NENS):
                        fed = dict()
                        for i, f in enumerate(
                            ["PCANA", "PCREF", "PCWLT", "ZROOT", "PZ", "OMGC"]
                        ):
                            fed[f] = list(FeddesParam_mat_ens[es, :, i])
                        FeddesParam_map_ensi.append(fed)
                    
                    self.update_soil(
                        FP_map=FeddesParam_map_ensi[ens_nb],
                        SPP_map=SPP_map,
                        zone3d=zone3d,
                        verbose=self.verbose,
                        filename=os.path.join(os.getcwd(), "input/soil"),
                        shellprint_update=shellprint_update,
                        backup=True,
                    )
                # Archie_p update
                # --------------------------------------------------------------
                elif key_root[0] in Archie_p_names:
                    self.console.print(
                        ":arrows_counterclockwise: [b]Update Archie parameters[/b]"
                    )

                    # self.Archie_parms = {'rFluid':rFluid, 'a':a, 'm':m, 'n':n, 'pert_sigma':pert_sigma}
                    idArchie = Archie_p_names.index(key_root[0])
                    if len(Archie_parms_mat_ens) == 0:
                        for p in Archie_p_names:
                            if len(self.Archie_parms[p]) != self.NENS:
                                self.Archie_parms[p] = list(
                                    self.Archie_parms[p] * np.ones(self.NENS)
                                )
                        Archie_parms_mat_ens = np.zeros(
                            [len(Archie_p_names), self.NENS]
                        )
                        Archie_parms_mat_ens[:] = np.nan
                    Archie_parms_mat_ens[idArchie, ens_nb] = dict_parm_pert[key][
                        update_key
                    ][ens_nb]

                    if np.isnan(Archie_parms_mat_ens[idArchie, :]).any() == False:
                        self.Archie_parms[key_root[0]] = list(
                            Archie_parms_mat_ens[idArchie, :]
                        )

                    # print(self.Archie_parms)

                else:
                    # variable perturbated to ommit: state, Archie param
                    var_pert_to_ommit = ["St. var."]

                    if not key in var_pert_to_ommit:
                        raise ValueError(
                            "key:" + str(key) + " to update is not existing!"
                        )

            shellprint_update = False

        pass



    def _get_data2assimilate(self, list_assimilated_obs, time_ass=None, match=False):
        """
        Loop over observation dictionnary and select corresponding data for a given assimilation time

        Parameters
        ----------
        list_assimilated_obs : list
            list of observation to assimilate.

        Returns
        -------
        data : list
            list of data values.
        """

        # find the method to map for this time step
        # --------------------------------------------------------
        obskey2map, obs2map = self._obs_key_select(list_assimilated_obs)

        def extract_data(sensor, time_ass, data):
            data2add = []
            if "ERT" in sensor:
                if "pygimli" in items_dict[time_ass][1][sensor]["data_format"]:
                    data2add.append(
                        np.array(items_dict[time_ass][1][sensor]["data"]["rhoa"])
                    )
                else:
                    data2add.append(
                        items_dict[time_ass][1][sensor]["data"]["resist"].to_numpy()
                    )
                data2add = np.hstack(data2add)
            else:
                data2add.append(items_dict[time_ass][1][sensor]["data"])

            return data2add

        if time_ass is None:
            time_ass = self.count_DA_cycle

        data = []
        # Loop trought observation dictionnary for a given assimilation time (count_DA_cycle)
        # -----------------------------------------------------------------
        items_dict = list(self.dict_obs.items())
        for sensor in items_dict[time_ass][1].keys():
            if "all" in list_assimilated_obs:
                data2add = extract_data(sensor, time_ass, data)
                data.append(data2add)
            else:
                if match:
                    if sensor in list_assimilated_obs:
                        data2add = extract_data(sensor, time_ass, data)
                        data.append(data2add)
                else:
                    str_sensor_rootname = re.sub(r'\d', '', sensor)
                    if str_sensor_rootname in list_assimilated_obs:
                        data2add = extract_data(sensor, time_ass, data)
                        data.append(data2add)

        return np.hstack(data), obskey2map

    def _obs_key_select(self, list_assimilated_obs):
        """find data to map from dictionnary of observations"""
        # Loop trought observation dictionnary for a given assimilation time (count_DA_cycle)
        # -----------------------------------------------------------------
        obskey2map = []  # data type 2 map
        obs2map = []  # data dict 2 map
        items_dict = list(self.dict_obs.items())
        for sensor in items_dict[self.count_DA_cycle][1].keys():
            if "all" in list_assimilated_obs:
                obskey2map.append(sensor)
                obs2map.append(items_dict[self.count_DA_cycle][1][sensor])
            else:
                str_sensor_rootname = re.sub(r'\d', '', sensor)
                if str_sensor_rootname in list_assimilated_obs:
                    obskey2map.append(sensor)
                    obs2map.append(items_dict[self.count_DA_cycle][1][sensor])
        return obskey2map, obs2map

    def map_states2Observations(
                                    self, 
                                    list_assimilated_obs="all", 
                                    parallel=False, 
                                    default_state="psi", 
                                    **kwargs
                                ):
        """
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

        Parameters
        ----------
        state : np.array([])
            [pressure heads, saturation] for each mesh nodes. The default is [None, None].
        list_assimilated_obs : TYPE, optional
            DESCRIPTION. The default is 'all'.
        parallel : Bool, optional
            DESCRIPTION. The default is False.
        **kwargs : TYPE
            DESCRIPTION.

        Returns
        -------
        Hx_ens : np.array (size: nb of observations x nb of ensembles)
            Ensemble of the simulated (predicted) observations.
        """
        #%%
        if parallel:
            path_fwd_CATHY_list = []  # list of ensemble path of cathy output
            for ens_nb in self.ens_valid:  # loop over ensemble files
                path_fwd_CATHY_list.append(
                    os.path.join(
                        self.workdir,
                        self.project_name,
                        "DA_Ensemble/cathy_" + str(ens_nb + 1),
                    )
                )

        write2shell_map = True
        if "write2shell_map" in kwargs:
            write2shell_map = kwargs["write2shell_map"]

        Hx_ens = []  # matrice of predicted observation for each ensemble realisation
        for ens_nb in self.ens_valid:  # loop over ensemble files
            # print('ens nb:' + str(ens_nb))

            path_fwd_CATHY = os.path.join(
                self.workdir, self.project_name, "DA_Ensemble/cathy_" + str(ens_nb + 1)
            )

            df_psi = self.read_outputs(
                filename="psi", path=os.path.join(path_fwd_CATHY, self.output_dirname)
            )
            df_sw, _ = self.read_outputs(
                filename="sw", path=os.path.join(path_fwd_CATHY, self.output_dirname)
            )

            # infer soil parameters properties
            # ---------------------------------
            # porosity = self.soil_SPP["SPP"][:, 4][0]
            SPP_ensi, _ = in_CT.read_soil(os.path.join(path_fwd_CATHY, 
                                                       self.input_dirname,'soil'
                                                       ),
                                          self.dem_parameters,
                                          self.cathyH["MAXVEG"]
                                          )
            # find data to map with dictionnary of observations
            # --------------------------------------------
            obskey2map, obs2map = self._obs_key_select(list_assimilated_obs)
            state = [df_psi.iloc[-1], df_sw.iloc[-1].values]            
            
            Hx_stacked = []  # stacked predicted observation
            # Loop over observations to map
            # ---------------------------------------------------------------------
            for i, obs_key in enumerate(obskey2map):

                # print('obs_key nb:' + obs_key)

                if "tensio" in obs_key:
                    obs2map_node = obs2map[i]["mesh_nodes"]
                    Hx_PH = mapper.tensio_mapper(state,obs2map_node)
                    Hx_stacked.append(Hx_PH)

                if "swc" in obs_key:
                    # case 2: sw assimilation (Hx_SW)
                    # --------------------------------------------------------------------
                    obs2map_node = obs2map[i]["mesh_nodes"]
                    Hx_SW = mapper.swc_mapper(state,
                                              obs2map_node,
                                              self.grid3d['mesh3d_nodes'],
                                              self.dem_parameters,
                                              SPP_ensi,
                                              self.console
                                              )
                    Hx_stacked.append(Hx_SW)


                if 'ET' in obs_key:
                    obs2map_node = obs2map[i]["mesh_nodes"]
                    Hx_ET = mapper.ET_mapper(obs2map_node, 
                                             self.console,
                                             path_fwd_CATHY,
                                             )
                    Hx_stacked.append(Hx_ET)

                if "stemflow" in obs_key:
                    # Atmact-vf(13) : Actual infiltration (+ve) or exfiltration (-ve) at atmospheric BC nodes as a volumetric flux [L^3/T]
                    # Atmact-v (14) : Actual infiltration (+ve) or exfiltration (-ve) volume [L^3]
                    # Atmact-r (15) : Actual infiltration (+ve) or exfiltration (-ve) rate [L/T]
                    # Atmact-d (16) : Actual infiltration (+ve) or exfiltration (-ve) depth [L]
                    print("Not yet implemented")

                if "scale" in obs_key:
                    # Atmact-vf(13) : Actual infiltration (+ve) or exfiltration (-ve) at atmospheric BC nodes as a volumetric flux [L^3/T]
                    # Atmact-v (14) : Actual infiltration (+ve) or exfiltration (-ve) volume [L^3]
                    # Atmact-r (15) : Actual infiltration (+ve) or exfiltration (-ve) rate [L/T]
                    # Atmact-d (16) : Actual infiltration (+ve) or exfiltration (-ve) depth [L]
                    print("Not yet implemented")

                    df_dtcoupling = self.read_outputs(
                        filename="dtcoupling",
                        path=os.path.join(path_fwd_CATHY, self.output_dirname),
                    )

                    # Hx_scale = state[1][obs2map[i]["mesh_nodes"]] * porosity
                    # Hx_stacked.append(Hx_scale)

                if "discharge" in obs_key:
                    # case 3: discharge
                    # need to read the hgsfdet file (Hx_Q)
                    # --------------------------------------------------------------------
                    # if key[0] in 'discharge':
                    # derivation of the dircharge, Q from file 'hgsfdet'
                    # Hx_Q = []
                    # Hx.vstack(Hx_Q)
                    print("Not yet implemented")

                if "ERT" in obs_key:
                    if parallel == False:
                        Hx_ERT, df_Archie, mesh2test = mapper._map_ERT(
                                                            state, 
                                                            path_fwd_CATHY, 
                                                            ens_nb, 
                                                            **kwargs
                                                            )
                        df_Archie_new = mapper._add_2_ensemble_Archie(self.df_Archie,
                                                                      df_Archie
                                                                      )
                        self.df_Archie = df_Archie_new

                        if "pygimli" in obs2map[i]["data_format"]:
                            Hx_stacked.append(Hx_ERT["rhoa"])
                        else:
                            Hx_stacked.append(Hx_ERT["resist"])
                        Hx_stacked = np.hstack(Hx_stacked)

            if len(Hx_stacked) > 0:
                Hx_ens.append(Hx_stacked)

        if isinstance(Hx_ens, float):
            pass
        else:
            if len(Hx_ens)>0:
                if np.shape(Hx_ens)!=(len(self.ens_valid),len(obskey2map)):
                    Hx_ens = np.hstack(Hx_ens)
                else:
                    Hx_ens = np.array(Hx_ens).T
               
        #%%
        # ---------------------------------------------------------------------
        # special case of ERT // during sequential assimilation
        # ---------------------------------------------------------------------
        for i, obs_key in enumerate(obskey2map):
            if "ERT" in obs_key:

                if write2shell_map:
                    self.console.print(":zap: Mapping to ERT prediction")
                    self.console.print(
                        r"""
                                    Fwd modelling ERT
                                    Noise level: {??}
                                    """,
                        style="green",
                    )

                    self.console.rule(
                        ":octopus: Parameter perturbation :octopus:", style="green"
                    )
                    self.console.print(
                        r"""
                                    Archie perturbation for DA analysis
                                    Archie rFluid: {}
                                    Pert. Sigma Archie: {}
                                    Nb of zones (?): {}
                                    """.format(
                            self.Archie_parms["rFluid_Archie"],
                            self.Archie_parms["pert_sigma_Archie"],
                            len(self.Archie_parms["rFluid_Archie"]),
                        ),
                        style="green",
                    )
                    self.console.rule("", style="green")

                if parallel:
                    Hx_ens_ERT, df_Archie, mesh2test = mapper._map_ERT_parallel(
                                                                    self.dict_obs,
                                                                    self.count_DA_cycle,
                                                                    self.project_name,
                                                                    self.Archie_parms,
                                                                    path_fwd_CATHY_list,
                                                                    list_assimilated_obs="all",
                                                                    default_state="psi",
                                                                    DA_cnb=self.count_DA_cycle,
                                                                    savefig=False,
                                                                    verbose=False,
                                                                    )
                    
                    # cell_array_names = mesh2test.array_names
                    # mesh2test['ER_converted1'].min()
                    # mesh2test['ER_converted1'].max()

                    # mesh2test.set_active_scalars('ER_converted1')
                    
                    df_Archie_new = mapper._add_2_ensemble_Archie(self.df_Archie,
                                                                  df_Archie
                                                                  )
                    self.df_Archie = df_Archie_new                
                    
                    
                    if len(Hx_ens) > 0:
                        Hx_ens = np.vstack([Hx_ens, Hx_ens_ERT])
                    else:
                        Hx_ens = Hx_ens_ERT
                        
        return Hx_ens  # meas_size * ens_size

    def _read_state_ensemble(self):
        # read grid to infer dimension of the ensemble X
        # --------------------------------------------------------------------
        # self.grid3d = {}
        if len(self.grid3d) == 0:
            self.grid3d = in_CT.read_grid3d(
                os.path.join(self.workdir, self.project_name)
            )

        M_rows = np.shape(self.grid3d["mesh3d_nodes"])[0]
        N_col = self.NENS
        # N_col = len(ens_valid)

        # Ensemble matrix X of M rows and N  + fill with psi values
        # --------------------------------------------------------------------
        psi = np.zeros([M_rows, N_col]) * 1e99
        sw = np.zeros([M_rows, N_col]) * 1e99
        for j in range(self.NENS):
            # for j in range(len(ens_valid)):
            try:
                df_psi = out_CT.read_psi(
                    os.path.join(
                        self.workdir,
                        self.project_name,
                        "DA_Ensemble/cathy_" + str(j + 1),
                        "output/psi",
                    )
                )
                psi[:, j] = df_psi.iloc[-1, :]
            except:
                pass

            try:
                df_sw, _ = out_CT.read_sw(
                    os.path.join(
                        self.workdir,
                        self.project_name,
                        "DA_Ensemble/cathy_" + str(j + 1),
                        "output/sw",
                    )
                )
                sw[:, j] = df_sw.iloc[-1, :]
            except:
                pass
        # check if there is still zeros
        if np.count_nonzero(psi == 1e99) != 0:
            print("!!!unconsistent filled X, missing values!!!")
            # sys.exit()

        return psi, sw, N_col, M_rows

    def update_pert_parm_dict(self, update_key, list_update_parm, analysis_param_valid):
        """update dict of perturbated parameters i.e. add a new entry with new params"""

        index_update = 0
        if self.count_DA_cycle > 0:
            for pp in list_update_parm:
                if "St. var." in pp:
                    index_update = index_update + 1
                    pass
                else:
                    update_key = "update_nb" + str(self.count_DA_cycle)
                    self.dict_parm_pert[pp][update_key] = analysis_param_valid[
                        index_update - 1
                    ]
                    index_update = index_update + 1

            # check after analysis
            # ----------------------------------------------------------------
            id_valid_after = self._check_after_analysis(update_key, list_update_parm)
            intersection_set = set.intersection(
                set(id_valid_after), set(self.ens_valid)
            )
            self.ens_valid = list(intersection_set)

        pass

    def _update_input_ensemble(
        self, ENS_times, list_parm2update=["St. var."], analysis=[], **kwargs
    ):
        """
        Select the time window of the hietograph.
        Write new variable perturbated/updated value into corresponding input file

        Parameters
        ----------
        parm_pert : dict
            dict of all perturbated variables holding values and metadata.

        Returns
        -------
        New file written/overwritten into the input dir.

        """

        self.console.print(":sponge: [b]update input ensemble[/b]")

        # resample atmbc file for the given DA window
        # ---------------------------------------------------------------------
        self.selec_atmbc_window(self.NENS, ENS_times)

        # ---------------------------------------------------------------------
        self.update_ENS_files(
                                self.dict_parm_pert,
                                list_parm2update,
                                cycle_nb=self.count_DA_cycle,
                                analysis=analysis,
                            )

        pass

    def selec_atmbc_window(self, NENS, ENS_times):
        """
        Select the time window of the hietograph
        == time between two assimilation observations

        Parameters
        ----------
        NENS : list
            # of valid ensemble.
        ENS_times : list
            Assimilation times (in sec).
        """

        if len(self.grid3d) == 0:
            self.run_processor(IPRT1=3, DAFLAG=0)
            self.grid3d = out_CT.read_grid3d(os.path.join(self.workdir,
                                                          self.project_name, 
                                                          'output', 'grid3d')
                                             )
        # read full simulation atmbc and filter time window
        # ----------------------------------------------------------------------
        df_atmbc, HSPATM, IETO = in_CT.read_atmbc(
            os.path.join(self.workdir, self.project_name, "input", "atmbc"),
            grid=self.grid3d,
        )
        
        # Lopp over ensemble
        # ------------------
        for ens_nb in range(self.NENS):
            
            cwd = os.path.join(self.workdir,
                               self.project_name,
                               "./DA_Ensemble/cathy_" + str(ens_nb + 1)
                               )
            # Case where atmbc are perturbated
            # -------------------------------
            if 'atmbc0' in self.dict_parm_pert.keys(): # case where atmbc are homogeneous
                if HSPATM!=0:
                    VALUE = np.array(self.dict_parm_pert['atmbc0']['time_variable_perturbation'])[:,ens_nb]
                    times = self.dict_parm_pert['atmbc0']['data2perturbate']['time']
                    df_atmbc = pd.DataFrame(np.c_[times,VALUE],
                                            columns=['time','value']
                                            )
                    
            if self.count_atmbc_cycle is not None:
                if HSPATM==0:
                    time_window_atmbc = [
                        df_atmbc.time.unique()[self.count_atmbc_cycle],
                        df_atmbc.time.unique()[self.count_atmbc_cycle+1]    
                        ]
                else:
                    time_window_atmbc = [
                        df_atmbc.iloc[self.count_atmbc_cycle]["time"],
                        df_atmbc.iloc[self.count_atmbc_cycle + 1]["time"],
                    ]
            else:
                time_window_atmbc = [ENS_times[0], ENS_times[1]]
    
            df_atmbc_window = df_atmbc[
                (df_atmbc["time"] >= time_window_atmbc[0])
                & (df_atmbc["time"] <= time_window_atmbc[1])
            ]

            diff_time = time_window_atmbc[1] - time_window_atmbc[0]
            self.update_parm(
                TIMPRTi=[0, diff_time],
                TMAX=diff_time,
                filename=os.path.join(cwd, "input/parm"),
                IPRT1=2,
                backup=True,
            )

            VALUE = []
            for t in df_atmbc_window["time"].unique():
                VALUE.append(df_atmbc_window[df_atmbc_window["time"] == t]["value"].values)

            self.update_atmbc(
                HSPATM=HSPATM,
                IETO=IETO,
                time=[0, diff_time],
                netValue=VALUE,
                filename=os.path.join(cwd, "input/atmbc"),
            )

        pass

    def _DA_df(
        self, state=[None, None], state_analysis=None, rejected_ens=[], **kwargs
    ):
        """
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

        """

        df_DA_ti_ni = pd.DataFrame()  # nested df for a given ensemble/given time
        df_DA_ti = pd.DataFrame()  # nested df for a given time

        cols_root = [
            "time",
            "Ensemble_nb",
            "psi_bef_update",
            "sw_bef_update_",
            "analysis",
            "OL",  # Open loop boolean
            "rejected",
        ]

        # Open loop kwargs
        # --------------------------
        t_ass = []
        if "t_ass" in kwargs:
            t_ass = kwargs["t_ass"]

        NENS = []
        if "ens_nb" in kwargs:
            ens_nb = kwargs["ens_nb"]
        # --------------------------

        if "openLoop" in kwargs:

            data_df_root = [
                t_ass * np.ones(len(state[0])),
                ens_nb * np.ones(len(state[0])),
                state[0],
                state[1],
                state[0],
                True * np.ones(len(state[0])),
                False * np.ones(len(state[0])),
            ]
            df_DA_ti = pd.DataFrame(np.transpose(data_df_root), columns=cols_root)

        else:

            for n in range(self.NENS):

                data_df_root = [
                    self.count_atmbc_cycle * np.ones(len(state[0])),
                    n * np.ones(len(state[0])),
                    state[0][:, n],
                    state[1][:, n],
                    state_analysis[:, n],
                    False * np.ones(len(state[0])),
                    int(rejected_ens[n]) ** np.ones(len(state[0])),
                ]
                df_DA_ti_ni = pd.DataFrame(
                    np.transpose(data_df_root), columns=cols_root
                )

                df_DA_ti = pd.concat([df_DA_ti, df_DA_ti_ni], axis=0)

        self.df_DA = pd.concat([self.df_DA, df_DA_ti], axis=0, ignore_index=True)

        pass

    def _add_to_perturbated_dict(self, var_per_2add):
        """
        Update dict of perturbated variable.
        """
        self.var_per_dict = self.var_per_dict | var_per_2add
        return self.var_per_dict

    def _create_subfolders_ensemble(self):
        """
        Create ensemble subfolders
        """

        if not os.path.exists(
            os.path.join(self.workdir, self.project_name, "DA_Ensemble/cathy_origin")
        ):
            os.makedirs(
                os.path.join(
                    self.workdir, self.project_name, "DA_Ensemble/cathy_origin"
                )
            )

            # copy input, output and vtk dir
            for dir2copy in enumerate(["input", "output", "vtk"]):
                shutil.copytree(
                    os.path.join(self.workdir, self.project_name, dir2copy[1]),
                    os.path.join(
                        self.workdir,
                        self.project_name,
                        "DA_Ensemble/cathy_origin",
                        dir2copy[1],
                    ),
                )

        try:
            # for ff in [self.processor_name,'cathy.fnames','prepro']
            # copy exe into cathy_origin folder
            shutil.copy(
                os.path.join(self.workdir, self.project_name, self.processor_name),
                os.path.join(
                    self.workdir,
                    self.project_name,
                    "DA_Ensemble/cathy_origin",
                    self.processor_name,
                ),
            )
            # copy cathy.fnames into cathy_origin folder
            shutil.copy(
                os.path.join(self.workdir, self.project_name, "cathy.fnames"),
                os.path.join(
                    self.workdir,
                    self.project_name,
                    "DA_Ensemble/cathy_origin/cathy.fnames",
                ),
            )
            # copy prepro into cathy_origin folder
            shutil.copytree(
                os.path.join(self.workdir, self.project_name, "prepro"),
                os.path.join(
                    self.workdir, self.project_name, "DA_Ensemble/cathy_origin/prepro"
                ),
            )

        except:
            self.console.print(":worried_face: [b]processor exe not found[/b]")

        path_origin = os.path.join(
            self.workdir, self.project_name, "DA_Ensemble/cathy_origin"
        )

        # copy origin folder to each ensemble subfolders
        for i in range(self.NENS):
            path_nudn_i = os.path.join(
                self.workdir, self.project_name, "DA_Ensemble/cathy_" + str(i + 1)
            )

            if not os.path.exists(path_nudn_i):
                shutil.copytree(path_origin, path_nudn_i)


    def set_Archie_parm(
        self,
        porosity=[],
        rFluid_Archie=[1.0],
        a_Archie=[1.0],
        m_Archie=[2.0],
        n_Archie=[2.0],
        pert_sigma_Archie=[0],
    ):
        """
        Dict of Archie parameters. Each type of soil is describe within a list
        Note that if pert_sigma is not None a normal noise is
        added during the translation of Saturation Water to ER
    
        Parameters
        ----------
        rFluid : TYPE, optional
            Resistivity of the pore fluid. The default is [1.0].
        a : TYPE, optional
            Tortuosity factor. The default is [1.0].
        m : TYPE, optional
            Cementation exponent. The default is [2.0].
            (usually in the range 1.3 -- 2.5 for sandstones)
        n : TYPE, optional
            Saturation exponent. The default is [2.0].
        pert_sigma_Archie : TYPE, optional
            Gaussian noise to add. The default is None.
    
        ..note:
            Field procedure to obtain tce he covariance structure of the model
            estimates is described in Tso et al () - 10.1029/2019WR024964
            "Fit a straight line for log 10 (S) and log 10 (ρ S ) using the least-squares criterion.
            The fitting routine returns the covariance structure of the model estimates, which can be used to de-
            termine the 68% confidence interval (1 standard deviation) of the model estimates.""
    
        """
        if len(porosity) == 0:
            porosity = self.soil_SPP["SPP"][:, 4][0]
        if not isinstance(porosity, list):
            porosity = [porosity]
        self.Archie_parms = {
            "porosity": porosity,
            "rFluid_Archie": rFluid_Archie,
            "a_Archie": a_Archie,
            "m_Archie": m_Archie,
            "n_Archie": n_Archie,
            "pert_sigma_Archie": pert_sigma_Archie,
        }
        
    def _performance_assessement(
                                self, 
                                list_assimilated_obs, 
                                data, 
                                prediction, 
                                t_obs, 
                                **kwargs
                                ):
        """
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
    
        """
    
        OL_bool = False
        if "openLoop" in kwargs:
            OL_bool = True
    
        if hasattr(self, "df_performance") is False:
            # initiate if simulation just started
            # -------------------------------------------------------------
            self.df_performance = pd.DataFrame()
            self.RMSE_avg_stacked = []
    
            # NAs default
            # -------------------------------------------------------------
            # RMSE_avg = np.nan # root mean square error (RMSE)
            # NMRMSE = np.nan # time-averaged normalized root mror: unexpected indent
    
        # average differences over the number of ensemble
        # ALL OBSERVATIONS
        # ------------------------------
    
        all_Obs_diff_mat = np.zeros(np.shape(prediction))
        for i in range(len(data)):
            for j in range(np.shape(prediction)[1]):  # Loop over ensemble collumns
                all_Obs_diff_mat[i, j] = abs(data[i] - prediction[i, j])
                # all_Obs_diff_mat[:,j] = abs(data[i]-prediction[i,j])
    
        all_Obs_diff_avg = np.nansum(all_Obs_diff_mat, axis=1) * (1 / len(prediction))
    
        # average differences over all the observations
        all_Obs_RMSE_avg_ti = np.sum(all_Obs_diff_avg, axis=0) * (1 / len(data))
    
        # # compute metrics for each observation variable
        # # ------------------------------------------
        # Loop trought observation dictionnary for a given assimilation time (count_DA_cycle)
        # -----------------------------------------------------------------
        obs2eval_key = []  # data type 2 eval
        obs2eval = []  # data dict 2 eval
    
        _, obs2eval_key = self._get_data2assimilate(list_assimilated_obs)
    
        start_line_obs = 0
        for s, name_sensor in enumerate(obs2eval_key):  # case where there is other observations than ERT
            obs2eval, _ = self._get_data2assimilate([name_sensor], match=True)
    
            # number of observation at a given time
            # for ERT number of observation = is the number of grid cells
            n_obs = len(obs2eval)
            prediction2eval = prediction[start_line_obs : start_line_obs + n_obs]
            obs2eval_diff_mat = np.zeros(np.shape(prediction2eval))
    
            for i in range(len(obs2eval)):
                for j in range(
                    np.shape(prediction2eval)[1]
                ):  # Loop over ensemble collumns
                    obs2eval_diff_mat[i, j] = abs(obs2eval[i] - prediction2eval[i, j])
    
            obs2eval_diff_avg = np.nansum(obs2eval_diff_mat, axis=1) * (
                1 / len(prediction2eval)
            )
    
            start_line_obs = start_line_obs + n_obs
    
            if "ERT" in name_sensor:
                RMSE_sensor_ti = np.sum(obs2eval_diff_avg, axis=0) * (1 / len(data))
                # Plot here scatter Transfer resistance scatter plot between the observed and simulated data
                # plt.scatter(obs2eval[i],prediction2eval[i,j])
            else:
                RMSE_sensor_ti = obs2eval_diff_avg
    
            # compute normalised RMSE
            # ---------------------------------------------------------
    
            self.RMSE_avg_stacked.append(all_Obs_RMSE_avg_ti)
            
            if 'ObsType' in self.df_performance:
                if name_sensor in self.df_performance['ObsType'].unique():
                    bool_sensor = self.df_performance['ObsType']==name_sensor
                    RMSE_sensor_stacked = (self.df_performance[bool_sensor]["RMSE" + name_sensor].sum() 
                                       + RMSE_sensor_ti
                                       )
                else:
                    RMSE_sensor_stacked = RMSE_sensor_ti
    
            # if hasattr(self, 'RMSE') is False:
            if t_obs == 0:
                NMRMSE_sensor_ti = RMSE_sensor_ti
                NMRMSE_avg_ti = all_Obs_RMSE_avg_ti
            else:
                NMRMSE_sensor_ti = (1 / (t_obs + 1)) * RMSE_sensor_stacked
                NMRMSE_avg_ti = (1 / (t_obs + 1)) * np.sum(self.RMSE_avg_stacked)
    
            # root names for the collumns name
            # -------------------------------------------------------------
            cols_root = [
                "time",
                "ObsType",
                "RMSE" + name_sensor,
                "RMSE_avg",
                "NMRMSE" + name_sensor,
                "NMRMSE_avg",
                "OL",
            ]
    
            # root data
            # -------------------------------------------------------------
            data_df_root = [
                [t_obs],
                [name_sensor],
                [RMSE_sensor_ti],
                [all_Obs_RMSE_avg_ti],
                [NMRMSE_sensor_ti],
                [NMRMSE_avg_ti],
                [OL_bool],
            ]
    
            df_performance_ti = pd.DataFrame(
                np.array(data_df_root, dtype=object).T, columns=cols_root
            )
    
            # concatenate with main RMSE dataframe
            # -------------------------------------------------------------
            self.df_performance = pd.concat(
                [self.df_performance, df_performance_ti], axis=0, ignore_index=True
            )
    
        return self.df_performance
    