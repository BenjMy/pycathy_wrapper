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

        The ensemble file update is controlled by `apply_ensemble_updates()`


        - Evaluate performance


"""
import glob

import multiprocessing
import concurrent.futures
import time

import os
import re
import shutil
import subprocess
import warnings
import copy

import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import pandas as pd
from rich.progress import track

from pyCATHY.cathy_utils import dictObs_2pd
from pyCATHY.cathy_tools import CATHY, subprocess_run_multi
from pyCATHY.DA import enkf, pf, mapper, localisation
from pyCATHY.importers import cathy_inputs as in_CT
from pyCATHY.importers import cathy_outputs as out_CT
# from pyCATHY.importers import sensors_measures as in_meas
import pyCATHY.meshtools as mt
from pyCATHY.plotters import cathy_plots as plt_CT

warnings.filterwarnings("ignore",
                        category=DeprecationWarning
                        )

# from functools import partial

REMOVE_CPU = 5

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


    if DA_type == "enkf_analysis_localized":
        L = kwargs.pop('localisationMatrix')
        Analysis = enkf.enkf_analysis_localized(data,
                               data_cov,
                               param,
                               ensembleX[id_state],
                               prediction,
                               L,  # localization matrix (state x obs)
                               Sakov=True,
                               )
        return Analysis

    if DA_type == "enkf_analysis_localized_with_inflation":
        
        Sakov = kwargs.pop('Sakov', True)
        if Sakov is True:
            warnings.warn("Sakov is True")
            
        L = kwargs.pop('localisationMatrix', None)
        if L is None:
            warnings.warn("L not provided. None as default.")
    
        inflate_states = kwargs.pop('damping', None)
        if inflate_states is None:
            inflate_states = 1.0
            warnings.warn("inflate_states not provided. Using default: 1.0 (no inflation).")
    
        inflate_params = kwargs.pop('inflate_params', None)
        if inflate_params is None:
            inflate_params = 1.0
            warnings.warn("inflate_params not provided. Using default: 1.0 (no inflation).")
    
        jitter_params = kwargs.pop('jitter_params', None)
        if jitter_params is None:
            jitter_params = 0.0
            warnings.warn("jitter_params not provided. Using default: 0.0 (no additive noise).")


        Analysis = enkf.enkf_analysis_localized_with_inflation(data,
                                                               data_cov,
                                                               ensembleX[id_state],
                                                               param,
                                                               prediction,
                                                               L,  # localization matrix (state x obs)
                                                               Sakov=Sakov,
                                                               inflate_states=inflate_states,
                                                               inflate_params=inflate_params,
                                                               jitter_params=jitter_params,
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



#%%
# DA class
# ---------

class DA(CATHY):
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
        self.localisationBool = False
        # self.localisationMatrix = []


        # backup all spatial ET iterations (self.backupfort777)
        # ET file are large!
        self.backupfort777 = False

        # Infer ensemble size NENS from perturbated parameter dictionnary
        # -------------------------------------------------------------------
        for name in self.dict_parm_pert:
            # print(name)
            # print(self.dict_parm_pert[name])
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

    def prepare_DA(
                    self,
                    dict_obs,
                    dict_parm_pert,
                    list_parm2update,
                    list_assimilated_obs,
                    parallel,
                    **kwargs
                    ):

        self.run_processor(DAFLAG=1,
                           runProcess=False,
                           **kwargs) # to create the processor exe
        callexe = "./" + self.processor_name

        self.damping = kwargs.get("damping", 1)
        self.localisation = kwargs.get("localisation", None)
        threshold_rejected = kwargs.get("threshold_rejected", 10)
        self.verbose = kwargs.get("verbose", False)
        # This is particulary useful for ERT data where Archie non-linear operator create a non gaussian data distribution
        self.use_log_transformed_obs = kwargs.get("use_log_transformed_obs", False)
        
        # dict_parm_pert.keys()
        # list_parm2update
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
        self.apply_ensemble_updates(
            dict_parm_pert,
            list_parm2update="all",
            cycle_nb=self.count_DA_cycle,
        )

        all_atmbc_times = np.copy(self.atmbc["time"])

        return (callexe,
                all_atmbc_times,
                threshold_rejected,
                ENS_times,
                update_key)

    def run_DA_smooth(
                        self,
                        DA_type,
                        dict_obs,
                        list_parm2update,
                        dict_parm_pert,
                        list_assimilated_obs,
                        parallel=True,
                        **kwargs
                    ):

        """

        Run Data Assimilation

        .. note::

            Steps are:
            1. DA init (create subfolders)
            2. run CATHY hydrological model stacked using ALL DA times
            3. check before analysis
            4. analysis
            5. check after analysis

        Parameters
        ----------
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

        (callexe,
         self.all_atmbc_times,
         threshold_rejected,
         self.ENS_times,
         update_key) = self.prepare_DA(
                                         dict_obs,
                                         dict_parm_pert,
                                         list_parm2update,
                                         list_assimilated_obs,
                                         parallel,
                                         **kwargs
                                         )

        # self._run_ensemble_hydrological_model(parallel, callexe)
        # os.chdir(os.path.join(self.workdir))

        # map states to observation = apply H operator to state variable
        # ----------------------------------------------------------------
        prediction = self.map_states2Observations(
                                                    list_assimilated_obs,
                                                    default_state="psi",
                                                    parallel=parallel,
                                                )
        # np.shape(prediction)
        # self.atmbc

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

        pass


    def get_mesh_mask_nodes_observed(self,list_assimilated_obs):
        mesh_mask_nodes_observed = np.zeros(int(self.grid3d['nnod3']), dtype=bool)
        obskey2map, obs2map = self._obs_key_select(list_assimilated_obs)
        mesh_nodes_observed = []
        if obskey2map[0] == 'ERT':
            raise ValueError("Cannot localise with apparent resistivity data.")
        elif obskey2map[0] == 'EM':
            nbofcoils =  len(obs2map[0]['coils'])
            for node_nb in obs2map[0]['data']['mesh_node']:
                mesh_mask_nodes_observed[int(node_nb)] = True
                mesh_nodes_observed.append(int(node_nb))
            mesh_mask_nodes_observed = [mesh_mask_nodes_observed]*nbofcoils
            mesh_nodes_observed = [mesh_nodes_observed]*nbofcoils
        else:
            for obsi in obs2map:
                mesh_mask_nodes_observed[obsi['mesh_nodes']] = True
                mesh_nodes_observed.append(obsi['mesh_nodes'])
        return mesh_nodes_observed, mesh_mask_nodes_observed

    def run_DA_sequential(
                            self,
                            DA_type,
                            dict_obs,
                            list_parm2update,
                            dict_parm_pert,
                            list_assimilated_obs,
                            parallel=True,
                            **kwargs,
                        ):

        """

        Run Data Assimilation

        .. note::

            Steps are:
            1. DA init (create subfolders)
            2. run CATHY hydrological model recursively using DA times
                <- Loop-->
                    3. check before analysis
                    4. analysis
                    5. check after analysis
                    6. update input files
                <- Loop-->

        Parameters
        ----------
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

        (callexe,
         self.all_atmbc_times,
         threshold_rejected,
         self.ENS_times,
         update_key) = self.prepare_DA(
                                         dict_obs,
                                         dict_parm_pert,
                                         list_parm2update,
                                         list_assimilated_obs,
                                         parallel,
                                         **kwargs
                                         )


        # -------------------------------------------------------------------
        # (time-windowed) update input files ensemble again
        # -------------------------------------------------------------------
        self._update_input_ensemble(list_parm2update)
        # -----------------------------------
        # Run hydrological model sequentially
        # = Loop over atmbc times (including assimilation observation times)
        # -----------------------------------
        for t_atmbc in self.all_atmbc_times:  # atmbc times MUST include assimilation observation times

            self._run_ensemble_hydrological_model(parallel, callexe)
            os.chdir(os.path.join(self.workdir))

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
            if t_atmbc in self.ENS_times:

                self.console.print(":world_map: [b] map states to observations [/b]")
                #%% map states to observation = apply H operator to state variable
                # ----------------------------------------------------------------
                prediction = self.map_states2Observations(
                                                            list_assimilated_obs,
                                                            default_state="psi",
                                                            parallel=parallel,
                                                        )
                # invalid = (prediction <= 0) | np.isnan(prediction)
                # if np.any(invalid):
                #     print(f"⚠️ Invalid values found in prediction: {np.sum(invalid)} entries invalid")

                # Check for negatives and NaNs
                invalid = (prediction <= 0) | np.isnan(prediction)
                if self.isET:
                    invalid = (prediction < 0) | np.isnan(prediction)

                num_invalid = np.sum(invalid)
                num_nan = np.sum(np.isnan(prediction))
                num_negative = np.sum(prediction < 0)
                
                if np.any(invalid):
                    print(f"⚠️ Invalid prediction values detected: {num_invalid} entries total")
                    print(f"   - NaN values: {num_nan}")
                    print(f"   - Negative or zero values: {num_negative}")
                    raise ValueError("prediction contains invalid values. Please clean or replace them before proceeding.")
                    
    
                #%% DA analysis preparation
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

                #%% DA localisation preparation
                # ----------------------------------------------------------------
                if self.localisation is not None:
                    self.localisationBool = any(item in self.localisation for item in ['veg_map',
                                                                                       'zones',
                                                                                       'nodes']
                                                )

                if self.localisationBool:
                    self.console.print(f":carpentry_saw: [b] Domain localisation: {self.localisation} [/b]")
                    
                    if hasattr(self, 'zone') is False: 
                        zone, _ = self.read_inputs('zone')
                        self.zone = zone
                    #%%
                    (mesh3d_nodes_valid,
                     mask_surface_nodes) = localisation.create_mask_localisation(self.localisation,
                                                                                self.veg_map,
                                                                                self.zone,
                                                                                self.hapin,
                                                                                self.grid3d,
                                                                                self.DEM
                                                                                )
                                                                                                    
                    #%%
                    ensemble_psi = []
                    ensemble_sw = []
                    data = []
                    analysis_with_localisation = []
                    analysis_param = []
                    # Loop over unique map nodes and perform the DA analysis
                    for loci, node_value in enumerate(np.unique(mask_surface_nodes)):
                        if node_value == -9999:
                            continue

                        # Initialize valid observation nodes as all true
                        obs_surface_nodes_bool = np.ones(int(self.grid3d['nnod']), dtype=bool)
                        valid_surface_nodes = mask_surface_nodes.ravel()
                        obs_surface_nodes_bool[valid_surface_nodes != node_value] = False

                        #%%
                        # check position of the sensors to make sure it is located in the right area
                        # --------------------------------------------------------------------------
                        (mesh_nodes_observed,
                        mesh_mask_nodes_observed) = self.get_mesh_mask_nodes_observed(list_assimilated_obs)
                        mesh_nodes_observed = np.hstack(mesh_nodes_observed)
                         
                        #%%                         
                        if 'EM' in list_assimilated_obs:
                            obs_surface_nodes_valid = []
                            for i in range(np.shape(mesh_mask_nodes_observed)[0]):
                                obs_surface_nodes_valid.append(obs_surface_nodes_bool)
                            obs_surface_nodes_valid = np.hstack(obs_surface_nodes_valid)
                        else:
                            obs_surface_nodes_valid = obs_surface_nodes_bool
                            
                        if len(mesh_nodes_observed) <= len(obs_surface_nodes_valid):
                        # else:
                            self.console.print(f"""[b]
                                               Observations nodes **might** not be contained within the surface mesh
                                               Cannot validate its positionning regarding to the localisation mask
                                               ['NOT IMPLEMENTED YET']
                                               [/b]"""
                                               )
                        # if None of the observation are within the localisation area:
                        #     continue without analysis
                        # ------------------------------------------------------------
                        if np.sum(obs_surface_nodes_valid)==0:
                            continue

                        selected_parm2update = list_parm2update
                        if any(item in self.localisation for item in ['veg_map', 'zones']):
                            if len(list_parm2update)>1:
                                selected_parm2update = [list_parm2update[0], list_parm2update[loci + 1]]
                            else:
                                selected_parm2update = [list_parm2update[0]]

                        self.console.print(f"""[b]
                                           SHORTCUT: ALL PARAMS ARE UPDATED but
                                           are not necessaraly associated
                                           with observation localisation [/b]
                                           """)
                        self.console.print(f"[b] Params to update (domain loc): {selected_parm2update} [/b]")
                        
                        #%%
                        # Perform DA analysis for each mask of localisation
                        # --------------------------------------------------
                        (
                            ensemble_psi_loci,
                            ensemble_sw_loci,
                            data_loci,
                            analysis_loci,
                            analysis_param_loci,
                        ) = self._DA_analysis(
                            prediction,
                            DA_type,
                            selected_parm2update,
                            list_assimilated_obs,
                            ens_valid=self.ens_valid,
                            obs_valid=obs_surface_nodes_valid,
                            obs_mesh_nodes_valid= mesh3d_nodes_valid[loci]
                        )
                        #%%

                        data.append(data_loci)
                        analysis_with_localisation.append(np.c_[mesh3d_nodes_valid[loci],analysis_loci])
                        analysis_param.append(analysis_param_loci)

                    self.console.print(f"""[b] ONLY TAKE THE first part of the data? 
                                       Need to refactor performance assesement for DA with localisation 
                                       self._performance_assessement_pq( [/b]
                                        """)
                    data = data[0]
                    ensemble_psi = ensemble_psi_loci
                    ensemble_sw = ensemble_sw_loci
                    analysis_with_localisation = np.vstack(analysis_with_localisation)
                    analysis_with_localisation = analysis_with_localisation[analysis_with_localisation[:, 0].argsort()]
                    analysis_param = np.hstack(analysis_param)

                    if len(analysis_with_localisation)!=len(ensemble_psi):
                        self.console.print(f"""
                                               [b] !! All the mesh nodes havent been analysed with DA !! [/b]
                                               Taking state values where missing analysis nodes
                                           """)
                        (
                            ensemble_psi_valid,
                            ensemble_sw_valid,
                            ens_size,
                            sim_size,
                        ) = self._read_state_ensemble()
                        analysis = ensemble_psi_valid
                        indice_nodes_with_analysis = [int(ii) for ii in analysis_with_localisation[:,0]]
                        analysis[indice_nodes_with_analysis]=analysis_with_localisation[:,1:]
                    else:
                        analysis=analysis_with_localisation[:,1:]

                else:
                    # No localisation, consider all nodes as valid
                    # ---------------------------------------------------------
                    obs_valid = np.ones(prediction.shape[0], dtype=bool)
                    mesh3d_nodes_valid =  np.ones(int(self.grid3d['nnod3']), dtype=bool)

                    # Perform DA analysis without localisation
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
                        obs_valid=obs_valid,
                        obs_mesh_nodes_valid= mesh3d_nodes_valid
                    )


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

                #%% check analysis quality
                # ----------------------------------------------------------------

                self.console.print(
                    ":face_with_monocle: [b]check analysis performance[/b]"
                )
                self._performance_assessement_pq(
                                                list_assimilated_obs,
                                                data,
                                                prediction_valid,
                                                t_obs=self.count_DA_cycle,
                                            )
                #%% the counter is incremented here
                # ----------------------------------------------------------------
                self.count_DA_cycle = self.count_DA_cycle + 1

                # # update parameter dictionnary
                # ----------------------------------------------------------------
                def check_ensemble(ensemble_psi_valid, ensemble_sw_valid):
                    if np.any(ensemble_psi_valid > 0):
                        print("!positive pressure heads observed!")

                check_ensemble(ensemble_psi_valid,
                               ensemble_sw_valid
                               )

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
                print("!shortcut here ensemble are not validated!")

            self.count_atmbc_cycle = self.count_atmbc_cycle + 1

            self.console.rule(
                ":round_pushpin: end of time step (s) "
                + str(int(t_atmbc))
                + "/"
                + str(int(self.all_atmbc_times[-1]))
                + " :round_pushpin:",
                style="yellow",
            )
            self.console.rule(
                ":round_pushpin: end of atmbc update # "
                + str(self.count_atmbc_cycle)
                + "/"
                + str(len(self.all_atmbc_times))
                + " :round_pushpin:",
                style="yellow",
            )

            # create dataframe _DA_var_pert_df holding the results of the DA update
            # ---------------------------------------------------------------------
            if len(ensemble_psi_valid_bef_update)==0:
                (
                    ensemble_psi_valid_bef_update,
                    ensemble_sw_valid_bef_update,
                    ens_size,
                    sim_size,
                ) = self._read_state_ensemble()
            # else:
            self.console.rule("Back up DA step")
            self._DA_df(
                state=[ensemble_psi_valid_bef_update,
                       ensemble_sw_valid_bef_update
                       ],
                state_analysis=analysis_valid,
                rejected_ens=rejected_ens,
            )

            self._DA_ET_xr()

            # export summary results of DA
            # ----------------------------------------------------------------
            meta_DA = {
                "listAssimilatedObs": list_assimilated_obs,
                "listUpdatedparm": list_parm2update,
            }

            # backup ET map for each iteration of ensemble numbers
            if self.backupfort777:
                self.backup_fort777()
            self.backup_results_DA(meta_DA)
            self.backup_simu()

            print(f'Parameters to update are: {list_parm2update}')
            # overwrite input files ensemble (perturbated variables)
            # ---------------------------------------------------------------------
            if (self.count_atmbc_cycle < (len(self.all_atmbc_times) - 1)):  # -1 cause all_atmbc_times include TMAX
                self._update_input_ensemble(
                                                list_parm2update,
                                                analysis=analysis_valid
                )
            else:
                self.console.rule(
                    ":red_circle: end of DA" ":red_circle:", style="yellow"
                )
                pass
        
            print(list_parm2update)

            self.console.rule(
                                ":red_circle: end of DA update "
                                + str(self.count_DA_cycle)
                                + "/"
                                + str(len(self.ENS_times))
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


    def _DA_analysis(
        self,
        prediction,
        DA_type="enkf_Evensen2009_Sakov",
        list_update_parm=["St. var."],
        list_assimilated_obs="all",
        ens_valid=[],
        obs_valid=[],
        obs_mesh_nodes_valid=[],
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

        #%% find the method to map for this time step
        # --------------------------------------------------------
        obskey2map, obs2map = self._obs_key_select(list_assimilated_obs)
        
        #%% prepare states
        # ---------------------------------------------------------------------
        ensemble_psi, ensemble_sw, ens_size, sim_size = self._read_state_ensemble()
        #%%
        # ---------------------------------------------------------------------
        # prepare data
        # ---------------------------------------------------------------------
        # select ONLY data to assimilate
        # ---------------------------------------------------------------------
        print(f'Assimilated Observations are: {list_assimilated_obs}')
        data, _ = self._get_data2assimilate(list_assimilated_obs)
        # data_test = np.hstack(data).T
        
        #%%
        self.console.print(
            r""":
                               - Data size: {}
                                 --> Observations --> {}
                           """.format(
                str(np.shape(data)), list_assimilated_obs
            )
        )
        # np.shape(data)
        # check size of data_cov
        # ---------------------------------------------------------------------
        expected = len(data)
        actual_shape = np.shape(self.stacked_data_cov[self.count_DA_cycle])
        if actual_shape[0] != expected:
            raise ValueError(
                f"Wrong stacked_data_cov shape: expected ({expected}, {expected}), "
                f"but got {actual_shape}. "
                "Consider re-computing data covariance."
            )

        data_cov = self.stacked_data_cov[self.count_DA_cycle]
        
        if hasattr(self,'localisationMatrix')>0:
            localisationMatrix = self.localisationMatrix[self.count_DA_cycle]
        else:
            localisationMatrix = []
            
        # # filter non-valid ensemble before Analysis
        # # ---------------------------------------------------------------------
        prediction_valid = prediction[:] # already been filtered out from bad ensemble
        ensemble_psi_valid = ensemble_psi[:, ens_valid]
        ensemble_sw_valid = ensemble_sw[:, ens_valid]

        # prepare parameters
        # ---------------------------------------------------------------------
        # When updating the states only, the elements of X are the
        # pressure heads at each node of the finite element grid, while
        # the state augmentation technique is used when also updat-
        # ing the parameters
        # list_update_parm = ['St. var.', 'ZROOT0']

        param = self.transform_parameters(list_update_parm, update_key)
        print("parm size: " + str(len(param)))

        param_valid = []
        if len(param) > 0:
            param_valid = param[:, ens_valid]

        # # filter non-valid observations before Analysis
        # # ---------------------------------------------------------------------
        # This is implemented particulary for DA with localisation

        data_valid = data[obs_valid]
        data_cov_valid = data_cov[obs_valid,:][:,obs_valid]
        prediction_valid = prediction[obs_valid, :]
        ensemble_psi_valid = ensemble_psi_valid[obs_mesh_nodes_valid, :]
        ensemble_sw_valid = ensemble_sw_valid[obs_mesh_nodes_valid, :]
        # ensemble_sw_valid.min()
        # ensemble_sw_valid.max()
        
        def check_cov_consistency(data_valid, 
                                  prediction_valid, 
                                  data_cov_valid, 
                                  use_log=False, 
                                  threshold_log=1e3, 
                                  threshold_linear=1e4):
            # Check data and prediction ranges
            data_range = data_valid.max() - data_valid.min()
            prediction_range = prediction_valid.max() - prediction_valid.min()
            cov_diag = np.diag(data_cov_valid)
        
            print(f"\n{'LOG' if use_log else 'LINEAR'} domain check:")
            print(f"  Observation range: {data_valid.min():.2e} – {data_valid.max():.2e}")
            print(f"  Prediction range:  {prediction_valid.min():.2e} – {prediction_valid.max():.2e}")
            print(f"  Cov diag range:    {cov_diag.min():.2e} – {cov_diag.max():.2e}")
        
            # Raise warning if ranges and covariance seem inconsistent
            # if use_log:
            #     if cov_diag.max() > threshold_log:
            #         raise ValueError("⚠️ Covariance matrix values are too large for log-space. Did you forget to transform R accordingly?")
            # else:
            #     if cov_diag.max() < threshold_linear:
            #         raise ValueError("⚠️ Covariance matrix values seem too small for linear domain. Are you missing a log-transform in data or model?")
        
        # Example: conditional log transform
        # use_log = True
        if self.use_log_transformed_obs:

            eps = 1e-3  
            # Before log-transform
            print("Before log-transform:")
            print(f"Observed ERa: min = {data_valid.min():.2f}, max = {data_valid.max():.2f}")
            print(f"Predicted ERa: min = {prediction_valid.min():.2f}, max = {prediction_valid.max():.2f}")
            
            # Example: predicted_era is a NumPy array
            invalid = (prediction_valid <= 0) | np.isnan(prediction_valid)
            if np.any(invalid):
                print(f'''
                      ⚠️ Invalid values found in prediction valid (before log-transform): 
                          {np.sum(prediction_valid <= 0)} entries ≤ 0 
                          {np.sum(np.isnan(prediction_valid))}  nan"
                        '''
                    )


            # Apply log10 transform
            data_valid_log = np.log10(data_valid + eps)
            prediction_valid_log = np.log10(prediction_valid + eps)
            # np.log1p(prediction_valid) 
            
            
            # After log-transform
            print("\nAfter log-transform:")
            print(f"Observed log10(ERa): min = {data_valid_log.min():.2f}, max = {data_valid_log.max():.2f}")
            print(f"Predicted log10(ERa): min = {prediction_valid_log.min():.2f}, max = {prediction_valid_log.max():.2f}")
    

            # Run check
            check_cov_consistency(data_valid, prediction_valid, 
                                  data_cov_valid, 
                                  use_log=self.use_log_transformed_obs
                                  )

            fig, ax = plt.subplots()
            ax.hist(np.log10(data_valid + 1e-3), bins=30, alpha=0.5, 
                    label='Observed',density=True)
            ax.hist(np.log10(prediction_valid + 1e-3).flatten(), bins=30, alpha=0.5, 
                    label='Predicted',density=True)
            ax.legend(); 
            plt.title("log10(ERa): Observed vs Predicted")
            
            savename=os.path.join(
                self.workdir,
                self.project_name,
                "DA_Ensemble",
                "log10_Obs_Pred_hist"  + str(self.count_DA_cycle),
            )
                        
            fig.savefig(savename + ".png", dpi=300)


            # run Analysis using log transform
            # ---------------------------------------------------------------------
            result_analysis = run_analysis(
                                            DA_type,
                                            data_valid_log,
                                            data_cov_valid,
                                            param_valid,
                                            list_update_parm,
                                            [ensemble_psi_valid, ensemble_sw_valid],
                                            prediction_valid_log,
                                            alpha=self.damping,
                                            localisationMatrix =localisationMatrix
                                        )
            
            
            # def run_analysis(
            #     DA_type,
            #     data,
            #     data_cov,
            #     param,
            #     list_update_parm,
            #     ensembleX,
            #     prediction,
            #     default_state="psi",
            #     **kwargs,
            # ):
                
                

        else:
            
                     
            savename=os.path.join(
                self.workdir,
                self.project_name,
                "DA_Ensemble",
                "Obs_Pred_hist" + str(self.count_DA_cycle),
            )
            fig, ax = plt.subplots()
            ax.boxplot(
                [data_valid.flatten(), prediction_valid.flatten()],
                labels=['Observed', 'Predicted'],
                showfliers=True,  # Show outliers
                notch=True,       # Optional: adds a notch for median confidence
                patch_artist=True # Optional: for filled boxes
            )
            ax.set_title("Observed vs Predicted")
            fig.savefig(savename + ".png", dpi=300)

            # run Analysis
            # ---------------------------------------------------------------------
            result_analysis = run_analysis(
                                            DA_type,
                                            data_valid,
                                            data_cov_valid,
                                            param_valid,
                                            list_update_parm,
                                            [ensemble_psi_valid, ensemble_sw_valid],
                                            prediction_valid,
                                            alpha=self.damping,
                                            localisationMatrix =localisationMatrix
                                        )
                
        #%%
        _ , obs2eval_key = self._get_data2assimilate(list_assimilated_obs)

        #%%


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

            self.console.print("[b]Plotting COV matrices[/b]")
            try:
                if np.shape(COV)[0]<1e9:

                    plt_CT.show_DA_process_ens(
                        ensemble_psi_valid,
                        data_valid,
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
                        label_sensor=str([*obskey2map]),
                    )
                else:
                    self.console.print("[b]Impossible to plot matrices[/b]")
            except:
                self.console.print("[b]Impossible to plot matrices[/b]")


        else:
            [analysis, analysis_param] = result_analysis

        # Back-transformation of the parameters
        # ---------------------------------------------------------------------
        self.transform_parameters(
            list_update_parm, param=np.hstack(analysis_param), back=True
        )
        # return ensemble_psi, Analysis, AnalysisParam, Data, data_cov
        # ---------------------------------------------------------------------
        return (ensemble_psi, ensemble_sw, data_valid, analysis, analysis_param)

    def _mark_invalid_ensemble(
        self, ens_valid, prediction, ensemble_psi, ensemble_sw, analysis, analysis_param
    ):
        """mark invalid ensemble - invalid ensemble are filled with NaN values"""

        ensemble_psi_valid = np.empty(ensemble_psi.shape)
        ensemble_psi_valid[:] = np.nan
        ensemble_psi_valid[:, ens_valid] = ensemble_psi[:, ens_valid]

        ensemble_sw_valid = np.empty(ensemble_psi.shape)
        ensemble_sw_valid[:] = np.nan
        ensemble_sw_valid[:, ens_valid] = ensemble_sw[:, ens_valid]

        analysis_valid = np.empty(ensemble_psi.shape)
        # analysis_valid = np.empty(analysis.shape)
        analysis_valid[:] = np.nan
        analysis_valid[:, ens_valid] = analysis

        prediction_valid = np.empty([prediction.shape[0], ensemble_psi.shape[1]])
        prediction_valid[:] = np.nan
        prediction_valid[:, ens_valid] = prediction

        analysis_param_valid = []
        if len(analysis_param[0]) > 0:
            analysis_param_valid = np.empty(
                [analysis_param.shape[1], ensemble_psi.shape[1]]
            )
            analysis_param_valid[:] = np.nan
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


    def _mapping_petro_init(self):
        """
        Initiate Archie and VGP petro/pedophysical parameters

        If Archie and VGP dictionnaries are not existing fill with default values
        """

        warnings_petro = []
        # check that porosity is a list
        # -------------------------------------
        
        if hasattr(self, 'soil_SPP') is False:
            self.read_inputs('soil')
            
        porosity = self.soil_SPP["SPP"]['POROS'].mean()
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


    def _update_soil_params_from_perturbations(
        self,
        dict_parm_pert,
        df_SPP_2fill,
        key_root,
        key,
        update_key,
        ens_nb,
        shellprint_update
    ):
        """
        Update the ensemble-specific soil parameter file (df_SPP_2fill) using 
        a dictionary of perturbed parameters.
    
        This method allows flexible updating of soil parameters for different
        ensemble members based on control strategies such as zones, layers, or
        vegetation-root mapping. The soil parameters are updated with the values
        from a perturbed parameter dictionary, typically generated for Ensemble
        Kalman Filter (EnKF) or Monte Carlo simulations.
    
        Parameters
        ----------
        dict_parm_pert : dict
            Dictionary containing perturbation information and sampled parameter values.
            It must include a 'pert_control_name' key (e.g., 'zone', 'layers', 'root_map'),
            and a nested dictionary for each `update_key` and ensemble member.
    
        df_SPP_2fill : list of pandas.DataFrame
            A list of DataFrames, each representing the soil parameter file for one ensemble.
            If the list is empty, it is initialized from soil input files.
    
        key_root : list or tuple of str and int
            Root identifier of the parameter being updated, e.g., ['ks', 0] or ['porosity', 1].
    
        key : str
            Full parameter key in the dictionary, used to locate the corresponding entry in
            `dict_parm_pert`.
    
        update_key : str
            Name of the key in `dict_parm_pert[key]` used to fetch perturbed values,
            typically something like 'samples' or 'transformed'.
    
        ens_nb : int
            Ensemble member index (zero-based) to update in `df_SPP_2fill`.
    
        shellprint_update : function
            A logging or print function (e.g., `print()` or custom shell logger), currently unused.
    
        Raises
        ------
        ValueError
            If the parameter name is unknown or if the required index is missing
            for zone/layer/root_map updates.
    
        Notes
        -----
        - This function assumes that `in_CT.read_soil(...)` loads soil input files
          from the CATHY model directory structure.
        - All updated columns must exist in the soil DataFrame.
        - The initialization block assumes perturbations are defined per layer or zone
          (e.g., ks0 = layer 0).
        """
        
        pert_control_name = dict_parm_pert[key]['pert_control_name']
        update_value = dict_parm_pert[key][update_key][ens_nb]
        param_name = key_root[0]
        index = int(key_root[1]) if len(key_root) > 1 else None
    
        # Define which columns to update per parameter type
        param_map = {
            'Ks': ['PERMX', 'PERMY', 'PERMZ'],
            'porosity': ['POROS'],
            'n_VG': ['VGNCELL'],
            'thetar_VG': ['VGRMCCELL'],
            'SATCELL_VG': ['VGPSATCELL'],
        }
        update_cols = param_map.get(param_name)
        if update_cols is None:
            raise ValueError(f"Unknown parameter '{param_name}'.")
            
        # Initialize if needed
        if len(df_SPP_2fill) == 0:
            for i in range(self.NENS):
                df_SPP, _ = in_CT.read_soil(
                    os.path.join(
                        self.workdir,
                        self.project_name,
                        f'DA_Ensemble/cathy_{i+1}',
                        'input/soil'
                    ),
                    self.dem_parameters,
                    self.cathyH["MAXVEG"]
                )
                df = df_SPP.copy()
                # Pre-fill the relevant columns with NaN
                df[update_cols] = np.nan
                df_SPP_2fill.append(df)
    
            print(
                "\nAssuming that the dictionary of perturbed parameters is built "
                "using 'pert_control_name', i.e. per layer or per zone like ks0, ks1, etc...\n"
            )
        
        # Apply the update depending on control type
        df = df_SPP_2fill[ens_nb]
        if pert_control_name == 'zone' and index is not None:
            df.loc[(index + 1, slice(None)), update_cols] = update_value
        elif pert_control_name == 'layers' and index is not None:
            df.loc[(slice(None), index), update_cols] = update_value
        elif pert_control_name == 'root_map' and index is not None:
            mask_vegi = (index + 1 == self.veg_map)
            zones_in_vegi = self.zone[mask_vegi]
            df.loc[zones_in_vegi, update_cols] = update_value
        else:
            df.loc[:, update_cols] = update_value

        #%% check if contains nan (if not update soil file)
        contains_nan = df_SPP_2fill[ens_nb].isna().any().any()
        if contains_nan==False:
            self.update_soil(
                SPP_map=df_SPP_2fill[ens_nb],
                FP_map=self.soil_FP["FP_map"],
                verbose=self.verbose,
                filename=os.path.join(os.getcwd(), "input/soil"),
                shellprint_update=shellprint_update,
                backup=True,
            )

            saveMeshPath = os.path.join(
                                        os.getcwd(),
                                        "vtk/",
                                        self.project_name + ".vtk"
                                        )
            if pert_control_name == 'layers':
                (
                POROS_mesh_cells,
                 POROS_mesh_nodes
                 ) = self.map_prop_2mesh_markers('POROS',
                                                df_SPP_2fill[ens_nb]['POROS'].groupby('layer').mean().to_list(),
                                                to_nodes=False,
                                                saveMeshPath=saveMeshPath,
                                                # zones_markers_3d=None,
                                                )
                self.map_prop2mesh({"POROS": POROS_mesh_nodes})


        return df_SPP_2fill

    def apply_ensemble_updates(self, dict_parm_pert, list_parm2update, **kwargs):
        """
        Applies updates to ensemble files by overwriting input files.
        
        Updates may include:
        - State variables (e.g., pressure head from data assimilation)
        - Model parameters: Ks, porosity, Feddes, Archie, initial conditions, atmbc
        - Boundary and initial conditions
        
        Parameters
        ----------
        dict_parm_pert : dict
            Dictionary of perturbation data per parameter.
        list_parm2update : list or "all"
            List of parameters to update. Use "all" to update everything.
        **kwargs : dict
            cycle_nb: int, optional
                Data assimilation cycle number (0 = prior).
            analysis: np.ndarray, optional
                State analysis results to update the pressure head.
        
        Returns
        -------
        None
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
                # dict_parm_pert.keys()
        # loop over ensemble files to update ic -->
        # this is always done since psi are state variable to update
        # ------------------------------------------------------------------
        shellprint_update = True
        for ens_nb in self.ens_valid:
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
        Feddes_withPertParam = []  # matrice of ensemble for Feddes parameters
        VG_parms_mat_ens = []  # matrice of ensemble for VG parameters
        VG_p_possible_names = ["n_VG", "thetar_VG", "alpha_VG", "SATCELL_VG"]
        VG_p_possible_names_positions_in_soil_table = [5, 6, 7, 7]

        # SPP_possible_names = ["porosity", "ks"]
        ks_het_ens = []  # matrice of ensemble for heterogeneous SPP parameters
        POROS_het_ens = []  # matrice of ensemble for heterogeneous SPP parameters
        VG_p_het_ens = []
        
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

        # Feddes
        # -------
        Feddes_p_names = [
            "PCANA",
            "PCREF",
            "PCWLT",
            "ZROOT",
            "PZ",
            "OMGC",
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
                        times = dict_parm_pert[key]['data2perturbate']['time']
                        VALUE = dict_parm_pert[key]['time_variable_perturbation'][ens_nb]
                        self.update_atmbc(
                                            HSPATM=self.atmbc['HSPATM'],
                                            IETO=self.atmbc['IETO'],
                                            time=list(times),
                                            netValue=VALUE,
                                            filename=os.path.join(os.getcwd(), "input/atmbc"),
                                            verbose=self.verbose,
                                            # show=True
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
                    
                    if dict_parm_pert[key]['pert_control_name'] == 'layers':
                        # print('test')
                    
                        
                    # if match_ic_withLayers:
                        DApath = os.getcwd()
                        saveMeshPath = os.path.join(DApath, "vtk",
                                                     self.project_name + '.vtk'
                                                     )
                        ic_values_by_layers = dict_parm_pert[key][f'ini_{key}_withLayers'][:,ens_nb]
                        ic_nodes = self.map_prop_2mesh_markers('ic',
                                                                ic_values_by_layers,
                                                                to_nodes=True,
                                                                saveMeshPath=saveMeshPath,
                                                                )
                        self.update_ic(
                            INDP=1,
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
                    ks_het_ens = self._update_soil_params_from_perturbations(
                                                       dict_parm_pert,
                                                       ks_het_ens,
                                                       key_root,
                                                       key,
                                                       update_key,
                                                       ens_nb,
                                                       shellprint_update
                                                       )


                elif key_root[0].casefold() in 'porosity':
                    POROS_het_ens = self._update_soil_params_from_perturbations(
                                                       dict_parm_pert,
                                                       POROS_het_ens,
                                                       key_root,
                                                       key,
                                                       update_key,
                                                       ens_nb,
                                                       shellprint_update
                                                       )

                # VG parameters update
                # --------------------------------------------------------------
                elif key_root[0] in VG_p_possible_names:
                    VG_p_het_ens = self._update_soil_params_from_perturbations(
                                                       dict_parm_pert,
                                                       VG_p_het_ens,
                                                       key_root,
                                                       key,
                                                       update_key,
                                                       ens_nb,
                                                       shellprint_update
                                                       )
                # FeddesParam update
                # --------------------------------------------------------------
                elif key_root[0] in Feddes_p_names:
                    # SPP_map = self.soil_SPP["SPP_map"]

                    df_SPP, _ = in_CT.read_soil(
                        os.path.join(
                            self.workdir,
                            self.project_name,
                            f'DA_Ensemble/cathy_{ens_nb+1}',
                            'input/soil'
                        ),
                        self.dem_parameters,
                        self.cathyH["MAXVEG"]
                    )
                    # df = df_SPP.copy()
  
                    if len(Feddes_withPertParam) == 0:
                        Feddes =  self.soil_FP["FP_map"]
                        Feddes_ENS = [Feddes]*self.NENS
                        Feddes_withPertParam = [copy.deepcopy(feddes) for feddes in Feddes_ENS]

                    # print(key,update_key)
                    new_Feddes_parm = dict_parm_pert[key][update_key]
                    Feddes_withPertParam[ens_nb][key_root[0]].iloc[int(key_root[1])] = new_Feddes_parm[ens_nb]
                    # solution to check to avoid warning: Feddes_withPertParam[ens_nb][key_root[0]].loc[int(key_root[1])] = new_Feddes_parm[ens_nb]

                    self.update_soil(
                        FP_map=Feddes_withPertParam[ens_nb],
                        SPP_map=df_SPP,
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

                    idArchie = Archie_p_names.index(key_root[0])

                    # pert_control_name = dict_parm_pert[key_root[0]]['pert_control_name']
                    # update_value = dict_parm_pert[key_root[0]][update_key][ens_nb]

                    if len(Archie_parms_mat_ens) == 0:
                        for p in Archie_p_names:
                            if len(self.Archie_parms[p]) != self.NENS:
                                self.Archie_parms[p] = list(
                                    self.Archie_parms[p] * np.ones(self.NENS)
                                )
                        Archie_parms_mat_ens = np.zeros([len(Archie_p_names), self.NENS])
                        Archie_parms_mat_ens[:] = np.nan
                   
                    Archie_parms_mat_ens[idArchie, ens_nb] = dict_parm_pert[key][update_key][ens_nb]

                    if np.isnan(Archie_parms_mat_ens[idArchie, :]).any() == False:
                        self.Archie_parms[key_root[0]] = list(Archie_parms_mat_ens[idArchie, :])
                        print(f'New Archie parm: self.Archie_parms[key_root[0]]')

                else:
                    # variable perturbated to ommit: state
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
        print(f'Assimilated Observations are: {list_assimilated_obs}')

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
            elif "EM" in sensor:
                data2add.append(
                    np.hstack(items_dict[time_ass][1][sensor]["data"][obs2map[0]['coils']].to_numpy())
                    # items_dict[time_ass][1][sensor]["data"][obs2map[0]['coils']].to_numpy()
                )
            else:
                data2add.append(items_dict[time_ass][1][sensor]["data"])

            return data2add

        if time_ass is None:
            # time_ass = self.count_DA_cycle + 1  # + 1 since we compare the predicted observation at ti +1
            print('!!!!Take care here not sure I am taking the right time!!!')
            time_ass = self.count_DA_cycle   # + 1 since we compare the predicted observation at ti +1

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

        return np.squeeze(data), obskey2map

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


    def _check_prediction_shape(self,Hx_ens,obskey2map):
        if isinstance(Hx_ens, list):
            try:
                Hx_ens = np.array(Hx_ens)
            except:
                pass
        if not isinstance(Hx_ens, float):
            Hx_ens_shape = np.shape(Hx_ens)
            expected_shape = (len(obskey2map),len(self.ens_valid))
            if Hx_ens_shape != expected_shape:
                if Hx_ens_shape == (len(self.ens_valid), 1, len(obskey2map)):
                    Hx_ens = np.vstack(Hx_ens).T
                elif Hx_ens_shape[0] == expected_shape[1]:
                    Hx_ens = np.array(Hx_ens).T

        return Hx_ens


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
        
        write2shell_map = True
        if "write2shell_map" in kwargs:
            write2shell_map = kwargs["write2shell_map"]

        Hx_ens = []  # matrice of predicted observation for each ensemble realisation
        for ens_nb in self.ens_valid:  # loop over ensemble files
        
            path_fwd_CATHY_i = os.path.join(
                self.workdir,
                self.project_name,
                f"DA_Ensemble/cathy_{ens_nb + 1}",
            )
                
            SPP_ensi, _ = in_CT.read_soil(
                os.path.join(path_fwd_CATHY_i, self.input_dirname, 'soil'),
                self.dem_parameters,
                self.cathyH["MAXVEG"]
            )
                
            POROS_mesh_cell_ensi, POROS_mesh_nodes_ensi = self.map_prop_2mesh_markers(
                          'POROS',
                          SPP_ensi['POROS'].to_list(),
                          to_nodes=False,
              )
              
            # print('ens nb:' + str(ens_nb))

            path_fwd_CATHY = os.path.join(
                self.workdir, self.project_name, "DA_Ensemble/cathy_" + str(ens_nb + 1)
            )
            # infer soil parameters properties
            # ---------------------------------
            # porosity = self.soil_SPP["SPP"][:, 4][0]
            # SPP_ensi, _ = in_CT.read_soil(os.path.join(path_fwd_CATHY,
            #                                            self.input_dirname,'soil'
            #                                            ),
            #                               self.dem_parameters,
            #                               self.cathyH["MAXVEG"]
            #                               )
            
            # find data to map with dictionnary of observations
            # --------------------------------------------
            obskey2map, obs2map = self._obs_key_select(list_assimilated_obs)

            if obskey2map != 'ET':
                df_psi = self.read_outputs(
                    filename="psi", path=os.path.join(path_fwd_CATHY, self.output_dirname)
                )
                df_sw, _ = self.read_outputs(
                    filename="sw", path=os.path.join(path_fwd_CATHY, self.output_dirname)
                )
                state = [df_psi.iloc[-1], df_sw.iloc[-1].values]

            Hx_stacked = []  # stacked predicted observation
            self.isET = any('ETact' in item for item in obskey2map)
            if self.isET:
                df_fort777 = out_CT.read_fort777(os.path.join(path_fwd_CATHY,
                                                              'fort.777'),
                                                  )
                df_fort777 = df_fort777.set_index('time_sec')
                t_ET = df_fort777.index.unique()
                df_fort777.loc[t_ET[-1]]['ACT. ETRA']
                # obs2map
                # dictObs_2pd(obs2map)

                obs2map_node=[]
                for obsi in obs2map:
                    if obsi['data_type']=='ETact':
                        obs2map_node.append(obsi['mesh_nodes'])
                ETobs_selected = df_fort777.loc[t_ET[-1]]['ACT. ETRA'].iloc[obs2map_node]
                Hx_stacked.append(ETobs_selected.values)

            # Loop over observations to map
            # ---------------------------------------------------------------------
            for i, obs_key in enumerate(obskey2map):

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

                if "stemflow" in obs_key:
                    # Atmact-vf(13) : Actual infiltration (+ve) or exfiltration (-ve) at atmospheric BC nodes as a volumetric flux [L^3/T]
                    # Atmact-v (14) : Actual infiltration (+ve) or exfiltration (-ve) volume [L^3]
                    # Atmact-r (15) : Actual infiltration (+ve) or exfiltration (-ve) rate [L/T]
                    # Atmact-d (16) : Actual infiltration (+ve) or exfiltration (-ve) depth [L]
                    print("Not yet implemented")

                if "scale" in obs_key:
                    print("Not yet implemented")
                    df_dtcoupling = self.read_outputs(
                        filename="dtcoupling",
                        path=os.path.join(path_fwd_CATHY, self.output_dirname),
                    )
                if "discharge" in obs_key:
                    # case 3: discharge
                    # need to read the hgsfdet file (Hx_Q)
                    # --------------------------------------------------------------------
                    print("Not yet implemented")

                if "EM" in obs_key:
                    # Electromagnetic surveys
                    # --------------------------------------------------------------------
                    Hx_EM, df_Archie = mapper._map_EM(
                                                       self.dict_obs,
                                                       self.project_name,
                                                       self.Archie_parms,
                                                       self.count_DA_cycle,
                                                       path_fwd_CATHY,
                                                       ens_nb,
                                                       POROS_mesh_nodes_ensi,
                                                       self.grid3d,
                                                       # savefig=True,
                                                       verbose=True,
                                                    )
                    Hx_stacked.append(Hx_EM)

                    
            if len(Hx_stacked) > 0:
                Hx_stacked_all_obs = np.hstack(Hx_stacked)
                Hx_ens.append(Hx_stacked_all_obs)

        Hx_ens = self._check_prediction_shape(Hx_ens,
                                              obskey2map
                                              )

        #%%
        # ---------------------------------------------------------------------
        # special case of ERT // during sequential assimilation
        # ---------------------------------------------------------------------
        
        if parallel:
            if "ERT" in obs_key:
                path_fwd_CATHY_list = []  # list of ensemble paths
                POROS_mesh_nodes_ensi_list = np.empty(self.NENS, dtype=object)  # preallocate object array
                # POROS_mesh_nodes_ensi_list = []
                print(self.ens_valid)
                
                for ens_nb in self.ens_valid:  # loop over ensemble files
                    # Build path to ensemble
                    path_fwd_CATHY_i = os.path.join(
                        self.workdir,
                        self.project_name,
                        f"DA_Ensemble/cathy_{ens_nb + 1}",
                    )
                    path_fwd_CATHY_list.append(path_fwd_CATHY_i)
            
                    # Read soil properties
                    SPP_ensi, _ = in_CT.read_soil(
                        os.path.join(path_fwd_CATHY_i, self.input_dirname, 'soil'),
                        self.dem_parameters,
                        self.cathyH["MAXVEG"]
                    )
            
                    # Map porosity to mesh nodes
                    POROS_mesh_cell_ensi, POROS_mesh_nodes_ensi = self.map_prop_2mesh_markers(
                        'POROS',
                        # SPP_ensi['POROS'].groupby('zone').mean().to_list(),
                        SPP_ensi['POROS'].to_list(),
                        to_nodes=False,
                    )
            
                    # Store porosity array for this ensemble
                    POROS_mesh_nodes_ensi_list[ens_nb] = POROS_mesh_nodes_ensi
                    # POROS_mesh_nodes_ensi_list.append(POROS_mesh_nodes_ensi)
                    
                    # if POROS_mesh_nodes_ensi is None or len(POROS_mesh_nodes_ensi) == 0:
                    #     print("Warning: first ensemble porosity is empty or None")
                    # else:
                    #     print("Nb of porosity used in Archie in the ensemble:", len(np.unique(POROS_mesh_nodes_ensi)))
        
                # Optional: check the first ensemble
                print('Nb of porosity used in Archie (first ensemble):', len(np.unique(POROS_mesh_nodes_ensi_list[0])))


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

                # self.mesh_pv_attributes['POROS']
                # self.mesh_pv_attributes["node_markers_layers"]
                # SPP_ensi['POROS'].groupby('layer').mean().to_list()
                # SPP_ensi['POROS'].groupby('zone').mean().to_list()
                print("Shape of POROS_mesh_nodes_ensi_list:", np.shape(POROS_mesh_nodes_ensi_list))
                print("Shape of path_fwd_CATHY_list:", np.shape(path_fwd_CATHY_list))

            
                if parallel:
                    Hx_ens_ERT, df_Archie, mesh2test = mapper._map_ERT_parallel(
                                                                    self.dict_obs,
                                                                    self.count_DA_cycle,
                                                                    self.project_name,
                                                                    self.Archie_parms,
                                                                    POROS_mesh_nodes_ensi_list,
                                                                    path_fwd_CATHY_list,
                                                                    list_assimilated_obs="all",
                                                                    default_state="psi",
                                                                    DA_cnb=self.count_DA_cycle,
                                                                    savefig=True,
                                                                    verbose=True,
                                                                    )
                    # mesh2test.array_names

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
        self, list_parm2update=["St. var."], analysis=[], **kwargs
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
        self.selec_atmbc_window()
        # ---------------------------------------------------------------------
        self.apply_ensemble_updates(
                                self.dict_parm_pert,
                                list_parm2update,
                                cycle_nb=self.count_DA_cycle,
                                analysis=analysis,
                            )

        pass

    def selec_atmbc_window(self):
        """
        Select the time window of the hietograph
        == time between two assimilation observations

        Parameters
        ----------
        NENS : list
            # of valid ensemble.
        """

        if len(self.grid3d) == 0:
            self.run_processor(IPRT1=3, DAFLAG=0, verbose=True)
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
                    VALUE = np.array(self.dict_parm_pert['atmbc0']['time_variable_perturbation'])[ens_nb]
                    times = self.dict_parm_pert['atmbc0']['data2perturbate']['time']
                    df_atmbc = pd.DataFrame(np.c_[times,VALUE],
                                            columns=['time','value']
                                            )

            if self.count_atmbc_cycle is not None:
                if HSPATM==0:
                    time_window_atmbc = [
                        int(df_atmbc.time.unique()[self.count_atmbc_cycle]),
                        int(df_atmbc.time.unique()[self.count_atmbc_cycle+1])
                        ]
                else:
                    time_window_atmbc = [
                        int(df_atmbc.time.iloc[self.count_atmbc_cycle]),
                        int(df_atmbc.time.iloc[self.count_atmbc_cycle + 1]),
                    ]
            else:
                time_window_atmbc = [self.ENS_times[0], self.ENS_times[1]]

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

    def _create_empty_ET_xr(self):
        x = np.unique(self.grid3d['mesh3d_nodes'][:,0])
        y = np.unique(self.grid3d['mesh3d_nodes'][:,1])
        nan_dataset = xr.Dataset(
                                {
                                    "time_sec": self.all_atmbc_times[self.count_atmbc_cycle],
                                    "SURFACE NODE": (("x", "y"), np.full((len(x), len(y)), np.nan)),
                                    "ACT. ETRA": (("x", "y"), np.full((len(x), len(y)), np.nan)),
                                },
                                coords={
                                    "x": x,
                                    "y": y,
                                }
                            )
        return nan_dataset

    def _DA_ET_xr(self):

        ET_ftime_ensi = []
        for nensi in range(self.NENS):
            if nensi in self.ens_valid:
                try:
                    df_fort777 = out_CT.read_fort777(
                        os.path.join(self.workdir,
                                     self.project_name,
                                     f'DA_Ensemble/cathy_{nensi+1}',
                                     'fort.777'
                                     )
                    )
                    df_fort777 = df_fort777.drop_duplicates()
                    df_fort777.columns = [col.lower() if col.lower() in ['x', 'y'] else col for col in df_fort777.columns]
                    df_fort777 = df_fort777.set_index(['time', 'x', 'y']).to_xarray()
                    df_fort777_selec_ti = df_fort777.isel(time=-1)
                    df_fort777_selec_ti = df_fort777_selec_ti.rio.set_spatial_dims('x', 'y')
                    df_fort777_selec_ti = df_fort777_selec_ti.reset_coords("time",
                                                                           drop=True
                                                                           )
                    df_fort777_selec_ti['time_sec'] = self.all_atmbc_times[self.count_atmbc_cycle]
                    ET_ftime_ensi.append(df_fort777_selec_ti)

                except:
                    print(f'Impossible to read ET - exclude scenario ensemble {nensi}')
                    nan_dataset = self._create_empty_ET_xr()
                    ET_ftime_ensi.append(nan_dataset)
            else:
                nan_dataset = self._create_empty_ET_xr()
                ET_ftime_ensi.append(nan_dataset)
        ET_DA_xr_time_i = xr.concat(ET_ftime_ensi,
                                dim="ensemble"
                                )

        # Concatenate along the "assimilation" dimension if `self.ET_DA_xr` exists, or create it otherwise
        if hasattr(self, 'ET_DA_xr'):
            self.ET_DA_xr = xr.concat([self.ET_DA_xr,
                                       ET_DA_xr_time_i.expand_dims(assimilation=[self.count_atmbc_cycle])],
                                      dim="assimilation"
                                      )
        else:
            self.ET_DA_xr = ET_DA_xr_time_i.expand_dims(assimilation=[self.count_atmbc_cycle])

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
        rFluid : TYPE, list
            Resistivity of the pore fluid. The default is [1.0].
        a : TYPE, list
            Tortuosity factor. The default is [1.0].
        m : TYPE, list
            Cementation exponent. The default is [2.0].
            (usually in the range 1.3 -- 2.5 for sandstones)
        n : TYPE, list
            Saturation exponent. The default is [2.0].
        pert_sigma_Archie : TYPE, list
            Gaussian noise to add. The default is None.

        ..note:
            Field procedure to obtain tce he covariance structure of the model
            estimates is described in Tso et al () - 10.1029/2019WR024964
            "Fit a straight line for log 10 (S) and log 10 (ρ S ) using the least-squares criterion.
            The fitting routine returns the covariance structure of the model estimates, which can be used to de-
            termine the 68% confidence interval (1 standard deviation) of the model estimates.""

        """
        if len(porosity) == 0:
            porosity = self.soil_SPP["SPP"]['POROS'].mean()

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
        pass

    def _compute_diff_matrix(self,data, prediction):
        """Computes the difference matrix between data and prediction."""
        diff_mat = np.zeros_like(prediction)
        for i in range(len(data)):
            for j in range(prediction.shape[1]):  # Loop over ensemble columns
                diff_mat[i, j] = abs(data[i] - prediction[i, j])
        return diff_mat

    def _compute_rmse(self,diff_matrix, num_predictions):
        """Computes RMSE from difference matrix."""
        diff_avg = np.nansum(diff_matrix, axis=1) / num_predictions
        return np.sum(diff_avg) / len(diff_avg)

    def _get_rmse_for_sensor(self, obs2eval, prediction, num_predictions):
        """Computes RMSE for a given sensor."""
        diff_mat = self._compute_diff_matrix(obs2eval, prediction)
        return self._compute_rmse(diff_mat, num_predictions)

    def _append_to_performance_df(self,df_performance, t_obs, name_sensor,
                                  rmse_sensor, rmse_avg, nrmse_sensor, nrmse_avg,
                                  ol_bool):
        """Appends performance metrics to the main DataFrame."""
        data = {
            "time": [t_obs],
            "ObsType": [name_sensor],
            "RMSE" + name_sensor: [rmse_sensor],
            "RMSE_avg": [rmse_avg],
            "NMRMSE" + name_sensor: [nrmse_sensor],
            "NMRMSE_avg": [nrmse_avg],
            "OL": [ol_bool]
        }
        return pd.concat([df_performance, pd.DataFrame(data)], ignore_index=True)

    def _save_performance_to_parquet(self,df_performance, file_path="performance_data.parquet"):
        """Saves the performance DataFrame to Parquet format, appending if the file exists."""
        if os.path.exists(file_path):
            df_existing = pd.read_parquet(file_path, engine="pyarrow")
            df_performance = pd.concat([df_existing, df_performance], ignore_index=True)

        df_performance.to_parquet(file_path, engine="pyarrow", compression="snappy", index=False)


    def _load_performance_from_parquet(self,file_path="performance_data.parquet"):
        """Loads performance data from a Parquet file if it exists."""
        if os.path.exists(file_path):
            return pd.read_parquet(file_path, engine="pyarrow")
        else:
            return pd.DataFrame()  # Return an empty DataFrame if no file exists

    def _performance_assessement_pq(self,
                                    list_assimilated_obs,
                                    data,
                                    prediction,
                                    t_obs,
                                    # file_path="performance_data.parquet",
                                    **kwargs
                                    ):
        """Calculates (N)RMSE for each observation in a data assimilation cycle."""

        ol_bool = kwargs.get("openLoop", False)
        file_path = self.project_name + '_df_performance.parquet'
        # Initialize df_performance if file exists
        if hasattr(self, 'df_performance') is False:
            self.df_performance = pd.DataFrame()
        else:
            self.df_performance = self._load_performance_from_parquet(os.path.join(
                                                                         self.workdir,
                                                                         self.project_name,
                                                                         file_path,
                                                                         )
                                                            )
        rmse_avg_stacked = []  # This list will track average RMSE values over time steps

        num_predictions = len(prediction)
        all_obs_diff_mat = self._compute_diff_matrix(data, prediction)
        all_obs_rmse_avg = self._compute_rmse(all_obs_diff_mat, num_predictions)

        # Example function to retrieve list of keys for assimilated observations
        obs2eval_key = self._get_data2assimilate(list_assimilated_obs)[1]
        start_line_obs = 0

        for name_sensor in obs2eval_key:
            obs2eval = self._get_data2assimilate([name_sensor], match=True)[0]
            n_obs = len(obs2eval)
            prediction2eval = prediction[start_line_obs:start_line_obs + n_obs]

            rmse_sensor = self._get_rmse_for_sensor(obs2eval, prediction2eval, num_predictions)
            rmse_avg_stacked.append(all_obs_rmse_avg)

            if "ObsType" in self.df_performance and name_sensor in self.df_performance["ObsType"].unique():
                prev_rmse = self.df_performance[self.df_performance["ObsType"] == name_sensor]["RMSE" + name_sensor].sum()
                rmse_sensor_stacked = prev_rmse + rmse_sensor
            else:
                rmse_sensor_stacked = rmse_sensor

            if t_obs == 0:
                nrmse_sensor = rmse_sensor
                nrmse_avg = all_obs_rmse_avg
            else:
                nrmse_sensor = rmse_sensor_stacked / (t_obs + 1)
                nrmse_avg = np.sum(rmse_avg_stacked) / (t_obs + 1)

            # t_obs = 3
            self.df_performance = self._append_to_performance_df(
                                                                self.df_performance,
                                                                t_obs,
                                                                name_sensor,
                                                                rmse_sensor,
                                                                all_obs_rmse_avg,
                                                                nrmse_sensor,
                                                                nrmse_avg,
                                                                ol_bool
                                                                )

            start_line_obs += n_obs

        # Save to Parquet and clear memory if too large
        # if len(self.df_performance) > 10000:  # Adjust batch size as needed
        self._save_performance_to_parquet(self.df_performance,
                                              os.path.join(
                                                self.workdir,
                                                self.project_name,
                                                file_path,
                                                )
                                     )
        self.df_performance = pd.DataFrame()  # Clear in-memory DataFrame
