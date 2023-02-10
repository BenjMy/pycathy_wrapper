# from pyCATHY.DA.cathy_DA import DA
import pandas as pd


def run_DA_sequential(
    self,
    callexe,
    parallel,
    DA_type,
    dict_obs,
    list_update_parm,
    dict_parm_pert,
    list_assimilated_obs,
    open_loop_run,
    threshold_rejected,
    verbose,
    **kwargs
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
    # from pyCATHY.DA.cathy_DA import DA
    # CATHYDA = DA()

    # Initiate
    # -------------------------------------------------------------------
    update_key = "ini_perturbation"

    # check if dict_obs is ordered in time
    # ------------------------------------
    if (
        all(i < j for i, j in zip(list(dict_obs.keys()), list(dict_obs.keys())[1:]))
    ) is False:
        raise ValueError("Observation List is not sorted.")

    # dict_obs.keys()
    # self.dict_obs.keys()
    if hasattr(self, "dict_obs") is False:
        self.dict_obs = dict_obs  # self.dict_obs is already assigned (in read observatio! change it to self.obs
    self.dict_parm_pert = dict_parm_pert
    self.df_DA = pd.DataFrame()

    # Infer ensemble size NENS from perturbated parameter dictionnary
    # -------------------------------------------------------------------
    for name in self.dict_parm_pert:
        NENS = len(self.dict_parm_pert[name]["ini_perturbation"])

    # Infer ensemble update times ENS_times from observation dictionnary
    # -------------------------------------------------------------------
    ENS_times = []
    for ti in self.dict_obs:
        ENS_times.append(float(ti))

    # data_measure_df = self.dictObs_2pd()

    # start DA cycle counter
    # -------------------------------------------------------------------
    self.count_DA_cycle = 0
    self.count_atmbc_cycle = 0
    # (the counter is incremented during the update analysis)

    # initiate DA
    # -------------------------------------------------------------------
    list_update_parm = self._DA_init(
        NENS=NENS,  # ensemble size
        ENS_times=ENS_times,  # assimilation times
        parm_pert=dict_parm_pert,
        update_parm_list=list_update_parm,
    )

    # initiate mapping petro
    # -------------------------------------------------------------------
    self._mapping_petro_init()

    # update the perturbated parameters
    # --------------------------------------------------------------------
    self.update_ENS_files(
        dict_parm_pert,
        update_parm_list="all",  # list_update_parm
        cycle_nb=self.count_DA_cycle,
    )

    all_atmbc_times = self.atmbc["time"]
    # -------------------------------------------------------------------
    if open_loop_run:
        self._DA_openLoop(ENS_times, list_assimilated_obs, parallel)
    # end of Open loop - start DA

    # -------------------------------------------------------------------
    # update input files ensemble again (time-windowed)
    # ---------------------------------------------------------------------
    self._update_input_ensemble(
        NENS, list(self.atmbc["time"]), dict_parm_pert, update_parm_list="all"
    )  # list_update_parm

    # -----------------------------------
    # Run hydrological model sequentially = Loop over atmbc times (including assimilation observation times)
    # -----------------------------------
    # self.sequential_DA()

    # TO CHANGE HERE
    for (
        t_atmbc
    ) in all_atmbc_times:  # atmbc times include assimilation observation times
        print(t_atmbc)
        # t_atmbc = self.atmbc['time'][-2]

        self._run_ensemble_hydrological_model(parallel, verbose, callexe)
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

        # define subloop here
        # if 'optimize_inflation' in DA_type:
        # threshold_convergence = 10
        # while self.df_performance < threshold_convergence:

        # self.atmbc['time']

        if t_atmbc in ENS_times:

            print("t=" + str(t_atmbc))

            # map states to observation = apply H operator to state variable
            # ----------------------------------------------------------------
            prediction = self.map_states2Observations(
                list_assimilated_obs,
                default_state="psi",
                parallel=parallel,
            )
            # print(len(prediction))
            # print(np.shape(prediction))
            # ValueError: operands could not be broadcast together with shapes (612,) (1836,)

            # data2test , obs2evaltest = self._get_data2assimilate(list_assimilated_obs)
            # # len(data2test)

            # x = np.linspace(prediction.min(),
            #                 prediction.max(),
            #                 100)
            # y = x
            # fig, ax = plt.subplots(2,1)
            # for ii in range(np.shape(prediction)[1]):
            #     ax[0].scatter(prediction[:,ii],data2test)
            #     ax[0].plot(x, y, '-r', label='y=x')
            # ax[0].set_title('nb of ens' + str(np.shape(prediction)[1]))

            # ax[1].plot(data2test)
            # plt.savefig('test' + str(DA_cnb))

            # plt.plot(data2test)
            # plt.plot(prediction[:,ii])
            # ax.scatter(prediction[:,ii],data2test)
            # ax.set_aspect('equal')

            # fig.tight_layout()

            # DA analysis
            # ----------------------------------------------------------------
            (
                ensemble_psi,
                ensemble_sw,
                data,
                analysis,
                analysis_param,
            ) = self._DA_analysis(
                prediction,
                DA_type,
                list_update_parm,
                list_assimilated_obs,
                ens_valid=self.ens_valid,
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

            # check analysis quality
            # ----------------------------------------------------------------

            self.console.print(":face_with_monocle: [b]check analysis performance[/b]")
            self._performance_assessement(
                list_assimilated_obs, data, prediction_valid, t_obs=self.count_DA_cycle
            )

            # the counter is incremented here
            # ----------------------------------------------------------------
            self.count_DA_cycle = self.count_DA_cycle + 1

            # # update parameter dictionnary
            # ----------------------------------------------------------------
            def check_ensemble(ensemble_psi_valid, ensemble_sw_valid):
                if np.any(ensemble_psi_valid > 0):
                    # raise ValueError('positive pressure heads observed')
                    print("!!!!!!positive pressure heads observed!!!!!")
                    psi_2replace = np.where(ensemble_psi_valid >= 0)
                    ensemble_psi_valid_new = ensemble_psi_valid
                    ensemble_psi_valid_new[psi_2replace] = -1e-3

            check_ensemble(ensemble_psi_valid, ensemble_sw_valid)

            self.update_pert_parm_dict(
                update_key, list_update_parm, analysis_param_valid
            )

        else:
            self.console.print(
                ":confused: No observation for this time - run hydrological model only"
            )
            print("!!!!!!!!! shoetcutttt here ensemble are anot validated!!!!!!!!!! S")
            (
                ensemble_psi_valid,
                ensemble_sw_valid,
                ens_size,
                sim_size,
            ) = self._read_state_ensemble()
            # analysis_valid = np.empty(ensemble_psi_valid.shape)
            # analysis_valid[:] = np.NaN
            analysis_valid = ensemble_psi_valid

        self.count_atmbc_cycle = self.count_atmbc_cycle + 1

        print(
            "------ end of time step (s) -------"
            + str(int(t_atmbc))
            + "/"
            + str(int(all_atmbc_times[-1]))
            + "------"
        )
        print(
            "------ end of atmbc update --------"
            + str(self.count_atmbc_cycle)
            + "/"
            + str(len(all_atmbc_times) - 1)
            + "------"
        )

        # create dataframe _DA_var_pert_df holding the results of the DA update
        # ---------------------------------------------------------------------
        self._DA_df(
            state=[ensemble_psi_valid, ensemble_sw_valid],
            state_analysis=analysis_valid,
            rejected_ens=rejected_ens,
        )

        # export summary results of DA
        # ----------------------------------------------------------------

        meta_DA = {
            "listAssimilatedObs": list_assimilated_obs,
            "listUpdatedparm": list_update_parm,
            # '':,
            # '':,
            # '':,
        }

        self.backup_results_DA(meta_DA)

        # test = self.Archie

        self.backup_simu()

        # overwrite input files ensemble (perturbated variables)
        # ---------------------------------------------------------------------

        if (
            self.count_atmbc_cycle < len(all_atmbc_times) - 1
        ):  # -1 cause all_atmbc_times include TMAX
            # len(all_atmbc_times)
            self._update_input_ensemble(
                self.NENS,
                all_atmbc_times,
                self.dict_parm_pert,
                update_parm_list=list_update_parm,
                analysis=analysis_valid,
            )
        else:
            print("------ end of DA ------")
            pass
        # print('------ end of update -------' + str(self.count_atmbc_cycle) + '/' + str(len(all_atmbc_times)-1) + '------')
        print(
            "------ end of DA update -------"
            + str(self.count_DA_cycle)
            + "/"
            + str(len(ENS_times))
            + "------"
        )
        print(
            "% of valid ensemble is: " + str((len(self.ens_valid) * 100) / (self.NENS))
        )
        # print(self.ens_valid)

        plt.close("all")
