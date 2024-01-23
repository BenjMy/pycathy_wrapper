"""Class managing data Assimilation mapping of state variable to observations
"""
import pyCATHY.meshtools as mt
import numpy as np
import os 
from pyCATHY.ERT import petro_Archie as Archie
from pyCATHY.importers import cathy_outputs as out_CT
import multiprocessing
from functools import partial
REMOVE_CPU = 1
from pyCATHY.cathy_tools import CATHY

#%%

def tensio_mapper(state,obs2map_node):
    # case: pressure head assimilation (Hx_PH)
    # -------------------------------------------------------------
    Hx_PH = state[0][obs2map_node]
    return Hx_PH


def swc_mapper(state,obs2map_node,grid3d,dem_parameters,
               SPP_ensi,console): 
    # find porosity associated to mesh node position
    # ----------------------------------------------
    xyz_swc = grid3d[obs2map_node]
    
    ltop, lbot = mt.get_layer_depths(dem_parameters)
    idlayer_swc = np.where(abs(xyz_swc[2])>ltop)[0]
    porosity_swc = SPP_ensi.xs((1,idlayer_swc[0]),
                               level=('zone', 'layer'))['POROS'].values
    console.print(f"Transform sat to SWC with porosity=' {str(porosity_swc)}")
    console.rule(
        ":warning: warning messages above :warning:", style="yellow"
    )
    console.print(
        r"""Current implementation does not support different porosity zones! 
            Only layers porosity is considered - Taking zone 1 as default
        """,
        style="yellow",
    )
    Hx_SW = state[1][obs2map_node] * porosity_swc
    return Hx_SW


def ET_mapper(obs2map_node,
              console,
              path_fwd_CATHY
              ):
    # Atmact-vf(13) : Actual infiltration (+ve) or exfiltration (-ve) at atmospheric BC nodes as a volumetric flux [L^3/T]
    # Atmact-v (14) : Actual infiltration (+ve) or exfiltration (-ve) volume [L^3]
    # Atmact-r (15) : Actual infiltration (+ve) or exfiltration (-ve) rate [L/T]
    # Atmact-d (16) : Actual infiltration (+ve) or exfiltration (-ve) depth [L]
    console.print(f"Read actual ET on a given node")

    df_fort777 = out_CT.read_fort777(os.path.join(path_fwd_CATHY,
                                                  'fort.777'),
                                      )
    df_fort777 = df_fort777.set_index('time_sec')                    
    t_ET = df_fort777.index.unique()
    # print('here I''m taking the wrong time no? ' + str(t_ET[-1]))
    Hx_ET = df_fort777.loc[t_ET[-1]].iloc[obs2map_node]['ACT. ETRA']
    return Hx_ET
                    
#%% ERT mapping

def _add_2_ensemble_Hx(Hx, Hx_2add):
    """
    Store in an array predicted value for all ensembles and all assimilation times
    """
    try:
        Hx.append(Hx_2add)
    except:
        Hx = list(Hx)
        Hx.append(Hx_2add)
    return Hx
    
def _map_ERT(state,
             dict_obs, 
             project_name,
             Archie_parms,
             count_DA_cycle,
             path_fwd_CATHY, 
             ens_nb, 
             **kwargs
             ):
    """
    Mapping of state variable to observation (predicted)
    ERT using pedophysical transformation H
    """

    savefig = False
    if "savefig" in kwargs:
        savefig = kwargs["savefig"]

    # search key value to identify time and method
    # --------------------------------------------
    tuple_list_obs = list(dict_obs.items())
    key_time = tuple_list_obs[count_DA_cycle]
    
    # Load ERT metadata information form obs dict
    # -------------------------------------------
    ERT_meta_dict = _parse_ERT_metadata(key_time)

    Hx_ERT, df_Archie = Archie.SW_2_ERa_DA(
                                            project_name,
                                            Archie_parms,
                                            Archie_parms["porosity"],
                                            ERT_meta_dict,
                                            path_fwd_CATHY,
                                            df_sw=state[1],  # kwargs
                                            DA_cnb=count_DA_cycle,  # kwargs
                                            Ens_nbi=ens_nb,  # kwargs
                                            savefig=savefig,  # kwargs
                                            noise_level=ERT_meta_dict["data_err"],  # kwargs
                                            dict_ERT=key_time[1]["ERT"],  #  kwargs
                                        )
    df_Archie["OL"] = np.ones(len(df_Archie["time"])) * False

    return Hx_ERT, df_Archie

def _map_ERT_parallel_DA(
                        dict_obs,
                        ENS_times,
                        ERT_meta_dict,
                        key_time,
                        path_fwd_CATHY_list,
                        DA_cnb,
                        project_name,
                        Archie_parms,
                    ):
    """
    Parallel mapping of ERT data using pedophysical transformation H
    """
    Hx_ERT_ens = []

    # freeze fixed arguments of Archie.SW_2_ERa_DA
    # -----------------------------------------------------------------
    ERTmapping_args = partial(
                                Archie.SW_2_ERa_DA,
                                project_name,
                                Archie_parms,
                                Archie_parms["porosity"],
                                ERT_meta_dict,
                                DA_cnb=DA_cnb,
                                savefig=False,
                                noise_level=ERT_meta_dict["data_err"],  # kwargs
                                dict_ERT=key_time[1]["ERT"],  #  kwargs
                            )

    # // run using ensemble subfolders path as a list
    # -----------------------------------------------------------------
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()-REMOVE_CPU) as pool:
        results_mapping = pool.map(ERTmapping_args, 
                                   path_fwd_CATHY_list
                                   )

    for ens_i in range(len(path_fwd_CATHY_list)):
        df_Archie = results_mapping[ens_i][1]
        df_Archie["OL"] = np.zeros(len(df_Archie))
        Hx_ERT_ens_i = results_mapping[ens_i][0]

        if "pygimli" in dict_obs[key_time[0]]["ERT"]["data_format"]:
            Hx_ERT_ens = _add_2_ensemble_Hx(Hx_ERT_ens, Hx_ERT_ens_i["rhoa"])
        else:
            Hx_ERT_ens = _add_2_ensemble_Hx(Hx_ERT_ens, Hx_ERT_ens_i["resist"])

    return Hx_ERT_ens, df_Archie

def _map_ERT_parallel_OL(
                            project_name,
                            Archie_parms,
                            ENS_times,
                            ERT_meta_dict,
                            key_time,
                            path_fwd_CATHY_list,
                        ):
    """
    Parallel mapping of ERT data using pedophysical transformation H
    case of the open Loop = nested loop with ensemble time
    """

    Hx_ERT_ens = []
    for t in range(len(ENS_times)):
        print("t_openLoop mapping:" + str(t))

        # freeze fixed arguments of Archie.SW_2_ERa
        # -----------------------------------------------------------------
        ERTmapping_args = partial(
                                    Archie.SW_2_ERa_DA,
                                    project_name,
                                    Archie_parms,
                                    Archie_parms["porosity"],
                                    ERT_meta_dict,
                                    time_ass=t,
                                    savefig=True,
                                    noise_level=ERT_meta_dict["data_err"],
                                    dict_ERT=key_time[1]["ERT"],
                                )
        #
        # -----------------------------------------------------------------
        with multiprocessing.Pool(processes=multiprocessing.cpu_count()-REMOVE_CPU) as pool:
            results_mapping_time_i = pool.map(ERTmapping_args, path_fwd_CATHY_list)
            # print(f"x= {path_fwd_CATHY_list}, PID = {os.getpid()}")

        for ens_i in range(len(path_fwd_CATHY_list)):
            Hx_ERT_time_i = results_mapping_time_i[ens_i][0]

            # print(np.shape(results_mapping_time_i))
            df_Archie = results_mapping_time_i[ens_i][1]

            # print(df_Archie)
            df_Archie["OL"] = np.ones(len(df_Archie))
            self._add_2_ensemble_Archie(df_Archie)

            if "pygimli" in ERT_meta_dict["data_format"]:
                Hx_ERT_ens = _add_2_ensemble_Hx(
                    Hx_ERT_ens, Hx_ERT_time_i["rhoa"]
                )
            else:
                Hx_ERT_ens = _add_2_ensemble_Hx(
                    Hx_ERT_ens, Hx_ERT_time_i["resist"]
                )

    # prediction_ERT = np.reshape(Hx_ERT_ens,[self.NENS,
    #                                         len(Hx_ERT_ens[0]),
    #                                         len(ENS_times)])  # (EnSize * data size * times)
    return Hx_ERT_ens

def _map_ERT_parallel(
                        dict_obs,
                        count_DA_cycle,
                        project_name,
                        Archie_parms,
                        path_fwd_CATHY_list,
                        list_assimilated_obs="all",
                        default_state="psi",
                        verbose=False,
                        **kwargs,
                    ):
    
                                                            
    """Mapping of state variable to observation (predicted) ERT using pedophysical transformation H,
    // run using ensemble subfolders path as a list
    """

    savefig = False
    if "savefig" in kwargs:
        savefig = kwargs["savefig"]
    ENS_times = []
    if "ENS_times" in kwargs:
        ENS_times = kwargs["ENS_times"]
    DA_cnb = []
    if "DA_cnb" in kwargs:
        DA_cnb = kwargs["DA_cnb"]

    # search key value to identify time and method
    tuple_list_obs = list(dict_obs.items())
    key_time = tuple_list_obs[count_DA_cycle]
    
    # Load ERT metadata information form obs dict
    # -------------------------------------------
    ERT_meta_dict = _parse_ERT_metadata(key_time)

    if len(ENS_times) > 0:  # case of the open Loop = nested loop with ensemble time
        Hx_ERT_ens = _map_ERT_parallel_OL(
                                            project_name,
                                            Archie_parms,
                                            ENS_times,
                                            ERT_meta_dict,
                                            key_time,
                                            path_fwd_CATHY_list,
                                        )
    else:
        Hx_ERT_ens, df_Archie = _map_ERT_parallel_DA(
                                                    dict_obs,
                                                    ENS_times,
                                                    ERT_meta_dict,
                                                    key_time,
                                                    path_fwd_CATHY_list,
                                                    DA_cnb,
                                                    project_name,
                                                    Archie_parms,
                                                )
    prediction_ERT = np.vstack(Hx_ERT_ens).T

    return prediction_ERT, df_Archie



def _parse_ERT_metadata(key_value):
    """
    Extract ERT metadata information form obs dict
    """
    ERT_meta_dict = {}
    ERT_meta_dict["forward_mesh_vtk_file"] = key_value[1]["ERT"][
        "forward_mesh_vtk_file"
    ]
    
    for kk in key_value[1]["ERT"].keys():
        if 'filename' in kk:
            ERT_meta_dict["pathERT"] = os.path.split(key_value[1]["ERT"]["filename"])[0]
        else:
            ERT_meta_dict[kk]= key_value[1]["ERT"][kk]

    return ERT_meta_dict


#%%

# class petro(CATHY):
    
#     def set_Archie_parm(
#         self,
#         porosity=[],
#         rFluid_Archie=[1.0],
#         a_Archie=[1.0],
#         m_Archie=[2.0],
#         n_Archie=[2.0],
#         pert_sigma_Archie=[0],
#     ):
#         """
#         Dict of Archie parameters. Each type of soil is describe within a list
#         Note that if pert_sigma is not None a normal noise is
#         added during the translation of Saturation Water to ER
    
#         Parameters
#         ----------
#         rFluid : TYPE, optional
#             Resistivity of the pore fluid. The default is [1.0].
#         a : TYPE, optional
#             Tortuosity factor. The default is [1.0].
#         m : TYPE, optional
#             Cementation exponent. The default is [2.0].
#             (usually in the range 1.3 -- 2.5 for sandstones)
#         n : TYPE, optional
#             Saturation exponent. The default is [2.0].
#         pert_sigma_Archie : TYPE, optional
#             Gaussian noise to add. The default is None.
    
#         ..note:
#             Field procedure to obtain tce he covariance structure of the model
#             estimates is described in Tso et al () - 10.1029/2019WR024964
#             "Fit a straight line for log 10 (S) and log 10 (ρ S ) using the least-squares criterion.
#             The fitting routine returns the covariance structure of the model estimates, which can be used to de-
#             termine the 68% confidence interval (1 standard deviation) of the model estimates.""
    
#         """
#         if len(porosity) == 0:
#             porosity = self.soil_SPP["SPP"][:, 4][0]
#         if not isinstance(porosity, list):
#             porosity = [porosity]
#         self.Archie_parms = {
#             "porosity": porosity,
#             "rFluid_Archie": rFluid_Archie,
#             "a_Archie": a_Archie,
#             "m_Archie": m_Archie,
#             "n_Archie": n_Archie,
#             "pert_sigma_Archie": pert_sigma_Archie,
#         }
    
#         pass
