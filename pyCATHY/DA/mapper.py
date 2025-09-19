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
import pandas as pd

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
    # console.print(f"Read actual ET on node {obs2map_node}")

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


def build_forward_profiles(sw_nodes: np.ndarray,
                           mesh3d_nodes: np.ndarray,
                           var: str = "EC"):
    """
    Build depths and conductivity arrays for all (x, y) columns from 3D node data.

    Parameters
    ----------
    sw_nodes : np.ndarray, shape (nnodes,)
        Array of SW or EC values at each node.
    mesh3d_nodes : np.ndarray, shape (nnodes, 3)
        Node coordinates, columns: [x, y, z].
    var : str
        Name of the variable (for labeling, not used in calculations).

    Returns
    -------
    depths : list of np.ndarray
        List of layer thicknesses for each (x, y) column.
    conds : list of np.ndarray
        List of conductivity arrays for each (x, y) column.
    xy_coords : list of tuple
        List of (x, y) coordinates corresponding to each profile.
    """

    depths = []
    conds = []
    xy_coords = []

    # extract x, y, z
    x_nodes = mesh3d_nodes[:, 0]
    y_nodes = mesh3d_nodes[:, 1]
    z_nodes = mesh3d_nodes[:, 2]

    # get unique (x, y) pairs
    xy_unique = np.unique(np.column_stack([x_nodes, y_nodes]), axis=0)

    for x, y in xy_unique:
        # mask for nodes at this (x, y)
        mask = (x_nodes == x) & (y_nodes == y)
        z_profile = z_nodes[mask]
        cond_profile = sw_nodes[mask]

        # sort so that surface is first (highest z)
        sort_idx = np.argsort(-z_profile)  # descending
        z_sorted = z_profile[sort_idx]
        cond_sorted = cond_profile[sort_idx]

        nz = len(cond_sorted)

        if nz == 1:
            depths.append(np.array([]).reshape(1, -1))
            conds.append(cond_sorted.reshape(1, -1))
        else:
            # interface depths = difference from surface to next layers
            z_surface = z_sorted[0]
            depth_interfaces = np.abs(z_surface - z_sorted[1:])  # (nz-1,)
            depths.append(depth_interfaces.reshape(1, -1))
            conds.append(cond_sorted.reshape(1, -1))

        xy_coords.append((float(x), float(y)))

    return depths, conds, xy_coords

    
  
def _map_EM(
             dict_obs, 
             project_name,
             Archie_parms,
             count_DA_cycle,
             path_fwd_CATHY, 
             ens_nb, 
             POROS_mesh_nodes_ensi,
             grid3d,
             **kwargs
             ):
   
    """
    Mapping of state variable to observation (predicted)
    EM using pedophysical transformation H
    """
    
    # search key value to identify time and method
    # --------------------------------------------
    tuple_list_obs = list(dict_obs.items())
    key_time = tuple_list_obs[count_DA_cycle]
    
    #%%
    # Load ERT metadata information form obs dict
    # -------------------------------------------
    EM_meta_dict = _parse_EM_metadata(key_time)
    EM_meta_dict.keys()
    
    ER_converted_ti, df_Archie, sw_nodes =   Archie.SW_2_ER0_DA(
                                            project_name,
                                            Archie_parms,
                                            POROS_mesh_nodes_ensi,
                                            EM_meta_dict,
                                            path_fwd_CATHY,
                                            DA_cnb=count_DA_cycle,  # kwargs
                                            Ens_nbi=ens_nb,  # kwargs
                                            # savefig=savefig,  # kwargs
                                            noise_level=EM_meta_dict["data_err"],  # kwargs
                                        )
    df_Archie["OL"] = np.ones(len(df_Archie["time"])) * False
    EC_converted_ti_mS_m = (1.0 / ER_converted_ti) * 1000
    
    # EC_converted_ti_mS_m stats:
    #   min = 21.94 mS/m
    #   max = 70.94 mS/m
      
    # print("ER_converted_ti stats:")
    # print(f"  min = {ER_converted_ti.min():.2f} Ohm.m")
    # print(f"  max = {ER_converted_ti.max():.2f} Ohm.m")
      
    print("EC_converted_ti_mS_m stats:")
    print(f"  min = {EC_converted_ti_mS_m.min():.2f} mS/m")
    print(f"  max = {EC_converted_ti_mS_m.max():.2f} mS/m")


    print('Build forward EM model')
    depths, conds, xy_coords = build_forward_profiles(EC_converted_ti_mS_m, 
                                                      grid3d['mesh3d_nodes'],
                                                      var="EC")
    
    print("depths stats:")
    print(f"  min = {np.array(depths).min():.2f} m")
    print(f"  max = {np.array(depths).max():.2f} m")
    

    from emagpy import Problem
    k = Problem()
    k.setModels(depths,conds)
    
    print('Forward EM model')
    # dfsFSeq = k.forward(forwardModel='FSeq', coils=coils, noise=5)
    EM_fwd_model_array = k.forward(forwardModel='FSlin',
                         coils=EM_meta_dict["coils"],
                         noise=EM_meta_dict["fwdEMnoise"]/100
                         # noise=0
                         )
    
    # EM_fwd_model_array = k.forward(forwardModel='CS',
    #                      coils=EM_meta_dict["coils"],
    #                      # noise=EM_meta_dict["fwdEMnoise"]
    #                      noise=0
    #                      )
    # depths[0]
    # conds[0]
    # EM_fwd_model_array[0]
    # np.shape(dfsFSlin)
    EM_fwd_model_2darray = np.vstack(EM_fwd_model_array)
    EM_fwd_model_1darray = np.hstack(EM_fwd_model_2darray)

    return EM_fwd_model_1darray, df_Archie
 


def _map_ERT_parallel_DA(
                        dict_obs,
                        ENS_times,
                        ERT_meta_dict,
                        key_time,
                        path_fwd_CATHY_list,
                        DA_cnb,
                        project_name,
                        Archie_parms,
                        POROS_nodes_ensi_list,
                    ):
    """
    Parallel mapping of ERT data using pedophysical transformation H
    """
    Hx_ERT_ens = []

    # freeze fixed arguments of Archie.SW_2_ERT_ERa_DA
    # -----------------------------------------------------------------
    # print("POROS_nodes_ensi shape:", np.shape(POROS_nodes_ensi))
    # print("fwd_path shape:", np.shape(path_fwd_CATHY_list))

    ERTmapping_args = partial(
                                Archie.SW_2_ERT_ERa_DA,
                                project_name,
                                Archie_parms,
                                POROS_nodes_ensi_list,
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

    Archie_df = []
    for ens_i in range(len(path_fwd_CATHY_list)):
        df_Archie2add = results_mapping[ens_i][1]
        df_Archie2add["OL"] = np.zeros(len(df_Archie2add))
        Hx_ERT_ens_i = results_mapping[ens_i][0]
        mesh2test = results_mapping[ens_i][2]

        
        Archie_df = _add_2_ensemble_Archie(Archie_df,df_Archie2add)

        if "pygimli" in dict_obs[key_time[0]]["ERT"]["data_format"]:
            Hx_ERT_ens = _add_2_ensemble_Hx(Hx_ERT_ens, Hx_ERT_ens_i["rhoa"])
        else:
            Hx_ERT_ens = _add_2_ensemble_Hx(Hx_ERT_ens, Hx_ERT_ens_i["resist"])

    return Hx_ERT_ens, Archie_df, mesh2test



def _add_2_ensemble_Archie(Archie_df, df_Archie_2add):
    """
    Append Archie relationship data for a specific ensemble and time step to the main ensemble DataFrame.

    This function concatenates a new DataFrame (`df_Archie_2add`) containing Archie-related
    simulation or assimilation results to an existing ensemble-level DataFrame (`Archie_df`).
    If `Archie_df` is empty, it initializes it with the appropriate columns.

    Parameters
    ----------
    Archie_df : pd.DataFrame
        DataFrame containing existing Archie relationship entries with columns:
        ["time", "ens_nb", "sw", "ER_converted", "OL"].

    df_Archie_2add : pd.DataFrame
        New data to append, structured with the same columns as `Archie_df`.

    Returns
    -------
    pd.DataFrame
        Updated `Archie_df` with the new entries from `df_Archie_2add` appended.
    """

    if len(Archie_df) == 0:
        Archie_df = pd.DataFrame(
            columns=["time", "ens_nb", "sw", "ER_converted", "OL"]
        )
    Archie_df = pd.concat([Archie_df, df_Archie_2add])
    
    return Archie_df
    
def _map_ERT_parallel(
                        dict_obs,
                        count_DA_cycle,
                        project_name,
                        Archie_parms,
                        POROS_nodes_ensi_list,
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
    Hx_ERT_ens, df_Archie, mesh2test = _map_ERT_parallel_DA(
                                                dict_obs,
                                                ENS_times,
                                                ERT_meta_dict,
                                                key_time,
                                                path_fwd_CATHY_list,
                                                DA_cnb,
                                                project_name,
                                                Archie_parms,
                                                POROS_nodes_ensi_list,
                                            )
    prediction_ERT = np.vstack(Hx_ERT_ens).T

    return prediction_ERT, df_Archie, mesh2test


def _parse_EM_metadata(key_value):
    """
    Extract EM metadata information form obs dict
    """
    EM_meta_dict = {}
    
    for kk in key_value[1]["EM"].keys():
        if 'filename' in kk:
            EM_meta_dict["pathEM"] = os.path.split(key_value[1]["EM"]["filename"])[0]
        else:
            EM_meta_dict[kk]= key_value[1]["EM"][kk]

    return EM_meta_dict

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


