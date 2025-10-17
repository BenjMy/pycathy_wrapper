"""Pedophysical transformations"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyvista as pv
import rich.console
from rich import print

from pyCATHY import meshtools as mt
from pyCATHY.ERT import simulate_ERT as simuERT
from pyCATHY.importers import cathy_outputs as out_CT
from pyCATHY.meshtools import CATHY_2_pg, CATHY_2_Resipy


def get_sw_ens_i(path_fwd_CATHY, **kwargs):
    """Return sw array for a given ensemble (if ensemble exist)"""
    # Get the state vector from kwargs if existing otherwise read output
    # ------------------------------------------------------------------------
    if "df_sw" in kwargs:
        df_sw = kwargs["df_sw"]
    else:
        df_sw, _ = out_CT.read_sw(os.path.join(path_fwd_CATHY, "output/sw"))
        # case of open Loop where df_sw is list of nb of assimilation time length
        # -------------------------------------------------------------------
        if "time_ass" in kwargs:
            time_ass = kwargs["time_ass"]
            df_sw = df_sw.iloc[time_ass, :]
            DA_cnb = time_ass
        # case of sequential assimilation
        # -------------------------------------------------------------------
        else:
            df_sw =  df_sw.iloc[-1, :]  # take the last time
            # (in case the sw file contains intermediate times for plotting observations)
    return df_sw.values


def get_Archie_ens_i(ArchieParms, Ens_nb):
    """Return Archie parameter for a given ensemble (if ensemble exist)"""
    ArchieParms2parse = {}
    if len(ArchieParms["rFluid_Archie"]) > 1:
        for p in [
            "porosity",
            "rFluid_Archie",
            "a_Archie",
            "m_Archie",
            "n_Archie",
            "pert_sigma_Archie",
        ]:
            ArchieParms2parse[p] = [ArchieParms[p][Ens_nb]]
    else:
        ArchieParms2parse = ArchieParms
    return ArchieParms2parse

def SW_2_ER0_DA(
                project_name, 
                ArchieParms, 
                porosity, 
                meta_dict, 
                path_fwd_CATHY, 
                **kwargs
                ):
    
    if "DA_cnb" in kwargs:
        DA_cnb = kwargs["DA_cnb"]
    if "Ens_nbi" in kwargs:
        Ens_nbi = kwargs["Ens_nbi"]
        
    # Get sw array for a given ensemble
    # ------------------------------------
    df_sw = get_sw_ens_i(path_fwd_CATHY, **kwargs)

    # Choose archie parameter for a given realisation (from the ensemble)
    # --------------------------------------------------------------------
    ArchieParms2parse = get_Archie_ens_i(ArchieParms, Ens_nbi)
    
    # Convert to SW to ER values
    # --------------------------------------------------------------------
    print("****")
 

    
    if isinstance(porosity[0], (list, np.ndarray)):
        porosity = porosity[Ens_nbi]
    else:
        porosity = porosity
        
    # print(f"porosity_ensi: {porosity_ensi}")
    # print(f"porosity_ensi: {porosity_ensi}")
    print(f"Ens_nbi indices: {Ens_nbi}")
    print(f"Porosity shape: {np.shape(porosity)}")
    print(f"Porosity statistics (selected indices): min={np.min(porosity[Ens_nbi]):.3f}, "
          f"mean={np.mean(porosity[Ens_nbi]):.3f}, max={np.max(porosity[Ens_nbi]):.3f}")
    print("****")

    
    ER_converted_ti = Archie_rho_DA(
                                    rFluid_Archie=ArchieParms2parse["rFluid_Archie"],
                                    sat=[df_sw],
                                    porosity=porosity,
                                    a_Archie=ArchieParms2parse["a_Archie"],
                                    m_Archie=ArchieParms2parse["m_Archie"],
                                    n_Archie=ArchieParms2parse["n_Archie"],
                                    pert_sigma_Archie=ArchieParms2parse["pert_sigma_Archie"],
                                )

    # build df Archie df
    df_Archie = pd.DataFrame(columns=["time", "ens_nb", "sw", "ER_converted"])
    df_Archie["time"] = DA_cnb * np.ones(len(ER_converted_ti))
    df_Archie["ens_nbi"] = Ens_nbi * np.ones(len(ER_converted_ti))
    df_Archie["sw"] = df_sw
    df_Archie["ER_converted"] = ER_converted_ti
    df_Archie["porosity"] = ArchieParms2parse["porosity"] * np.ones(
        len(ER_converted_ti)
    )
    return ER_converted_ti, df_Archie, df_sw

# def ER0_to_EM_ECa_DA(ER_converted_ti,EM_meta_dict):
    
#     # ------------------------------------------------------------------------
#     print('----fwd ER data-----')
#     from emagpy import Problem
#     k= Problem()
#     k.createSurvey(EM_meta_dict['filename'])
    
    
    
def SW_2_ERT_ERa_DA(
                project_name, 
                ArchieParms, 
                porosity, 
                ERT_meta_dict, 
                path_fwd_CATHY, 
                **kwargs
                ):

    """
    Data Assimilation
    -----------------
    Map saturation water from CATHY model to apparent Electrical Resistivities
    1. Import/read CATHY sw file result
    2. Convert SW to ER uisng Archie law
    3. Fwd ERT model to obtain apparent resistivities (predicted observation for DA)

    Parameters
    ----------
    project_name : str
        name of the current project.
    ArchieParms : dict
        Archie params.
    porosity : np.array([])
        medium porosity.
    pathERT : str
        path of the ERT forward mesh.
    meshERT : str
        filename of the ERT fwd mesh.
    elecs : TYPE
        electrode positions.
    sequenceERT : TYPE
        ERT sequence.
    path_fwd_CATHY: list of folder containing CATHY simulation outputs
        Use for // simulations
    **kwargs :
        df_sw : pd df
            saturation mesh values to convert
        savefig : TYPE, optional
            DESCRIPTION. The default is True.
    Returns
    -------
    df_ERT_predicted : pd df
        Dataframe of predicted ER.
    df_Archie : pd df
        Dataframe of the current (given ensemble and given assimilation time) Archie relationship.

    """
    
    path_CATHY = os.path.join(path_fwd_CATHY, "vtk/")

    # Some flag for DA assimilation
    # ------------------------------------------------------------------------
    DA_cnb = None
    Ens_nbi = int(os.path.split(path_fwd_CATHY)[-1].split("_", 1)[1]) - 1

    savefig = False
    if "DA_cnb" in kwargs:
        DA_cnb = kwargs["DA_cnb"]
    if "Ens_nbi" in kwargs:
        Ens_nbi = kwargs["Ens_nbi"]
    if "savefig" in kwargs:
        savefig = kwargs["savefig"]
                
    ER_converted_ti, df_Archie, df_sw = SW_2_ER0_DA(
                                            project_name,
                                            ArchieParms, 
                                            porosity, 
                                            ERT_meta_dict, 
                                            path_fwd_CATHY, 
                                            Ens_nbi=Ens_nbi,
                                            **kwargs
                                            )
    
    # Read the input mesh using pyvista
    # ------------------------------------
    if DA_cnb is not None:
        mesh_CATHY_ref = pv.read(os.path.join(path_fwd_CATHY, "vtk/100.vtk"))
    else:
        mesh_CATHY_ref = pv.read(os.path.join("vtk/100.vtk"))
        
    # add attribute converted to CATHY mesh
    # ------------------------------------------------------------------------
    mesh_CATHY_new_attr, active_attr = mt.add_attribute_2mesh(
        ER_converted_ti,
        mesh_CATHY_ref,
        "ER_converted" + str(DA_cnb),
        overwrite=True,
        path=path_CATHY,
    )

    if "pygimli" in ERT_meta_dict["data_format"]:
        # copy attribute to pg mesh
        # ------------------------------------------------------------------------
        mesh_geophy_new_attr, scalar_new = CATHY_2_pg(
            mesh_CATHY_new_attr,
            ERT_meta_dict,
            scalar="ER_converted" + str(DA_cnb),
            show=False,
            path=path_CATHY,
            **kwargs
        )
    elif "resipy" in ERT_meta_dict["data_format"]:
        # copy attribute to resipy mesh
        # ------------------------------------------------------------------------
        mesh_geophy_new_attr, scalar_new = CATHY_2_Resipy(
            mesh_CATHY_new_attr,
            ERT_meta_dict,
            scalar="ER_converted" + str(DA_cnb),
            show=False,
            path=path_CATHY,
            **kwargs
        )
    else:
        raise ValueError("Mesh format not recognized")
        

    res0 = mesh_geophy_new_attr.get_array(scalar_new)

    # ------------------------------------------------------------------------
    print('----fwd ER data-----')
    # print(scalar_new)
    # print("Shape:", res0.shape)
    # print("ERT_meta_dict:", ERT_meta_dict)
    # print("sequenceERT:", ERT_meta_dict["sequenceERT"])
    
    # res0 = list(res0)
    # res0 = [int(round(x)) for x in res0]

    # --- Ens_nbi check ---
    if 'Ens_nbi' in locals():
        print(f"Ens_nbi: {Ens_nbi}")
    else:
        print("Warning: Ens_nbi variable not found in scope.")
    
    # --- res0 checks ---
    res0 = np.array(res0, dtype=float)
    
    if res0.size == 0:
        raise ValueError("res0 is empty. Please provide a valid initial resistivity array.")
    
    # Check for NaNs
    if np.any(np.isnan(res0)):
        # Optional: check porosity for this Ens_nbi index if exists
        if 'porosity' in locals() and 'Ens_nbi' in locals():
            try:
                poro_value = porosity[Ens_nbi]
                if np.any(np.isnan(poro_value)):
                    print("Porosity contains NaN values!")
            except IndexError:
                print("Ens_nbi index out of bounds for porosity.")
        raise ValueError("res0 contains NaN values. Please clean or replace them before proceeding.")
    
    # Check for negative values
    if np.any(res0 <= 0):
        raise ValueError("res0 contains negative values. Resistivity must be positive.")
        
    if "pygimli" in ERT_meta_dict["data_format"]:       
        ERT_predicted = simuERT.create_ERT_survey_pg(
            os.path.join(ERT_meta_dict["pathERT"], project_name, "predicted"),
            sequence=ERT_meta_dict["sequenceERT"],
            mesh=ERT_meta_dict["forward_mesh_vtk_file"],
            res0=res0,
            DAcnb=DA_cnb,
            pathfig=path_CATHY,
            **kwargs
        )
    elif "resipy" in ERT_meta_dict["data_format"]:

        ERT_predicted = simuERT.create_ERT_survey_Resipy(
            os.path.join(ERT_meta_dict["pathERT"], project_name, "predicted"),
            sequence=ERT_meta_dict["sequenceERT"],
            mesh=ERT_meta_dict["forward_mesh_vtk_file"],
            res0=res0,
            **kwargs
        )
    # ERT_predicted to dict
    # ------------------------
    d = {
        "a": ERT_predicted["a"],
        "b": ERT_predicted["b"],
        "k": ERT_predicted["k"],
        "m": ERT_predicted["m"],
        "n": ERT_predicted["n"],
        "rhoa": ERT_predicted["rhoa"],
        "valid": ERT_predicted["valid"],
    }
    df_ERT_predicted = pd.DataFrame(data=d)


    # savefig = True
    if savefig:
        print('backup figures')
        plotter = pv.Plotter(shape=(3, 1), off_screen=True)  # notebook = True
        plotter.subplot(0, 0)
        mesh_CATHY_df, name_new_attr_CATHY = mt.add_attribute_2mesh(
            df_sw[-1], mesh_CATHY_ref, "saturation_df", overwrite=True
        )
        mesh_CATHY_df.set_active_scalars("saturation_df")
        my_colormap = "Blues"
        _ = plotter.add_mesh(mesh_CATHY_df,
                             cmap=my_colormap,
                             show_edges=False,
                             # clim=[
                             #     0.3,
                             #     0.7,
                             # ],
                             )
        # plotter.update_scalar_bar_range([0, 1])  # max saturation is 1
        plotter.show_grid()

        plotter.subplot(1, 0)
        mesh_CATHY_new_attr.set_active_scalars(active_attr)
        _ = plotter.add_mesh(
            mesh_CATHY_new_attr,
            show_edges=False,
            clim=[
                min(mesh_CATHY_new_attr[active_attr]),
                max(mesh_CATHY_new_attr[active_attr]),
            ],
        )
        plotter.show_grid()
        plotter.view_xz()

        plotter.subplot(2, 0)
        mesh_geophy_new_attr.set_active_scalars(scalar_new)
        _ = plotter.add_mesh(
            mesh_geophy_new_attr,
            show_edges=False,
            # clim=[
            #     min(mesh_CATHY_new_attr[active_attr]),
            #     max(mesh_CATHY_new_attr[active_attr]),
            # ],
        )

        # if "pygimli" in ERT_meta_dict["data_format"]:
        plotter.view_xy()
        # plotter.view_xz()

        # else:
        #     plotter.view_xz()

        plotter.show_grid()
        plotname = "suplot_ER" + str(DA_cnb)
        plotter.save_graphic(
            path_CATHY + plotname + str(".svg"), title="", 
            raster=True, 
            painter=True
        )

        plotter.close()
    print('end of ER prediction')

    return df_ERT_predicted, df_Archie, mesh_CATHY_new_attr


def Archie_rho_DA(
    rFluid_Archie=[],
    sat=[],
    porosity=[],
    a_Archie=[1.0],
    m_Archie=[2.0],
    n_Archie=[2.0],
    pert_sigma_Archie=[0],
):
    #     This should be moved to DA Class
    """
    Compute ER values at each mesh nodes.
    If pert_sigma, add a nornal noise of sigma standart deviation and 0 mean (# See eq. 4.4 thesis Isabelle p.95)

    Parameters
    ----------
    rFluid : list
        Conductivity of the pore fluid for each soil zone
    sat : list
        Water saturation for each soil zone
    porosity : list
        Porosity for each soil zone
    a : list, optional
        Tortuosity factor. The default is 1.0.
    m : list, optional
        Cementation exponent. The default is 2.0.(usually in the range 1.3 -- 2.5 for sandstones)
    n : list, optional
        Saturation exponent. The default is 2.0.
    pert_sigma_Archie: list, optional
        Normal noise to add to the data

    Returns
    -------
    rho: Converted Apparent Electrical Resistivity (Ohm.m).

    """
    verbose = True
    console = rich.console.Console(stderr=True, quiet=not verbose)
    
    print('------------- Archie transformation -------------')
    print(f'Shape Saturation predicted: {np.shape(sat)}')
    print(f'Shape Porosity nodes: {np.shape(porosity)}')
    print(f'Shape rFluid_Archie: {np.shape(rFluid_Archie)}')

    # Loop over soil type (as many soil type than a_Archie list length)
    # -----------------------------------------------
    for i in range(len(rFluid_Archie)):  # loop over Archie heterogeneity i.e. soil zones
        # print('!!! shortcut sat[:] is not valid for heteregeneous soil!')
        rho = (
            rFluid_Archie[i]
            * a_Archie[i]
            * porosity[i] ** (-m_Archie[i])
            * sat[i] ** (-n_Archie[i])
        )
        sigma = 1 / rho
        
        print(f'min RHO converted SW nodes: {np.min(rho)}')
        print(f'max RHO converted SW nodes: {np.max(rho)}')
        print(f'nb of porosities in the soil: {len(np.unique(porosity))}')
        print(f'sigma: {np.shape(sigma)}')
        # sat[0]
        try: # CASE WITH DATA ASSIMILATION
            for ti in range(np.shape(sigma)[0]):  # Loop over assimilation times
                for meas_nb in range(np.shape(sigma)[1]):  # Loop over mesh nodes
                    # See eq. 4.4 thesis Isabelle p.95
                    # ------------------------------------
                    noise = np.random.normal(
                        0,
                        (pert_sigma_Archie[i] * (sigma[ti, meas_nb])) ** 2,
                        len(sigma[ti, :]),
                    )
                sigma[ti, :] = sigma[ti, :] + noise
                
            if i == 0: # to avoid displaying this for all the soil layers
                console.rule(
                    ":octopus: Parameter perturbation :octopus:", style="green"
                )
                console.print(
                    """
                                Archie perturbation for fwd model \n
                                Archie rFluid: {} \n
                                Archie pert_sigma: {}  \n
                                Nb of zones (need to check this): {}
                                """.format(
                        rFluid_Archie, pert_sigma_Archie, len(rFluid_Archie)
                    ),
                    style="green",
                )
                console.rule("", style="green")
        except:
            if i == 0:
                print('Add noise to Archie converted resistivity - See eq. 4.4 thesis Isabelle p.95')
            for meas_nb in range(len(sigma)):  # Loop over mesh nodes
                # See eq. 4.4 thesis Isabelle p.95
                # ------------------------------------
                noise = np.random.normal(
                    0, (pert_sigma_Archie[i] * (sigma[meas_nb])) ** 2, 1
                )  # See eq. 4.4 thesis Isabelle p.95
                sigma[meas_nb] = sigma[meas_nb] + noise

    # print('max res after Archie')
    # print(np.max(1 / sigma))
    return 1 / sigma
