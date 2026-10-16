"""Pedophysical transformations"""

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyvista as pv
import logging

from pyCATHY import meshtools as mt
from pyCATHY.ERT import ERT_checks
from pyCATHY.ERT import simulate_ERT as simuERT
from pyCATHY.importers import cathy_outputs as out_CT
from pyCATHY.meshtools import CATHY_2_pg, CATHY_2_Resipy

logger = logging.getLogger(__name__)


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
    
    logger.debug("Archie parameters: %s", ArchieParms2parse)
    return ArchieParms2parse

def SW_2_ER0_DA(
                project_name, 
                ArchieParms, 
                porosity, 
                meta_dict, 
                path_fwd_CATHY, 
                **kwargs
                ):
    """Convert the saturation of one ensemble member to nodal resistivity (Archie).

    Parameters
    ----------
    project_name : str
    ArchieParms : dict
        Archie parameters; per-member lists are indexed with the member number.
    porosity : array-like
        Nodal porosity (one array per member, or a single array shared by all).
    meta_dict : dict
        ERT metadata (unused here, kept for API symmetry).
    path_fwd_CATHY : str
        Folder of the member simulation ("<prefix>_<n>").
    **kwargs
        ``DA_cnb`` (cycle number), ``Ens_nbi`` (0-based member index),
        ``df_sw`` / ``time_ass`` (see `get_sw_ens_i`).

    Returns
    -------
    ER_converted_ti : np.ndarray
        Resistivity (Ohm.m) at the mesh nodes.
    df_Archie : pandas.DataFrame
        time, ens_nbi, sw, ER_converted, porosity.
    df_sw : np.ndarray
        Saturation used.
    """

    # DA cycle number (None outside a DA run) and member index (0-based; by
    # default read from the ensemble folder name "<prefix>_<n>")
    DA_cnb = kwargs.get("DA_cnb")
    Ens_nbi = kwargs.get("Ens_nbi")
    if Ens_nbi is None:
        Ens_nbi = int(os.path.split(path_fwd_CATHY)[-1].split("_", 1)[1]) - 1

    # Get sw array for a given ensemble
    # ------------------------------------
    df_sw = get_sw_ens_i(path_fwd_CATHY, **kwargs)

    # Choose archie parameter for a given realisation (from the ensemble)
    # --------------------------------------------------------------------
    ArchieParms2parse = get_Archie_ens_i(ArchieParms, Ens_nbi)
    
    # Convert to SW to ER values
    # --------------------------------------------------------------------
 

    
    if isinstance(porosity[0], (list, np.ndarray)):
        porosity = porosity[Ens_nbi]
    else:
        porosity = porosity
        
    # print(f"porosity_ensi: {porosity_ensi}")
    # print(f"porosity_ensi: {porosity_ensi}")
    # (statistics over all the nodes: the former print indexed the nodal porosity
    # with the member number, hence always showed a single value)
    logger.debug(
        "member %d: porosity %s, min=%.3f mean=%.3f max=%.3f",
        Ens_nbi + 1,
        np.shape(porosity),
        np.min(porosity),
        np.mean(porosity),
        np.max(porosity),
    )

    
    ER_converted_ti = Archie_rho_DA(
                                    rFluid_Archie=ArchieParms2parse["rFluid_Archie"],
                                    sat=[df_sw],
                                    porosity=porosity,
                                    a_Archie=ArchieParms2parse["a_Archie"],
                                    m_Archie=ArchieParms2parse["m_Archie"],
                                    n_Archie=ArchieParms2parse["n_Archie"],
                                    pert_sigma_Archie=ArchieParms2parse["pert_sigma_Archie"],
                                )

    ERT_checks.check_archie_inputs(
        df_sw,
        porosity,
        ArchieParms2parse,
        name=f"Archie member {Ens_nbi + 1}",
    ).warn()

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
                                            **{**kwargs, "Ens_nbi": Ens_nbi}  # no duplicate keyword
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
    logger.debug("forward ERT: member %d, DA cycle %s", Ens_nbi + 1, DA_cnb)
    # print(scalar_new)
    # print("Shape:", res0.shape)
    # print("ERT_meta_dict:", ERT_meta_dict)
    # print("sequenceERT:", ERT_meta_dict["sequenceERT"])
    
    # res0 = list(res0)
    # res0 = [int(round(x)) for x in res0]

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
                    logger.error("porosity of member %d contains NaN values", Ens_nbi + 1)
            except IndexError:
                logger.error("member index %d out of bounds for porosity", Ens_nbi)
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

    # checks: forward run (raw data) and forward vs observed quadrupoles / k
    member = f"member {Ens_nbi + 1}"
    check_report = simuERT.get_last_check_report() if "pygimli" in ERT_meta_dict["data_format"] else None
    if check_report is None:
        check_report = ERT_checks.check_forward_data(df_ERT_predicted, res0=res0, name="fwd ERT data")
        check_report.warn()
    check_report.name = f"fwd ERT {member}"
    vs_obs = ERT_checks.check_forward_vs_obs(
        df_ERT_predicted, ERT_meta_dict.get("data"), name=f"fwd vs obs {member}"
    ).warn()
    df_ERT_predicted.attrs["ert_check"] = ERT_checks.ERTCheckReport(
        check_report.name, list(check_report.issues), dict(check_report.stats)
    ).merge(vs_obs).to_dict()


    # savefig = False
    if savefig:
        logger.debug("saving ERT backup figures")
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
    ert_check = df_ERT_predicted.attrs["ert_check"]
    logger.info(
        "ERT forward, member %d (DA cycle %s): res0 %.1f-%.1f Ohm.m -> rhoa %.3g to %.3g Ohm.m | checks: %s",
        Ens_nbi + 1,
        DA_cnb,
        res0.min(),
        res0.max(),
        ert_check["stats"].get("rhoa_min", np.nan),
        ert_check["stats"].get("rhoa_max", np.nan),
        ERT_checks.issues_summary(ert_check["issues"]),
    )

    return df_ERT_predicted, df_Archie, mesh_CATHY_new_attr



# def Archie_rho_DA(
#     rFluid_Archie,
#     sat,
#     porosity,
#     a_Archie=None,
#     m_Archie=None,
#     n_Archie=None,
#     pert_sigma_Archie=None,
#     verbose=True
# ):
#     """
#     Compute apparent electrical resistivity (rho) at each mesh node using Archie’s law.
#     Optionally adds Gaussian noise for data assimilation.

#     Returns a list of resistivity arrays per soil zone. If porosity[i] is None, the corresponding
#     entry in rho_list is None.
#     """
#     console = rich.console.Console(stderr=True, quiet=not verbose)

#     n_soil_units = len(rFluid_Archie)
#     if a_Archie is None:
#         a_Archie = [1.0] * n_soil_units
#     if m_Archie is None:
#         m_Archie = [2.0] * n_soil_units
#     if n_Archie is None:
#         n_Archie = [2.0] * n_soil_units
#     if pert_sigma_Archie is None:
#         pert_sigma_Archie = [0.0] * n_soil_units

#     if verbose:
#         console.rule("------------- Archie transformation -------------")
#         console.print(f"Number of soil zones: {n_soil_units}")
#         for i in range(n_soil_units):
#             console.print(
#                 f"Zone {i}: sat shape {sat[i].shape if sat[i] is not None else None}, "
#                 f"porosity shape {porosity[i].shape if porosity[i] is not None else None}, "
#                 f"rFluid {rFluid_Archie[i]}"
#             )

#     # rho_list = []

#     for i in range(n_soil_units):
#         if porosity[i] is None:
#             # rho_list.append(None)
#             sigma = None
#             continue

#         # Archie law
#         rho = rFluid_Archie[i] * a_Archie[i] * porosity[i] ** (-m_Archie[i]) * sat[i] ** (-n_Archie[i])
#         sigma = 1 / rho

#         # Add Gaussian noise for DA
#         if pert_sigma_Archie[i] > 0:
#             try:
#                 if sigma.ndim == 2:
#                     for zonei in range(sigma.shape[0]):
#                         noise = np.random.normal(
#                             0,
#                             pert_sigma_Archie[i] * sigma[zonei, :],
#                             sigma.shape[1]
#                         )
#                         sigma[zonei, :] += noise
#                 else:
#                     noise = np.random.normal(
#                         0,
#                         pert_sigma_Archie[i] * sigma,
#                         sigma.shape
#                     )
#                     sigma += noise
#             except Exception:
#                 noise = np.random.normal(
#                     0,
#                     pert_sigma_Archie[i] * sigma,
#                     sigma.shape
#                 )
#                 sigma += noise

#             if i == 0 and verbose:
#                 console.rule(":octopus: Parameter perturbation :octopus:", style="green")
#                 console.print(
#                     f"Archie perturbation for fwd model\n"
#                     f"rFluid_Archie: {rFluid_Archie}\n"
#                     f"pert_sigma_Archie: {pert_sigma_Archie}\n"
#                     f"Number of zones: {n_soil_units}",
#                     style="green",
#                 )
#                 console.rule("", style="green")

#         # rho_list.append(1 / sigma)

#     # if verbose:
#         # max_rho = max(np.max(rho) for rho in rho_list if rho is not None)
#         # console.print(f"Max resistivity after Archie: {max_rho}")

#     return 1 / sigma

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
    logger.debug(
        "Archie transformation: sat %s, porosity %s, rFluid %s",
        np.shape(sat),
        np.shape(porosity),
        np.shape(rFluid_Archie),
    )

    # Loop over soil type (as many soil type than a_Archie list length)
    # -----------------------------------------------
    for i in range(len(rFluid_Archie)):  # loop over Archie heterogeneity i.e. soil zones
        # print('!!! shortcut sat[:] is not valid for heteregeneous soil!')
        if porosity[i] is not None: #Case for an ensemble non valid
            rho = (
                rFluid_Archie[i]
                * a_Archie[i]
                * porosity[i] ** (-m_Archie[i])
                * sat[i] ** (-n_Archie[i])
            )
            sigma = 1 / rho

            logger.debug(
                "Archie zone %d: rho from sw nodes min=%.4g max=%.4g (%d distinct porosities)",
                i,
                np.min(rho),
                np.max(rho),
                np.unique(np.asarray(porosity).flatten()).size,
            )

            # Gaussian noise on the conductivity with a standard deviation proportional
            # to the conductivity: std = pert_sigma_Archie[i] * sigma (relative noise, e.g.
            # 0.05 = 5 %). It used to be (pert * sigma)**2 (a variance passed as numpy's
            # scale, i.e. ~1e-6 relative noise for pert=0.05), and the DA branch applied
            # the draw of the last node to the whole row. The perturbed conductivity is
            # kept positive (>= 1 % of the noise-free value).
            pert = pert_sigma_Archie[i]
            if i == 0:  # to avoid displaying this for all the soil layers
                logger.debug(
                    "Archie perturbation for fwd model: rFluid=%s, pert_sigma=%s, nb of zones=%d",
                    rFluid_Archie,
                    pert_sigma_Archie,
                    len(rFluid_Archie),
                )
            if pert:
                sigma_ref = np.asarray(sigma, dtype=float)
                # same code for 1-D (one value per node) and 2-D (assimilation time x node)
                noise = np.random.normal(0.0, 1.0, size=sigma_ref.shape) * pert * np.abs(sigma_ref)
                sigma = np.maximum(sigma_ref + noise, 0.01 * sigma_ref)

    logger.debug("max resistivity after Archie: %.4g", np.max(1 / sigma))
    return 1 / sigma
