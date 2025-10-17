#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plotting functions for pyCATHY DA results(2D and 3D)

Plotting:
    - State evolution (SW anf PSI)
    - DA Performance
    - Parameters evolution with assimilation
    - Archie Mapping performance
    - Comparison between scenarios
"""

import os

import numpy as np
from pyCATHY.cathy_utils import (
    change_x2date,
    convert_time_units,
    label_units,
    transform2_time_delta,
)
from pyCATHY.importers import cathy_inputs as in_CT
from pyCATHY.importers import cathy_outputs as out_CT
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.stats as stats

# mpl.style.use('default')
mpl.rcParams["grid.color"] = "k"
mpl.rcParams["grid.linestyle"] = ":"
mpl.rcParams["grid.linewidth"] = 0.25
# mpl.rcParams['font.family'] = 'Avenir'
plt.rcParams["font.size"] = 12
plt.rcParams["axes.linewidth"] = 0.75

# nice_fonts = {
#     # "font.family": "serif",
#     "font.serif": "Times New Roman",
# }
# matplotlib.rcParams.update(nice_fonts)


#%% ---------------------------------------------------------------------------
# -----------------Plot PERFORMANCE DA ----------------------------------------
# -----------------------------------------------------------------------------

def plot_perf(id2plot,obs_df,sensor_match):


    fig, ax = plt.subplots(3,1,sharex=True)
    for sensor in sensor_match:
        prepare_Noble_plots.DA_RMS(results['df_performance'],
                                      sensor,
                                      atmbc_times=list(obs_df.loc[sensor].index),
                                      start_date=simu_args.startD,
                                      ax=ax
                                    )
        prepare_Noble_plots.plot_secondaxis(ax[0],
                                                df_atmbc,
                                                start_date=simu_args.startD,
                                                )

        # savename = (os.path.join(simu_Noble.workdir,simu_Noble.project_name,
        #                          'perf_' + sensor + '.png')
        #          )
    savename = (os.path.join(simu_Noble.workdir,simu_Noble.project_name,
                             'perf_' + sensor + '.png')
             )
    fig.savefig(savename, dpi=450)



#%% ---------------------------------------------------------------------------
# -----------------Plot STATE DA ----------------------------------------------
# -----------------------------------------------------------------------------

def DA_plot_time_dynamic(
    DA, state="psi", nodes_of_interest=[], savefig=False, **kwargs
):
    """Plot result of Data Assimilation: state estimation evolution over the time"""
    prep_DA = prepare_DA_plot_time_dynamic(
        DA, state=state, nodes_of_interest=nodes_of_interest, **kwargs
    )
    if "ax" in kwargs:
        ax = kwargs["ax"]
    else:
        fig = plt.figure(figsize=(6, 3), dpi=350)
        ax = fig.add_subplot()

    alpha = 0.5
    colors_minmax = 'blue'
    if "colors_minmax" in kwargs:
        colors_minmax = kwargs["colors_minmax"]


    keytime = "time"
    xlabel = "time (h)"
    if "start_date" in kwargs:
        keytime = "time_date"
        xlabel = "date"

    ylabel = r"pressure head $\psi$ (m)"
    if "sw" in state:
        ylabel = "water saturation (-)"

    # Plot
    # --------------------------------------------------------------------------

    if len(prep_DA["isENS"]) > 0:
        prep_DA["ens_mean_isENS_time"].pivot(
            index=keytime,
            columns=["idnode"],
            # columns=['idnode'],
            values=["mean(ENS)"],
        ).plot(ax=ax, style=[".-"], color=colors_minmax,
               )

        prep_DA["ens_max_isENS_time"].pivot(
            index=keytime,
            columns=["idnode"],
            # columns=['idnode'],
            values=["max(ENS)"],
        ).plot(ax=ax, style=[".-"], color=colors_minmax,
               )

        prep_DA["ens_min_isENS_time"].pivot(
            index=keytime,
            columns=["idnode"],
            # columns=['idnode'],
            values=["min(ENS)"],
        ).plot(
            ax=ax,
            style=[".-"],
            color=colors_minmax,
            xlabel="(assimilation) time - (h)",
            ylabel="pressure head $\psi$ (m)",
        )

        lgd = ax.fill_between(
            prep_DA["ens_max_isENS_time"][keytime],
            prep_DA["ens_min_isENS_time"]["min(ENS)"],
            prep_DA["ens_max_isENS_time"]["max(ENS)"],
            alpha=0.2,
            color=colors_minmax,
            label="minmax DA",
        )
        # lgd= ax.fill_between(prep_DA['ens_max_isENS_time'][keytime],
        #                 prep_DA['ens_min_isENS_time'][[keytime,'min(ENS)']],
        #                 prep_DA['ens_max_isENS_time'][[keytime,'max(ENS)']],
        #                 alpha=0.2,
        #                 color='blue',
        #                 label='minmax DA')

    # if len(prep_DA['isOL'])>0:
    if "ens_mean_isOL_time" in prep_DA.keys():
        prep_DA["ens_mean_isOL_time"].pivot(
            index=keytime,
            # columns=["Ensemble_nb",'idnode'],
            columns=["idnode"],
            values=["mean(ENS)_OL"],
        ).plot(
            ax=ax,
            style=["-"],
            color="grey",
            label=False,
            ylabel="pressure head $\psi$ (m)",
        )  # water saturation (-)
        prep_DA["ens_min_isOL_time"].pivot(
            index=keytime,
            columns=["idnode"],
            # columns=['idnode'],
            values=["min(ENS)_OL"],
        ).plot(ax=ax, style=["--"], color="grey", label=False)

        prep_DA["ens_max_isOL_time"].pivot(
            index=keytime,
            columns=["idnode"],
            # columns=['idnode'],
            values=["max(ENS)_OL"],
        ).plot(ax=ax, style=["--"], color="grey", label=False)

        ax.fill_between(
            prep_DA["ens_mean_isOL_time"][keytime],
            prep_DA["ens_min_isOL_time"]["min(ENS)_OL"],
            prep_DA["ens_max_isOL_time"]["max(ENS)_OL"],
            alpha=0.2,
            color="grey",
            label="minmax OL",
        )

    ax.set_xlabel(xlabel)
    # ax.set_ylabel('pressure head (m)')
    ax.set_ylabel(ylabel)

    savename = "showDA_dynamic"
    if "savename" in kwargs:
        savename = kwargs["savename"]

    if savefig == True:
        plt.savefig(savename + ".png", dpi=300)
    pass


def DA_plot_ET_dynamic(ET_DA,
                        nodePos=None,
                        nodeIndice=None,
                        observations=None,
                        ax=None,
                        unit='m/s',
                        **kwargs
                        ):

    meanETacolor = 'red'
    if 'color' in kwargs:
        meanETacolor = kwargs.pop('color')
    alphaENS = 0.1
    if 'alphaENS' in kwargs:
        alphaENS = kwargs.pop('alphaENS')
        
    if nodePos is not None:
        # Select data for the specific node
        ET_DA_node = ET_DA.sel(x=nodePos[0],
                               y=nodePos[1],
                               method="nearest"
                               )
        ET_DA_act_etra = ET_DA_node["ACT. ETRA"]
    else:
        ET_DA_act_etra = ET_DA.mean(dim=['x','y'])

    if unit=='mm/day':
        ET_DA_act_etra = ET_DA_act_etra*(1e3*86400)

    # Plot each ensemble member in grey
    ET_DA_act_etra.plot(
        ax=ax,
        x="assimilation",
        hue="ensemble",
        color="grey",
        alpha=alphaENS,
        add_legend=False
    )

    # Plot the mean across ensembles in red
    ET_DA_mean = ET_DA_act_etra.mean(dim="ensemble")
    ET_DA_mean.plot(ax=ax, x="assimilation", 
                    color=meanETacolor,
                    alpha=0.5,
                    linewidth=1,
                    label="Mean pred."
                    )

    if observations is not None:
        if nodePos is not None:
            obs2plot_selecnode = observations.xs(f'ETact{nodeIndice}')[['data','data_err','datetime']]
            obs2plot = obs2plot_selecnode.iloc[:len(ET_DA_mean)][['data','datetime']]
        else:
            import copy
            obs2plot = copy.copy(observations)

        if unit=='mm/day':
            obs2plot['data'] = obs2plot['data']*(1e3*86400)

        ax.scatter(
            obs2plot.datetime[:],
            obs2plot.data[0:len(obs2plot.datetime)],
            label="Observed",
            color='darkgreen',
            s=6
        )
    # ax.set_title('')
    if nodePos is not None:
        ax.set_title(f"Node at ({nodePos[0]}, {nodePos[1]})")

    ax.set_ylabel(f'ETa - {unit}')

# Function to calculate R² and p-value
def calculate_r2_p_value(modelled_data, observed_data):
    corr_coeff, p_value = stats.pearsonr(modelled_data, observed_data)
    r2 = corr_coeff ** 2  # R² value
    return r2, p_value

# Function to annotate the plot with R² and p-value
def annotate_r2_p_value(axi, r2, p_value):
    # annotation_text = f"R² = {r2:.2f}\np-value = {p_value:.2e}"
    annotation_text = f"R² = {r2:.2f}"
    axi.annotate(annotation_text,
                 xy=(0.05, 0.95),
                 xycoords='axes fraction',
                 fontsize=12,
                 ha='left',
                 va='top',
                 bbox=dict(facecolor='white', alpha=0.6, edgecolor='none', boxstyle="round,pad=0.5")
                 )


def DA_plot_ET_performance(ET_DA,
                            observations,
                            axi,
                            nodeposi=None,
                            nodei=None,
                            unit='m/s'
                            ):

    if nodeposi is not None:
    # Select data for the specific node
        ET_DA_node = ET_DA.sel(x=nodeposi[0], y=nodeposi[1], method="nearest")
        obs2plot_selecnode = observations.xs(f'ETact{nodei}')[['data','data_err']]
        # Extract the "ACT. ETRA" variable as a DataArray
        ET_DA_act_etra = ET_DA_node["ACT. ETRA"]
    else:
        ET_DA_act_etra = ET_DA.mean(dim=['x','y'])["ACT. ETRA"]
        obs2plot_selecnode = observations[['data','data_err']].groupby(level=1).mean()

    ET_DA_mean = ET_DA_act_etra.mean(dim="ensemble")
    if unit=='mm/day':
        ET_DA_act_etra = ET_DA_act_etra*(1e3*86400)
        ET_DA_mean = ET_DA_act_etra.mean(dim="ensemble")

    # Plot data for each ensemble member
    for ensi in range(len(ET_DA_act_etra.ensemble)):
        ET_DA_act_etra_ensi = ET_DA_act_etra.isel(ensemble=ensi)

        obs2plot = obs2plot_selecnode.iloc[0:len(ET_DA_act_etra_ensi)].data
        if unit=='mm/day':
            obs2plot = obs2plot*(1e3*86400)

        print(len(obs2plot))
        print(len(ET_DA_act_etra_ensi.values[:]))

        axi.scatter(
            ET_DA_act_etra_ensi.values[:],
            obs2plot.values,
            label=f"Ensemble {ensi}" if ensi == 0 else "",  # Label only once
            color='grey',
            alpha=0.1,
        )

    axi.scatter(
        ET_DA_mean.values[:],
        obs2plot,
        alpha=0.5,
        color='red',
    )

    # Add 1:1 line
    min_val = min(np.nanmin(ET_DA_act_etra.values), np.nanmin(obs2plot_selecnode.data))
    max_val = max(np.nanmax(ET_DA_act_etra.values), np.nanmax(obs2plot_selecnode.data))
    axi.plot([min_val, max_val], [min_val, max_val], "k--", label="1:1 Line")
    axi.set_xlim([min_val, max_val])
    axi.set_ylim([min_val, max_val])


    # Calculate R² and p-value using the separate function
    r2, p_value = calculate_r2_p_value(ET_DA_mean.values, obs2plot)

    # Example usage
    # rmse, nrmse = calculate_rmse_nrmse(ET_DA_mean.values, 
    #                                    obs2plot, 
    #                                    normalization="mean"
    #                                    )
    # print(f'RMSE:{rmse}')
    # print(f'nrmse:{nrmse}')


    # Annotate the plot with R² and p-value using the separate function
    annotate_r2_p_value(axi, r2, p_value)

    # Customize subplot
    if nodeposi is not None:
        axi.set_title(f"Node at ({nodeposi[0]}, {nodeposi[1]})")
    axi.set_xlabel("Modelled ETa")
    axi.set_ylabel("Observed ETa")
    axi.legend()
    axi.set_aspect('equal')
    

#%% ---------------------------------------------------------------------------
# -----------------Plot PARMS DA ----------------------------------------------
# -----------------------------------------------------------------------------


#%% ---------------------------------------------------------------------------
# -----------------Plot ARCHIE DA ---------------------------------------------
# -----------------------------------------------------------------------------
