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


from pyCATHY.cathy_utils import (
    change_x2date,
    convert_time_units,
    label_units,
    transform2_time_delta,
)
from pyCATHY.importers import cathy_inputs as in_CT
from pyCATHY.importers import cathy_outputs as out_CT

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



#%% ---------------------------------------------------------------------------
# -----------------Plot PARMS DA ----------------------------------------------
# -----------------------------------------------------------------------------


#%% ---------------------------------------------------------------------------
# -----------------Plot ARCHIE DA ---------------------------------------------
# -----------------------------------------------------------------------------
