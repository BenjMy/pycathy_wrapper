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
import math

import numpy as np
from pyCATHY.cathy_utils import (
    change_x2date,
    convert_time_units,
    label_units,
    transform2_time_delta,
    get_UNIT_divisor
)
from pyCATHY.importers import cathy_inputs as in_CT
from pyCATHY.importers import cathy_outputs as out_CT
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import scipy.stats as stats
import pandas as pd
import matplotlib.gridspec as gridspec
import logging

logger = logging.getLogger(__name__)

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

#%% ---------------------------------------------------------------------------
# -----------------Figure helpers (date axis, shared legend) ------------------
# -----------------------------------------------------------------------------

LEGEND_RENAME = {
    "Mean pred.": "ETa mean pred.",
    "mean DA": "sw mean pred. (DA)",
    "minmax DA": "sw min-max (DA)",
    "mean OL": "sw mean (open loop)",
    "minmax OL": "sw min-max (open loop)",
    "solution": "sw solution",
}


def format_date_axis(axes, fmt=None, minticks=4, maxticks=7, rotation=30, minor_days=False):
    """Readable date ticks: few labels, short format, rotated.

    Parameters
    ----------
    axes : matplotlib Axes or iterable of Axes
    fmt : str, optional
        ``strftime`` format of the tick labels (e.g. ``'%d %b'``, ``'%b %Y'``).
        Default ``None``: matplotlib's ``ConciseDateFormatter``, which adapts to the
        span (hours for sub-daily data, days, months...) and never repeats a label.
        A fixed daily format such as ``'%d %b'`` gives duplicated labels on a 1-2 day span.
    minticks, maxticks : int
        Bounds on the number of major ticks chosen by ``AutoDateLocator``.
    rotation : float
        Rotation of the tick labels (degrees).
    minor_days : bool
        Draw a minor tick every day (only useful for spans of weeks/months).
    """
    if hasattr(axes, "xaxis"):
        axes = [axes]
    for a in np.ravel(axes):
        loc = mdates.AutoDateLocator(minticks=minticks, maxticks=maxticks)
        a.xaxis.set_major_locator(loc)
        a.xaxis.set_major_formatter(
            mdates.ConciseDateFormatter(loc) if fmt is None else mdates.DateFormatter(fmt)
        )
        if minor_days:
            a.xaxis.set_minor_locator(mdates.DayLocator())
        # tick_params (not setp on the labels) so the rotation survives redraws
        a.tick_params(axis="x", which="major", labelrotation=rotation)
        for lab in a.get_xticklabels():
            lab.set_ha("right" if rotation else "center")


def add_fig_legend(fig, axes, rename=None, ncol=3, fontsize=9):
    """One legend below the figure gathering the entries of all ``axes``.

    Twin axes must be passed explicitly. Labels starting with ``_`` are skipped
    and repeated labels (several nodes) appear once. Call it *after* the last
    ``tight_layout``: it re-runs it leaving room for the legend.

    Parameters
    ----------
    rename : dict, optional
        ``{old_label: new_label}``; defaults to ``LEGEND_RENAME``.
    """
    rename = LEGEND_RENAME if rename is None else rename
    seen = {}
    for a in np.ravel(axes):
        for h, l in zip(*a.get_legend_handles_labels()):
            if l.startswith("_"):
                continue
            seen.setdefault(rename.get(l, l), h)
    labels = list(seen)
    if not labels:
        return None
    ncol = max(1, min(ncol, len(labels)))
    rows = math.ceil(len(labels) / ncol)
    h_frac = (0.24 * rows + 0.12) / fig.get_figheight()
    fig.tight_layout(rect=[0, h_frac, 1, 1])
    return fig.legend(
        [seen[l] for l in labels], labels, loc="lower center",
        ncol=ncol, frameon=False, fontsize=fontsize,
    )


def _wide(df, keytime, col):
    """Long format (keytime, idnode, col) -> wide frame (time x node), sorted by time."""
    return df.pivot(index=keytime, columns="idnode", values=col).sort_index()


def _plot_stat_band(
    ax, df_mean, df_min, df_max, cols, keytime, color, label, mean_fmt, bound_fmt, fill_alpha=0.2
):
    """Mean line, min/max lines and min-max band of an ensemble, one set per node.

    Plain matplotlib calls on numpy arrays (no ``DataFrame.plot``): the x-axis stays
    a regular matplotlib date axis, so ``fill_between`` and any later ``ax.plot``
    with datetimes line up (pandas' period axis made ``fill_between`` fail with a
    DTypePromotionError for regularly spaced dates).
    """
    c_mean, c_min, c_max = cols
    mx = _wide(df_max, keytime, c_max)
    mn = _wide(df_min, keytime, c_min).reindex_like(mx)
    mean = _wide(df_mean, keytime, c_mean).reindex_like(mx)
    x = mx.index.to_numpy()
    multi = mx.shape[1] > 1
    for node in mx.columns:
        suffix = f" (node {node})" if multi else ""
        y_mean, y_min, y_max = (d[node].to_numpy() for d in (mean, mn, mx))
        ax.plot(x, y_mean, mean_fmt, color=color, label=f"mean {label}{suffix}")
        ax.plot(x, y_min, bound_fmt, color=color, label="_nolegend_")
        ax.plot(x, y_max, bound_fmt, color=color, label="_nolegend_")
        ax.fill_between(
            x, y_min, y_max, alpha=fill_alpha, color=color, label=f"minmax {label}{suffix}"
        )


def DA_plot_time_dynamic(
    DA, state="psi", nodes_of_interest=None, savefig=False, legend=True, **kwargs
):
    """Plot result of Data Assimilation: state estimation evolution over the time.

    Parameters
    ----------
    DA : pd.DataFrame
        DA dataframe (``results['df_DA']``).
    state : str
        ``"psi"`` (pressure head) or ``"sw"`` (water saturation).
    nodes_of_interest : int or list
        Node ids to plot.
    savefig : bool
        Save the figure to ``kwargs['savename'] + '.png'``.
    legend : bool
        Draw the legend (as ``DataFrame.plot`` did).
    **kwargs
        ``ax``, ``start_date`` (+ ``atmbc_times``) to plot against dates, ``keytime``,
        ``colors_minmax``, ``savename``; the rest goes to ``prepare_DA_plot_time_dynamic``.
        Extra options (not forwarded):

        - ``solution``: reference curve drawn in red and labelled, either a
          ``pd.Series`` indexed by date or an ``(x, y)`` tuple.
        - ``solution_label``: its legend label (default ``"solution"``).
        - ``date_format``: apply :func:`format_date_axis` when plotting against
          dates (default ``True``; ``False`` to keep matplotlib's default).
    """
    solution = kwargs.pop("solution", None)
    solution_label = kwargs.pop("solution_label", "solution")
    date_format = kwargs.pop("date_format", True)

    prepare = globals().get("prepare_DA_plot_time_dynamic")
    if prepare is None:  # lazy import: avoids loading the heavy plotting module at import
        from pyCATHY.plotters.cathy_plots import prepare_DA_plot_time_dynamic as prepare

    if nodes_of_interest is None:
        nodes_of_interest = []
    prep_DA = prepare(DA, state=state, nodes_of_interest=nodes_of_interest, **kwargs)

    ax = kwargs.get("ax")
    if ax is None:
        fig = plt.figure(figsize=(6, 3), dpi=350)
        ax = fig.add_subplot()

    colors_minmax = kwargs.get("colors_minmax", "blue")

    keytime, xlabel = "time", "time (h)"
    if kwargs.get("start_date") is not None:
        keytime, xlabel = "time_date", "date"
    if "keytime" in kwargs:
        keytime, xlabel = kwargs["keytime"], "assimilation_times"

    ylabel = "water saturation (-)" if "sw" in state else r"pressure head $\psi$ (m)"

    # Plot
    # --------------------------------------------------------------------------
    if len(prep_DA["isENS"]) > 0 and "ens_mean_isENS_time" in prep_DA:
        _plot_stat_band(
            ax,
            prep_DA["ens_mean_isENS_time"],
            prep_DA["ens_min_isENS_time"],
            prep_DA["ens_max_isENS_time"],
            ("mean(ENS)", "min(ENS)", "max(ENS)"),
            keytime,
            colors_minmax,
            "DA",
            mean_fmt=".-",
            bound_fmt=".-",
        )

    if "ens_mean_isOL_time" in prep_DA:
        _plot_stat_band(
            ax,
            prep_DA["ens_mean_isOL_time"],
            prep_DA["ens_min_isOL_time"],
            prep_DA["ens_max_isOL_time"],
            ("mean(ENS)_OL", "min(ENS)_OL", "max(ENS)_OL"),
            keytime,
            "grey",
            "OL",
            mean_fmt="-",
            bound_fmt="--",
        )

    if solution is not None:
        if isinstance(solution, pd.Series):
            x_sol, y_sol = solution.index, solution.values
        else:
            x_sol, y_sol = solution
        ax.plot(np.asarray(x_sol), np.asarray(y_sol), color="r", marker=".",
                markersize=5, alpha=0.8, label=solution_label)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if date_format and keytime == "time_date":
        format_date_axis(ax)
    if legend:
        ax.legend()

    if savefig:
        savename = kwargs.get("savename", "showDA_dynamic")
        ax.figure.savefig(savename + ".png", dpi=300)


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
        ET_DA_act_etra = ET_DA.mean(dim=['x','y'])["ACT. ETRA"]

    if unit=='mm/day':
        ET_DA_act_etra = ET_DA_act_etra*(1e3*86400)

    # Plot each ensemble member in grey (one matplotlib call instead of xarray's
    # hue machinery)
    members = ET_DA_act_etra.transpose(..., "ensemble")
    member_lines = ax.plot(
        members["assimilation"].values,
        members.values,
        color="grey",
        alpha=alphaENS,
    )
    if member_lines:  # a single legend entry for all the members
        member_lines[0].set_label("ensemble members")

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

    if unit == 'mm/day':
        ET_DA_act_etra = ET_DA_act_etra * (1e3 * 86400)
    ET_DA_mean = ET_DA_act_etra.mean(dim="ensemble")

    # observations matching the modelled time steps (same for every member)
    obs2plot = obs2plot_selecnode.iloc[0:len(ET_DA_mean)].data
    if unit == 'mm/day':
        obs2plot = obs2plot * (1e3 * 86400)
    logger.debug(len(obs2plot))

    # All ensemble members in a single scatter call
    members = ET_DA_act_etra.transpose("ensemble", ...).values  # (n_ens, n_times)
    axi.scatter(
        members.ravel(),
        np.tile(obs2plot.values, members.shape[0]),
        label="Ensemble 0",
        color='grey',
        alpha=0.1,
    )

    axi.scatter(
        ET_DA_mean.values[:],
        obs2plot,
        alpha=0.5,
        color='red',
    )

    # Add 1:1 line (same units for model and observations)
    min_val = min(np.nanmin(members), np.nanmin(obs2plot))
    max_val = max(np.nanmax(members), np.nanmax(obs2plot))
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
    



def plot_stats_ensemble_ET(ds, var="ACT. ETRA", output_fig="ensemble_ET_stats.png",
                           time_unit="hours", origin=None):
    """
    Plot ensemble ET statistics:
      1. Ensemble mean ET over time (spatial mean ± std)
      2. Spatial map of ensemble mean at last time step
      3. Spatial map of ensemble std at last time step

    Parameters
    ----------
    ds : xr.Dataset
    var : str
        Variable to plot. Default is 'ACT. ETRA'.
    output_fig : str
        Output figure path.
    time_unit : str
        Time unit for x-axis when no origin is given.
        One of 'seconds', 'minutes', 'hours', 'days'. Default 'hours'.
    origin : str or datetime-like, optional
        If provided (e.g. '2023-06-01 06:00'), timedeltas are converted to
        absolute datetimes and displayed as calendar dates. Overrides time_unit.
    """
    if time_unit not in get_UNIT_divisor():
        raise ValueError(f"time_unit must be one of {list(get_UNIT_divisor().keys())}")

    da = ds[var]
    time_vals = da.time.values  # timedelta64[ns]
    origin_ts = pd.Timestamp(origin) if origin is not None else None

    use_datetime = origin_ts is not None
    if use_datetime:
        time_axis = [origin_ts + pd.to_timedelta(t) for t in time_vals]
        xlabel = "Date / Time"
        last_time_label = time_axis[-1].strftime("%Y-%m-%d %H:%M")
    else:
        time_axis = time_vals.astype("int64") / get_UNIT_divisor()[time_unit]
        xlabel = f"Time ({time_unit})"
        last_time_label = f"{time_axis[-1]:.2f} {time_unit}"

    # Spatial mean across X, Y for each ensemble member and time
    da_spatial_mean = da.mean(dim=["X", "Y"])      # (ensemble, time)
    ens_mean = da_spatial_mean.mean("ensemble")     # (time,)
    ens_std  = da_spatial_mean.std("ensemble")      # (time,)

    # Last time step for spatial maps
    last_time = time_vals[-1]
    da_last  = da.sel(time=last_time)               # (ensemble, Y, X)
    map_mean = da_last.mean("ensemble")             # (Y, X)
    map_std  = da_last.std("ensemble")              # (Y, X)

    fig = plt.figure(figsize=(14, 10))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.4, wspace=0.35)

    # --- 1. Time series: ensemble mean ± std ---
    ax1 = fig.add_subplot(gs[0, :])
    ax1.plot(time_axis, ens_mean.values, "k-o", label="Ensemble mean", linewidth=2)
    ax1.fill_between(
        time_axis,
        ens_mean.values - ens_std.values,
        ens_mean.values + ens_std.values,
        alpha=0.3, color="steelblue", label="±1 std"
    )
    if use_datetime:
        fig.autofmt_xdate(rotation=30)
    else:
        plt.setp(ax1.get_xticklabels(), rotation=30, ha="right", fontsize=9)
    ax1.set_xlabel(xlabel)
    ax1.set_ylabel(var)
    ax1.set_title(
        f"Ensemble mean ± std of spatial-mean '{var}' over time\n"
        f"(n={ds.sizes['ensemble']} members)"
    )
    ax1.legend()
    ax1.grid(True, linestyle="--", alpha=0.5)

    # --- 2. Spatial map: ensemble mean at last time step ---
    ax2 = fig.add_subplot(gs[1, 0])
    im2 = ax2.pcolormesh(
        ds.X.values, ds.Y.values, map_mean.values,
        cmap="YlGn", shading="auto"
    )
    plt.colorbar(im2, ax=ax2, label=var)
    ax2.set_title(f"Ensemble mean\n(t = {last_time_label})")
    ax2.set_xlabel("X"); ax2.set_ylabel("Y")

    # --- 3. Spatial map: ensemble std at last time step ---
    ax3 = fig.add_subplot(gs[1, 1])
    im3 = ax3.pcolormesh(
        ds.X.values, ds.Y.values, map_std.values,
        cmap="OrRd", shading="auto"
    )
    plt.colorbar(im3, ax=ax3, label=f"std({var})")
    ax3.set_title(f"Ensemble spread (std)\n(t = {last_time_label})")
    ax3.set_xlabel("X"); ax3.set_ylabel("Y")

    plt.suptitle(f"Ensemble ET Statistics  —  '{var}'", fontsize=13, fontweight="bold")
    plt.savefig(output_fig, dpi=150, bbox_inches="tight")
    plt.show()
    logger.info(f"Figure saved to: {output_fig}")


# --- Usage examples ---
# Relative time in hours (default):
#   plot_stats_ensemble_ET(ds, var="ACT. ETRA", time_unit="hours")
#
# Relative time in days:
#   plot_stats_ensemble_ET(ds, var="ACT. ETRA", time_unit="days")
#
# Absolute datetime (origin = simulation start):
#   plot_stats_ensemble_ET(ds, var="ACT. ETRA", origin="2023-06-01 06:00")



#%% ---------------------------------------------------------------------------
# -----------------Plot PARMS DA ----------------------------------------------
# -----------------------------------------------------------------------------


def DA_plot_parm_dynamic_scatter(
    parm="ks", dict_parm_pert=None, list_assimilation_times=None, savefig=False, **kwargs
):
    """Boxplot of a parameter ensemble after each update (one box per update).

    Only the ``update_nb*`` entries of ``dict_parm_pert[parm]`` are plotted (the
    initial perturbation is *not* a column), so ``list_assimilation_times`` must
    hold the **dates (or numbers) of the updates**, one per update.

    Parameters
    ----------
    parm : str
        Key of ``dict_parm_pert`` (e.g. ``'ZROOT0'``).
    dict_parm_pert : dict
        ``results['dict_parm_pert']``.
    list_assimilation_times : sequence
        One entry per update. If it has one more entry than there are updates
        the first one is assumed to be the initial time and is dropped (warning).
    **kwargs
        ``ax``, ``color``, ``log`` (log y-scale), ``date_fmt`` (strftime format of
        the box labels, default adapts to the span), ``savename``.

    Returns
    -------
    ax, df : the axis and the DataFrame (ensemble member x update) that was plotted.
    """
    import seaborn  # lazy: only needed here

    if dict_parm_pert is None:
        dict_parm_pert = {}
    ax = kwargs.get("ax")
    if ax is None:
        fig = plt.figure(figsize=(6, 3), dpi=350)
        ax = fig.add_subplot()

    # one column per update (keys like 'update_nb0', 'update_nb1', ...)
    upd_keys = [k for k in dict_parm_pert[parm].keys() if "upd" in k]
    df = pd.DataFrame({k: np.hstack(dict_parm_pert[parm][k]) for k in upd_keys})
    df.index.name = "Ensemble_nb"

    times = pd.Index(list_assimilation_times) if list_assimilation_times is not None else pd.RangeIndex(1, df.shape[1] + 1)
    if len(times) == df.shape[1] + 1:
        logger.warning(
            f"{parm}: {len(times)} times for {df.shape[1]} updates -> dropping the first (initial) time"
        )
        times = times[1:]
    if len(times) != df.shape[1]:
        raise ValueError(
            f"{parm}: {df.shape[1]} updates in dict_parm_pert but {len(times)} times given; "
            "pass one time per update"
        )

    # too many updates: keep one every 15
    if df.shape[1] > 30:
        nii = np.arange(0, df.shape[1], 15)
        df = df.iloc[:, nii]
        times = times[nii]

    is_date = pd.api.types.is_datetime64_any_dtype(times)
    df.columns = range(df.shape[1])  # positional, labels are set below
    boxplot = seaborn.boxplot(df, color=kwargs.get("color", "grey"), ax=ax)

    if is_date:
        fmt = kwargs.get("date_fmt")
        if fmt is None:
            sub_daily = (times.normalize() != times).any() or (times.max() - times.min()) < pd.Timedelta(days=3)
            fmt = "%d %b %H:%M" if sub_daily else "%d %b"
        labels = [t.strftime(fmt) for t in times]
        ax.set_xticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=30, ha="right")
        ax.set_xlabel("assimilation date")
    else:
        ax.set_xticks(range(len(times)))
        ax.set_xticklabels([str(t) for t in times])
        ax.set_xlabel("assimilation #")
    ax.set_ylabel(parm)

    if kwargs.get("log"):
        boxplot.set_yscale("log")
    if savefig:
        ax.figure.savefig(kwargs.get("savename", "showDA_parm") + ".png", dpi=300)
    df.columns = times
    return ax, df



#%% ---------------------------------------------------------------------------
# -----------------Plot ARCHIE DA ---------------------------------------------
# -----------------------------------------------------------------------------
