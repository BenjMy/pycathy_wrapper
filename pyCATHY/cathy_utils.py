#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utilities for pyCATHY
- Unit conversion and labelling
- Others utilities
"""

import pandas as pd

from datetime import datetime, timedelta
# from rosetta import rosetta, SoilData

#%% ---------------------------------------------------------------------------
#############    Unit conversion and label ####################################
## ---------------------------------------------------------------------------


def kPa2m(kpa):
    return kpa * 0.101971621


def kPa2cm(kpa):
    return kpa * 0.101971621 * 1e2


def transform2_time_delta(t, x_units):
    """
    Time to time delta
    """
    delta_t = pd.to_timedelta(t, unit=x_units)
    return delta_t


def convert_time_units(t, x_units):
    """
    convert time units
    """
    xlabel = " (s)"
    if x_units == "days":
        xlabel = " (days)"
        t_new = [x / (24 * 60 * 60) for x in t]
        t_new = round(t_new[0], 2)
    if x_units == "hours":
        xlabel = " (h)"
        t_new = [x / (60 * 60) for x in t]
        t_new = round(t_new[0], 1)
    return xlabel, t_new


# Atmpot-vf (9) : Potential atmospheric forcing (rain +ve / evap -ve) as a volumetric flux [L^3/T]
# Atmpot-v (10) : Potential atmospheric forcing volume [L^3] (See parm input file for units)
# Atmpot-r (11) : Potential atmospheric forcing rate [L/T]
# Atmpot-d (12) : Potential atmospheric forcing depth [L]
# Atmact-vf(13) : Actual infiltration (+ve) or exfiltration (-ve) at atmospheric BC nodes as a volumetric flux [L^3/T]
# Atmact-v (14) : Actual infiltration (+ve) or exfiltration (-ve) volume [L^3]
# Atmact-r (15) : Actual infiltration (+ve) or exfiltration (-ve) rate [L/T]
# Atmact-d (16) : Actual infiltration (+ve) or exfiltration (-ve) depth [L]


def label_units(units, **kwargs):
    """
    label units
    """
    if units == "SW":
        label = "Soil Water Content \n ($m^{3}/m^{3}$)"
    elif units == "PH":
        label = "Pressure head (m)"
    elif units == "CKRW":
        label = "Relative hydraulic conductivity"
    elif units == "QTRANIE":
        label = "Root‐zone water uptake"
    elif units == "Atmpot-vf":
        label = "Potential atmospheric forcing (rain +ve / evap -ve) \n ($[L^{3}/T]$)"
    elif units == "Atmpot-v":
        label = "Potential atmospheric forcing volume \n ($[L^{3}/T]$)"
    elif units == "Atmpot-v":
        label = "Potential atmospheric forcing rate \n ($[L/T]$)"
    elif units == "Atmpot-d":
        label = "Potential atmospheric forcing depth \n ($[L]$)"
    elif units == "Atmact-vf":
        label = "Actual infiltration (+ve) or exfiltration (-ve) \n at atmospheric BC nodes \n as a volumetric flux  \n ($[L^3/T]$)"
    elif units == "Atmact-v ":
        label = "Actual infiltration (+ve) or exfiltration (-ve) volume \n ($[L^3]$)"
    elif units == "Atmact-r":
        label = "Actual infiltration (+ve) or exfiltration (-ve) rate \n ($[L/T]$)"
    elif units == "Atmact-d ":
        label = "Actual infiltration (+ve) or exfiltration (-ve) depth \n ($[L]$)"

    return label


def change_x2date(time_in_sec, start_date):
    """change x axis in sec to datetime"""
    date0 = pd.to_datetime(start_date, format="%Y%m%d")
    date_label = [date0]
    date_label_str = [date0.strftime("%Y-%m-%d %H:%M:%S")]
    for d in time_in_sec[1:]:
        date_label.append(date0 + timedelta(seconds=int(d)))
        date_label_str.append(
            (date0 + timedelta(seconds=int(d))).strftime("%Y-%m-%d %H:%M:%S")
        )
    dates = pd.to_datetime(date_label_str, format="%Y-%m-%d %H:%M:%S")
    return dates


#%% ---------------------------------------------------------------------------
#############  Other utils functions             #############################
## ---------------------------------------------------------------------------


def dictObs_2pd(obs2plot):
    """transform Observation dictionnary to panda dataframe for a given observation"""

    df_obs = pd.DataFrame.from_dict(obs2plot).stack().to_frame()
    df_obs = pd.DataFrame(df_obs[0].values.T.tolist(), index=df_obs.index)
    df_obs.index.names = ["sensorNameidx", "assimilation time"]
    return df_obs
