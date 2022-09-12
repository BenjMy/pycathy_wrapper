#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utilities for pyCATHY
- Unit conversion and labelling
- Others utilities
"""

import pandas as pd
#from datetime import datetime, timedelta
#from rosetta import rosetta, SoilData

#%% ---------------------------------------------------------------------------
#############    Unit conversion and label ####################################
## ---------------------------------------------------------------------------


def kPa2m(kpa):
    return kpa*0.101971621
def kPa2cm(kpa):
    return kpa*0.101971621*1e2

def transform2_time_delta(t,x_units):
    '''
    Time to time delta
    '''
    delta_t = pd.to_timedelta(t,unit=x_units) 
    return delta_t

def convert_time_units(t, x_units):
    '''
    convert time units
    '''
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


def label_units(units,**kwargs):
    '''
    label units
    '''
    if units == "SW":
        label = "Soil Water Content \n ($m^{3}/m^{3}$)"
    elif units == "PH":
        label = "Pressure head (m)"
    elif units == "CKRW":
        label = "Relative hydraulic conductivity"
    elif units == "QTRANIE":
        label = "Root‐zone water uptake"
    else:
        label = ''                
        
    return label


def change_x2date(time_in_sec,start_date):
    ''' change x axis in sec to datetime '''
    date0 = pd.to_datetime(start_date, format='%Y%m%d')
    date_label = [date0] 
    date_label_str = [date0.strftime("%Y-%m-%d %H:%M:%S")]
    for d in time_in_sec[1:]:
        date_label.append(date0 + timedelta(seconds=int(d)))
        date_label_str.append((date0 + timedelta(seconds=int(d))).strftime("%Y-%m-%d %H:%M:%S"))   
    dates = pd.to_datetime(date_label_str,format="%Y-%m-%d %H:%M:%S")       
    return dates


#%% ---------------------------------------------------------------------------
#############  Other utils functions             #############################
## ---------------------------------------------------------------------------

def dictObs_2pd(obs2plot):
    ''' transform Observation dictionnary to panda dataframe for a given observation'''
    
    df_obs = pd.DataFrame.from_dict(obs2plot).stack().to_frame()
    df_obs = pd.DataFrame(df_obs[0].values.T.tolist(), 
                            index=df_obs.index)
    df_obs.index.names= ['sensorNameidx','assimilation time']
    return df_obs

