#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utilities for pyCATHY
- Unit conversion and labelling
- Others utilities
"""

import pandas as pd
from shapely.geometry import mapping

from datetime import datetime, timedelta
# from rosetta import rosetta, SoilData

#%% ---------------------------------------------------------------------------
#############    Unit conversion and label ####################################
## ---------------------------------------------------------------------------

def MPa2m(mpa):
    return mpa * 101.99773339984

def kPa2m(kpa):
    return kpa * 0.101971621

def m2kPa(m):
    return m/0.101971621

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
        label = "ETp (rain +ve / evap -ve) \n ($[L^{3}/T]$)"
    elif units == "Atmpot-v":
        label = "ETp volume \n ($[L^{3}]$)"
    elif units == "Atmpot-r":
        label = "ETp rate \n ($[L/T]$)"
    elif units == "Atmpot-d":
        label = "ETp depth \n ($[L]$)"
    elif units == "Atmact-vf":
        label = "Actual infiltration (+ve) or exfiltration (-ve) \n at atmospheric BC nodes \n as a volumetric flux  \n ($[L^3/T]$)"
    elif units == "Atmact-v":
        label = "Actual infiltration (+ve) or exfiltration (-ve) volume \n ($[L^3]$)"
    elif units == "Atmact-r":
        label = "Actual infiltration (+ve) or exfiltration (-ve) rate \n ($[L/T]$)"
    elif units == "Atmact-d":
        label = "Actual infiltration (+ve) or exfiltration (-ve) depth \n ($[L]$)"

    return label


def change_x2date(time_in_sec, start_date,
                  formatIn="%Y%m%d",
                  formatOut="%Y-%m-%d %H:%M:%S"
                  ):
    """change x axis in sec to datetime"""
    
    date0 = pd.to_datetime(start_date, format=formatIn)
    # date0 = pd.to_datetime(start_date)
    # date_label = [date0]
    date_label_str = [date0.strftime("%Y-%m-%d %H:%M:%S")]
    dates = []
    for d in time_in_sec[1:]:
        # date_label.append(date0 + timedelta(seconds=int(d)))
        date_label_str.append(
            (date0 + timedelta(seconds=int(d))).strftime("%Y-%m-%d %H:%M:%S")
        )
        
        # dates.append((date0 + timedelta(seconds=int(d))).strftime("%Y-%m-%d %H:%M:%S"), 
        #              format=formatOut)


    dates = pd.to_datetime(date_label_str, format="%Y-%m-%d %H:%M:%S")
    
    
    # date0 = pd.to_datetime(start_date)
    # Vectorized operation to add seconds directly to the start date
    # dates = date0 + pd.to_timedelta(time_in_sec, unit='s')


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


def backup_simulog_DA(args,filename='DAlog.csv'):
    '''
    Get simulation arguments and save it to a csv file to keep track of the Data
    Assimilation done.
    Returns the index if the simulation is already existing
    '''
    results_df = pd.read_csv(filename,index_col=0)
    now = datetime.now()
    results_df_cols = vars(args).keys()
    results_df_new = pd.DataFrame([vars(args)])
    cols2check = list(vars(args).keys())
    
    values = results_df_new[cols2check].values
    matching_index = results_df.index[(results_df[cols2check] == values).all(axis=1)].tolist()
    if matching_index:
        now = datetime.now()
        results_df.loc[matching_index, 'datetime'] = now
        matching_index = matching_index[0]
    else:
        results_df_new['datetime']=now
        results_df = pd.concat([results_df,results_df_new],ignore_index=True)
        matching_index = len(results_df)-1
    results_df.to_csv(filename)
    return results_df, matching_index


def clip_ET_withLandCover(
    LCnames,
    gdf_AOI,
    ETxr,
    ETname='ACT. ETRA',
    crs_ET=None,
):
    """
    Clips an ET xarray dataset using land cover masks and returns the modified dataset.

    Parameters
    ----------
    LCnames : list of str
        Names of the land cover classes to use as masks, corresponding to the 'POI/AOI' column in `gdf_AOI`.
    gdf_AOI : geopandas.GeoDataFrame
        GeoDataFrame containing the area of interest (AOI) polygons with a 'POI/AOI' column for land cover types.
    ETxr : xarray.Dataset
        The xarray dataset containing the evapotranspiration (ET) data.
    ETname : str, optional
        Name of the ET variable in `ETxr`. Default is 'ACT. ETRA'.
    crs_ET : str, optional
        Coordinate reference system of the ET dataset. If not provided, it is inferred from the dataset.

    Returns
    -------
    xarray.Dataset
        The modified `ETxr` dataset with added layers for each land cover mask under the names "<LCname>_CLCmask".
        
    Notes
    -----
    - Each land cover mask is created using the geometry in `gdf_AOI` for the specified `LCnames`.
    - The `clip` method uses the geometries to create masks applied to the ET dataset.
    - The function returns the modified dataset without generating any plots.

    Examples
    --------
    >>> ETxr_updated = clip_ET_withLandCover(
    ...     LCnames=['Lake', 'Intensive Irrigation'],
    ...     gdf_AOI=gdf,
    ...     ETxr=et_xr,
    ...     ETname='ET',
    ...     crs_ET='EPSG:4326'
    ... )
    """
    for lcn in LCnames:  # axs is still included for flexibility
        CLC_mask = gdf_AOI.set_index('POI/AOI').loc[lcn].geometry
        ETxr = ETxr.rio.write_crs(crs_ET)
        mask_ETA = ETxr[ETname].rio.clip(
            CLC_mask.apply(mapping), 
            crs_ET, 
            drop=False
        )

        ETxr[lcn + '_CLCmask'] = mask_ETA

    return ETxr

        
        