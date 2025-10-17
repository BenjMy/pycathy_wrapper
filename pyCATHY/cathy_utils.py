#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Utilities for pyCATHY
- Unit conversion and labelling
- Others utilities
"""

import pandas as pd
from shapely.geometry import mapping
import json
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

        

# utils observations
# ------------------


def dictObs_2pd(dict_obs):
    """dict of observation to dataframe of observation"""
    df_obs = pd.DataFrame.from_dict(dict_obs).stack().to_frame()
    df_obs = pd.DataFrame(df_obs[0].values.T.tolist(), index=df_obs.index)
    df_obs.index.names = ["sensorNameidx", "assimilation time"]
    return df_obs


def resynchronise_times(data_measure, atmbc_df, time_offset):
    """resynchronise dict as old key is elapsed time in second from the first observation,
    while new key is from the first atmbc time
    """
    data_measure_sync = dict(data_measure)
    try:
        new_keys = []
        for d in range(len(data_measure_sync.keys())):
            items_dict = list(data_measure_sync.items())
            # new_key = (
            #     atmbc_df["diff"]
            #     .dt.total_seconds()
            #     .unique()[d]
            # )

            old_key = list(data_measure_sync.keys())[d]
            new_key = list(data_measure_sync.keys())[d] - time_offset

            for sensor in list(items_dict[d][1].keys()):
                data_measure_sync[old_key][sensor]["assimilation_times"] = new_key
            new_keys.append(new_key)
        data_measure_sync = dict(zip(new_keys, data_measure_sync.values()))

    except:
        print("datetime first atmbc point")
        print(atmbc_df["datetime"][0])
        print("datetime first measurement")
        print(data_measure[0][sensor]["datetime"])
        print("cant synchronise times - continue without")
    return data_measure_sync



#%% ---------------------------------------------------------------------------
#############  Data Assimilation utils functions          #####################
## ---------------------------------------------------------------------------

import json
COMMON_PARAM_VALUES = {
    "ic":       {"nom": -5,    "mean": -5,   "sigma": 1.75, "bounds": 'None', "type": 'None', "sampling": "normal", "transf": 'None', "update": "St. var."},
    "Ks":       {"nom": 4e-4,  "mean": 1,    "sigma": 0.5,  "bounds": [0,1], "type": "multiplicative", "sampling": "normal", "transf": 'None', "update": "Ks"},
    "porosity": {"nom": 0.35,  "mean": 0.35, "sigma": 0.05, "bounds": [0,1], "type": 'None', "sampling": "normal", "transf": 'None', "update": "porosity"}
}

def create_scenario_file_single_control(
    scenario_name,
    parameters=None,      # List of parameters to include
    use_suggested=True,   # If False, use placeholders
    control_type=None,    # 'layers', 'zone', 'root_map', or None
    nlay=0,
    nzones=0,
    nveg=0,
    filetype="json",
    filename=None
):
    import json

    if filename is None:
        filename = f"{scenario_name}.{filetype}"

    if parameters is None:
        parameters = list(COMMON_PARAM_VALUES.keys())

    def generate_list(value):
        if control_type == "layers":
            return [value] * nlay
        elif control_type == "zone":
            return [value] * nzones
        elif control_type == "root_map":
            return [[value] * nlay for _ in range(nveg)]
        else:
            return value

    per_name, per_type, per_nom, per_mean, per_sigma = [], [], [], [], []
    per_bounds, sampling_type, transf_type, listUpdateParm, pert_control_name = [], [], [], [], []

    for pname in parameters:
        if use_suggested and pname in COMMON_PARAM_VALUES:
            vals = COMMON_PARAM_VALUES[pname]
        else:
            # Placeholder values
            vals = {"nom": None, "mean": None, "sigma": None, "bounds": None, "type": None,
                    "sampling": None, "transf": None, "update": pname}

        per_name.append(pname)
        per_type.append(generate_list(vals["type"]))
        per_nom.append(generate_list(vals["nom"]))
        per_mean.append(generate_list(vals["mean"]))
        per_sigma.append(generate_list(vals["sigma"]))
        per_bounds.append(generate_list(vals["bounds"]))
        sampling_type.append(generate_list(vals["sampling"]))
        transf_type.append(generate_list(vals["transf"]))
        listUpdateParm.append(vals["update"])
        pert_control_name.append(control_type if control_type else None)

    scenario_data = {
        scenario_name: {
            "per_name": per_name,
            "per_type": per_type,
            "per_nom": per_nom,
            "per_mean": per_mean,
            "per_sigma": per_sigma,
            "per_bounds": per_bounds,
            "sampling_type": sampling_type,
            "transf_type": transf_type,
            "listUpdateParm": listUpdateParm,
            "pert_control_name": pert_control_name
        },
        "_explanation": {
            "per_name": "Parameter name",
            "per_type": "Perturbation type: None or multiplicative",
            "per_nom": "Nominal value",
            "per_mean": "Mean value for sampling",
            "per_sigma": "Standard deviation",
            "per_bounds": "Bounds for parameter",
            "sampling_type": "Sampling method: normal, uniform",
            "transf_type": "Transformation type: None, log, etc.",
            "listUpdateParm": "Which parameters are updated",
            "pert_control_name": "Control type: 'layers', 'zone', 'root_map' or None"
        }
    }

    with open(filename, "w") as f:
        json.dump(scenario_data, f, indent=4)

    print(f"Scenario file '{filename}' created successfully with parameters {parameters} and control '{control_type}'.")
    return 

def create_scenario_file(
    scenario_name,
    param_names=None,
    control_type=None,
    filename=None,
    filetype="json",
    use_common_values=True
):
    """
    Create a scenario file using optionally a set of common suggested values.
    
    Parameters
    ----------
    scenario_name : str
        Name of the scenario.
    param_names : list of str, optional
        List of parameter names to include. Defaults to all keys in COMMON_PARAM_VALUES.
    control_type : str or None
        Optional control type for all parameters ('layers', 'zones', 'root_map', or None).
    filename : str, optional
        Output file name.
    filetype : str
        'json' (default) or 'toml'.
    use_common_values : bool
        If True, fill parameter values from COMMON_PARAM_VALUES table.
    """
    if param_names is None:
        param_names = list(COMMON_PARAM_VALUES.keys())
    n = len(param_names)
    
    # Build placeholders
    placeholders = {
        "per_type": ["<type>"] * n,
        "per_nom": ["<value>"] * n,
        "per_mean": ["<mean>"] * n,
        "per_sigma": ["<sigma>"] * n,
        "per_bounds": ["<bounds>"] * n,
        "sampling_type": ["<sampling>"] * n,
        "transf_type": ["<transf>"] * n,
        "listUpdateParm": ["<param>"] * n
    }
    
    # Fill placeholders from common values
    if use_common_values:
        for i, pname in enumerate(param_names):
            vals = COMMON_PARAM_VALUES.get(pname, {})
            placeholders["per_type"][i] = vals.get("type", None)
            placeholders["per_nom"][i] = vals.get("nom", None)
            placeholders["per_mean"][i] = vals.get("mean", None)
            placeholders["per_sigma"][i] = vals.get("sigma", None)
            placeholders["per_bounds"][i] = vals.get("bounds", None)
            placeholders["sampling_type"][i] = vals.get("sampling", None)
            placeholders["transf_type"][i] = vals.get("transf", None)
            placeholders["listUpdateParm"][i] = vals.get("update", pname)
    
    pert_control = [control_type]*n if control_type else [None]*n
    
    scenario_data = {
        scenario_name: {
            "per_name": param_names,
            "per_type": placeholders["per_type"],
            "per_nom": placeholders["per_nom"],
            "per_mean": placeholders["per_mean"],
            "per_sigma": placeholders["per_sigma"],
            "per_bounds": placeholders["per_bounds"],
            "sampling_type": placeholders["sampling_type"],
            "transf_type": placeholders["transf_type"],
            "listUpdateParm": placeholders["listUpdateParm"],
            "pert_control_name": pert_control
        }
    }
    
    if filetype.lower() == "json":
        scenario_data["_explanation"] = {
            "per_type": "Type of parameter: None, multiplicative, etc.",
            "per_name": "Parameter name",
            "per_nom": "Nominal value",
            "per_mean": "Mean for sampling",
            "per_sigma": "Standard deviation",
            "per_bounds": "Bounds: [min,max] or None",
            "sampling_type": "Sampling method",
            "transf_type": "Transformation: None, log, etc.",
            "listUpdateParm": "Parameters updated in model",
            "pert_control_name": "Control type: 'layers', 'zones', 'root_map', or None"
        }
        if filename is None:
            filename = f"{scenario_name}.json"
        with open(filename, "w") as f:
            json.dump(scenario_data, f, indent=4)
    
    elif filetype.lower() == "toml":
        try:
            import tomli_w  # For writing
        except ImportError:
            print("TOML library (tomli_w) not installed. Use: pip install tomli_w")
            print("Falling back to JSON format.")
            # Fallback to JSON
            scenario_data["_explanation"] = {
                "per_type": "Type of parameter: None, multiplicative, etc.",
                "per_name": "Parameter name",
                "per_nom": "Nominal value",
                "per_mean": "Mean for sampling",
                "per_sigma": "Standard deviation",
                "per_bounds": "Bounds: [min,max] or None",
                "sampling_type": "Sampling method",
                "transf_type": "Transformation: None, log, etc.",
                "listUpdateParm": "Parameters updated in model",
                "pert_control_name": "Control type: 'layers', 'zones', 'root_map', or None"
            }
            filename = f"{scenario_name}.json" if filename is None else filename.replace('.toml', '.json')
            with open(filename, "w") as f:
                json.dump(scenario_data, f, indent=4)
            return
        
        if filename is None:
            filename = f"{scenario_name}.toml"
        
        # TOML doesn't support None values, convert to empty string or remove
        def convert_none_values(obj):
            """Recursively convert None values to empty strings for TOML compatibility."""
            if isinstance(obj, dict):
                return {k: convert_none_values(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [convert_none_values(item) for item in obj]
            elif obj is None:
                return ""  # Convert None to empty string
            else:
                return obj
        
        # Convert the data for TOML
        toml_data = convert_none_values(scenario_data)
        
        # FIXED: Use binary mode "wb" for tomli_w
        with open(filename, "wb") as f:
            tomli_w.dump(toml_data, f)
    
    else:
        raise ValueError(f"Unsupported filetype: {filetype}. Use 'json' or 'toml'.")
    
    print(f"Scenario file '{filename}' created successfully for parameters: {param_names}")
    return filename

# def create_scenario_file(
#     scenario_name,
#     param_names=None,
#     control_type=None,
#     filename=None,
#     filetype="json",
#     use_common_values=True
# ):
#     """
#     Create a scenario file using optionally a set of common suggested values.
    
#     Parameters
#     ----------
#     scenario_name : str
#         Name of the scenario.
#     param_names : list of str, optional
#         List of parameter names to include. Defaults to all keys in COMMON_PARAM_VALUES.
#     control_type : str or None
#         Optional control type for all parameters ('layers', 'zones', 'root_map', or None).
#     filename : str, optional
#         Output file name.
#     filetype : str
#         'json' (default) or 'toml'.
#     use_common_values : bool
#         If True, fill parameter values from COMMON_PARAM_VALUES table.
#     """
#     if param_names is None:
#         param_names = list(COMMON_PARAM_VALUES.keys())

#     n = len(param_names)

#     # Build placeholders
#     placeholders = {
#         "per_type": ["<type>"] * n,
#         "per_nom": ["<value>"] * n,
#         "per_mean": ["<mean>"] * n,
#         "per_sigma": ["<sigma>"] * n,
#         "per_bounds": ["<bounds>"] * n,
#         "sampling_type": ["<sampling>"] * n,
#         "transf_type": ["<transf>"] * n,
#         "listUpdateParm": ["<param>"] * n
#     }

#     # Fill placeholders from common values
#     if use_common_values:
#         for i, pname in enumerate(param_names):
#             vals = COMMON_PARAM_VALUES.get(pname, {})
#             placeholders["per_type"][i] = vals.get("type", None)
#             placeholders["per_nom"][i] = vals.get("nom", None)
#             placeholders["per_mean"][i] = vals.get("mean", None)
#             placeholders["per_sigma"][i] = vals.get("sigma", None)
#             placeholders["per_bounds"][i] = vals.get("bounds", None)
#             placeholders["sampling_type"][i] = vals.get("sampling", None)
#             placeholders["transf_type"][i] = vals.get("transf", None)
#             placeholders["listUpdateParm"][i] = vals.get("update", pname)

#     pert_control = [control_type]*n if control_type else [None]*n

#     scenario_data = {
#         scenario_name: {
#             "per_name": param_names,
#             "per_type": placeholders["per_type"],
#             "per_nom": placeholders["per_nom"],
#             "per_mean": placeholders["per_mean"],
#             "per_sigma": placeholders["per_sigma"],
#             "per_bounds": placeholders["per_bounds"],
#             "sampling_type": placeholders["sampling_type"],
#             "transf_type": placeholders["transf_type"],
#             "listUpdateParm": placeholders["listUpdateParm"],
#             "pert_control_name": pert_control
#         }
#     }

#     if filetype.lower() == "json":
#         scenario_data["_explanation"] = {
#             "per_type": "Type of parameter: None, multiplicative, etc.",
#             "per_name": "Parameter name",
#             "per_nom": "Nominal value",
#             "per_mean": "Mean for sampling",
#             "per_sigma": "Standard deviation",
#             "per_bounds": "Bounds: [min,max] or None",
#             "sampling_type": "Sampling method",
#             "transf_type": "Transformation: None, log, etc.",
#             "listUpdateParm": "Parameters updated in model",
#             "pert_control_name": "Control type: 'layers', 'zones', 'root_map', or None"
#         }
#         if filename is None:
#             filename = f"{scenario_name}.json"
#         with open(filename, "w") as f:
#             json.dump(scenario_data, f, indent=4)
#     else:
#         try:
#             # import tomllib
#             import tomli_w  # For writing
#         except ImportError:
#             print("TOML not installed. Use JSON instead.")
#             return
#         if filename is None:
#             filename = f"{scenario_name}.toml"
#         with open(filename, "w") as f:
#             tomli_w.dump(scenario_data, f)

#     print(f"Scenario file '{filename}' created successfully for parameters: {param_names}")
#     return 



    # Save file
    # if filetype.lower() == "toml":
    #     try:
    #         import tomllib
    #     except ImportError:
    #         raise ImportError("Module 'toml' not installed. Use JSON instead.")
    #     with open(filename, "w") as f:
    #         tomllib.dump(scenario_data, f)
    # else:
    #     with open(filename, "w") as f:
    #         json.dump(scenario_data, f, indent=4)

    print(f"Scenario file '{filename}' created successfully with control type '{control_type}'.")



def read_scenario_file(filename, filetype=None):
    """
    Read a scenario file created by `create_scenario_file_single_control`.

    Parameters
    ----------
    filename : str
        Path to the scenario file.
    filetype : str or None
        "json" or "toml". If None, inferred from file extension.

    Returns
    -------
    dict
        Dictionary containing the scenario data.
    """
    if filetype is None:
        if filename.endswith(".toml"):
            filetype = "toml"
        elif filename.endswith(".json"):
            filetype = "json"
        else:
            raise ValueError("Cannot infer file type. Specify filetype='json' or 'toml'.")

    if filetype.lower() == "toml":
        try:
            import tomllib
        except ImportError:
            raise ImportError("Module 'toml' not installed. Use JSON instead.")
        with open(filename, "rb") as f:
            data = tomllib.load(f)
    else:
        with open(filename, "r") as f:
            data = json.load(f)

    return data

def scenario_dict_to_df_list(scenario_dict, scenario_name):
    """
    Convert a scenario dictionary to a DataFrame where each parameter
    is a row and its attributes are stored as lists.
    """
    scenario = scenario_dict[scenario_name]
    params = scenario["per_name"]
    
    data = {}
    for i, pname in enumerate(params):
        data[pname] = {
            "nom": scenario["per_nom"][i] if not isinstance(scenario["per_nom"][i], list) else scenario["per_nom"][i],
            "mean": scenario["per_mean"][i] if not isinstance(scenario["per_mean"][i], list) else scenario["per_mean"][i],
            "sigma": scenario["per_sigma"][i] if not isinstance(scenario["per_sigma"][i], list) else scenario["per_sigma"][i],
            "bounds": scenario["per_bounds"][i] if not isinstance(scenario["per_bounds"][i], list) else scenario["per_bounds"][i],
            "type": scenario["per_type"][i],
            "sampling": scenario["sampling_type"][i],
            "transf": scenario["transf_type"][i],
            "control": scenario.get("pert_control_name", [None]*len(params))[i]
        }
    
    df = pd.DataFrame.from_dict(data, orient="index")
    return df.T

def scenario_dict_to_multiindex_df(scenario_dict, scenario_name):
    """
    Convert a scenario dictionary to a MultiIndex DataFrame:
    - Level 0: per_name
    - Level 1: attributes ('nom', 'mean', 'sigma', etc.)
    Rows represent each instance (layer, zone, root_map element).
    """
    scenario = scenario_dict[scenario_name]
    
    # Determine number of rows
    n_rows = max(len(x) if isinstance(x, list) else 1 for x in scenario["per_nom"])
    
    # Collect data
    data = {}
    for i, pname in enumerate(scenario["per_name"]):
        attrs = {
            "nom": scenario["per_nom"][i],
            "mean": scenario["per_mean"][i],
            "sigma": scenario["per_sigma"][i],
            "bounds": scenario["per_bounds"][i],
            "type": scenario["per_type"][i],
            "sampling": scenario["sampling_type"][i],
            "transf": scenario["transf_type"][i],
            "control": scenario.get("pert_control_name", [None]*len(scenario["per_name"]))[i]
        }
        
        # Expand lists to n_rows
        for key, val in attrs.items():
            if isinstance(val, list):
                val_expanded = val + [None]*(n_rows - len(val))
            else:
                val_expanded = [val]*n_rows
            data[(pname, key)] = val_expanded
    
    # Create DataFrame with MultiIndex columns
    df = pd.DataFrame(data)
    df.columns = pd.MultiIndex.from_tuples(df.columns)
    
    return df.T
