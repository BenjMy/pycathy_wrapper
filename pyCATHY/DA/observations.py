""" Managing Data Assimilation process observation. 
    Read observation, prepare covariances
    Check consistency
    Prepare for DA class
"""

import re
from collections import OrderedDict

import numpy as np

from pyCATHY import cathy_utils as utils_CT
from pyCATHY.importers import sensors_measures as in_meas
from pyCATHY.DA.cathy_DA import dictObs_2pd


def read_observations(dict_obs, obs_2_add, data_type, data_err, show=False, 
                      **kwargs):
    """
    read measures (real observations) from file
    and prepare for DA analysis (covariance matrice and perturbation)
    Uses `pyCATHY.importers`, `sensors_measures` to read files with common standards
    Need to be call for each time/ each observation

    Parameters
    ----------
    dict_obs : dict
        existing observation dictionnary (or empty one)
    obs_2_add : str or dataframe
        filename (or data) of the observation dataset.
    data_type : str
        key tring to identify what type of measure is to read.
    data_err : float
        % of error for a given measurement dataset.
    show : Bool
        plot measure graph. The default is False.

    Returns
    -------
    dict_obs : dict
        dict merging all observations + metadatas.
        First key level is data assimilation time
        Loop over data assimilation times means:
        for ti in self.dict_obs:
            print(ti)


    OrderedDict([(0.0,
                  {'swc': {'filename': None,
                           'data_type': 'swc',
                           'units': '$m^{3}/m^{3}$',
                           'data': 0.3391464625,
                           'data_err': 0.01,
                           'mesh_nodes': [159],
                           'assimilation_times': 0.0,
                           'data_cov': [],
                           'dataPert': [],
                           'sensor_name': 'swc'},
                   'swc1': {'filename': None,
                            'data_type': 'swc',
                            'units': '$m^{3}/m^{3}$',
                            'data': 0.3391464625,
                            'data_err': 0.01,
                            'mesh_nodes': [159],
                            'assimilation_times': 0.0,
                            'data_cov': [],
                            'dataPert': [],
                            'sensor_name': 'swc'},
                   }
                 3600,
                  {'swc': {'filename': None,
                           'data_type': 'swc',
                           'units': '$m^{3}/m^{3}$',
                           'data': 0.3391464625,
                           'data_err': 0.01,
                           'mesh_nodes': [159],
                           'assimilation_times': 0.0,
                           'data_cov': [],
                           'dataPert': [],
                           'sensor_name': 'swc'},
                   'swc1': {'filename': None,
                            'data_type': 'swc',
                            'units': '$m^{3}/m^{3}$',
                            'data': 0.3391464625,
                            'data_err': 0.01,
                            'mesh_nodes': [159],
                            'assimilation_times': 0.0,
                            'data_cov': [],
                            'dataPert': [],
                            'sensor_name': 'swc'},
                   .
                   .
                   .
                 )])

    """

    dict_obs_2add = {}
    # specify mesh node position for point measurements
    # ---------------------------------------------------------------------
    mesh_nodes = []
    if "mesh_nodes" in kwargs:
        mesh_nodes = kwargs["mesh_nodes"]

    # specify assimilation time if not contained in the file
    # ---------------------------------------------------------------------
    tA = []
    if "tA" in kwargs:
        tA = kwargs["tA"]
    datetime = None
    if "datetime" in kwargs:
        datetime = kwargs["datetime"]    
    # discharge type read
    # ---------------------------------------------------------------------
    if data_type == "discharge":
        df = in_meas.read_discharge(obs_2_add)
        units = "$m^{3}/s$"

    # weight type read
    # ---------------------------------------------------------------------
    if data_type == "scale":
        # infer ET from weight data
        if isinstance(obs_2_add, str):
            df = in_meas.read_scale(obs_2_add)
        else:
            df = obs_2_add
            filename = None

        units = "$mm/s$"
        obs_cov_type = None

    # discharge type read
    # ---------------------------------------------------------------------
    elif data_type == "swc":
        if isinstance(obs_2_add, str):
            df = in_meas.read_swc(obs_2_add)
        else:
            df = obs_2_add
            filename = None
        units = "$m^{3}/m^{3}$"
        obs_cov_type = None
        # point sensors --> covariance between the sensors is formed later

        # if 'date' in df.keys:
        #     date_yyyymmdd_hhmmss = df['date']
        #     dict_obs_2add.update(
        #          date_yyyymmdd_hhmmss = date_yyyymmdd_hhmmss
        #          )

    # tensiometer type read
    # ---------------------------------------------------------------------
    elif "tensio" in data_type:
        if isinstance(obs_2_add, str):
            df = in_meas.read_tensiometers(obs_2_add)
        else:
            df = obs_2_add
            filename = None

        units = 'm'
        if "units" in kwargs:
            if kwargs['units'] == "$kPa$":
                # convert in m
                # -------------------
                df = utils_CT.kPa2m(df)
            
        obs_cov_type = None
        
        
    # Actual ET type read
    # ---------------------------------------------------------------------
    elif data_type == "ETact":
        df = obs_2_add
        filename = None
        units = 'm'
        obs_cov_type = None

    # EM (Electromagnetic survey) type read
    # ---------------------------------------------------------------------
    elif data_type == "EM":
        obs_cov_type = "reciprocal_err"
        if "obs_cov_type" in kwargs:
            obs_cov_type = kwargs["obs_cov_type"]
        data_format = "resipy"
        if "data_format" in kwargs['meta']:
            data_format = kwargs['meta']["data_format"]
            dict_obs_2add.update(data_format=data_format)
            
        df, k = in_meas.read_EM(obs_2_add, data_format)
        dict_obs_2add.update(coils=k.coils)
        dict_obs_2add.update(k_emagpy=k)
        filename = obs_2_add
        units = "$mS/m$"
        # df = df

    # ERT type read
    # ---------------------------------------------------------------------
    elif data_type == "ERT":

        obs_cov_type = "reciprocal_err"
        if "obs_cov_type" in kwargs:
            obs_cov_type = kwargs["obs_cov_type"]

        data_format = "resipy"
        if "data_format" in kwargs['meta']:
            data_format = kwargs['meta']["data_format"]
            dict_obs_2add.update(data_format=data_format)
            
        elecs = []
        if "elecs" in kwargs:
            elecs = kwargs.pop('elecs')

        df, dict_ERT = in_meas.read_ERT(obs_2_add, data_format)
        filename = obs_2_add
        dict_obs_2add.update(sequenceERT=filename)

        units = "$\Omega$"

        if "elecs" in dict_ERT.keys():
            elecs = dict_ERT["elecs"]
            
        if obs_cov_type == "reciprocal_err":
           data_err =  df['rec_err'].values

        fwdNoiseLevel = 5
        if "fwdNoiseLevel" in kwargs:
            fwdNoiseLevel = kwargs.pop('fwdNoiseLevel')
            
        sequenceERT = None
        if "sequenceERT" in kwargs:
            sequenceERT = kwargs.pop('sequenceERT')
            # import pygimli as pg
            # shm = pg.load(sequenceERT)
            
        dict_obs_2add.update(elecs=elecs)
        dict_obs_2add.update(fwdNoiseLevel=fwdNoiseLevel)

    # no file specified (raise error)
    # ---------------------------------------------------------------------
    else:
        print("no file specified")

    dict_obs_2add.update(
        filename=filename,
        data_type=data_type,
        units=units,  # units
        data=df,
        data_err=data_err,
        mesh_nodes=mesh_nodes,
        assimilation_times=tA,
        datetime=datetime,
    )

    # add optionnal data metadata to the dictionnary
    # ---------------------------------------------------------------------
    if "meta" in kwargs:
        meta = kwargs["meta"]
        dict_obs_2add.update(meta)

        # print(dict_obs_2add.keys())
        # print(dict_obs_2add)

    # data covariance and perturbation (if required)
    # ---------------------------------------------------------------------

    data_cov, dataPert, stacked_data_cov = prepare_observations(
        dict_obs_2add, perturbate=False, obs_cov_type=obs_cov_type
    )
    dict_obs_2add.update(data_cov=data_cov, dataPert=dataPert)

    dict_obs_2add = OrderedDict(dict_obs_2add)

    # check if assimilation already existing andincrement sensor name if so
    # ---------------------------------------------------------------------
    sensor_name = data_type

    if tA in dict_obs.keys():
        print("already existing assimilation time")
        for k in dict_obs[tA]:
            k
        if data_type in k:
            match = re.match(r"([a-z]+)([0-9]+)", k, re.I)
            if match:
                items = match.groups()
                it = int(items[1]) + 1
            else:
                it = 1
            sensor_name = sensor_name + str(it)
        dict_obs[tA][sensor_name] = {}
        dict_obs_2add.update(sensor_name=sensor_name)
    else:
        dict_obs[tA] = {}
        dict_obs[tA][sensor_name] = {}
        dict_obs_2add.update(sensor_name=sensor_name)

    for key in dict_obs_2add.keys():
        dict_obs[tA][sensor_name][key] = dict_obs_2add[key]

    return dict_obs



def make_data_cov(simu, dict_obs, list_assimilated_obs="all", **kwargs):
    """
    Build covariance matrices for assimilated observations.

   Combine covariance between different observations

   cov_swc = [ 1.  0   0
               0   1.  0
               0   0   1
           ]

   cov_tensio = [  err_tensio      0               0
                   0               err_tensio.     0
                   0               0               err_tensio
               ]

   cov_ert =   [   1.  0   0
                   0   1.  0
                   0   0   1
               ]

   COV_swc_tensio = [  1.  0   0   0   0   0
                       0   1.  0   0   0   0
                       0   0   1   0   0   0
                       0   0   0   1   0   0
                       0   0   0   0   1   0
                       0   0   0   0   0   1
                   ]
   
    Parameters
    ----------
    simu : object
        Simulation object that will store stacked covariance.
    dict_obs : dict
        Dictionary of observations, passed to dictObs_2pd().
    list_assimilated_obs : list or "all", optional
        Which sensors to assimilate. Default = "all".
    kwargs : dict
        Optional arguments (e.g., coils for EM sensors).

    Returns
    -------
    data_cov : np.ndarray
        Covariance matrix (per time step).
    data_pert : np.ndarray
        Perturbed data (from prepare_observations).
    stacked_data_cov : list of np.ndarray
        Time-dependent covariance matrices.
    """

    data_measure_df = dictObs_2pd(dict_obs)
    assimilation_times = np.unique(data_measure_df.index.get_level_values(1))

    def extract_sensor_data(sensor, tA):
        """Helper: extract data and covariance for one sensor at time tA."""
        entry = data_measure_df.xs((sensor, tA))

        if "ERT" in sensor:
            values = np.array(entry["data"]["rhoa"])
        elif "EM" in sensor:
            coils = kwargs.get("coils", None)
            if coils is None:
                raise ValueError("Missing `coils` argument for EM data.")
            values = np.hstack(np.array(entry["data"][coils]))
        else:
            values = [np.array(entry["data"])]

        cov = np.ones(len(values)) * entry["data_err"]**2 # **2 as we need to take the variance
        return values, cov

    # Store covariance per assimilation time
    data_cov_list = []

    for tA in assimilation_times:
        sensors_at_t = data_measure_df.xs(tA, level=1).index

        # Determine which sensors to use
        if list_assimilated_obs == "all":
            selected_sensors = sensors_at_t
        else:
            selected_sensors = [s for s in sensors_at_t if any(l in s for l in list_assimilated_obs)]

        # Collect data and covariance
        all_values, all_cov = [], []
        for sensor in selected_sensors:
            values, cov = extract_sensor_data(sensor, tA)
            all_values.append(values)
            all_cov.append(cov)

        # Flatten and build diagonal covariance
        data_flat = np.hstack(all_values)
        cov_flat = np.hstack(all_cov)

        cov_diag = np.zeros((len(cov_flat), len(cov_flat)))
        np.fill_diagonal(cov_diag, cov_flat)

        data_cov_list.append(cov_diag)

    # Prepare final covariance structures
    data_cov, data_pert, stacked_data_cov = prepare_observations(stacked_data_cov=data_cov_list)
    simu.stacked_data_cov = stacked_data_cov

    return data_cov, data_pert, stacked_data_cov




# def make_data_cov(simu, dict_obs, list_assimilated_obs="all", **kwargs):
#     """
#     Combine covariance between different observations

#     cov_swc = [ 1.  0   0
#                 0   1.  0
#                 0   0   1
#             ]

#     cov_tensio = [  err_tensio      0               0
#                     0               err_tensio.     0
#                     0               0               err_tensio
#                 ]

#     cov_ert =   [   1.  0   0
#                     0   1.  0
#                     0   0   1
#                 ]

#     COV_swc_tensio = [  1.  0   0   0   0   0
#                         0   1.  0   0   0   0
#                         0   0   1   0   0   0
#                         0   0   0   1   0   0
#                         0   0   0   0   1   0
#                         0   0   0   0   0   1
#                     ]

#     Parameters
#     ----------
#     simu : TYPE
#         DESCRIPTION.
#     dict_obs : TYPE
#         DESCRIPTION.

#     Returns
#     -------
#     None.

#     """

#     # merge all dataset and compute covariance matrice
#     # for exemple merging ERT data with SWC data would consist in
#     #
#     #  [R,Theta] --> cov[R,Theta]
#     #
#     # initialize the matrice [meas*meas] size
#     # Here we consider that the covariance doesnt change with time

#     items_dict = list(dict_obs.items())
#     data_measure_df = dictObs_2pd(dict_obs) 
#     nb_assimilation_times = len(np.unique(data_measure_df.index.get_level_values(1)))
#     assimilation_times = np.unique(data_measure_df.index.get_level_values(1))

#     # covariance change with time
#     # ----------------------------------------------------------------
#     data_cov_list = []
#     # for nbt in range(nb_assimilation_times):
#     for tA in assimilation_times:
#         data = []
#         data_ERT = []
#         data_EM = []
#         data_cov = []

#         # self.get_selected_data()
#         # There can be multiple sensors measuring at the same time
#         # -----------------------------------------------------------------
#         if "all" in list_assimilated_obs:
#             for sensor in data_measure_df.xs(tA,level=1).index:
#                 print(sensor)
#                 if "ERT" in sensor:
#                     data_ERT.append(
#                         np.array(data_measure_df.xs((sensor,tA))["data"]["rhoa"])

#                     )
#                     data_cov.append(
#                         np.ones(
#                             len(np.array(data_measure_df.xs((sensor,tA))["data"]["rhoa"]))
#                         )
#                         * data_measure_df.xs((sensor,tA))["data_err"]
#                     )
#                 elif "EM" in sensor:
#                    print('EMs')
#                    coils = kwargs.pop('coils') 
#                    data_EM.append(
#                        np.array(data_measure_df.xs((sensor,tA))["data"][coils])

#                    )
#                    data_cov.append(
#                        np.ones(
#                            len(np.array(data_measure_df.xs((sensor,tA))["data"][coils]))
#                        )
#                        * data_measure_df.xs((sensor,tA))["data_err"]
#                    )
#                 else:
#                     data.append(data_measure_df.xs((sensor,tA))["data"])

#             # data_cov.append(np.ones(len(data)) * items_dict[nbt][1][sensor]["data_err"])
#             data_cov.append(np.ones(len(data)) * data_measure_df.xs((sensor,tA))["data_err"])
#             data = data + data_ERT

#         else:
#             for l in list_assimilated_obs:
#                 # for sensor in items_dict[nbt][1].keys():
#                 for sensor in data_measure_df.xs(tA,level=1).index:
#                     if l in sensor:
#                         if "ERT" in l:
#                             data_ERT.append(
#                                 np.array(data_measure_df.xs((sensor,tA))["data"]["rhoa"])
#                             )
#                             data_cov.append(
#                                 np.ones(
#                                     len(np.array(data_measure_df.xs((sensor,tA))["data"]["rhoa"]))
                     
#                                 )
#                                 * data_measure_df.xs((sensor,tA))["data_err"]
#                             )
#                         elif "EM" in sensor:
#                            print('EMs')
#                            coils = kwargs.pop('coils') 
#                            data_EM.append(
#                                np.array(data_measure_df.xs((sensor,tA))["data"][coils])
        
#                            )
#                            data_cov.append(
#                                np.ones(
#                                    len(np.array(data_measure_df.xs((sensor,tA))["data"][coils]))
#                                )
#                                * data_measure_df.xs((sensor,tA))["data_err"]
#                            )
#                         else:
#                             data.append(data_measure_df.xs((sensor,tA))["data"])

#             # data_cov.append(np.ones(len(data)) * items_dict[nbt][1][sensor]["data_err"])
#             data_cov.append(np.ones(len(data)) * data_measure_df.xs((sensor,tA))["data_err"])
#             data = data + data_ERT

#         data = np.hstack(data)
#         data_cov = np.hstack(data_cov)

#         data_cov_diag = np.zeros([len(data_cov), len(data_cov)])
#         np.fill_diagonal(data_cov_diag, data_cov)
#         data_cov_list.append(data_cov_diag)

#     # # covariance doesnt change with time
#     # # ----------------------------------------------------------------
#     # for nbt in range(nb_assimilation_times):
#     #     data_cov_list.append(data_cov)

#     data_cov, data_pert, stacked_data_cov = prepare_observations(
#         stacked_data_cov=data_cov_list
#     )

#     simu.stacked_data_cov = stacked_data_cov

#     return data_cov, data_pert, stacked_data_cov


def prepare_observations(
    dict_obs=[], perturbate=False, obs_cov_type=None, Bishop=False, **kwargs
):
    """
    prepare observations before DA


    1. Measurement perturbation (if required)
    2. Compute matrice covariance
    3. normalise data (if necessary)

    Parameters
    ----------
    dict_obs : dict
        dict containing measure data + metadata.

    Returns
    -------
    TYPE
        data_cov, DataPert, stacked_data_cov

    """

    # merged data covariance matrice of all the observation data and all times
    stacked_data_cov = []

    # compute individual observation data covariance
    # ----------------------------------------------------
    data_cov = []
    if obs_cov_type is not None:
        data_cov = init_obs_cov_matrice(dict_obs, obs_cov_type)

    # compute individual observation perturbation covariance
    # ----------------------------------------------------
    data_pert = []
    if perturbate == True:
        data_pert = perturbate_obs(dict_obs)

    # stacked data covariance (multiple observations at the same time)
    # ----------------------------------------------------
    if "stacked_data_cov" in kwargs:
        stacked_data_cov = kwargs["stacked_data_cov"]
    else:
        stacked_data_cov.append(data_cov)

    # (Bishop et al., 2001) --> not require the perturbation of observations
    # taking advantage of the high time resolution of the collected data.
    # ---------------------------------------------------------------------
    # The EnKF algorithm implemented here is actually an en-
    # semble transform Kalman filter (Bishop et al., 2001) that
    # does not require the perturbation of observations. On the
    # other hand, the measurement error covariance matrix, R,
    # must be assumed to be known a priori.

    # unit conversion (?)
    # ------------------------------------------------------------------
    # When assimilating multiple variables, proper normaliza-
    # tion of the measurement error covariance matrices, anoma-
    # lies of the simulated data, and innovation vectors were per-
    # formed, using values of 0.6 m, 0.58, and 4.17 × 10 −5 m 3 s −1
    # for pressure head, water content and subsurface outflow,
    # respectively. The normalization ensures that in multivari-
    # ate assimilation scenarios the covariance matrices in the
    # Kalman gain are not ill-conditioned (Evensen, 2003; Cam-
    # porese et al., 2009b).

    # normalisation
    # ------------------------------------------------------------------

    return data_cov, data_pert, stacked_data_cov


def perturbate_obs(dict_obs):
    """
    Not yet implemented
    """
    # if (self.data_cov.ndim == 0):
    #     DataPerturbation=np.sqrt(self.data_cov)*rn.randn(1, self.EnSize)
    # elif (self.data_cov.ndim == 2):
    #     # Compute SVD of Data Covariance to generate noise
    #     U, s, V=np.linalg.svd(self.data_cov, full_matrices=False)
    #     DataPerturbation=np.dot(np.dot(U, np.diag(np.sqrt(s))),
    #                               rn.randn(self.Data.shape[1], self.EnSize))
    # else:
    #     print('Data Covariance should be matrix or scalar', '\n')

    # DataArray=np.tile(
    #     self.Data[i, :], (self.EnSize, 1)).transpose() + DataPerturbation
    pass


def init_obs_cov_matrice(dict_obs, obs_cov_type):
    """
    compute measurement cov matrice and add it to dict_obs

    Returns
    -------
    None.

    """

    # Attribute to define the data-error covariance matrix.
    # In this example the data noise is assumed to be a scalar
    # that is constant at each observation. We use an array so
    # that the 'shape' of the data covariance is understood by
    # numpy in a linear algebra sense.

    if obs_cov_type == "reciprocal_err":
        # print(dict_obs)
        # data_cov = np.diag(1 / abs(dict_obs["data"]["rec_err"].to_numpy()))
        data_cov = np.diag(abs(dict_obs["data"]["rec_err"].to_numpy()))
    elif obs_cov_type == "data_err":
        # resipy format
        try:
            err_arr = np.tile(dict_obs["data_err"], len(dict_obs["data"]))
            data_cov = np.diag(err_arr)
        # pygimli format
        except:
            err_arr = np.tile(dict_obs["data_err"], len(dict_obs["data"]["a"]))
            data_cov = np.diag(err_arr)
    return data_cov
