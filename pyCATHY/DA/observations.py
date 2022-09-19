""" Managing Data Assimilation process observation. 
    Read observation, prepare covariances
    Check consistency
    Prepare for DA class
"""

from pyCATHY.importers import sensors_measures as in_meas
from pyCATHY import cathy_utils as utils_CT
from collections import OrderedDict
import re
import numpy as np

def read_observations(dict_obs, obs_2_add, data_type, data_err, show=False, **kwargs):
    '''
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
        % of error for a given measurement dataset. The default is 0.05.
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

    '''

    dict_obs_2add = {}
    # specify mesh node position for point measurements
    # ---------------------------------------------------------------------
    mesh_nodes = []
    if 'mesh_nodes' in kwargs:
        mesh_nodes = kwargs['mesh_nodes']


    # specify assimilation time if not contained in the file
    # ---------------------------------------------------------------------
    tA = []
    if 'tA' in kwargs:
        tA = kwargs['tA']
    datetime =  None
    if 'datetime' in kwargs:
        datetime = kwargs['datetime']


    # discharge type read
    # ---------------------------------------------------------------------
    if data_type == 'discharge':
        df = in_meas.read_discharge(obs_2_add)
        units = '$m^{3}/s$'


    # weight type read
    # ---------------------------------------------------------------------
    if data_type == 'scale':
        # infer ET from weight data
        if isinstance(obs_2_add, str):
            df = in_meas.read_scale(obs_2_add)
        else:
            df = obs_2_add
            filename = None

        units = '$mm/s$'
        obs_cov_type = None



    # discharge type read
    # ---------------------------------------------------------------------
    elif data_type == 'swc':
        if isinstance(obs_2_add, str):
            df = in_meas.read_swc(obs_2_add)
        else:
            df = obs_2_add
            filename = None
        units = '$m^{3}/m^{3}$'
        obs_cov_type = None
        # point sensors --> covariance between the sensors is formed later

        # if 'date' in df.keys:
        #     date_yyyymmdd_hhmmss = df['date']
        #     dict_obs_2add.update(
        #          date_yyyymmdd_hhmmss = date_yyyymmdd_hhmmss
        #          )


    # tensiometer type read
    # ---------------------------------------------------------------------
    elif data_type == 'tensiometer':
        if isinstance(obs_2_add, str):
            df = in_meas.read_tensiometers(obs_2_add)
        else:
            df = obs_2_add
            filename = None

        units = '$kPa$'
        # convert in m
        # -------------------

        df = utils_CT.kPa2m(df)
        obs_cov_type = None

    # ERT type read
    # ---------------------------------------------------------------------
    elif data_type == 'ERT':
           
        obs_cov_type = 'reciprocal_err'
        if 'obs_cov_type' in kwargs:
            obs_cov_type = kwargs['obs_cov_type']


        data_format = 'resipy'
        if 'data_format' in kwargs:
            data_format = kwargs['data_format']
            dict_obs_2add.update(
                 data_format = data_format
                 )

        elecs = []
        if 'elecs' in kwargs:
            elecs = kwargs['elecs']

        df, dict_ERT = in_meas.read_ERT(obs_2_add,data_format)
        filename = obs_2_add

        units = '$\Omega$'

        if 'elecs' in dict_ERT.keys():
            elecs = dict_ERT['elecs']


        dict_obs_2add.update(
                 elecs = elecs
                 )


    # no file specified (raise error)
    # ---------------------------------------------------------------------
    else:
        print('no file specified')


    dict_obs_2add.update(
                        filename = filename,
                        data_type =  data_type,
                        units= units,  # units
                        data=  df,
                        data_err=  data_err,
                        mesh_nodes =  mesh_nodes,
                        assimilation_times=  tA,
                        datetime =  datetime
                        )


    # add optionnal data metadata to the dictionnary
    # ---------------------------------------------------------------------
    if 'meta' in kwargs:
        meta = kwargs['meta']
        dict_obs_2add.update(meta)



    # data covariance and perturbation (if required)
    # ---------------------------------------------------------------------

    data_cov, dataPert, stacked_data_cov = prepare_observations(
                                                                dict_obs_2add,
                                                                perturbate = False,
                                                                obs_cov_type = obs_cov_type
                                                                )
    dict_obs_2add.update(
                         data_cov = data_cov,
                         dataPert = dataPert
                         )


    dict_obs_2add = OrderedDict(dict_obs_2add)

    # check if assimilation already existing andincrement sensor name if so
    # ---------------------------------------------------------------------
    sensor_name =  data_type

    if tA in dict_obs.keys():
        print('already existing assimilation time')
        for k in dict_obs[tA]:
            k
        if data_type in k:
            match = re.match(r"([a-z]+)([0-9]+)", k, re.I)
            if match:
                items = match.groups()
                it = int(items[1]) +1
            else:
                it = 1
            sensor_name = sensor_name + str(it)
        dict_obs[tA][sensor_name] = {}
        dict_obs_2add.update(sensor_name = sensor_name)
    else:
        dict_obs[tA] = {}
        dict_obs[tA][sensor_name] = {}
        dict_obs_2add.update(sensor_name = sensor_name)

    for key in dict_obs_2add.keys():
        dict_obs[tA][sensor_name][key] = dict_obs_2add[key]


    return dict_obs



def prepare_observations(
                         dict_obs = [],
                         perturbate = False,
                         obs_cov_type = None,
                         Bishop=False, **kwargs
                         ):
    '''
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

    '''

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
    if 'stacked_data_cov' in kwargs:
        stacked_data_cov = kwargs['stacked_data_cov']
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
    '''
    Not yet implemented
    '''
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
    '''
    THIS SHOULD BE MOVED TO DA CLASS

    compute measurement cov matrice and add it to dict_obs

    Returns
    -------
    None.

    '''

    # Attribute to define the data-error covariance matrix.
    # In this example the data noise is assumed to be a scalar
    # that is constant at each observation. We use an array so
    # that the 'shape' of the data covariance is understood by
    # numpy in a linear algebra sense.

    if obs_cov_type == 'reciprocal_err':
        # print(dict_obs)
        data_cov = np.diag(1/abs(dict_obs['data']['recipError'].to_numpy()))
    elif obs_cov_type == 'data_err':
        # resipy format
        try:
            err_arr = np.tile(dict_obs['data_err'],len(dict_obs['data']))
            data_cov = np.diag(err_arr)
        # pygimli format
        except:
            err_arr = np.tile(dict_obs['data_err'],len(dict_obs['data']['a']))
            data_cov = np.diag(err_arr)
    return data_cov