""" Managing Data Assimilation process observation. 
    Read observation, prepare covariances
    Check consistency
    Prepare for DA class
"""

import logging
import re
from collections import OrderedDict

import numpy as np

from pyCATHY import cathy_utils as utils_CT
from pyCATHY.ERT import ERT_checks
from pyCATHY.importers import sensors_measures as in_meas
from pyCATHY.cathy_utils import dictObs_2pd

logger = logging.getLogger(__name__)


# the DA loop (DA._obs_key_select) and the covariance (make_data_cov) MUST select the same sensors, in the same order:
# one rule, defined here
def sensor_root_name(sensor):
    """Name of a sensor without its digits: ``swc12`` -> ``swc`` (the key used to select what is assimilated)."""
    return re.sub(r"\d", "", str(sensor))


def select_sensors(sensor_names, list_assimilated_obs="all"):
    """Sensors (given order kept) that are assimilated.

    ``list_assimilated_obs`` is ``"all"``, a sensor type (``"swc"``) or a list of them (``["swc", "ERT"]``); a sensor is
    selected when its name without digits equals one of them. A plain string is ONE type, not a sequence of letters.
    """
    wanted = [list_assimilated_obs] if isinstance(list_assimilated_obs, str) else list(list_assimilated_obs)
    if "all" in wanted:
        return list(sensor_names)
    return [s for s in sensor_names if sensor_root_name(s) in wanted]


def unique_sensor_name(existing, data_type):
    """``data_type`` if free, else ``data_type<k>`` with k one above the highest number already used."""
    existing = {str(k) for k in existing}
    if data_type not in existing:
        return data_type
    pattern = re.compile(rf"^{re.escape(data_type)}(\d+)$")
    used = [int(m.group(1)) for m in map(pattern.match, existing) if m]
    return f"{data_type}{max(used, default=0) + 1}"


def _read_or_pass(obs_2_add, reader_name, data_type):
    """Read a file with ``importers.sensors_measures.<reader_name>``; data given directly (value, dataframe) pass through."""
    if not isinstance(obs_2_add, str):
        return obs_2_add
    reader = getattr(in_meas, reader_name, None)
    if reader is None:
        raise NotImplementedError(
            f"read_observations: no file reader for data_type {data_type!r} "
            f"(`importers.sensors_measures.{reader_name}` does not exist): pass the data itself (value or dataframe)"
        )
    return reader(obs_2_add)


SUPPORTED_DATA_TYPES = ("discharge", "scale", "swc", "tensio*", "ETact", "EM", "ERT")




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
    **kwargs
        ``tA`` (required): assimilation time, ONE number (s). ``mesh_nodes``, ``datetime``, ``meta`` (dict merged into the
        entry, e.g. ``data_format``, ``data_err_unit``), ``obs_cov_type``, ``elecs``, ``units`` ('$kPa$' for tensiometers).

    Notes
    -----
    A second observation of the same ``data_type`` at the same ``tA`` is stored as ``<type>1``, ``<type>2`` ...
    (first free number); an observation is never replaced. The DA loop selects sensors by their name WITHOUT digits
    (``swc12`` -> ``swc``), so a ``data_type`` should not contain digits.
    Raises ``ValueError`` for an unknown ``data_type`` or a missing ``tA``, ``NotImplementedError`` when a file is given
    for a type that has no reader (pass the data itself then).

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

    mesh_nodes = kwargs.get("mesh_nodes", [])
    datetime = kwargs.get("datetime", None)
    tA = kwargs.get("tA", None)
    if tA is None or isinstance(tA, (list, tuple, set, dict, np.ndarray)):
        raise ValueError(f"read_observations: `tA` (the assimilation time, ONE number in s) is required, got {tA!r}")
    meta = kwargs.get("meta") or {}
    if re.search(r"\d", str(data_type)):
        logger.warning("read_observations: data_type %r contains digits; the DA selects sensors by their name without "
                       "digits, so it cannot be picked by name in list_assimilated_obs", data_type)

    dict_obs_2add = {}
    filename = obs_2_add if isinstance(obs_2_add, str) else None
    obs_cov_type = None          # point sensors: the covariance between the sensors is formed later (make_data_cov)

    # one chain: every branch defines df, units, filename, obs_cov_type
    # ---------------------------------------------------------------------
    if data_type == "discharge":
        df = _read_or_pass(obs_2_add, "read_discharge", data_type)
        units = "$m^{3}/s$"

    elif data_type == "scale":                     # infer ET from weight data
        df = _read_or_pass(obs_2_add, "read_scale", data_type)
        units = "$mm/s$"

    elif data_type == "swc":
        df = _read_or_pass(obs_2_add, "read_swc", data_type)
        units = "$m^{3}/m^{3}$"

    elif "tensio" in data_type:
        df = _read_or_pass(obs_2_add, "read_tensiometers", data_type)
        units = "m"
        if kwargs.get("units") == "$kPa$":
            df = utils_CT.kPa2m(df)                # convert in m

    elif data_type == "ETact":
        df = obs_2_add
        units = "m"

    elif data_type in ("EM", "ERT"):
        obs_cov_type = kwargs.get("obs_cov_type", "reciprocal_err")
        data_format = meta.get("data_format", "resipy")
        dict_obs_2add.update(data_format=data_format)      # the DA reads it to know which column is the datum
        filename = obs_2_add

        if data_type == "EM":
            df, k = in_meas.read_EM(obs_2_add, data_format)
            dict_obs_2add.update(coils=k.coils, k_emagpy=k)
            units = "$mS/m$"
        else:
            elecs = kwargs.get("elecs", [])
            df, dict_ERT = in_meas.read_ERT(obs_2_add, data_format)
            ERT_checks.check_observations(df, name=f"ERT obs {obs_2_add}").warn()
            dict_obs_2add.update(sequenceERT=filename)
            units = r"$\Omega$"
            if "elecs" in dict_ERT:
                elecs = dict_ERT["elecs"]
            if obs_cov_type == "reciprocal_err":
                data_err = df["rec_err"].values
            dict_obs_2add.update(elecs=elecs, fwdNoiseLevel=kwargs.get("fwdNoiseLevel", 5))

    else:
        raise ValueError(f"read_observations: unknown data_type {data_type!r}; supported: {', '.join(SUPPORTED_DATA_TYPES)}")

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

    # optional data metadata (data_format, data_err_unit, ...)
    # ---------------------------------------------------------------------
    dict_obs_2add.update(meta)

    # data covariance and perturbation (if required)
    # ---------------------------------------------------------------------
    data_cov, dataPert, stacked_data_cov = prepare_observations(
        dict_obs_2add, perturbate=False, obs_cov_type=obs_cov_type
    )
    dict_obs_2add.update(data_cov=data_cov, dataPert=dataPert)

    # a second sensor of the same type at the same time gets the next free name; nothing is ever replaced
    # ---------------------------------------------------------------------
    sensors_at_tA = dict_obs.setdefault(tA, {})
    sensor_name = unique_sensor_name(sensors_at_tA.keys(), data_type)
    if sensor_name != data_type:
        logger.debug("assimilation time %s already has a %s: stored as %s", tA, data_type, sensor_name)
    dict_obs_2add.update(sensor_name=sensor_name)
    sensors_at_tA[sensor_name] = dict(dict_obs_2add)

    return dict_obs



def obs_error_variance(values, data_err, unit=None):
    """Variance of the observation errors (diagonal of R).

    ``unit`` selects how ``data_err`` is interpreted:

    - ``None`` / ``"absolute"`` (default, historical): ``data_err`` is a standard
      deviation in the units of the data (Ohm.m for rhoa), ``var = data_err**2``.
    - ``"percent"`` (alias ``"relative"``): ``data_err`` is a relative error in %,
      ``var = (data_err/100 * |d|)**2``. A fixed 5 Ohm.m is 5 % of 100 Ohm.m but
      33 % of 15 Ohm.m, so percent is usually what an "error in %" means.

    ``data_err`` may be a scalar or one value per datum.
    """
    values = np.asarray(values, dtype=float)
    err = np.asarray(data_err, dtype=float)
    if unit in (None, "absolute"):
        return np.ones(values.shape) * err**2
    if unit in ("percent", "relative"):
        return (err / 100.0 * np.abs(values)) ** 2
    raise ValueError(f"unknown data_err_unit {unit!r}: use 'absolute' or 'percent'")


def _sensor_values(sensor, entry, coils=None):
    """Data of one sensor as a flat array, in the SAME layout as ``DA._get_data2assimilate``."""
    data = entry["data"]
    if "ERT" in sensor:
        fmt = str(entry.get("data_format") or "")
        # the DA assimilates `resist` for resipy and `rhoa` for pygimli data (rhoa when the format is not recorded)
        for col in (("rhoa", "resist") if (not fmt or "pygimli" in fmt) else ("resist", "rhoa")):
            try:
                return np.asarray(data[col]).ravel()
            except KeyError:
                continue
        raise KeyError(f"ERT sensor {sensor!r}: neither 'rhoa' nor 'resist' in its data")
    if "EM" in sensor:
        coils = coils if coils is not None else entry.get("coils")
        if coils is None:
            raise ValueError("Missing `coils` argument for EM data (not in the observation entry either).")
        return np.hstack(np.asarray(data[coils]))
    return np.atleast_1d(np.asarray(data)).ravel()


def make_data_cov(simu, dict_obs, list_assimilated_obs="all", **kwargs):
    """
    Build the (diagonal) data-error covariance of every assimilation time.

    For each time, in chronological order (the DA refuses an unsorted ``dict_obs``), the sensors selected by
    ``list_assimilated_obs`` are stacked in the order they have in ``dict_obs[tA]``: this is the order of the data vector
    and of the predictions of ``DA._obs_key_select``, so R(i, i) belongs to datum i.

    The variance of a sensor is ``data_err**2`` (``data_err`` = standard deviation, scalar or one value per datum), or
    ``(data_err/100 * |d|)**2`` when the entry has ``data_err_unit = 'percent'`` (see :func:`obs_error_variance`).

    Parameters
    ----------
    simu : object
        Simulation object that will store ``stacked_data_cov``.
    dict_obs : dict
        ``{tA: {sensor_name: entry}}`` built by :func:`read_observations`.
    list_assimilated_obs : "all", str or list, optional
        Which sensors to assimilate: ``"all"``, a type (``"swc"``) or a list (``["swc", "ERT"]``); a sensor is selected when
        its name without digits is in the list (same rule as the DA loop, :func:`select_sensors`).
    kwargs : dict
        ``coils`` (EM; default: the ``coils`` stored in the entry). ``diagonal=True`` stores each cycle's covariance as its
        diagonal (a vector of m values) instead of a dense (m, m) matrix: ``cycles x m^2`` values become ``cycles x m``. The
        analysis rebuilds the dense matrix of one cycle at a time, so results are identical. Default False (dense).

    Returns
    -------
    data_cov, data_pert, stacked_data_cov : the first two are empty (kept for the historical signature); the third is the
        list of covariances, one per assimilation time, also stored in ``simu.stacked_data_cov``.

    Raises
    ------
    ValueError : no sensor selected at a time, ``data_err`` of the wrong length, EM without coils.
    """
    diagonal = bool(kwargs.get("diagonal", False))
    data_cov_list = []

    for tA in sorted(dict_obs):
        sensors_at_t = dict_obs[tA]
        selected = select_sensors(list(sensors_at_t), list_assimilated_obs)
        if not selected:
            raise ValueError(f"make_data_cov: no sensor of {list_assimilated_obs!r} at tA={tA} "
                             f"(sensors there: {list(sensors_at_t)})")

        data_parts, cov_parts = [], []
        for sensor in selected:
            entry = sensors_at_t[sensor]
            values = _sensor_values(sensor, entry, kwargs.get("coils"))
            err = np.asarray(entry["data_err"], dtype=float)
            if err.ndim > 0 and err.size not in (1, len(values)):
                raise ValueError(f"make_data_cov: sensor {sensor!r} at tA={tA} has {err.size} data_err values "
                                 f"for {len(values)} data")
            unit = entry.get("data_err_unit")
            if isinstance(unit, float) and np.isnan(unit):
                unit = None
            if unit is None or unit == "absolute":
                cov = np.ones(len(values)) * err ** 2            # variance (not std) of the errors
            else:
                cov = obs_error_variance(values, err, unit)
            data_parts.append(values)
            cov_parts.append(cov)

        data_flat = np.hstack(data_parts)
        cov_flat = np.hstack(cov_parts)
        if any("ERT" in s for s in selected):
            ERT_checks.check_data_cov(cov_flat, values=data_flat, name=f"ERT data cov tA={tA}").warn()

        if diagonal:
            # keep only the diagonal (m values): the dense matrix of a cycle is rebuilt when
            # that cycle is analysed (see DA._data_cov_of_cycle)
            data_cov_list.append(np.asarray(cov_flat, dtype=float))
        else:
            data_cov_list.append(np.diag(np.asarray(cov_flat, dtype=float)))

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
    dict_obs=None, perturbate=False, obs_cov_type=None, Bishop=False, **kwargs
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

    if dict_obs is None:
        dict_obs = {}

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
        # variance = rec_err**2, rec_err being an absolute error in the units of the data
        # (Ohm.m for rhoa), as in make_data_cov (the DA path); it used to be rec_err itself
        data_cov = np.diag(dict_obs["data"]["rec_err"].to_numpy() ** 2)
    elif obs_cov_type == "data_err":
        # resipy format
        try:
            err_arr = np.tile(dict_obs["data_err"], len(dict_obs["data"]))
            data_cov = np.diag(err_arr**2)
        # pygimli format
        except Exception:  # pygimli container: no len() and no dict access
            err_arr = np.tile(dict_obs["data_err"], len(dict_obs["data"]["a"]))
            data_cov = np.diag(err_arr**2)
    return data_cov
