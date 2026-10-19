"""Performance table of the assimilation (RMSE per sensor): computation and parquet storage.
Moved out of ``cathy_DA`` without change; every name is still importable from
``pyCATHY.DA.cathy_DA``."""
from __future__ import annotations

import os

import numpy as np
import pandas as pd


def _compute_diff_matrix(data, prediction):
    """Absolute difference between data (n_obs,) and prediction (n_obs, n_ens)."""
    data = np.atleast_1d(np.asarray(data, dtype=float))  # () -> (1,)
    prediction = np.asarray(prediction, dtype=float)
    if prediction.ndim == 1:  # (n_ens,) -> (1, n_ens)
        prediction = prediction.reshape(len(data), -1)
    if data.shape[0] != prediction.shape[0]:  # never broadcast one datum over all rows
        raise ValueError(
            f"{data.shape[0]} observations but {prediction.shape[0]} predictions"
        )
    return np.abs(data[:, None] - prediction)  # broadcast, no loops


def _compute_rmse(diff_matrix):
    """Mean absolute difference: over ensemble members (axis 1), then over observations."""
    return float(np.mean(np.nanmean(diff_matrix, axis=1)))


def _get_rmse_for_sensor(obs2eval, prediction):
    """Computes RMSE for a given sensor."""
    diff_mat = _compute_diff_matrix(obs2eval, prediction)
    return _compute_rmse(diff_mat)


def _append_to_performance_df(
    df_performance,
    t_obs,
    name_sensor,
    rmse_sensor,
    rmse_avg,
    nrmse_sensor,
    nrmse_avg,
    ol_bool,
):
    """Appends performance metrics to the main DataFrame.

    The wide layout and the ``NMRMSE`` column names are legacy: they are read by
    ``pyCATHY.plotters.cathy_plots.DA_RMS`` and by existing parquet files.
    Change them only together with that plotter.
    """
    data = {
        "time": [t_obs],
        "ObsType": [name_sensor],
        "RMSE" + name_sensor: [rmse_sensor],
        "RMSE_avg": [rmse_avg],
        "NMRMSE" + name_sensor: [nrmse_sensor],
        "NMRMSE_avg": [nrmse_avg],
        "OL": [ol_bool],
    }
    return pd.concat([df_performance, pd.DataFrame(data)], ignore_index=True)


def _save_performance_to_parquet(df_performance, file_path="performance_data.parquet"):
    """Saves the performance DataFrame to Parquet format, appending if the file exists."""
    if os.path.exists(file_path):
        df_existing = pd.read_parquet(file_path, engine="pyarrow")
        df_performance = pd.concat([df_existing, df_performance], ignore_index=True)

    df_performance.to_parquet(
        file_path, engine="pyarrow", compression="snappy", index=False
    )


def _load_performance_from_parquet(file_path="performance_data.parquet"):
    """Loads performance data from a Parquet file if it exists."""
    if os.path.exists(file_path):
        return pd.read_parquet(file_path, engine="pyarrow")
    else:
        return pd.DataFrame()  # Return an empty DataFrame if no file exists
