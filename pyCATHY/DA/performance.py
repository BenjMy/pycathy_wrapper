"""Class managing data Assimilation performance assessement
"""
import pandas as pd
import numpy as np
from pyCATHY.cathy_tools import CATHY


def _append_to_performance_df(self,df_performance, t_obs, name_sensor,
                              rmse_sensor, rmse_avg, nrmse_sensor, nrmse_avg,
                              ol_bool):
    """Appends performance metrics to the main DataFrame."""
    data = {
        "time": [t_obs],
        "ObsType": [name_sensor],
        "RMSE" + name_sensor: [rmse_sensor],
        "RMSE_avg": [rmse_avg],
        "NMRMSE" + name_sensor: [nrmse_sensor],
        "NMRMSE_avg": [nrmse_avg],
        "OL": [ol_bool]
    }
    return pd.concat([df_performance, pd.DataFrame(data)], ignore_index=True)
        
#%%
def _save_performance_to_parquet(self,df_performance, file_path="performance_data.parquet"):
    """Saves the performance DataFrame to Parquet format, appending if the file exists."""
    # if os.path.exists(file_path):
        # df_performance.to_parquet(file_path, engine="pyarrow",
        #                           compression="snappy",
        #                           index=False,
        #                           append=True
        #                           )

    if os.path.exists(file_path):
        df_existing = pd.read_parquet(file_path, engine="pyarrow")
        df_performance = pd.concat([df_existing, df_performance], ignore_index=True)

    df_performance.to_parquet(file_path, engine="pyarrow", compression="snappy", index=False)


def _load_performance_from_parquet(self,file_path="performance_data.parquet"):
    """Loads performance data from a Parquet file if it exists."""
    if os.path.exists(file_path):
        return pd.read_parquet(file_path, engine="pyarrow")
    else:
        return pd.DataFrame()  # Return an empty DataFrame if no file exists

def _performance_assessement_pq(self,
                                list_assimilated_obs,
                                data,
                                prediction,
                                t_obs,
                                # file_path="performance_data.parquet",
                                **kwargs
                                ):
    """Calculates (N)RMSE for each observation in a data assimilation cycle."""

    ol_bool = kwargs.get("openLoop", False)
    file_path = self.project_name + '_df_performance.parquet'
    # Initialize df_performance if file exists
    if hasattr(self, 'df_performance') is False:
        self.df_performance = pd.DataFrame()
    else:
        self.df_performance = self._load_performance_from_parquet(os.path.join(
                                                                     self.workdir,
                                                                     self.project_name,
                                                                     file_path,
                                                                     )
                                                        )
    rmse_avg_stacked = []  # This list will track average RMSE values over time steps

    num_predictions = len(prediction)
    all_obs_diff_mat = self._compute_diff_matrix(data, prediction)
    all_obs_rmse_avg = self._compute_rmse(all_obs_diff_mat, num_predictions)

    # Example function to retrieve list of keys for assimilated observations
    obs2eval_key = self._get_data2assimilate(list_assimilated_obs)[1]
    start_line_obs = 0

    for name_sensor in obs2eval_key:
        obs2eval = self._get_data2assimilate([name_sensor], match=True)[0]
        n_obs = len(obs2eval)
        prediction2eval = prediction[start_line_obs:start_line_obs + n_obs]

        rmse_sensor = self._get_rmse_for_sensor(obs2eval, prediction2eval, num_predictions)
        rmse_avg_stacked.append(all_obs_rmse_avg)

        if "ObsType" in self.df_performance and name_sensor in self.df_performance["ObsType"].unique():
            prev_rmse = self.df_performance[self.df_performance["ObsType"] == name_sensor]["RMSE" + name_sensor].sum()
            rmse_sensor_stacked = prev_rmse + rmse_sensor
        else:
            rmse_sensor_stacked = rmse_sensor

        if t_obs == 0:
            nrmse_sensor = rmse_sensor
            nrmse_avg = all_obs_rmse_avg
        else:
            nrmse_sensor = rmse_sensor_stacked / (t_obs + 1)
            nrmse_avg = np.sum(rmse_avg_stacked) / (t_obs + 1)

        # t_obs = 3
        self.df_performance = self._append_to_performance_df(
                                                            self.df_performance,
                                                            t_obs,
                                                            name_sensor,
                                                            rmse_sensor,
                                                            all_obs_rmse_avg,
                                                            nrmse_sensor,
                                                            nrmse_avg,
                                                            ol_bool
                                                            )

        start_line_obs += n_obs

    # Save to Parquet and clear memory if too large
    # if len(self.df_performance) > 10000:  # Adjust batch size as needed
    self._save_performance_to_parquet(self.df_performance,
                                          os.path.join(
                                            self.workdir,
                                            self.project_name,
                                            file_path,
                                            )
                                 )
    self.df_performance = pd.DataFrame()  # Clear in-memory DataFrame

    # return df_performance


def _performance_assessement(
                            self,
                            list_assimilated_obs,
                            data,
                            prediction,
                            t_obs,
                            **kwargs
                            ):
    """
    (Normalized) root mean square errors (NRMSEs)
    RMSE is compute separately for each observation assimilated


    ----------
    list_assimilated_obs : TYPE
        DESCRIPTION.
    Data : TYPE
        Refers to the measured data.
    Observation : TYPE
        Refers to the simulated data.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """

    OL_bool = False
    if "openLoop" in kwargs:
        OL_bool = True

    if hasattr(self, "df_performance") is False:
        # initiate if simulation just started
        # -------------------------------------------------------------
        self.df_performance = pd.DataFrame()
        self.RMSE_avg_stacked = []

        # NAs default
        # -------------------------------------------------------------
        # RMSE_avg = np.nan # root mean square error (RMSE)
        # NMRMSE = np.nan # time-averaged normalized root mror: unexpected indent

    # average differences over the number of ensemble
    # ALL OBSERVATIONS
    # ------------------------------

    all_Obs_diff_mat = np.zeros(np.shape(prediction))
    for i in range(len(data)):
        for j in range(np.shape(prediction)[1]):  # Loop over ensemble collumns
            all_Obs_diff_mat[i, j] = abs(data[i] - prediction[i, j])
            # all_Obs_diff_mat[:,j] = abs(data[i]-prediction[i,j])

    all_Obs_diff_avg = np.nansum(all_Obs_diff_mat, axis=1) * (1 / len(prediction))

    # average differences over all the observations
    all_Obs_RMSE_avg_ti = np.sum(all_Obs_diff_avg, axis=0) * (1 / len(data))

    # # compute metrics for each observation variable
    # # ------------------------------------------
    # Loop trought observation dictionnary for a given assimilation time (count_DA_cycle)
    # -----------------------------------------------------------------
    obs2eval_key = []  # data type 2 eval
    obs2eval = []  # data dict 2 eval

    _, obs2eval_key = self._get_data2assimilate(list_assimilated_obs)

    start_line_obs = 0
    for s, name_sensor in enumerate(obs2eval_key):  # case where there is other observations than ERT
        obs2eval, _ = self._get_data2assimilate([name_sensor], match=True)

        # number of observation at a given time
        # for ERT number of observation = is the number of grid cells
        n_obs = len(obs2eval)
        prediction2eval = prediction[start_line_obs : start_line_obs + n_obs]
        obs2eval_diff_mat = np.zeros(np.shape(prediction2eval))

        for i in range(len(obs2eval)):
            for j in range(
                np.shape(prediction2eval)[1]
            ):  # Loop over ensemble collumns
                obs2eval_diff_mat[i, j] = abs(obs2eval[i] - prediction2eval[i, j])

        obs2eval_diff_avg = np.nansum(obs2eval_diff_mat, axis=1) * (
            1 / len(prediction2eval)
        )

        start_line_obs = start_line_obs + n_obs

        if "ERT" in name_sensor:
            RMSE_sensor_ti = np.sum(obs2eval_diff_avg, axis=0) * (1 / len(data))
            # Plot here scatter Transfer resistance scatter plot between the observed and simulated data
            # plt.scatter(obs2eval[i],prediction2eval[i,j])
        else:
            RMSE_sensor_ti = obs2eval_diff_avg

        # compute normalised RMSE
        # ---------------------------------------------------------

        self.RMSE_avg_stacked.append(all_Obs_RMSE_avg_ti)

        if 'ObsType' in self.df_performance:
            if name_sensor in self.df_performance['ObsType'].unique():
                bool_sensor = self.df_performance['ObsType']==name_sensor
                RMSE_sensor_stacked = (self.df_performance[bool_sensor]["RMSE" + name_sensor].sum()
                                   + RMSE_sensor_ti
                                   )
            else:
                RMSE_sensor_stacked = RMSE_sensor_ti

        # if hasattr(self, 'RMSE') is False:
        if t_obs == 0:
            NMRMSE_sensor_ti = RMSE_sensor_ti
            NMRMSE_avg_ti = all_Obs_RMSE_avg_ti
        else:
            NMRMSE_sensor_ti = (1 / (t_obs + 1)) * RMSE_sensor_stacked
            NMRMSE_avg_ti = (1 / (t_obs + 1)) * np.sum(self.RMSE_avg_stacked)

        # root names for the collumns name
        # -------------------------------------------------------------
        cols_root = [
            "time",
            "ObsType",
            "RMSE" + name_sensor,
            "RMSE_avg",
            "NMRMSE" + name_sensor,
            "NMRMSE_avg",
            "OL",
        ]

        # root data
        # -------------------------------------------------------------
        data_df_root = [
            [t_obs],
            [name_sensor],
            [RMSE_sensor_ti],
            [all_Obs_RMSE_avg_ti],
            [NMRMSE_sensor_ti],
            [NMRMSE_avg_ti],
            [OL_bool],
        ]

        df_performance_ti = pd.DataFrame(
            np.array(data_df_root, dtype=object).T, columns=cols_root
        )

        # concatenate with main RMSE dataframe
        # -------------------------------------------------------------
        self.df_performance = pd.concat(
            [self.df_performance, df_performance_ti], axis=0, ignore_index=True
        )

    return self.df_performance
