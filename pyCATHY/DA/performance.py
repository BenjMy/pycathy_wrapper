"""Class managing data Assimilation performance assessement
"""
import pandas as pd
import numpy as np
from pyCATHY.cathy_tools import CATHY

#%%

class perf(CATHY):
    
    def _performance_assessement(
        self, list_assimilated_obs, data, prediction, t_obs, **kwargs
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