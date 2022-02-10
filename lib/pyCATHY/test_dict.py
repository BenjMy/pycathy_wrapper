#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 08:39:17 2022

@author: ben
"""

dict_obs_tmp = {}
dict_obs_2add = {'filename': 'gdg2',
                'data_type': 'gdg',
                'units': 'fdgg', # units
                }


self_dict_obs= {}

for l in ['0','0']:
    dict_obs_tmp[l]= {}
    for m in ['ERT','Tensio']:
        dict_obs_tmp[l][m] = {}
        for key in dict_obs_2add.keys():
            dict_obs_tmp[l][m][key] = dict_obs_2add[key]
            

self_dict_obs = self_dict_obs | dict_obs_tmp
