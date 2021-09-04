#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  3 13:40:34 2021

@author: ben
"""

import os 
import numpy as np
import matplotlib.pyplot as plt

length=1;
width=1;

# Read the Header
# str_hd_dem = {'north':0,'south':0,'east':0,'west':0,'rows':0,'cols':0}
str_hd_dem = {}
with open('/home/ben/Documents/CATHY/pyCATHY/synthetic/weiltopo/prepro/dem') as f: # open the file for reading
    count = 0
    for line in f: # iterate over each line
        if count < 6:
            str_hd, value_hd = line.split() # split it by whitespace
            str_hd_dem[str_hd.replace(':', '')]=value_hd
        count += 1
            

dem_file = open(os.path.join('/home/ben/Documents/CATHY/pyCATHY/synthetic/weiltopo/prepro/dem'), 'r')
dem_mat = np.loadtxt(dem_file,skiprows=6)
dem_file.close()

x=np.zeros(int(str_hd_dem['rows']))
y=np.zeros(int(str_hd_dem['cols']))

for a in range(int(str_hd_dem['rows'])):
    x[a]=float(str_hd_dem['west'])+length*a;
    
for a in range(int(str_hd_dem['cols'])):
    y[a]=float(str_hd_dem['south'])+width*a;
    
# x=x-width/2
# y=y-length/2


fig = plt.figure()
ax = plt.axes(projection='3d')
# Make data.
X, Y = np.meshgrid(x, y)
# Plot the surface.
surf = ax.plot_surface(X, Y, dem_mat.T, cmap='viridis')
# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
ax.set(xlabel='Easting (m)', ylabel='Northing (m)', zlabel='Elevation (m)')

plt.show()

