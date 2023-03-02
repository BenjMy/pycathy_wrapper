#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 16:55:33 2023

@author: ben
"""

import numpy as np
from scipy.spatial import KDTree

def get_top_bottom_surface_nodes(x, y, z, topography, eps=1e-6):
    x_min, x_max = np.min(x), np.max(x)
    y_min, y_max = np.min(y), np.max(y)
    z_min, z_max = np.min(z), np.max(z)
    
    # outer_mask = np.logical_or.reduce((x == x_min, x == x_max, y == y_min, y == y_max))
    # nodes = np.column_stack((x, y, z))
    nodes = np.array((x, y, z))

    # Find the top and bottom surface nodes using the topography
    top_mask = z == z_max - topography
    bottom_mask = z == z_min
    
    np.shape(top_mask)
    np.shape(nodes[0,:,:,:])
    
    top_nodes = nodes[0,:,:,:][top_mask]
    bottom_nodes = nodes[0,:,:,:][bottom_mask]
    
    # Build a KDTree to find the nearest neighbors of each node
    tree = KDTree(np.column_stack((x, y, z)))
    top_node_indices = np.arange(len(x))[top_mask]
    bottom_node_indices = np.arange(len(x))[bottom_mask]
    
    top_surface_nodes = []
    for i in top_node_indices:
        distances, indices = tree.query(top_nodes[i], k=4)
        outer_neighbors = np.isin(indices, top_node_indices)
        if np.count_nonzero(outer_neighbors) <= 3:
            top_surface_nodes.append(top_nodes[i])
    
    bottom_surface_nodes = []
    for i in bottom_node_indices:
        distances, indices = tree.query(bottom_nodes[i], k=4)
        outer_neighbors = np.isin(indices, bottom_node_indices)
        if np.count_nonzero(outer_neighbors) <= 3:
            bottom_surface_nodes.append(bottom_nodes[i])
    
    return np.array(top_surface_nodes), np.array(bottom_surface_nodes)

#%%
# Define the node coordinates
x = np.linspace(0, 1, 10)
y = np.linspace(0, 1, 10)
z = np.linspace(0, 1, 5)
x, y, z = np.meshgrid(x, y, z, indexing='ij')

def gaussian_filter(kernel_size, sigma=0.6, muu=0):
    x, y = np.meshgrid(np.linspace(-1, 1, kernel_size),
                       np.linspace(-1, 1, kernel_size))
    dst = np.sqrt(x**2+y**2)
    gauss = np.exp(-((dst-muu)**2 / (2.0 * sigma**2))) 
    return gauss
 
 
kernel_size=10
topography = gaussian_filter(kernel_size)

ztopo = np.zeros(np.shape(x))
for i in range(5):
    ztopo[:,:,i]=(z[:,:,i]+topography)
    
fig, ax = plt.subplots(1)
cmap = ax.imshow(ztopo[0])
plt.colorbar(cmap)

%matplotlib auto
import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x,y,ztopo)

# Find the outer top and bottom surface nodes of the
get_outer_top_bottom_surface_nodes(x, y, z, topography, eps=1e-6)