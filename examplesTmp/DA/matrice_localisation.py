#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 12:39:12 2024

@author: z0272571a
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Initial covariance matrix (5x5) with 1e-2 on the diagonal
# cov_matrix = np.eye(36) * 1e-2

#%%

import numpy as np

# Example time series data (time steps x spatial locations)
# Suppose we have 5 time steps and 3 spatial locations
measurements = np.array([
    [1.1, 2.0, 3.2],
    [0.9, 2.1, 3.0],
    [1.0, 1.9, 3.1],
    [1.2, 2.2, 3.3],
    [1.0, 2.0, 3.1]
])

# If ground truth values are known, you can use them. Here we use the mean as an estimate
ground_truth = np.mean(measurements, axis=0)

# Calculate the errors (subtract ground truth from each measurement)
errors = measurements - ground_truth

# Each row of 'errors' is an error vector for a time step
# We need to stack these error vectors to compute the covariance matrix

# Compute the covariance matrix of the errors
error_covariance_matrix = np.cov(errors, rowvar=False)

print("Measurement Error Covariance Matrix:")
print(error_covariance_matrix)

fig, axs = plt.subplots()

im1 = axs.imshow(error_covariance_matrix, 
                  cmap='viridis', 
                  interpolation='none',
                  )
plt.colorbar(im1, ax=axs)


#%%
# Observation and its covariance matrice 

import numpy as np

# Generate a 10x10 matrix filled with random points from 0 to 1
random_matrix = np.random.rand(10, 10)

# Scale to the range [0.94, 1]
scaled_random_matrix = 0.1 + random_matrix * 0.4

print("10x10 Matrix filled with random points from 0.94 to 1:")
print(scaled_random_matrix)
fig, axs = plt.subplots(2)
im1 = axs[0].imshow(scaled_random_matrix, 
                  cmap='viridis', 
                  interpolation='none',
                  )
plt.colorbar(im1, ax=axs[0])

im1 = axs[1].imshow(np.cov(scaled_random_matrix), 
                  cmap='viridis', 
                  interpolation='none',
                  )
plt.colorbar(im1, ax=axs[1])



#%%

# Function to create a Gaussian covariance matrix
def gaussian_covariance_matrix(size, sigma):
    """Create a Gaussian covariance matrix."""
    cov_matrix = np.zeros((size, size))
    for i in range(size):
        for j in range(size):
            cov_matrix[i, j] = np.exp(-((i - j)**2) / (2 * sigma**2))
    return cov_matrix

# Define size and sigma for the Gaussian covariance matrix
size = 16
sigma = 3.0e-9

# Create the Gaussian covariance matrix
cov_matrix = gaussian_covariance_matrix(size, sigma)


# Define zones
# zone1 = [0, 1]
# zone2 = [2, 3, 4]

zone1 = np.arange(0,int(len(cov_matrix)/2))
zone2 = np.arange(int(len(cov_matrix)/2),len(cov_matrix))

# Localization function (Gaussian decay)
def localization_function(distance, length_scale):
    return np.exp(-(distance**2) / (2 * length_scale**2))

# localization_function(2, 1)

# Length scale for localization
length_scale = 5

# Apply localization to the covariance matrix
localized_cov_matrix = np.zeros_like(cov_matrix)
localized_matrix = np.zeros_like(cov_matrix)

# Populate the localized covariance matrix
for i in range(cov_matrix.shape[0]):
    for j in range(cov_matrix.shape[1]):
        if (i in zone1 and j in zone1) or (i in zone2 and j in zone2):
            distance = abs(i - j)
            print(cov_matrix[i, j],localization_function(distance, length_scale))
            # print(localization_function(distance, length_scale))
            localized_cov_matrix[i, j] = cov_matrix[i, j] * localization_function(distance, length_scale)
            localized_matrix[i, j] = localization_function(distance, length_scale)
            # print(cov_matrix[i, j]-localized_cov_matrix[i, j])

fig, axs = plt.subplots(1, 3, figsize=(12, 6))

# Original covariance matrix
im1 = axs[0].imshow(cov_matrix, cmap='viridis', interpolation='none')
axs[0].set_title('Original Covariance Matrix')
fig.colorbar(im1, ax=axs[0])

# Localized covariance matrix
im2 = axs[1].imshow(localized_cov_matrix, cmap='viridis', interpolation='none')
axs[1].set_title('Localization Matrix')
fig.colorbar(im2, ax=axs[1])


# Localized covariance matrix
im2 = axs[2].imshow(localized_cov_matrix, cmap='viridis', interpolation='none')
axs[2].set_title('Localized Covariance Matrix')
fig.colorbar(im2, ax=axs[2])

# Add zones to the plots
for ax in axs:
    # Zone 1
    rect1 = patches.Rectangle((zone1[0] - 0.5, zone1[0] - 0.5), len(zone1), len(zone1), linewidth=2, edgecolor='red', facecolor='none')
    ax.add_patch(rect1)
    # Zone 2
    rect2 = patches.Rectangle((zone2[0] - 0.5, zone2[0] - 0.5), len(zone2), len(zone2), linewidth=2, edgecolor='blue', facecolor='none')
    ax.add_patch(rect2)

plt.show()