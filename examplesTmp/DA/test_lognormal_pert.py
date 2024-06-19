#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 10:05:13 2023

@author: ben
"""

import numpy as np
import matplotlib.pyplot as plt

# Define the mean and standard deviation of the lognormal distribution for sand and clay
mean_sand = 1.0  # Mean for sand conductivity
std_dev_sand = 0.1  # Standard deviation for sand conductivity

mean_clay = 0.5  # Mean for clay conductivity
std_dev_clay = 0.05  # Standard deviation for clay conductivity

# Number of realizations or perturbations
num_realizations = 100

# Generate random multiplicative perturbations for sand and clay
perturbations_sand = np.random.lognormal(mean=np.log(mean_sand), sigma=std_dev_sand, size=num_realizations)
perturbations_clay = np.random.lognormal(mean=np.log(mean_clay), sigma=std_dev_clay, size=num_realizations)

# Apply the perturbations to the original conductivities
original_sand_conductivity = 10.0  # Example value for sand conductivity
original_clay_conductivity = 5.0  # Example value for clay conductivity

perturbed_sand_conductivities = original_sand_conductivity * perturbations_sand
perturbed_clay_conductivities = original_clay_conductivity * perturbations_clay

# Print or use the perturbed conductivities as needed
print("Perturbed Sand Conductivities:")
print(perturbed_sand_conductivities)

print("Perturbed Clay Conductivities:")
print(perturbed_clay_conductivities)


# Create a histogram plot for perturbed sand conductivities
plt.figure(figsize=(10, 5))
plt.hist(perturbed_sand_conductivities, bins=20, alpha=0.5, color='blue', label='Sand Conductivity')
plt.hist(perturbed_clay_conductivities, bins=20, alpha=0.5, color='red', label='Clay Conductivity')
plt.xlabel('Perturbed Conductivity')
plt.ylabel('Frequency')
plt.title('Perturbed Saturated Hydraulic Conductivities')
plt.legend(loc='upper right')
plt.grid(True)

