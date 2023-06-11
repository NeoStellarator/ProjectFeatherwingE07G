# -*- coding: utf-8 -*-
"""
Created on Wed May  3 14:21:28 2023

@author: aleja
"""
import os
import numpy as np
import matplotlib.pyplot as plt

"""
#Can manually enter data or run scrip below
loads = [1.012, 2.626, 10.009, 13.647, 17.2752] # Example data, replace with your own
scaling_factors = [8685.77, 8717.38,  8785.493, 8793.14, 8798.04]
cell_values = [8790,22891.85,87934.209,120000,151988]
"""

folder_path = './scale'
loads = []
scaling_factors = []
cell_values = [] 
for filename in os.listdir(folder_path):
    filepath = os.path.join(folder_path, filename)
    if os.path.isfile(filepath):
        with open(filepath, 'r') as file:
            last_line = file.readlines()[-1]
            values = last_line.split(" ")
            loads.append(float(values[9]))
            scaling_factors.append(float(values[7]))
            cell_values.append(float(values[3]))


def get_scaling_factor(load, loads, scaling_factors):
    
    """
    Given a load in kg, returns the corresponding scaling factor based on the
    input data in `loads` and `scaling_factors`.
    """    
    # Use linear interpolation to estimate the scaling factor
    scaling_factor = np.interp(load, loads, scaling_factors)
    return scaling_factor


def plot_scaling_factor(loads, scaling_factors):
    
    z = np.polyfit(loads, scaling_factors, 3)
    p = np.poly1d(z)
    x = np.linspace(loads[0], loads[-1], 100)

    # Create the plot
    plt.plot(loads, scaling_factors, 'o', label='Data')
    plt.plot(x, p(x), label='Approximation')
    plt.xlabel('Load')
    plt.ylabel('Scaling factor')    
    plt.legend()
    plt.annotate(f'y = {p[0]:.2f}x^3 + {p[1]:.2f}x^2 + {p[2]:.2f}x + {p[3]:.2f}', 
             xy=(0.4, 0.55), xycoords='axes fraction', fontsize=10.5)
    plt.show()
    return p


def plot_cell_values(loads, cell_values):
    
    z = np.polyfit(loads, cell_values, 1)
    p = np.poly1d(z)
    x = np.linspace(loads[0], loads[-1], 100)
    
    loads_mean = np.mean(loads)
    cell_values_pred = np.polyval(p, loads)
    TSS = np.sum((loads_mean - cell_values_pred)**2)
    RSS = np.sum((cell_values - cell_values_pred)**2)
    r_sq = 1 - (RSS/TSS)
    # Create the plot
    plt.plot(loads, cell_values, 'o', label='Data')
    plt.plot(x, p(x), label='Approximation')
    plt.xlabel('Load [Kg]')
    plt.ylabel('Cell Load Value [-]')    
    plt.legend()
    plt.annotate(f'y = {p[1]:.2f}x + {p[0]:.2f}    RSQ {r_sq:.9f}', 
             xy=(0.4, 0.15), xycoords='axes fraction', fontsize=10.5, color = 'black')
    plt.show()
    return p, r_sq
