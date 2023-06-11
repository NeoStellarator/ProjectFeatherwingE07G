# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 18:11:34 2023

@author: aleja
"""

import pandas as pd
import matplotlib.pyplot as plt

# Read in the CSV file
df = pd.read_csv('pfwtestdata2.csv')

# Extract the columns
disp = df['disp_calib']
load = list(df['load_calib'])
timetest = df['time']

# Create a figure with higher DPI
fig, axs = plt.subplots(1, 2, figsize=(10, 3.5), dpi=800)

# Plot load against time
axs[0].plot(timetest, load, 'r-', label='Load')
axs[0].set_xlabel('Time [s]')
axs[0].set_ylabel('Load [N]')
axs[0].set_title('Load vs. Time')
axs[0].grid(True)
axs[0].legend()

# Plot displacement against time
axs[1].plot(timetest, disp, 'b-', label='Displacement')
axs[1].set_xlabel('Time [s]')
axs[1].set_ylabel('Displacement [mm]')
axs[1].set_title('Displacement vs. Time')
axs[1].grid(True)
axs[1].legend()

# Adjust spacing between subplots
plt.tight_layout()

plt.style.use("seaborn")

# Show the figure
plt.show()