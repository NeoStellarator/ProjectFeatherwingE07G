# -*- coding: utf-8 -*-
"""
Created on Sun May 21 19:21:02 2023

@author: aleja
""" 

import pandas as pd
import matplotlib.pyplot as plt

# Read in the Excel file
df = pd.read_csv('pfwtestdata2.csv')

# Extract the two columns you want to plot
disp = df['disp_calib']
load = list(df['load_calib'])
timetest = df['time']
#load[0:148] = [load[i]+(346-timetest[i]*320/72.5) for i in range(len(load[0:148]))]
df['load_calib'] = load
# Plot the data
plt.scatter(disp[0:145], load[0:145], s=0.8, c='r', marker='o', label='Datapoints before slip')
plt.scatter(disp[145:-1], load[145:-1], s=0.8, c='g', marker='o', label='Datapoints after  slip')
plt.xlabel('Displacement [mm]')
plt.ylabel('Load [N]')
plt.title('Force vs. Displacement')

# Showing in the figure the highest load value
maxload = max(load)
maxdisp = disp[load.index(maxload)]
plt.scatter(maxdisp, maxload, s=10, c='b', marker='o', label='Max Load')
# Writing the value of the highest load, with two decimals, and the value of its displacement, with a small offset

# Specify the x and y position of the label for the maximum load
label_x = maxdisp  # X position of the label
label_y = maxload  # Y position of the label

plt.grid('True')
plt.text(100, 780, f'{round(maxload, 2)} N', fontsize=14)

# Add an arrow from the text label to the point
arrow_properties = dict(facecolor='black', arrowstyle='->')
plt.annotate('', xy=(89, 847), xytext=(100, 790 + 10),
             arrowprops=arrow_properties)

plt.style.use("seaborn")

plt.legend(loc='lower right')
plt.show()
