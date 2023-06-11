# -*- coding: utf-8 -*-
"""
Created on Wed May  3 19:57:48 2023

@author: Alfon
"""

import mechanicstools as mech
from matplotlib import pyplot as plt
import numpy as np

# Main Parameters
design = 2
F = 1500 # [N] applied load at point S of whiffletree
root = ['plots\it', str(design), str(F)+'N_']
root = "_".join(root)

# Beam Properties
beam_length = 2250e-3 # [m]

if design == 1:
    I = 1.6715019702667e-6 # [m4] (see App A for section calculations)
elif design == 2:
    I = 1.3408472769333e-6 # [m4] (see App A for section calculations)
    
# Static Analysis

forces, moments = mech.find_loading(F)

# Plot internal force and bending moment diagrams
z_tab, V_tab, M_tab = mech.IFD_plotter(forces, moments, beam_length, showplot=True, saveplot=False, filename=root+'00_IFD.png')

# Realizing that the largest internal forces and moments are at the root
M0 = M_tab[0]
V0 = V_tab[0]

# Plot bending and shear stresses
y1, sigma_tab = mech.bending_dist(M0, I, iteration=design)
y2, tau_tab = mech.shear_dist(V0, I, iteration=design)

# Obtaining the maximum shear stress at all points
size_y = np.size(y1)

abs_shear_stress = np.zeros((size_y))

for idx, y_val in np.ndenumerate(y1):
    abs_shear_stress[idx] = mech.plane_transform(sigma_tab[idx], tau_tab[idx])[3]

# max shear
max_stress = np.max(abs_shear_stress)
print(max_stress)


# ----------------------- Plotting and Saving Graphs -----------------------

# ---- 00 Change Stress Units to MPa ----
sigma_tab /= 1e6
tau_tab /= 1e6
abs_shear_stress /= 1e6

# ---- 01 Normal and Shear Stress Distribution ----
plt.style.use("seaborn-v0_8")
plt.vlines(0, sigma_tab[0],sigma_tab[-1], linestyle='--', lw=0.8)
plt.hlines(0, y1[0],y1[-1], linestyle='--', lw=0.8)
plt.plot(y1, sigma_tab,color='blue', label='Normal stress '+r'$\sigma$'+' [MPa]')
plt.plot(y1, tau_tab,color='green', label='Shear stress '+r'$\tau$'+' [MPa]')
plt.legend()
plt.xlabel(r'$y$' + " [m]")
plt.ylabel("Stress [MPa]")

# plt.savefig(root+"01_Normal_shear.png", dpi=500)
plt.show()

# ---- 02 Normal, Shear and Maximum Shear Stress Distribution ----

plt.style.use("seaborn-v0_8")
plt.vlines(0, sigma_tab[0],55, linestyle='--', lw=0.8)
plt.hlines(0, y1[0],y1[-1], linestyle='--', lw=0.8)
plt.plot(y1, sigma_tab, color='blue', alpha=0.2, label='Normal stress '+r'$\sigma$'+' [MPa]')
plt.plot(y1, tau_tab, color='green', alpha=0.2, label='Shear stress '+r'$\tau$'+' [MPa]')
plt.plot(y1, abs_shear_stress,color='red', label='Maximum shear stress (magnitude) '+r'$\tau_{max}$'+' [MPa]')
plt.legend(facecolor='white', edgecolor='grey')
plt.xlabel(r'$y$' + " [m]")
plt.ylabel("Stress [MPa]")

# plt.savefig(root+"02_Normal_shear_maxshear.png", dpi=500)
plt.show()

# ---- 03 Maximum Shear Stress Distribution Only ----
# plt.title('Maximum shear stress (magnitude)')
plt.style.use("seaborn-v0_8")
plt.vlines(0, 0,abs_shear_stress[0], linestyle='--', lw=0.8)
plt.hlines(0, y1[0],y1[-1], linestyle='--', lw=0.8)
plt.plot(y1, abs_shear_stress,color='red')
plt.xlabel(r'$y$' + " [m]")
plt.ylabel(r'$\tau_{max}$'+' [MPa]')
# plt.savefig(root+"03_max_shear_only.png", dpi=500)
plt.show()

# ---- 04 Normal Stress Distribution Only ----

y3, sigma_tab2 = mech.bending_dist(M0, I,consistent=False)
sigma_tab2 /= 1e6 # change units to MPa
# plt.title("Normal Stress Distribution")
plt.style.use("seaborn-v0_8")
plt.xlabel(r'$y$' + " [m]")
plt.ylabel(r'$\sigma$' + " [MPa]")
plt.vlines(0, np.min(sigma_tab2), np.max(sigma_tab2), linestyle='--', lw=0.8)
plt.hlines(0, y3[0],y3[-1], linestyle='--', lw=0.8)
plt.plot(y3, sigma_tab2, color='blue')
   
# plt.savefig(root+"04_normal_only.png", dpi=500)
plt.show()

# ---- 05 Shear Stress Distribution Only ----

# plt.title("Shear Stress Distribution")
plt.style.use("seaborn-v0_8")
plt.xlabel(r'$y$' + " [m]")
plt.ylabel(r'$\tau$' + " [MPa]")
plt.vlines(0, np.min(tau_tab), 0, linestyle='--', lw=0.8)
plt.hlines(0, y2[0],y2[-1], linestyle='--', lw=0.8)
plt.plot(y2, tau_tab, color='green')
   
# plt.savefig(root+"05_shear_only.png", dpi=500)
plt.show()
