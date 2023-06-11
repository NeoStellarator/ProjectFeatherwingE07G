# -*- coding: utf-8 -*-
"""
Created on Wed May  3 18:40:11 2023

@author: Alfon
"""

import math
import numpy as np, pandas as pd
from matplotlib import pyplot as plt

# -------- STATIC ANALYSIS --------

def whiffle_forces(F):
    '''
    Calculate the forces acting on the whiffletree by solving equilibrium equations.

    Parameters
    ----------
    F : FLOAT
        Applied load (in N) at point S in whiffletree.

    Returns
    -------
    LIST
        Forces acting on the whiffletree and their associated location based
        on the coordiante system specified in the Design Chapter. Format of elements
        in the list is [z, F], with z (in m) and F (in N). Sign of the force in
        accordance to the coordinate system (see FBD in Design Chapter).

    '''
    
    # Parameters that define whiffletree setup (see reader [1])
    a = 1.120 # [m]
    off = 30e-3 # [m]
    l1 = 0.830 # [m]
    l2 = 0.370 # [m]
    
    # a system of equations is built from sum of forces and sum of moments
    
    A = np.array([[1, 1],[l1, -l2]])
    b = np.array([[F],[0]])
    
    F1, F2 = np.linalg.solve(A,b)
    
    # the following distances define the positions of F1, F2 relative
    # to the coordinate system of the beam
    z1 = (a-off)-l1
    z2 = (a-off)+l2
    
    # sign of forces is changed in accordance to coordinate system of the wing 
    # spar beam. Note that only vertical forces are considered.
    return [[z1, -float(F1)], [z2, -float(F2)]]

def reactions_beam(forces):
    '''
    Determines the reaction forces on the beam from a given load case. The beam
    is clamped at the wing root, so the relevantreaction forces and moments are
    returned.

    Parameters
    ----------
    forces : LIST
        Forces applied on the beam and their location relative to the coordinate
        system at the root of the beam (see FBD in Design Chapter)

    Returns
    -------
    A : FLOAT 
        Vertical reaction force (in N) at the root. Positive in the positive direction
        of the coordinate sytem.
    M : TYPE
        Reaction moment (in Nm) at the root. Positive in the positive direction of the
        coordinate system
    '''
    
    A = 0
    M = 0
    
    for z, f in forces:
        A -= f    # negative sign, to find the reaction force in the positive direction
        M += z*f  # positive sign, to find the reaction force in the positive direction
    
    return A, M


def find_loading(F):
    '''
    For a given applied load at point S in the whiffletree, the loading on the 
    beam is returned: all forces and moments (and their locations) acting on it.
    It is different from whiffle_forces(F), as it includes the reaction forces.

    Parameters
    ----------
    F : FLOAT
        Applied load (in N) at point S in the whiffletree.

    Returns
    -------
    forces : FLOAT
        All forces acting on the whiffletree and their associated location based
        on the coordiante system specified in the Design Chapter. Format of elements
        in the list is [z, F], with z (in m) and F (in N). Sign of the force in
        accordance to the coordinate system (see FBD in Design Chapter).
    moments : FLOAT
        All moments acting on the whiffletree and their associated location based
        on the coordiante system specified in the Design Chapter. Format of elements
        in the list is [z, M], with z (in m) and M (in Nm). Sign of the moement in
        accordance to the coordinate system (see FBD in Design Chapter).

    '''
    # returns a list of forces and moments acting on the beam for a given applied load F
    forces = whiffle_forces(F)
    F0, M0 = reactions_beam(forces)

    forces.insert(0, [0, F0])
    moments = [[0, M0]]
    
    return forces, moments

def IFD_plotter(forces, moments, beam_length, showplot=False, saveplot=False, filename="IFD_plot.png"):
    '''
    Internal Force (and Moment) Diagram for loaded cantilevered beam.

    Parameters
    ----------
    forces : LIST
        A two dimensional array with each row assosiated to a point vertical force. 
        The first column defines the location of the force (as a function of z,
        with z=0 at the rigid support, in metres) and the second column defines
        the magnitude of the force in Newtons. A force is positive in accordance 
        to the established coordinate system (see Design Chapter).
        
    moments : LIST
        A two dimensional array with each row assosiated to a point moment. 
        The first column defines the location of the moment (as a function of z,
        with z=0 at the rigid support, in metres) and the second column defines
        the magnitude of the moment in Newtons. A moment is positive in accordance 
        to the established coordinate system (see Design Chapter).
        
    beam_length : FLOAT
        The length of the beam is specified (in metres).
        
    showplot : BOOL, optional
        True if the internal force (and moment) diagrams want to be shown.
        The default is False.
        
    saveplot : BOOL, optional
        True if the internal force (and moment) diagrams want to be saved.
        The default is False.
        
    filename : STRING, optional
        If the saveplot option is selected, the filename should be introduced.
        The default is "IFD_plot".

    Returns
    -------
    z_tab : LIST
        Datapoints in the z-direction (in metres)
    V_tab : LIST
        Internal shear forces (in Newtons) corresponding to z-values in z_tab. 
        Signs in accordance with sign convention (see Design Chapter).
    M_tab : LIST
        Internal bending moments (in Newtons) corresponding to z-values in z_tab.
        Signs in accordance with sign convention (see Design Chapter).

    '''
    
    z_tab = [] # [m]
    V_tab = [] # [N]
    M_tab = [] # [Nm]
    
    z = 0 # [m]
    
    dz = 0.0001 # [m]
    
    while z <= beam_length:
        V = 0
        M = 0
          
        for f in forces:
            if f[0] <= z:
                V -= f[1]
                M -= (z-f[0])*f[1]
        
        for m in moments:
            if m[0] <= z:
                M -= m[1]
                
        z_tab.append(z)
        V_tab.append(V)
        M_tab.append(M)
        
        z += dz
    
    if showplot or saveplot:
        # Plotting
        
        plt.figure(figsize=(12, 3))
        
        plt.subplot(1,2,1)
        plt.style.use("seaborn-v0_8")
        plt.title("Internal Shear Force Diagram")
        plt.xlabel(r'$z$'+"  [ m ]")
        plt.ylabel(r'$V$'+"  [ N ]")
        # plt.hlines(0, z_tab[0], z_tab[-1], colors="black", lw=0.5)
        plt.plot(z_tab, V_tab, color="g")
        plt.vlines(0, 0, V_tab[0], color="g")
        
        plt.subplot(1,2,2)
        plt.style.use("seaborn-v0_8")
        plt.title("Internal Bending Moment Diagram")
        plt.xlabel(r'$z$'+"  [ m ]")
        plt.ylabel(r'$M$'+"  [ Nm ]")
        # plt.hlines(0, z_tab[0], z_tab[-1], colors="black", lw=0.5)
        plt.plot(z_tab, M_tab, color="darkblue")
        plt.vlines(0, 0, M_tab[0], color="darkblue")
        
        plt.tight_layout()
        
        if saveplot:
            plt.savefig(filename, dpi=500)
        if showplot:
             plt.show()  
        
    return z_tab, V_tab, M_tab

# -------- STRESS DISTRIBUTIONS --------

def bending_dist(M, I, consistent=True, iteration=2):
    
    '''
    Generate the bending stress distribution for the selected shape (see Design
    Chapter), between the outer surfaces of the stringers
    
    Parameters
    ----------
    M : FLOAT
        Internal bending moment at a given point (in Newton metre). Sign in
        accordance to sign convention (see Design Chapter).
    I : FLOAT
        Area moment of inertia (in m^4)
    consistent : BOOL, optional
        If the bending stress and shear stress distribution are to be compared, 
        to ease calculations the y_max is altered to match that taken for the 
        shear stress distribution.  
    iteration : INT, optional
        Specify the iteration (1 or 2) of the design being analyzed.
    
    Returns
    -------
    y : ARRAY
        Datapoints in the y-direction (in metres)
    sigma_tab : ARRAY
        Normal stresses (in Pascals) corresponding to y-values in y_tab
    '''
    
    if not consistent:
        y_max = 75e-3 # [m] dist from N.A. to the outermost surface
    elif iteration == 1:
        y_max = 72.7e-3 # [m]  dist from N.A. to inner surface of stringer (to be consistent with shear distribution, see shear_dist)
    elif iteration == 2:
        y_max = 73.5e-3 # [m] dist from N.A. to inner surface of stringer (to be consistent with shear distribution, see shear_dist)
       
    
    dy = 0.0001  # [m] resolution of 0.1mm
    y = np.arange(-y_max, y_max+dy, dy)
    
    sigma_tab = (M/I)*y
    
    return y, sigma_tab


def shear_dist(V, I, iteration=2):
    '''
    Generate the shear stress distribution for the selected shape, between the 
    interior surfaces of the stringers
    
    Parameters
    ----------
    V : FLOAT
        Internal shear force at a given point (in Newtons). Sign in accordance 
        to sign convention (see Design Chapter).
    I : FLOAT
        Area moment of inertia (in m^4)
    iteration : INT, optional
        Specify the iteration (1 or 2) of the design being analyzed.
    
    Returns
    -------
    y : ARRAY
        Datapoints in the y-direction (in metres)
    tau_tab : ARRAY
        Shear stresses (in Pascals) corresponding to y-values in y_tab
    '''
    
    dy = 0.0001
    
    if iteration == 1:
        y_max = 72.7e-3 # [m] - only study the shear distribution from the inner surface of the stringers
        
        y = np.arange(y_max, -(y_max+dy), -dy)
        size_y = np.shape(y)[0]
        
        # Define thickness as functions of y
        t = np.zeros((size_y))
        t = np.where(np.abs(y)>=54.2e-3, 3.8e-3, 0.8e-3) # specify the point where the thickness changes (y>54.2mm)
        
        # Initialize Q
        Q = np.zeros((size_y)) # [m3]
        Q_val = np.array(6682.34e-9) # value calculated manually for the segment above the inner surface of stringer 
    elif iteration == 2:
        y_max = 73.5e-3 # [m] - only study the shear distribution from the inner surface of the stringers
        
        y = np.arange(y_max, -(y_max+dy), -dy)
        size_y = np.shape(y)[0]
        
        # Define thickness as functions of y
        t = np.zeros((size_y))
        t = np.where(np.abs(y)>=55e-3, 3.8e-3, 0.8e-3) # specify the point where the thickness changes (y>55mm)
        
        # Initialize Q
        Q = np.zeros((size_y)) # [m3]
        Q_val = np.array(4496.356e-9) # value calculated manually for the segment above the inner surface of stringer  

    
    dy = np.array(dy)
    
    for idx, y_val in np.ndenumerate(y):
        Q_val += t[idx]*dy*(y_val-dy/2)
        Q[idx] = Q_val
    
    tau_tab = (V*Q)/(I*t)
    
    y = np.flip(y)
    
    return y, tau_tab



# -------- PLANE STRESS TRANSFORMATIONS --------

def plane_transform(sig, tau):
    '''
    Determine principal stresses, orientation of the principal planes, maximum
    shear stresses and maximum shear stress plane from a given state of stress.
    Only normal (bending) stresses in one direction are needed.

    Parameters
    ----------
    sig : FLOAT
        Normal stress (in Pa) at the point (in z-direction). Positive for tensile case.
    tau : FLOAT
        Shear stress (in Pa) at the point. Positive if it causes the state of stress
        to rotate in the counterclockwise direction.

    Returns
    -------
    sig_max : FLOAT
        Maximum nomal stress (in Pa) of the state of stress. (Principal planes orientation)
    sig_min : FLOAT
        Minimum nomal stress (in Pa) of the state of stress. (Principal planes orientation)
    theta : FLOAT
        Counterclockwise angle (in rad) of the orientation of the principal planes.
    r : FLOAT
        Maximum shear stress (in Pa) of the state of stress. (Maximum shear stress plane)
    phi : FLOAT
        Counterclockwise angle (in rad) of the orientation of the maximum shear stress plane.

    '''
    
    sig_ave = sig / 2
    r = math.sqrt(sig_ave ** 2 + tau ** 2)
    if abs(sig_ave + r) > abs(sig_ave - r):
        sig_max = sig_ave + r
        sig_min = sig_ave - r  # sig_min isn't actually the lowest normal stress
    else:
        sig_max = sig_ave - r
        sig_min = sig_ave + r  # sig_min isn't actually the lowest normal stress
    try:
        if sig > 0:
            theta = -math.atan(tau / abs(sig - sig_ave)) / 2
            phi = (math.pi / 2 - math.atan(tau / abs(sig - sig_ave))) / 2
        else:
            theta = math.atan(tau / abs(sig - sig_ave)) / 2
            phi = -(math.pi / 2 - math.atan(tau / abs(sig - sig_ave))) / 2
    except ZeroDivisionError:
        if sig > 0:
            theta = -math.pi/2 / 2
            phi = (math.pi / 2 - math.pi/2) / 2
        else:
            theta = math.pi/2 / 2
            phi = -(math.pi / 2 - math.pi/2) / 2

    # print("τ max: ", r)
    # print("σ max: ", sig_max)
    # print("σ in other direction:", sig_min)
    # print("turn element", theta, "rad anti-clockwise for σ max (=", (theta * 180 / math.pi), "°)")
    # print("turn element", phi, "rad anti-clockwise for τ max (=", (phi * 180 / math.pi), "°)")
    
    return sig_max, sig_min, theta, r, phi


# -------- MISCELLANEOUS --------

def savedata(header, dat, filename="Data.xlsx"):
    '''
    Save data of a two dimensional array (as nested lists) into an excel file.

    Parameters
    ----------
    header : LIST
        Names of columns specified in list of strings.
    dat : LIST
        Data to be saved (two dimensinal array / nested lists)
    filename : STRING, optional
        Name of excel file to save the data. The default is "Data.xlsx".

    Returns
    -------
    None.

    '''
    
    dat_array = np.array(dat)
    dat_array = dat_array.transpose()
    df = pd.DataFrame(dat_array, columns=header)
    df.to_excel(filename)
    
    return None
