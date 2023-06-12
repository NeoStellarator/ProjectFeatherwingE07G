# ProjectFeatherwingE07G
Repository containing all the python code that was developed for the TU Delft AE Bsc Year 1 Q3 project, from design to testing.

During this engineering project, the analysis of structural integrity played a critical role in ensuring reliability and safety. To aid in the design and testing of structures, a set of Python codes has been developed. These codes provide functionalities for stress analysis, data visualization, and measurement processing. 

The following codes have played a crucial role in the design process. The first one contains the main functions, and the second one utilizes the functions to perform the stress analysis relevant for the spar design.


**mechanicstools.py**
 - Performs static analysis on the wing spar with the whiffletree setup.
 - Generates internal force and bending moment diagrams for a given loading on the spar.
 - Generates bending and shear stress distributions for a given internal force and moment on a spar cross-section.
 - Performs plane stress transformations on the state of stress with a normal stress component due to a bending distribution and a shear stress component due to a transverse shear distribution, to determine the principal stresses, absolute maximum shear stress and the orientation of the state of stress in either case.


**Featherwing Stress Analysis.py**
- Performs stress analysis specifically on the root of the spar, performing stress transformations at each point to find the maximum shear stress distribution.
- Generate plots for the stress distributions

These codes have been extensively used during the design phase to ensure the structural integrity and performance of the analyzed components.

The following codes have been utilized for testing purposes:

**scalingfactorw2.py**
 - Calculates the scaling factor based on input data of loads and scaling factors.
 - Provides a function to obtain the scaling factor for a given load.
 - Includes a plotting function to visualize the relationship between load and scaling factor.

**graphingpfw.py**
 - Generates a scatter plot of force vs. displacement.
 - Reads data from a CSV file.
 - Plots the datapoints.
 - Identifies and marks the maximum load value on the plot.

**IAC\_DAQ\_MCP2221.py**
 - Performs measurements using load cell and time-of-flight sensors.
 - Includes functionalities for calibration, zero calibration, and scale calibration.
 - Saves the measured data into text files and Excel files for further analysis.

**timegraphs2.py**
 - Reads data from a CSV file containing load and displacement measurements over time.
 - Creates a figure with two subplots: load vs. time and displacement vs. time.

These testing codes have played a crucial role in evaluating the performance of the analyzed components and collecting valuable data for analysis and improvement.

The combination of design and testing codes has significantly facilitated the development and assessment of structures, ensuring their reliability and safety.
