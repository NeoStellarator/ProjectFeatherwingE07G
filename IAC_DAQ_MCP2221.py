# Just in case the environment variables were not properly set

"""
Checklist
    1.Development Mode ON/OFF
    3.Calibration --> Remember to change MODE
        3.1.Zero Calibration --> Place no load on cell, and record returned value. 
                                 Then substract this value
        3.2.Scale Calibration --> Load Sensor with a known weight, and record the returned value. 
                                  Obtain scaling factor from this.
    4.Board must be connected, no working modules if else.
"""

import os
os.environ["BLINKA_MCP2221"] = "1"
os.environ["BLINKA_MCP2221_RESET_DELAY"] = "-1"

import board
import busio
import time
import pandas as pd

import adafruit_vl53l0x
from cedargrove_nau7802 import NAU7802

timestr = time.strftime("%Y%m%d-%H%M%S")

#Calibration Values
calibration_mode = "none" #"none", "Zero", "Scale"


# Load cell
loadCelSensor = NAU7802(board.I2C(), address=0x2a, active_channels=1)

# Time of flight sensor
i2c = busio.I2C(board.SCL, board.SDA)
tofSensor = adafruit_vl53l0x.VL53L0X(i2c)

print("Starting measurements. \n")

# Perform measurements
try:
    data = []
    if calibration_mode == "zero":
        naming_add = "fullsetup"
        num_samples = 0
        total_load = 0
        timetest = 0
        while True:
            # Get sensor readings
            load_cell = loadCelSensor.read()
            time_of_flight = tofSensor.range
    
        
            # Sleep
            time.sleep(0.5) #Change if wants 
            timetest = time.time()
            # Zero calibration
            num_samples += 1
            total_load += load_cell
            average_load = total_load / num_samples #Make the average of the measured load
            # Output sensor data
            print(f"Load cell: {load_cell}, Distance: {time_of_flight}, Iteration {num_samples}, Zero Val {average_load}")
            
            data.append({'load_cell': load_cell, 'Zero Val': average_load, 'Time': time})
            df = pd.DataFrame(data)
            df.to_excel("Zero_Val " + str(timestr) + ".xlsx", index=False)
            # Save calibration value
            with open("zerocalibration" + str(timestr) + naming_add +".txt", "a") as f:
                f.write(f"Zero_value {average_load}, Cell {load_cell}, TOF {time_of_flight}, time_since_start {timetest}\n")            
    
    
    elif calibration_mode == 'scale':
        known_weight = 21.237 #kg, replace for each
        num_samples = 0
        zero_value = 86604.28 #DETERMINED FROM BEFORE
        total_load = 0
        scale_load_tot = 0
        naming_add = "known_weight_kg"
        while True:
            # Get sensor readings
            load_cell = loadCelSensor.read()
            time_of_flight = tofSensor.range
                        
            # Sleep
            time.sleep(0.5)
            num_samples += 1
            
            scaled_load = (load_cell - zero_value)
            total_load += scaled_load
            average_load_calib = (total_load/ num_samples)
            scale_factor = average_load_calib / known_weight
            
            #print(scale_value)
            print(f"Load cell: {load_cell}, scaleloadtot {scale_load_tot}, numsample {num_samples}, Scale Factor {scale_factor}\n")
            data.append({'load_cell': load_cell, 'Av. Load': average_load_calib, 'Scale Factor': scale_factor, 'TOF': time_of_flight, 'Known Weight': known_weight})
            
            # Save calibration value
            with open("./scale/scalecalibration_date_" + str(timestr) + '_weight_' + str(known_weight) +".txt", "a") as f:
                f.write(f"Scaled_value {scaled_load} Average_Load_Calibrated {average_load_calib} Load_Cell {load_cell} Scale_Factor {scale_factor} Known_Weight {known_weight} \n")
            
            # Create a dataframe from the data list
            df = pd.DataFrame(data)
            # Save the dataframe to an excel file
            df.to_excel("TestData " + str(timestr) + ".xlsx", index=False)           
            
    else:
        while True:
            # Get sensor readings
            load_cell = loadCelSensor.read()
            time_of_flight = tofSensor.range
    
            # Output sensor data
            print("Load cell: {:.0f}, Distance: {:.0f}".format(load_cell, time_of_flight))
            
            # Sleep
            time.sleep(0.5)
            
            #Save data into text file
            with open("TestData" + str(timestr) +".txt", 'a') as f:
                f.write("Load cell: {:.0f}, Distance: {:.0f}".format(load_cell, time_of_flight) + '\n')
                
            # Append the values to the data list
            data.append({'load_cell': load_cell, 'time_of_flight': time_of_flight})
    
            # Create a dataframe from the data list
            df = pd.DataFrame(data)
        
            # Save the dataframe to an excel file
            df.to_excel("TestData " + str(timestr) + ".xlsx", index=False)


            #All data processing will be done afterwards, important to make sure
            #that raw data has been measured and stored adequately.
            #Apply zero-load for test set-up and scaling factor! (determine before testing)
            
        

# Exit
except KeyboardInterrupt:
    print("\nexiting...\n")
    #Nothing to add here!
