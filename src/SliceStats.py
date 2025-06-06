# -*- coding: utf-8 -*-
import sys
import os
import random
import math
import time
import csv
from tkinter import Tk
from tkinter.filedialog import askopenfilename
import numpy as np
import pandas as pd
from skimage import measure 

############################################################################################################
#   Eastern Michigan University
#   Backues Lab  
#   Author: Payton Dunning and Steven Backues
#   Last Date Modified: Dec. 18th 2024
#
#   A script for analyzing the contents of an Autophagic Vacuole Simulation (AVS) project formatted 
#   Compucell 3D (CC3D) simulation. The script takes in a PIF file (.piff), that must contain "Body" and
#   "Wall" cells, and takes a slice through the simulation comparable to a TEM image of a cell. The bodies
#   within this slice are then analyzed to determine their relative areas. This area data is then recorded
#   for later compilation and statistical analysis.
############################################################################################################
''' Todo Dec. 2024:
1. Fix lines 68-82: right now they are setting the slice size based on the body cluster, not the wall. Make them read the vacuole size from the combined csv made by vacuole_gen
2. Get values of mu and sigma for body size and number from Vacuole_gen csv and add those to the output csv (line 399))
3. Make it automatically read the CC3D PIFF-dumped file when mass-runs = true
4. Fix error handling for vacuole slice limit (line 135)
5. Add option to take serial slices?'''

# paramsFile is used to keep track of several variables used by multiple scripts.
paramsFile = './attributes/Model_Parameters.txt'   # For the linux server
# paramsFile = 'src/attributes/Model_Parameters.txt'   # For Windows

def main(fileSelectOpt, MassRunCheck, inputPiff):
    initialTime = time.asctime(time.localtime(time.time()))

    if MassRunCheck:
        # Directly use the expected PIFF file from CC3D
        inputName = "./Output/10_24Simulation000.piff"
        print(f"Running SliceStats.py with: {inputName}")
        # Ensure the file exists before proceeding
        while not os.path.exists(inputName):
            time.sleep(2)  # Wait 2 seconds before checking again
    else:
        if fileSelectOpt:
            print("Please Select piff file...")
            print("(The file selection screen may appear BEHIND your current application)")
            Tk().withdraw()
            filename = askopenfilename()
            inputName = filename
            print("Given File Name: %s" % filename)
        else:
            print(">>Enter INPUT PIF file path+name:")
            inputName = input()

    # Scan file to get x-coordinate range from Body entries
    with open(inputName, "r") as f:
        min_x = float('inf')
        max_x = float('-inf')
        for line in f:
            data = line.split()
            if len(data) >= 4 and data[1] == "Body":
                x1 = float(data[2])
                x2 = float(data[3])
                min_x = min(min_x, x1, x2)
                max_x = max(max_x, x1, x2)

    print("Grabbing AVS Model Parameters...\n")
    modelParams = load_parameters_from_file(paramsFile)
    print(modelParams)
    scaleFactor = int(modelParams['Scale_Factor'])  # Keep scale factor from params file
    unScaledSliceThickness = int(modelParams['unScaledSliceThickness'])
    unScaledVacMin = int(float(modelParams['unScaledVacMin']))
    unScaledminBodyRadius = int(modelParams['unScaledminBodyRadius'])
    # Scale vacuole coordinates
    centerX = float(modelParams.get("Vacuole_x", 0)) / scaleFactor
    centerY = float(modelParams.get("Vacuole_y", 0)) / scaleFactor
    centerZ = float(modelParams.get("Vacuole_z", 0)) / scaleFactor
    
    # Convert to integers (for array indexing)
    centerX = int(round(centerX))
    centerY = int(round(centerY))
    centerZ = int(round(centerZ))
    wallRadius = int(float(modelParams["Wall_Radius_mu"]))
 
    print(f"Using Vacuole Center as Slice Reference: X={centerX}, Y={centerY}, Z={centerZ}")

    print("Current Model Parameters:\n")
    print("\tScale_Factor: %d\n" % scaleFactor)
    print("\tunScaledSliceThickness: %d\n" % unScaledSliceThickness)
    print("\tunScaledVacMin: %d\n" % unScaledVacMin)
    print("\tunScaledminBodyRadius: %d\n" % unScaledminBodyRadius)
    
    if not MassRunCheck:
        print(">>Would you like to use these parameters?[y/n]")
        paramSelect = input()
        if paramSelect == "n":
            print(">>Please enter new values for parameters:\n")
            print("(The Wall radius parameter value should be a post-scaling value)")
            print("\n>>Enter new scaling factor: ")
            scaleFactor = int(input())
            print("\n>>Enter the given wall's radius", end='')
            wallRadius = int(input())
            print("\n>>Enter the given wall's central x-coordinate:", end='')
            centerX = int(input())
      
    vacMin = (unScaledVacMin / scaleFactor)
    print("Default slice recognition limit (radius) = %d units" % vacMin)
    
    if not MassRunCheck:
        print(">>Would you like to use this default minimum vacuole slice threshold?[y/n]")
        minDInput = input()
        if minDInput.lower() == "n":
            print("\n>>Enter new minimum vacuole threshold (scaled): ")
            vacMin = int(input())
    
    # Adjust recognition limit to match file scale
    wallRecDiff = (wallRadius**2) - (vacMin**2)
    diamRangeVar = math.sqrt(wallRecDiff) if wallRecDiff > 0 else wallRadius
    if wallRecDiff <= 0:
        print("\n!!!Using full coordinate range due to small wall radius")
        
    # Adjust minX and maxX using vacuole center and wall radius
    minX = max(int(centerX - diamRangeVar), int(min_x))
    maxX = min(int(centerX + diamRangeVar), int(max_x))
    
    print(f"Adjusted valid slice range: {minX} to {maxX}")
    print(f"Computed centerX: {centerX}, min_x: {min_x}, max_x: {max_x}")

    sliceCoord = -1
    if MassRunCheck:
        sliceCoord = random.randint(minX, maxX)
    else:
        print(">>Finally, select an option for determining where a slice will be taken:")
        print("\t[0 for slice to be taken at centerX coordinate]")
        print("\t[1 for slice to be taken at a randomly selected coordinate]")
        print("\t[2 for slice to be taken at a user specified coordinate]")
        sliceChoice = int(input())
        if sliceChoice == 0:
            sliceCoord = centerX
        elif sliceChoice == 1:
            sliceCoord = random.randint(minX, maxX)
        elif sliceChoice == 2:
            print("\n>>Enter the x coordinate you'd like the slice to be taken at:")
            sliceCoord = int(input())
        else:
            print("\nInput was found to be invalid. Please enter 0, 1, or 2 for your slice selection method.")
            return
            
    print(f"Taking slice at X coordinate: {sliceCoord}")
    
    # Adjust recognition limit based on body size
    minBodyRadius = (unScaledminBodyRadius / scaleFactor)
    recogLimit = math.pi * (minBodyRadius**2)
    print(f"Recognition limit: {recogLimit}")
        
    lineCollection = take_slice(inputName, sliceCoord, unScaledSliceThickness, scaleFactor)
    
    if len(lineCollection) == 0:
        print("No bodies found in slice")
        add_empty_line(initialTime)
        return
        
    overalldfsk = build_projection(lineCollection, wallRadius, recogLimit)
    
    if overalldfsk is None or overalldfsk.empty:
        print("No bodies met size criteria")
        add_empty_line(initialTime)
    else:    
        overalldfsk_new = split_duplicates(overalldfsk, recogLimit)
        to_nm(overalldfsk_new, scaleFactor, initialTime)

def take_slice(inputName, sliceCoord, unScaledSliceThickness, scaleFactor): 
    '''Sorts the pixels within the PIFF file into wallText and bodyText.
       The lines within bodyText that fall within the slice are sorted into lineCollection.'''  
    wallText = []   
    bodyText = []     
    bodyWholeVol = []   
    bodyTotalVolumeNums = []
    
    print(f"Opening file: {inputName}")
    print(f"Looking for slice at coordinate: {sliceCoord}")
    
    inStream = open(inputName, "r")
    
    # Track x-coordinate ranges
    min_x = float('inf')
    max_x = float('-inf')
    total_lines = 0
    body_lines = 0
    
    for line in inStream:
        total_lines += 1
        data = line.split()
        if len(data) < 7:
            continue
        bodyID = int(data[0])
        if bodyID in bodyTotalVolumeNums:
            index = bodyTotalVolumeNums.index(bodyID)
            bodyWholeVol[index] += 1
        else:
            bodyTotalVolumeNums.append(bodyID)
            bodyWholeVol.append(1)
        if data[1] == "Wall":
            wallText.append(line)
        elif data[1] == "Body":
            body_lines += 1
            bodyText.append(line)
            x1 = float(data[2])
            x2 = float(data[3])
            min_x = min(min_x, x1, x2)
            max_x = max(max_x, x1, x2)
            if body_lines <= 5:
                print(f"Sample body coordinates: Body {bodyID} at x1={x1}, x2={x2}")
    inStream.close()
    
    print(f"Total lines read: {total_lines}")
    print(f"Body lines found: {body_lines}")
    print(f"X-coordinate range in file: {min_x} to {max_x}")
           
    bodySliceNums = []   
    bodySliceVolCounts = []  
    lineCollection = []   
    
    sliceThickness = unScaledSliceThickness / scaleFactor
    HalfSliceThickness = round((sliceThickness - 1) / 2)
    
    print(f"Slice thickness: {sliceThickness}")
    print(f"Half slice thickness: {HalfSliceThickness}")
    print(f"Looking for x coordinates between {sliceCoord - HalfSliceThickness} and {sliceCoord + HalfSliceThickness}")

    bodies_in_slice = 0
    for bodyEntry in bodyText:
        bodyLine = bodyEntry.split()
        x1 = float(bodyLine[2])
        x2 = float(bodyLine[3])
        if (min(x1, x2) <= (sliceCoord + HalfSliceThickness) and 
            max(x1, x2) >= (sliceCoord - HalfSliceThickness)):
            bodyID = int(bodyLine[0])
            bodies_in_slice += 1
            try:
                index = bodySliceNums.index(bodyID)
                bodySliceVolCounts[index] += 1
            except:
                bodySliceNums.append(bodyID)
                lineCollection.append([])
                bodySliceVolCounts.append(1)
            posi = bodySliceNums.index(bodyID)
            lineCollection[posi].append(bodyEntry)
    
    print(f"Bodies found in slice range: {bodies_in_slice}")
    print(f"Unique bodies collected: {len(lineCollection)}")
    return lineCollection

def build_projection(lineCollection, wallRadius, recogLimit):
    '''Creates a projection of each body slice and analyzes its properties'''
    overalldfsk = pd.DataFrame()      
    min_y = float('inf')
    max_y = float('-inf')
    min_z = float('inf')
    max_z = float('-inf')
    
    # Find coordinate ranges from slice data
    for array in lineCollection:
        for line in array:
            lineData = line.split()
            y1, y2 = float(lineData[4]), float(lineData[5])
            z1, z2 = float(lineData[6]), float(lineData[7])
            min_y = min(min_y, y1, y2)
            max_y = max(max_y, y1, y2)
            min_z = min(min_z, z1, z2)
            max_z = max(max_z, z1, z2)
    
    y_size = int(max_y - min_y + 3)  # +3 for padding
    z_size = int(max_z - min_z + 3)  # +3 for padding
    
    print(f"Array dimensions: {y_size} x {z_size}")
    print(f"Coordinate ranges: y={min_y} to {max_y}, z={min_z} to {max_z}")

    for array in lineCollection:  
        index2 = 0
        currentPixels = []
        projectionData = []
        projectionData.append(array[index2].split())
        lineData = array[index2].split()
        y_coord = int(float(lineData[4]) - min_y + 1)
        z_coord = int(float(lineData[6]) - min_z + 1)
        currentPixels.append([y_coord, z_coord])
        for line in array:
            lineData = array[index2].split()
            index3 = 0
            found = 0
            while index3 < len(projectionData) and found < 1:
                if (lineData[4] == projectionData[index3][4] and 
                    lineData[6] == projectionData[index3][6]):
                    found += 1
                else:
                    index3 += 1 
            if found < 1:
                projectionData.append(lineData)
                y_coord = int(float(lineData[4]) - min_y + 1)
                z_coord = int(float(lineData[6]) - min_z + 1)
                currentPixels.append([y_coord, z_coord])
            index2 += 1
        currentArea = len(projectionData)
        if currentArea >= recogLimit:    
            # Create array with proper dimensions
            imageArray = np.zeros((y_size, z_size), dtype=int)
            for pix in currentPixels:
                imageArray[pix[0], pix[1]] = 1
            all_labels = measure.label(imageArray)
            propertylist = ['label', 'bbox', 'area', 'centroid', 'convex_area',
                            'eccentricity', 'euler_number', 'filled_area',
                            'major_axis_length', 'minor_axis_length', 'perimeter']
            if np.sum(all_labels) > 0:
                props2 = measure.regionprops_table(all_labels, properties=propertylist)
                df_skimage = pd.DataFrame(props2)  
                df_skimage['imgnum'] = int(array[0].split()[0])
                overalldfsk = pd.concat([overalldfsk, df_skimage], ignore_index=True)
    return overalldfsk 

def split_duplicates(overalldfsk, recogLimit):
    '''Handles APBs that had more than 1 big-enough region by creating a pivot table
       and renaming bodies with greater than one area with new numbers.'''
    if not overalldfsk.empty:  
        big_enough = overalldfsk['area'] >= recogLimit
        overalldfsk_big_enough = overalldfsk[big_enough]
        pvt_df = overalldfsk_big_enough.pivot_table(
            values='area',
            index='imgnum',
            aggfunc=["count", "mean", "std", "max"]
        )
        pvt_df.columns = list(map("_".join, pvt_df.columns))
        pvt_df_split = pvt_df[pvt_df['count_area'] >= 2]
        if len(pvt_df_split) >= 1:
            pvt_df_single = pvt_df[pvt_df['count_area'] < 2]
            singles = pvt_df_single.index
            splits = pvt_df_split.index
            new_bod_nums_req = len(pvt_df_split)
            new_bod_nums = range(1000, 1001 + new_bod_nums_req, 1)
            overalldfsk_single = overalldfsk_big_enough[overalldfsk_big_enough["imgnum"].isin(singles)]
            overalldfsk_splits = overalldfsk_big_enough[overalldfsk_big_enough["imgnum"].isin(splits)]
            overalldfsk_splits.loc[:, "imgnum"] = new_bod_nums
            overalldfsk_new = pd.concat([overalldfsk_single, overalldfsk_splits], ignore_index=True)
        else:
            overalldfsk_new = overalldfsk_big_enough
    return overalldfsk_new
    
def to_nm(overalldfsk_new, scaleFactor, initialTime):
    '''Adjusts the area and perimeter by the scale factor to get actual nm values then exports the statistics we want'''
    overalldfsk_new = overalldfsk_new.copy()
    overalldfsk_new.loc[:, "area_scaled"] = scaleFactor**2 * overalldfsk_new["area"]
    overalldfsk_new.loc[:, "perimeter_scaled"] = scaleFactor * overalldfsk_new["perimeter"]
    overalldfsk_new.loc[:, "AR"] = overalldfsk_new["major_axis_length"] / overalldfsk_new["minor_axis_length"]
    overalldfsk_new.loc[:, "circularity"] = 4 * math.pi * overalldfsk_new["area_scaled"] / (overalldfsk_new["perimeter_scaled"]**2)
    overalldfsk_new.loc[:, "time"] = initialTime
    overalldfsk_new.rename(columns={"imgnum": "body_number"}, inplace=True)
    finalOutput = overalldfsk_new[["time", "body_number", "area_scaled", "perimeter_scaled", "circularity", "AR"]]
    print(finalOutput)
    finalOutput.to_csv("sliceData/sliceMeasurements.csv", mode='a')  

def add_empty_line(initialTime):
    '''If the slice is empty, adds a line of NAs to the frame so that AVSStats can count it as an image with no bodies'''
    data_NA = pd.DataFrame(np.nan, index=range(1), columns=["time", "body_number", "area_scaled", "perimeter_scaled", "circularity", "AR"])
    data_NA['time'] = initialTime
    print(data_NA)
    data_NA.to_csv("sliceData/sliceMeasurements.csv", mode='a')

def load_parameters_from_file(file_path):
    """
    Reads parameters from a specified file and returns them as a dictionary.
    """
    try:
        parameters = {}
        with open(file_path, 'r') as file:
            for line in file:
                line = line.split('#')[0].strip()
                if not line:
                    continue
                parts = line.split('=', 1)
                if len(parts) == 2:
                    key = parts[0].strip()
                    value = parts[1].strip().strip('"')
                    parameters[key] = value
        # Load vacuole parameters from latest run folder, if available
        runs_dir = "./runs/"
        subfolders = [f.path for f in os.scandir(runs_dir) if f.is_dir()]
        if not subfolders:
            print("Warning: No subfolders found in ./runs/. Cannot load vacuole parameters.")
            return parameters
        latest_folder = max(subfolders, key=os.path.getmtime)
        vacuole_csv_path = os.path.join(latest_folder, "Vacuole_Data_combined.csv")
        if os.path.exists(vacuole_csv_path):
            df = pd.read_csv(vacuole_csv_path)
            if not df.empty:
                latest_row = df.iloc[-1]
                parameters["Body_Radius_Mu"] = float(latest_row["Body_Radius_Mu"])
                parameters["Body_Radius_Sigma"] = float(latest_row["Body_Radius_Sigma"])
                parameters["Vacuole_x"] = float(latest_row["Vacuole_x"])
                print(f"Loaded Body_Radius_Mu: {parameters['Body_Radius_Mu']}, Body_Radius_Sigma: {parameters['Body_Radius_Sigma']}")
            else:
                print(f"Warning: {vacuole_csv_path} exists but is empty.")
        else:
            print(f"Warning: {vacuole_csv_path} not found.")
        return parameters
    except Exception as e:
        print(f"Error reading parameters file: {e}")
        return None
        
if __name__ == "__main__":
    main(fileSelectOpt=True, MassRunCheck=True, inputPiff="src/output.piff")
