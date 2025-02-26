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
#   within this slice are then analyzed to determine their relative areas. This area data is the recorded
#   and reported for later compilation and statistical analysis.
############################################################################################################
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.
############################################################################################################
''' Todo Dec. 2024:
1. Fix lines 68-82: right now they are setting the slice size based on the body cluster, not the wall. Make them read the vacuole size from the combined csv made by vacuole_gen
2. Get values of mu and sigma for body size and number from Vacuole_gen csv and add those to the output csv (line 399))
3. Make it automatically read the CC3D PIFF-dumped file when mass-runs = true (And in general, mass runs need to not ask for any user input).  
4. Fix error handling for vacuole slice limit (line 135) - should ouput nothing into the csv file.  But could print a message saying that this vacuole is being skipped because too small.  
5. Add option to take serial slices?'''

#paramsFile is used to keep track of several variables used by multiple scipts.
paramsFile = './attributes/Model_Parameters.txt'   # For the linux server
#paramsFile = 'src/attributes/Model_Parameters.txt'   #For Windows

#MassRunCheck is a boolean variable. It is false when SliceStats is run alone, and true when run as part of the AVS cycle.
##### Note that mass runs have not been tested recently, and will probably need to be updated to the new Vacuolegen (2024)
def main(fileSelectOpt, MassRunCheck, inputPiff):
    initialTime = time.asctime(time.localtime(time.time()))

    if(MassRunCheck == True):
        Tk().withdraw()
        filename = askopenfilename()    #This is only temporary until we figure out how to have it automatically grab from the current run folder  
        inputName = filename 
    else:
        if(fileSelectOpt):
            print("Please Select piff file...")
            print("(The file selection screen may appear BEHIND your current application)")
            Tk().withdraw()
            filename = askopenfilename()
            inputName = filename
            print("Given File Name: %s" %(filename))
        else: 
            print(">>Enter INPUT PIF file path+name:")
            inputName = input()

    # First scan file to get actual coordinate range
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
    
    # Use actual file dimensions for parameters
    centerX = int((min_x + max_x) / 2)   #Needs to be in pixels
    wallRadius = int((max_x - min_x) / 2)   #Needs to be in pixels
    
    print("Grabbing AVS Model Parameters...\n")
    modelParams = load_parameters_from_file(paramsFile)
    print(modelParams)
    scaleFactor = int(modelParams['Scale_Factor'])  # Keep scale factor from params file
    unScaledSliceThickness = int(modelParams['unScaledSliceThickness'])
    unScaledVacMin = int(modelParams['unScaledVacMin'])
    unScaledminBodyRadius = int(modelParams['unScaledminBodyRadius'])

    
    print("Current Model Parameters:\n")
    print("\tScale_Factor: %d\n" %(scaleFactor))
    print("\t:unScaledSliceThickness %d\n" %(unScaledSliceThickness))
    print("\t:unScaledVacMin %d\n" %(unScaledVacMin))
    print("\t:unScaledminBodyRadius %d\n" %(unScaledminBodyRadius))
    
    if(MassRunCheck == False):
        print(">>Would you like to use these parameters?[y/n]")
        paramSelect = input()
        
        if(paramSelect == "n"):
            print(">>Please enter new values for parameters:\n")
            print("(The Wall radius parameter value should be a post-scaling value)")
            
            print("\n>>Enter new scaling factor: ")
            scaleFactor = int(input())
            
            print("\n>>Enter the given wall's radius", end='')
            wallRadius = int(input())
            
            print("\n>>Enter the given wall's central x-coordinate:", end='')
            centerX = int(input()) 
    
    unScaledVacMin = 300.0    
    vacMin = (unScaledVacMin / scaleFactor)
    print("Default slice recognition limit (radius) = %dunits" %(vacMin))
    
    if(MassRunCheck == False):
        print(">>Would you like to use this default minimum vacuole slice threshold?[y/n]")
        minDInput = input()
    
        if(minDInput == "n" or minDInput == "N"):
            print("\n>>Enter new minimum vacuole threshold (scaled): ")
            vacMin = int(input())
    
    # Adjust recognition limit to match file scale
    wallRecDiff = (wallRadius**2)-(vacMin**2)
    diamRangeVar = 0
    
    if(wallRecDiff > 0):
        diamRangeVar = math.sqrt(wallRecDiff)
    else:
        print("\n!!!Using full coordinate range due to small wall radius")  # What this should do instead is skip that vacuole and output nothing (not even a time stamp)
        diamRangeVar = wallRadius
        
    minX = max(int(centerX - diamRangeVar), int(min_x))
    maxX = min(int(centerX + diamRangeVar), int(max_x))
    
    print(f"Adjusted valid slice range: {minX} to {maxX}")
    
    sliceCoord = -1
    
    if(MassRunCheck == True):
        sliceCoord = random.randint(minX, maxX)
    else:
        print(">>Finally, select an option for determining where a slice will be taken:")
        print("\t[0 for slice to be taken at centerX coordinate]")
        print("\t[1 for slice to be taken at a randomly selected coordinate]")
        print("\t[2 for slice to be taken at a user specified coordinate]")
            
        sliceChoice = int(input())
        sliceCoord = 0
        
        if(sliceChoice == 0):
            sliceCoord = centerX
            
        elif(sliceChoice == 1):
            sliceCoord = random.randint(minX, maxX)
            
        elif(sliceChoice == 2):
            print("\n>>Enter the x coordinate you'd like the slice to be taken at:")
            sliceCoord = int(input())
        else:
            print("\nInput was found to be invalid. Please enter in 0, 1, or 2 for your choice of slice selection method.")
            return
            
    print(f"Taking slice at X coordinate: {sliceCoord}")
    
    # Adjust recognition limit
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
    '''Sorts the pixels within the PIFF file into wallText and bodyText. The lines within bodyText are parsed through, 
    and those pixels that fall within the slice are sorted into lineCollection.'''  

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
        idCheck = bodyID in bodyTotalVolumeNums
        
        if(len(bodyTotalVolumeNums) > 0 and (idCheck == True)):
            index = bodyTotalVolumeNums.index(bodyID)
            bodyWholeVol[index] += 1
        else:
            bodyTotalVolumeNums.append(bodyID)
            bodyWholeVol.append(1)

        if(data[1] == "Wall"):
            wallText.append(line)
        elif(data[1] == "Body"):
            body_lines += 1
            bodyText.append(line)
            # Track x-coordinate range
            x1 = float(data[2])
            x2 = float(data[3])
            min_x = min(min_x, x1, x2)
            max_x = max(max_x, x1, x2)
            
            # Print first few body coordinates for debugging
            if body_lines <= 5:
                print(f"Sample body coordinates: Body {bodyID} at x1={x1}, x2={x2}")

    inStream.close()
    
    print(f"Total lines read: {total_lines}")
    print(f"Body lines found: {body_lines}")
    print(f"X-coordinate range in file: {min_x} to {max_x}")
           
    bodySliceNums = []   
    bodySliceVolCounts = []  
    lineCollection = []   
    index = 0    
    
 
    sliceThickness = unScaledSliceThickness / scaleFactor
    HalfSliceThickness = round((sliceThickness - 1)/2)
    
    print(f"Slice thickness: {sliceThickness}")
    print(f"Half slice thickness: {HalfSliceThickness}")
    print(f"Looking for x coordinates between {sliceCoord-HalfSliceThickness} and {sliceCoord+HalfSliceThickness}")

    bodies_in_slice = 0
    
    for bodyEntry in bodyText:
        bodyLine = bodyEntry.split()
        x1 = float(bodyLine[2])
        x2 = float(bodyLine[3])
        
        if (min(x1, x2) <= (sliceCoord+HalfSliceThickness) and 
            max(x1, x2) >= (sliceCoord-HalfSliceThickness)):
            
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

    bodyAreas = []   
    bodyImages = []  
    overalldfsk = pd.DataFrame()      
    
    # Get the actual coordinate ranges from the data
    min_y = float('inf')
    max_y = float('-inf')
    min_z = float('inf')
    max_z = float('-inf')
    
    # First pass to find coordinate ranges
    for array in lineCollection:
        for line in array:
            lineData = line.split()
            y1, y2 = float(lineData[4]), float(lineData[5])
            z1, z2 = float(lineData[6]), float(lineData[7])
            min_y = min(min_y, y1, y2)
            max_y = max(max_y, y1, y2)
            min_z = min(min_z, z1, z2)
            max_z = max(max_z, z1, z2)
    
    # Calculate array dimensions based on actual coordinates
    y_size = int(max_y - min_y + 3)  # +3 for padding
    z_size = int(max_z - min_z + 3)  # +3 for padding
    
    print(f"Array dimensions: {y_size} x {z_size}")
    print(f"Coordinate ranges: y={min_y} to {max_y}, z={min_z} to {max_z}")

    index1 = 0
    for array in lineCollection:  
        index2 = 0
        currentArea = 0
        currentPixels = []
        projectionData = []
        projectionData.append(array[index2].split())
        lineData = array[index2].split()
        
        # Adjust coordinates to array indices
        y_coord = int(float(lineData[4]) - min_y + 1)
        z_coord = int(float(lineData[6]) - min_z + 1)
        
        currentPixels.append([y_coord, z_coord])

        for line in array:
            lineData = array[index2].split()
            currentBody = int(lineData[0])
            index3 = 0
            found = 0
            while index3 < len(projectionData) and found < 1:
                if (lineData[4] == projectionData[index3][4] and 
                    lineData[6] == projectionData[index3][6]):
                    found += 1
                    index3 += 1
                elif index3 < len(projectionData):    
                    index3 += 1 
            if found < 1:
                projectionData.append(lineData)
                y_coord = int(float(lineData[4]) - min_y + 1)
                z_coord = int(float(lineData[6]) - min_z + 1)
                currentPixels.append([y_coord, z_coord])
            index2 += 1

        currentArea = len(projectionData)

        if(currentArea >= recogLimit):    
            bodyAreas.append([currentBody, currentArea])    
            
            # Create array with proper dimensions
            imageArray = np.zeros((y_size, z_size), dtype=int)
            
            for pix in currentPixels:
                imageArray[pix[0], pix[1]] = 1
                
            bodyImages.append(imageArray)
            all_labels = measure.label(imageArray)
            
            propertylist=['label', 'bbox', 'area', 'centroid', 'convex_area',
                         'eccentricity', 'euler_number', 'filled_area',
                         'major_axis_length', 'minor_axis_length', 'perimeter']
            
            if(np.sum(all_labels) > 0):
                props2 = measure.regionprops_table(all_labels, properties=propertylist)
                df_skimage = pd.DataFrame(props2)  
                df_skimage['imgnum'] = currentBody
                overalldfsk = pd.concat([overalldfsk, df_skimage], ignore_index=True)
        index1 += 1
        
    return overalldfsk 

def split_duplicates(overalldfsk, recogLimit):
    '''Handles APBs that had more than 1 big-enough region by creating a pivot table
    and renaming bodies with greater than one area with new numbers.'''
    
    if(overalldfsk.empty == False):  
        # Filter for big enough areas
        big_enough = overalldfsk['area'] >= recogLimit
        overalldfsk_big_enough = overalldfsk[big_enough]

        # Create pivot table with just the essential aggregations we need
        pvt_df = overalldfsk_big_enough.pivot_table(
            values='area',
            index='imgnum',
            aggfunc=["count", "mean", "std", "max"]
        )
        
        # Rename columns for clarity
        pvt_df.columns = list(map("_".join, pvt_df.columns))
        
        # Find bodies that are split in two or more parts
        pvt_df_split = pvt_df[pvt_df['count_area'] >= 2]
        
        if (len(pvt_df_split) >= 1):
            # Get single bodies
            pvt_df_single = pvt_df[pvt_df['count_area'] < 2]
            singles = pvt_df_single.index
            splits = pvt_df_split.index

            # Generate new body numbers for split bodies
            new_bod_nums_req = len(pvt_df_split)
            new_bod_nums = range(1000, 1001+new_bod_nums_req, 1)
            
            # Filter and combine the data
            overalldfsk_single = overalldfsk_big_enough[overalldfsk_big_enough["imgnum"].isin(singles)]
            overalldfsk_splits = overalldfsk_big_enough[overalldfsk_big_enough["imgnum"].isin(splits)]
            overalldfsk_splits.loc[:,"imgnum"] = new_bod_nums
            overalldfsk_new = pd.concat([overalldfsk_single, overalldfsk_splits], ignore_index=True)
        else:
            overalldfsk_new = overalldfsk_big_enough
            
    return overalldfsk_new
    
def to_nm(overalldfsk_new, scaleFactor, initialTime):
    '''Adjusts the area and perimeter by the scale factor to get actual nm values then exports the statistics we want'''   #This is where we will need to add the parameters like mu and sigma
    overalldfsk_new["area_scaled"] = scaleFactor**2*overalldfsk_new["area"]
    overalldfsk_new["perimeter_scaled"] = scaleFactor*overalldfsk_new["perimeter"]
    # Now to do some more calculations to get some of the data I want, Aspect Ratio (AR) and Circularity 
    overalldfsk_new["AR"]=overalldfsk_new["major_axis_length"] / overalldfsk_new["minor_axis_length"]  #Adds Aspect ratio column
    overalldfsk_new["circularity"]= 4*math.pi*overalldfsk_new["area_scaled"] / (overalldfsk_new["perimeter_scaled"]**2)  #Adds circularity column
    overalldfsk_new["time"] = initialTime
    overalldfsk_new.rename(columns = {"imgnum":"body_number"}, inplace=True)
    # Now to export just what we want, in a nice format
    finalOutput = overalldfsk_new[["time", "body_number", "area_scaled", "perimeter_scaled", "circularity", "AR"]]   #### Will need to include more columns with parameters from which the bodies were made
    print (finalOutput)
    finalOutput.to_csv ("sliceData/sliceMeasurements.csv", mode='a')  

def add_empty_line(initialTime):
    '''If the slice is empty, adds a line of NAs to the frame so that AVSStats can count it as an image with no bodies'''
    data_NA = pd.DataFrame(np.nan, index = range(1), columns = ["time", "body_number", "area_scaled", "perimeter_scaled", "circularity", "AR"])
    data_NA['time'] = initialTime
    print(data_NA)
    data_NA.to_csv ("sliceData/sliceMeasurements.csv", mode='a')

def load_parameters_from_file(file_path):
    """
    Reads parameters from a specified file and returns them as a dictionary.

    Parameters:
        file_path (str): Path to the file containing parameters.

    Returns:
        dict: Dictionary with parameter names as keys and their values.
    """
    try:
        parameters = {}
        with open(file_path, 'r') as file:
            for line in file:
                # Remove comments and trim whitespace
                line = line.split('#')[0].strip()
                
                # Skip empty lines
                if not line:
                    continue
                
                # Split on the first '=' only
                parts = line.split('=', 1)
                if len(parts) == 2:
                    key = parts[0].strip()
                    value = parts[1].strip().strip('"')  # Remove quotes if present
                    parameters[key] = value
        
        return parameters
    except Exception as e:
        print(f"Error reading parameters file: {e}")
        return None


main(fileSelectOpt = True, MassRunCheck = True, inputPiff = "src/output.piff")
