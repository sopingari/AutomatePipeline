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
#   Last Date Modified: July 20th, 2021
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
''' Todo Nov. 2024:
1. update paramsfile grabs to the new parameters file, and include things including vacMin (line 116), slice thickness (line 216), and recognition limit (line 267)
2. Get values of mu and sigma for body size and number (and vacuole size? from Vacuole_gen csv and add those to the output csv (line 399))
3. Fix error handling for vacuole slice limit (line 135)
4. Verify that the mass runs are being handled appropriately
5. Add option to take serial slices?'''

#paramsFile is used to keep track of several variables used by multiple scipts.
print('Please choose the old parameters file (Old_Model_Parameters.txt)')   #### NEEDS TO BE UPDATED!
paramsFile = askopenfilename() 

#MassRunCheck is a boolean variable. It is false when SliceStats is run alone, and true when run as part of the AVS cycle.
##### Note that mass runs have not been tested recently, and will probably need to be updated to the new Vacuolegen (2024)
def main(fileSelectOpt, MassRunCheck, inputPiff):
    
    #For a given run of SliceStats_M, all body measurment output lines to the output
    #file will have the same date and time stamp.
    initialTime = time.asctime(time.localtime(time.time()))
    
    print("Now running SliceStats.py")
    #The master version of SliceStats will just use a predefined output file for easier use of this script.
    #outputName = "As of 20241112/Master/sliceData/sliceCoords.txt"  ## We shouldn't need this data, but leaving this here in case we do
    inputName = "C:/Users/sbackues/Documents/Program testing/AVS/Starting Bodies/scale 8 piffs/SphereGenRun20bod4_9_2_0Scale8.piff"
    
    if(MassRunCheck == True):
        inputName = inputPiff
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

    
    #### This will need to be updated for the new parameters file (2024)
    print("Grabbing AVS Model Parameters...\n")
    modelParams = grabParams()
    scaleFactor = modelParams[0]
    wallRadius = modelParams[1]
    centerX = modelParams[2]
    
    
    print("Current Model Parameters:\n")
    print("\tScale_Factor: %d\n" %(scaleFactor))
    print("\tWall_Radius: %d\n" %(wallRadius))
    print("\tWall_Diameter: %d\n" %(wallRadius*2))
    print("\tWall_X_Coordinate: %d\n" %(centerX))
    
    if(MassRunCheck == False):
        print(">>Would you like to use these parameters?[y/n]")
        paramSelect = input()
        
        if(paramSelect == "n"):
            print(">>Please enter new values for parameters:\n")
            print("(The Wall radius parameter value should be a post-scaling value)")
            
            print("\n>>Enter new scaling factor: ")
            scaleFactor = int(input())
            
            #The known Diameter of the simulation's Wall sphere.
            print("\n>>Enter the given wall's radius", end='')
            wallRadius = int(input())
            
            #The X value representing the X-coordinate of the Wall sphere's center.
            print("\n>>Enter the given wall's central x-coordinate:", end='')
            centerX = int(input()) 
    
    wallD = (wallRadius*2)
    
    #The minimum vacuole radius needed to perform the slice and analyze the body areas.
    #Used to define the usable range of coordinates a slice can be taken at.
    #Essentially used as a threshold within to take slices.
    unScaledVacMin = 300.0    #### This value should actually be read from the params file. 
    vacMin = (unScaledVacMin / scaleFactor)
    print("Default slice recognition limit (radius) = %dunits" %(vacMin))
    
    if(MassRunCheck == False):
        
        print(">>Would you like to use this default minimum vacuole slice threshold?[y/n]")
        minDInput = input()
    
        if(minDInput == "n" or minDInput == "N"):
            print("\n>>Enter new minimum vacuole threshold (scaled): ")
            vacMin = int(input())
    
    #Useable range of x-coordinates for the main slice.
    wallRecDiff = (wallRadius**2)-(vacMin**2)
    diamRangeVar = 0
    
    #(wallRadius**2)-(vacMin**2) must be checked to be non-negative before attempted to find its square root.
    if(wallRecDiff > 0):
        diamRangeVar = math.sqrt(wallRecDiff)
        
    if(diamRangeVar <= 0):
        sys.exit("\n!!!This model does not support a vacuole slice threshold of %s" %(vacMin))   #### Need better error handling here that doesn't crash the system
        
    #The starting and ending X-coordinates viable for a slice to be taken at.
    minX = int(centerX - diamRangeVar)
    maxX = int(centerX + diamRangeVar)
    
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
            sys.ext("\nInput was found to be invalid. Please enter in 0, 1, or 2 for your choice of slice selection method.")
        
    #The minimum recognition limit for the sub-slices of the bodies. This should be updated to pull from the paramaters text.  
    unScaledminBodyRadius = 50.0
    minBodyRadius = (unScaledminBodyRadius / scaleFactor)
    recogLimit = math.pi * (minBodyRadius**2)
        
    lineCollection = take_slice(inputName, sliceCoord, scaleFactor)   #Takes the slice, storing it as a list of bodies, each of which is a list of pixels, in lineCollection
    #print('line collection', lineCollection)
    overalldfsk = build_projection(lineCollection, wallRadius, recogLimit)    #Builds a 2D projection out of each body in the slice, turns it into a numpy array, then returns a dataframe with statistics on each body
    #print('overalldfsk', overalldfsk)
    if overalldfsk.empty == True:
        add_empty_line(initialTime)
    else:    
        overalldfsk_new = split_duplicates(overalldfsk, recogLimit)
        to_nm(overalldfsk_new, scaleFactor, initialTime)

def take_slice(inputName, sliceCoord, scaleFactor): 
    '''Sorts the pixes within the PIFF file into wallText and bodyText. The lines within bodyText are parsed through, 
    and those pixels that fall within the slice are sorted into the list of lists, lineCollection.
    lineCollection has one sublist for each body that falls into the slice, and that sublist contains all of the pixels of that 
    body that fall within the slice.'''  

    wallText = []   #Stores all lines in input PIFF file that contain Wall data.
    bodyText = []     #Stores all lines in input PIFF file that contain Body data.
    bodyWholeVol = []   #Keeps track of the entire volume of every body in the piff file, not just ones that make it into a slice.  I don't think we use this later, though.
    bodyTotalVolumeNums = []
    
    inStream = open(inputName, "r")
    '''Fills wallText and bodyText with relevent data from the input file.'''
    for line in inStream:
        data = line.split()
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
            
        if(data[1] == "Body"):
            bodyText.append(line)

    inStream.close()
           
    bodySliceNums = []   #A list of the bodies that make it into the slice. 
    bodySliceVolCounts = []  #  A list that keeps a count of the number of points for each body that fall within the slice.  This effectively keeps track of the volumes of each body's slice. 
    
    lineCollection = []   #A list of lists. Each sub-list contains all of the lines associated with a body.
    index = 0    
    unScaledSliceThickness = 70  # Changes the thickness of the slice depending on the scale.   #### This value should actually be read from the params file.  
    sliceThickness = unScaledSliceThickness / scaleFactor
    HalfSliceThickness = round((sliceThickness - 1)/2)

    '''Here is where we are building up lineCollection from bodyText, sorting into it the pixels that fall into the slice, sorted into sublists by body number'''        
    for bodyEntry in bodyText:
        bodyLine = bodyEntry.split()
        xValue = int(bodyLine[2])
        
        if(xValue <= (sliceCoord+HalfSliceThickness) and xValue >= (sliceCoord-HalfSliceThickness)):
            bodyID = int(bodyLine[0])
            
            try:
                index = bodySliceNums.index(bodyID)
                bodySliceVolCounts[index] += 1

                
            except:
                bodySliceNums.append(bodyID)
                lineCollection.append([])
                bodySliceVolCounts.append(1)
                
            posi = bodySliceNums.index(bodyID)
            lineCollection[posi].append(bodyEntry)
            
    index1 = 0
    
    ## The lines that made it into lineCollection may be recordered to the primary output file - uncomment this if so desired
    #outStream = open(outputName, "w")
    # for array in lineCollection:
    #     index2 = 0
    #     for line in lineCollection[index1]:
    #         outStream.write(lineCollection[index1][index2])
    #         index2 += 1
    #     index1 += 1
    # outStream.close()
    
    return lineCollection    #Do I also need to return bodySliceNums and bodySliceVolCounts?  Do we ever use these? I don't see them later in the code anywhere. If not, should I delete them?
    
def build_projection(lineCollection, wallRadius, recogLimit):
    '''The slice needs to be compacted to 2 dimensions, like a TEM image.  This iterates through lineCollection and builds a projection out of it.
     Each body, contained in each sub-list, is individually analyzed to determine every yz pixel that it contains, regardless of which 
     x coordinate it has (as long as the x coordinate was within the slice).  All of these yz pixels are collected together into a single list. 
     So, if that yz pixel is part of a body with an x coordinate, it will be part of that body in the final list.
     This creates a projection of the body slice - the shadow it would cast if a light were shined through it, and the pixels were completely opaque.
     The projection of each body is then turned into a numpy array (a binary image) for easier analysis.
     Last, a dataframe is where each line containts statistical measurements of a body'''

    bodyAreas = []   #Will be a list of the areas of all of the bodies.  
    bodyImages = []  #Will be a list of the numpy arrays (binary images) for each body
    overalldfsk = pd.DataFrame() #Will be a dataframe with statistics on all of the bodies (one body per line).     
    
    index1 = 0

    '''We build up the projection by going through each line in the lineCollection, checking if it is already in the projection, and, if not, adding it. 
    The outer loop is for each line in lineCollection. The inner loop goes through each line in the projectionData string and compares it to the current line from the lineCollection
    If it finds a match with identical Y and Z coordiantes it sets "found" to 1 and stops looking.  
    If it goes through the entire projectionData without finding a match, it also stops looking.
    Then it checks why it stopped looking, and if it wasn because it didn't find it (found = 0), at adds that line to the projectionData'''    
   
    ArDim = 2*(wallRadius+1) # Dimension of the array will be just slightly larger than of the simulation, to make sure that no pixels are right on the edge (needed later)

    for array in lineCollection:  #This loops over each body in lineCollection
        index2 = 0
        currentArea = 0
        currentPixels = []
        projectionData = []
        projectionData.append(array[index2].split())  #gets ProjectionData started by adding the very first line
        lineData = array[index2].split()
        Pixels = [int(lineData[4]), int(lineData[6])] #list of the just the yz pixels, for making the binary image later
        Shift = 1
        shiftedPixels = [x + Shift for x in Pixels]   #So that no pixels are right on the edge, later
        currentPixels.append(shiftedPixels) 

        for line in array:
            lineData = array[index2].split()
            currentBody = int(lineData[0])
            index3 = 0
            found = 0
            while index3 < len(projectionData) and found <1:
              
                if lineData[4] == projectionData[index3][4] and lineData[6] == projectionData[index3][6]:   #checks if both y and z match
                    found += 1
                    index3 += 1
                elif index3 < len(projectionData):    
                     index3 += 1 
            if found <1:
                projectionData.append(lineData)
                Pixels = [int(lineData[4]), int(lineData[6])]  #list of the just the yz pixels, for making the binary image later
                shiftedPixels = [x + Shift for x in Pixels]   #So that no pixels are right on the edge, later
                currentPixels.append(shiftedPixels)    
            index2 += 1

        currentArea = len(projectionData)

        if(currentArea >= recogLimit):    #Only keeps bodies that are above the recognition limit
            bodyAreas.append([currentBody, currentArea])    
            '''Now to make the imageArray for each body - a Numpy array the size of the simulation, with "1's" at every pixel location, and "0's" 
        everywhere there isn't a pixel'''
            imageArray = np.zeros ((ArDim, ArDim), dtype=int)  #Initiate a large array full of 0's
            for pix in currentPixels:
                imageArray[pix[0],pix[1]] = 1    #Put a "1" wherever the body is.  
            bodyImages.append(imageArray)
            all_labels = measure.label(imageArray)  #Labels connected regions of an integer array.  Since each array has only one body, there's no danger of confusion. 
            propertylist=['label', 'bbox', 'area', 'centroid', 'convex_area','eccentricity','euler_number','filled_area','major_axis_length','minor_axis_length','perimeter']
            # There aremore things on this property list than we currently are using, though it probably doesn't add much computation time to have them in there.
            # Area and perimeter units are pixels, not nanometers (nm); we'll translate numbers in the dataframe to nm all at once later.'''
            if( np.sum(all_labels) > 0):
                props2 = measure.regionprops_table(all_labels,properties=propertylist)
                df_skimage = pd.DataFrame(props2)  
                df_skimage['imgnum'] = currentBody
                overalldfsk = pd.concat([overalldfsk, df_skimage],ignore_index=True)
        index1 += 1
        
    return overalldfsk     #Not returning bodyAreas or bodyImages, because those aren't being used later.  

def split_duplicates(overalldfsk, recogLimit):
    '''Looks for APBs that had more than 1 big-enough region.We'll do what's called a "pivot table".
    Any imagenum with count_area >= 2 has multiple regions in it (and they are not just tiny ones, since we filtered those out already)
    Then we'll rename those bodies with greater than one area with new numbers - a new body number for each area.
    This is necessary because an irregularly shaped body could look like two separate regions in a slice, and if this were a real
    TEM image we would assume that those regions came from separate bodies and count them as such.'''
    if(overalldfsk.empty == False):  
        """ Check that each body slice is big enough, and filter it out if not."""   
        big_enough = overalldfsk['area'] >= recogLimit
        too_small = np.invert(big_enough)
        overalldfsk_big_enough = overalldfsk[big_enough]
        print(too_small)   #Should be empty!
        #assert too_small.empty == True          #It should be, since we only made numpy arrays for bodies that were over the recognition limit.  But it's not?  Check up later.        
        '''To quantify how spread-out the areas are for any bodynumber,  we'll use the statistical range (max-minus-min), which Python calls 'ptp'=peak-to-peak,
        We'll also take the StdDev, though that gives NaN when there's only 1 region for a bodynumber.''' #Why are we quantifying how spread out the areas are? 

        pvt_df=overalldfsk_big_enough.pivot_table(values='area',index='imgnum',aggfunc=["count",np.mean,np.std,np.amax,np.ptp])
        pvt_df.columns = list(map("_".join, pvt_df.columns)) #renames the columns with simpler names
        #print(pvt_df.columns)
        #print (pvt_df)  #this is a pandas dataframe
        pvt_df_split = pvt_df[pvt_df['count_area'] >= 2] #subsetting just those bodies that are split in two
        
        if (len(pvt_df_split) >= 1):
            pvt_df_single =pvt_df[pvt_df['count_area'] < 2] #subsetting just those bodies that are whole
            singles = pvt_df_single.index
            splits = pvt_df_split.index

            #now to rename the bodies that are split in parts, giving each part a new body number
            new_bod_nums_req = len(pvt_df_split) # how many new body numbers we need
            new_bod_nums = range(1000, 1001+new_bod_nums_req, 1)
            overalldfsk_single = overalldfsk_big_enough[overalldfsk_big_enough["imgnum"].isin(singles)] # filters the original list by just the single bodies
            overalldfsk_splits = overalldfsk_big_enough[overalldfsk_big_enough["imgnum"].isin(splits)] # filters the original list by just the split bodies
            overalldfsk_splits.loc[:,"imgnum"] = new_bod_nums  #giving the splits data frame the new body numbers
            overalldfsk_new = pd.concat([overalldfsk_single, overalldfsk_splits], ignore_index=True) #this has the data on all of the bodies, with unique body numbers
            print (overalldfsk_new)
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
       
        
    '''The following are all alternate output formats that were in the code as of Nov. 2024, but I commented them out because I think they are all vestigal.
        They work, but don't give us anything that we need.'''
        # numpy_array = finalOutput.to_numpy()
        # np.savetxt("test_file.txt", numpy_array, fmt = "%s")
        # finalOutputArray = finalOutput.to_numpy()
                
        # print("\nFinal Output Entries to be written to files:")
        # headerLine = "Time/Date , Body_Number , Body_Area , Body_Volume , Perimeter , Circularity , AR , Wall_Radius"
        # print(headerLine)
        
        # for ele in finalOutputArray:
        #     bodyNum = int(ele[1])
        #     bodyVolume = bodyWholeVol[bodyTotalVolumeNums.index(bodyNum)]
        #     scaledVolume = bodyVolume * scaleFactor
        #     scaledWallRadius = int(wallRadius) * scaleFactor
        #     outputLine = "%s , %s , %s , %s , %s , %s , %s, %s" %(ele[0], ele[1], ele[2], scaledVolume, ele[3], ele[4], ele[5], scaledWallRadius)
            
        #     print(outputLine)
            
        #     with open('As of 20241112/Master/sliceData/sliceDataOutput.csv', mode='a+') as csvOutputFile:
        #         sliceWriter = csv.writer(csvOutputFile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                
        #         if(os.path.getsize('As of 20241112/Master/sliceData/sliceDataOutput.csv')==0):
        #             sliceWriter.writerow(["Time/Data", "Body_Number", "Body_Area", "Body_Volume", "Perimeter", "Circularity", "AR", "Wall_Radius"])
        #         sliceWriter.writerow([ele[0], ele[1], ele[2], scaledVolume, ele[3], ele[4], ele[5], scaledWallRadius])   
            
        #     outStream2 = open("As of 20241112/Master/sliceData/sliceDefault.txt", "a+")
        #     if(os.path.getsize("As of 20241112/Master/sliceData/sliceDefault.txt")==0):
        #         outStream2.write("%s\n" %(headerLine))
        #     outStream2.write(outputLine) # Outputs each area value in the result array seperated by commas.
            
        # outStream2.close()
    
    # else:
    #     print("\n---Dataframe is empty, no bodies caught in slice.---")   #### This needs to be fixed
    


def grabParams():   # This needs to be updated to use the new parameters file
    params = []
    paramsInStream = open(paramsFile, "r")
    inStreamLines = paramsInStream.readlines()
    paramsInStream.close()
    
    for line in inStreamLines:
        splitLine = line.split()
        params.append(int(splitLine[1].strip()))
        
    return params

main(fileSelectOpt = True, MassRunCheck = False, inputPiff = None)
