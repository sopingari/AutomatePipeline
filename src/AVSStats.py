# -*- coding: utf-8 -*-
import numpy as np
import scipy.stats as stats
import statsmodels.api as sm
import matplotlib.pyplot as plt
from ridgeplot import ridgeplot
import pandas as pd 
from statsmodels.graphics.gofplots import qqplot_2samples
from tkinter import Tk
from tkinter.filedialog import askopenfilename
import os as os
import seaborn as sns
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from numpy.linalg import solve
from mpl_toolkits.mplot3d import Axes3D  

def main(fileSelectOpt = True, manual = True):
    if manual == True:  
        print(">>Please select an option: ")
        print("[1]: Estimate body size from body slice areas")
        print("[2]: Estimate body number from number of bodies per image")
        print("NOTE: You must generate an estimate of body size and vacuole size before you can estimate body number") 
        programMode = input()
        for _ in range(100):     #So that the user can run multiple tests
            print(">>Please select an option: ")
            print("[1]: Load your data") 
            print("[2]: Calculate statistics on your data")   
            print("[3]: Perform a KS (Kolmogorov-Smirnov) test and an ES (Epps-Singleton) test")  
            print("[4]: Generate a single graph to visualize the differences between your real and simulated data")
            print("[5]: Generate a set of graphs to visualize the differences between your real and simulated data for different mu and sigma combinations ")  
            print("[6]: Choose to analyze body size or body number")
            print("[0]: Exit Script")
            
            userSelection = input()
            if(userSelection == "1"):
                if programMode == "1":
                    real_slices, sim_slices, directory = loadDataArea(fileSelectOpt)
                elif programMode == "2":
                    real_slices, sim_slices, directory = loadDataNumber(fileSelectOpt)
                else: 
                    print("Please choose either option 1 or 2 by typing that number")
                    print("[1]: Estimate body size from body slice areas")
                    print("[2]: Estimate body number from number of bodies per image")
                    programMode = input() 
            if(userSelection == "2"):
                if programMode == "1":
                    print(sim_slices.head())
                    #try:
                    findAverage_size(real = real_slices, sim = sim_slices, directory = directory)
                    # except:
                        # loadDataMessage()
                elif programMode == "2":
                    print(sim_slices.head())
                    # try:
                    findAverage_num(real = real_slices, sim = sim_slices, directory = directory)
                    # except:
                        # loadDataMessage()
                else: 
                    print("Please choose either option 1 or 2 by typing that number")
                    print("[1]: Estimate body size from body slice areas")
                    print("[2]: Estimate body number from number of bodies per image")
                    programMode = input() 

            elif(userSelection == "3"):
                if programMode == "1":
                    try:
                        KS_results, ES_results = multi_compare_area(real = real_slices, sim = sim_slices, directory = directory)
                    except:
                        loadDataMessage()
                elif programMode == "2":
                    try:
                        KS_results, ES_results = multi_compare_number(real = real_slices, sim = sim_slices, directory = directory)
                    except:
                        loadDataMessage()

                else: 
                    print("Please choose either option 1 or 2 by typing that number")
                    print("[1]: Estimate body size from body slice areas")
                    print("[2]: Estimate body number from number of bodies per image")
                    programMode = input() 

            elif(userSelection == "4"):
                make_graph(real = real_slices, sim = sim_slices, directory = directory, programMode = programMode)

            elif(userSelection == "5"):
                make_graph_multi(real = real_slices, sim = sim_slices, directory = directory, programMode = programMode)

            elif(userSelection == "6"):
                print('NOTE!! You will need to reload your data after this for it to be valid')
                print(">>Please select an option: ")    
                print("[1]: Estimate body size from body slice areas")
                print("[2]: Estimate body number from number of bodies per image")
                programMode = input()
               
            elif(userSelection == "0"):
                raise SystemExit
            else:
                print("Please choose an option 0 through 7 by typing that number")

def loadDataArea(fileSelectOpt):
    if fileSelectOpt  == True:
        print(">>Select the csv file that contains your real data.  It must have the area of the body slices as the second column:")
        print("(The file selection screen may appear BEHIND your current application)")
        Tk().withdraw()
        inputFile = askopenfilename()
        directory = os.path.dirname(inputFile)  #Getting the directory of the input file, to output to later. 
        real_slices = pullData(inputFile, head = None) 
        real_slices.columns = ['image', 'area']
        real_slices = real_slices['area']
        print("Which simulated data to you want to compare to your real data?")
        print(">>Please select an option: ")    
        print("[1]: Analyze the simulated data generated by the latest run")
        print("[2]: Analyze a different set of simulated data that you have stored in a csv file")
        whichSim = input()
        if whichSim == "1":
            latest_run = os.listdir("./runs")[-1]
            print("Latest run:", latest_run)
            sim_slices = pullData(os.path.join("./runs", latest_run, "sliceMeasurements.csv"))
            print(sim_slices.head())
        elif whichSim == "2":    
            print(">>Now select the csv file that contains your simulated data.  It must have the area of the body slices in a column labeled 'area_scaled':")
            print("(The file selection screen may appear BEHIND your current application)")
            Tk().withdraw()
            inputFile = askopenfilename()
            sim_slices = pullData(inputFile)
            print(sim_slices.head())
        sim_slices = sim_slices[sim_slices.time != 'time'].dropna(axis = 0)  #Removing non-number rows (left-over headers) and NaN's.   
        sim_slices = sim_slices[['area_scaled', 'size_mu', 'size_sigma']].astype(float) #So we can filter by mu or sigma later
        print(sim_slices.head()) 
        print ('Your body area data has been loaded and is ready to use')
        return real_slices, sim_slices, directory

def loadDataNumber(fileSelectOpt):
    if fileSelectOpt  == True:
        print(">>Select the csv file that contains your real data.  It must have the number of body slices per image in the second column:")
        print("(The file selection screen may appear BEHIND your current application)")
        Tk().withdraw()
        inputFile = askopenfilename()
        directory = os.path.dirname(inputFile)  #Getting the directory of the input file, to output to later. 
        real_slices = pullData(inputFile, head = None) 
        real_slices.columns = ['image', 'number']
        # Adding zeros for the empty images. Verified to work.  
        noBodies = int(input("How many additional vacuole images were not analyzed because they did not contain any bodies?"))
        if noBodies > 0:
            real_slices = real_slices['number'].astype(int)
            print(real_slices)  
            addZeros = ([0]*noBodies)
            addZeros = pd.DataFrame(addZeros, columns=['number'], dtype=int)
            print(addZeros)
            real_slices = pd.concat([real_slices, addZeros], ignore_index = True) 
            print(real_slices)
            real_slices.columns = ['number']
        real_body_number = real_slices['number'].astype(int) 
        print(real_body_number.head())  
        print(real_body_number.value_counts())  #For verification
        print("Which simulated data to you want to compare to your real data?")
        print(">>Please select an option: ")    
        print("[1]: Analyze the simulated data generated by the latest run")
        print("[2]: Analyze a different set of simulated data that you have stored in a csv file")
        whichSim = input()
        if whichSim == "1":
            latest_run = os.listdir("./runs")[-1]
            print("Latest run:", latest_run)
            sim_slices = pullData(os.path.join("./runs", latest_run, "sliceMeasurements.csv"))
            print(sim_slices.head())
        elif whichSim == "2":
            print(">>Now select the csv file that contains your simulated data.  In this case the 'body_number' column contains a unique identifying number for each body sliced:")
            print("(The file selection screen may appear BEHIND your current application)")
            Tk().withdraw()
            inputFile = askopenfilename()
            sim_slices = pullData(inputFile)
        size_mus = sorted(sim_slices['size_mu'].value_counts().index.tolist())
        print("size mus", size_mus)
        size_sigmas = sorted(sim_slices['size_sigma'].value_counts().index.tolist())
        print("size sigmas", size_sigmas)
        number_mus = sorted(sim_slices['number_mu'].value_counts().index.tolist())
        print("number mus", number_mus)
        number_sigmas = sorted(sim_slices['number_sigma'].value_counts().index.tolist())
        print("number sigmas", number_sigmas)
        sim_body_numbers = pd.DataFrame()  #Creating an empty dataframe to hold the final body number data
        for size_mu in size_mus: 
            split_data = sim_slices.loc[sim_slices['size_mu'] == size_mu]  #splitting up the data
            for size_sigma in size_sigmas:
                split_data2 = split_data.loc[split_data['size_sigma'] == size_sigma]
                for number_mu in number_mus:
                    split_data3 = split_data2.loc[split_data2['number_mu'] == number_mu]
                    for number_sigma in number_sigmas:
                        split_slices = split_data3.loc[split_data3['number_sigma'] == number_sigma]
                        split_slices = split_slices[sim_slices.time != 'time']  #Removing non-number rows (left-over headers).  This works 
                        split_slices_noNaN = split_slices.dropna( axis = 0)   #Removing rows that had "NaN" because there were no bodies captured in that slice
                        empty_slice_num = split_slices.shape[0] - split_slices_noNaN.shape[0]   #Calculating the number of rows that had "NaN" because there were no bodies captured in that slice
                        empty_slices = [0]*empty_slice_num      #Creating a list of 0's to represent the empty slices
                        sim_body_number = pd.DataFrame({'number': split_slices_noNaN['time'].value_counts().to_list()})  #Each unique timestamp is a slice
                        sim_body_number = pd.concat([sim_body_number, pd.DataFrame({'number' : empty_slices })], ignore_index = True)   #Adding in the rows for the empty slices
                        size_mu_list = [float(size_mu)]*len(sim_body_number)  #Creating a list of the size_mu value to add to the dataframe
                        size_sigma_list = [float(size_sigma)]*len(sim_body_number)
                        number_mu_list = [float(number_mu)]*len(sim_body_number)  
                        number_sigma_list = [float(number_sigma)]*len(sim_body_number)
                        sim_body_number['size_mu'] = size_mu_list #Adding the size_mu value to the dataframe
                        sim_body_number['size_sigma'] = size_sigma_list  #Adding the size_sigma value to the dataframe
                        sim_body_number['number_mu'] = number_mu_list  #Adding the number_mu value  to the dataframe
                        sim_body_number['number_sigma'] = number_sigma_list  #Adding the number_sigma value to the dataframe
                        sim_body_numbers = pd.concat([sim_body_numbers, sim_body_number], ignore_index = True)  # putting it all together
        print(sim_body_numbers.head())  #For verification
        print ('Your body number data has been loaded and is ready to use')
        return real_body_number, sim_body_numbers, directory

def loadDataMessage():
    print('Remember that you must load your data before performing a test, and it must be in the proper format')
    print('The real data must have either body areas or body number in the second column, depending on which you are estimating.')
    print('The simulated data (from slicestats) should have columns labeled "area_scaled" and "body_number.')

def pullData(dataFile, head = 0):
    inStream = open(dataFile, "r")
    slices = pd.read_csv(inStream, header = head)
    inStream.close()
    return slices

def findAverage_size(real, sim, directory):
    data = real   
    # Summarizes the real data from the real slices
    print ("\nHere are the statistics for your real data:") 
    print(f"You have slices from {len(data)} bodies.")
    real_Average = data.mean()
    print("Average Body Size = %d" %(data.mean()))
    print("Largest Body Size = %d" %(data.max()))
    print("Smallest Body Size = %d" %(data.min()))
    real_stdDev = data.std()
    print("Standard Deviation of data set = %d" %(data.std()))

    # Outputting the real body area statistics to a csv file
    columns_real = ['length', 'average', 'largest', 'smallest', 'stdDev']
    real_results = pd.DataFrame([[len(data), data.mean(), data.max(), data.min(), data.std()]], columns = columns_real)
    with open(os.path.join(directory, 'real_body_size_statistics.csv'), 'w') as f:  
        real_results.to_csv(f, index = False) 
    print("Saved real body area statistics to 'real_body_size_statistics.csv' in the same directory as your original real body data")

    # Summarizes the simulated data from the simulated slices
    columns_sim = ['mu', 'sigma', 'length', 'average', 'largest', 'smallest', 'stdDev']
    multi_results = pd.DataFrame(columns = columns_sim)
    mus = sorted(sim['size_mu'].value_counts().index.tolist())      # Extracts all of the different values of mu, sorted
    sigmas = sorted(sim['size_sigma'].value_counts().index.tolist())      # Extracts all of the different values of sigma

    for mu in mus:
        split_data = sim.loc[sim['size_mu'] == mu]
        for sigma in sigmas:
            splitter_data = split_data.loc[split_data['size_sigma'] == sigma]
            data = splitter_data['area_scaled']
            Length = len(data)
            Average = data.mean()
            Largest = data.max()
            Smallest = data.min()
            stdDev = data.std()
            results = pd.DataFrame([[mu, sigma, Length, Average, Largest, Smallest, stdDev]], columns = columns_sim)
            multi_results = pd.concat([multi_results, results], ignore_index= True)
   
    #Doing a linear regression for mu and sigma vs mean and standard deviation for the simulated data, 
    # then finding a best fit mu and sigma from the mean and standard deviation of the real data
    if len(mus) >= 2 and len(sigmas) >= 2:
        try:
            best_mu, best_sigma, best_mean, best_SD = linear_regression(multi_results, real_Average, real_stdDev, directory)
            df_best_fit = pd.DataFrame({'mu' : [best_mu], 'sigma' : [best_sigma], 'length' : ["best fit"], 'average' : [best_mean], 'largest' : ["N/A"], 'smallest' : ["N/A"], 'stdDev' : [best_SD] })
            multi_results_full = check_in_range(multi_results, df_best_fit, best_mu, best_sigma, sim)
        except:
            print("No best fit mu and sigma could be found.")
            multi_results_full = multi_results
    elif len(sigmas) == 1:
        best_mu, best_mean = mu_linear_regression(multi_results, real_Average, directory)
        best_sigma, best_SD = np.nan, np.nan
        df_best_fit = pd.DataFrame({'mu' : [best_mu], 'sigma' : [best_sigma], 'length' : ["best fit"], 'average' : [best_mean], 'largest' : ["N/A"], 'smallest' : ["N/A"], 'stdDev' : [best_SD] })
        multi_results_full = check_in_range(multi_results, df_best_fit, best_mu, best_sigma, sim)
    elif len(mus) == 1: 
        best_sigma, best_SD = sigma_linear_regression(multi_results, real_stdDev, directory)
        best_mu, best_mean = np.nan, np.nan
        df_best_fit = pd.DataFrame({'mu' : [best_mu], 'sigma' : [best_sigma], 'length' : ["best fit"], 'average' : [best_mean], 'largest' : ["N/A"], 'smallest' : ["N/A"], 'stdDev' : [best_SD] })
        multi_results_full = check_in_range(multi_results, df_best_fit, best_mu, best_sigma, sim)

    with open(os.path.join(directory, 'sim_body_size_statistics.csv'), 'w') as f:  
        multi_results_full.to_csv(f, index = False) 
    print("Saved simulated body area statistics to 'sim_body_size_statistics.csv' in the same directory as your original real body data")

def findAverage_num(real, sim, directory):
    data = real   
    # Summarizes the real data from the real slices
    print ("\nHere are the statistics for your real data:") 
    print(f"You have slices from {len(data)} images.")
    real_Average = data.mean()
    print("Average Body Number per Slice = %d" %(data.mean()))
    print("Largest Body Number per Slice = %d" %(data.max()))
    print("Smallest Body Number per Slice = %d" %(data.min()))
    real_stdDev = data.std()
    print("Standard Deviation of data set = %d" %(data.std()))

    # Outputting the real body number statistics to a csv file
    columns_real = ['length', 'average', 'largest', 'smallest', 'stdDev']
    real_results = pd.DataFrame([[len(data), data.mean(), data.max(), data.min(), data.std()]], columns = columns_real)
    with open(os.path.join(directory, 'real_body_number_statistics.csv'), 'w') as f:  
        real_results.to_csv(f, index = False) 
    print("Saved real body number statistics to 'real_body_number_statistics.csv' in the same directory as your original real body data")

    # Summarizes the simulated data from the simulated slices
    sim_slices = sim
    
    size_mus = sorted(sim_slices['size_mu'].value_counts().index.tolist())  #extracts all of the different values of mu, sorted
    size_sigmas = sorted(sim_slices['size_sigma'].value_counts().index.tolist())
    number_mus = sorted(sim_slices['number_mu'].value_counts().index.tolist())
    number_sigmas = sorted(sim_slices['number_sigma'].value_counts().index.tolist())
    columns = ['size_mu', 'size_sigma', 'number_mu', 'number_sigma', 'length', 'average', 'largest', 'smallest', 'stdDev']
    multi_results = pd.DataFrame(columns = columns)
    df_best_fit_all = pd.DataFrame(columns = columns)
    for size_mu in size_mus: 
        split_data = sim_slices.loc[sim_slices['size_mu'] == size_mu]  #splitting up the data
        for size_sigma in size_sigmas:
            split_data2 = split_data.loc[split_data['size_sigma'] == size_sigma]
            for number_mu in number_mus:
                split_data3 = split_data2.loc[split_data2['number_mu'] == number_mu]
                for number_sigma in number_sigmas:
                    split_slices = split_data3.loc[split_data3['number_sigma'] == number_sigma]
                    data = split_slices['number']
                    Length = len(data)
                    Average = data.mean()
                    Largest = data.max()
                    Smallest = data.min()
                    stdDev = data.std()
                    results = pd.DataFrame([[size_mu, size_sigma, number_mu, number_sigma, Length, Average, Largest, Smallest, stdDev]], columns = columns)
                    multi_results = pd.concat([multi_results, results], ignore_index= True)
            if len(multi_results) > 0:
                #Doing a linear regression for number mu and number sigma vs mean and standard deviation for the simulated data, for each combination of size mu and sigma
                # then finding a best fit number mu and number sigma from the mean and standard deviation of the real data
                mus = multi_results['number_mu'].unique()
                sigmas = multi_results['number_sigma'].unique()
                if len(mus) >= 2 and len(sigmas) >= 2:
                    fittable_results = multi_results[['number_mu', 'number_sigma', 'average', 'stdDev']]
                    fittable_results.columns = ['mu', 'sigma', 'average', 'stdDev']
                    try:
                        best_mu, best_sigma, best_mean, best_SD = linear_regression(fittable_results, real_Average, real_stdDev, directory)
                        df_best_fit = pd.DataFrame({'size_mu' : [size_mu], 'size_sigma' : [size_sigma], 'number_mu' : [best_mu], 'number_sigma' : [best_sigma], 'length' : ["best fit"], 'average' : [best_mean], 'largest' : ["N/A"], 'smallest' : ["N/A"], 'stdDev' : [best_SD] })
                        include =check_in_range_number(fittable_results, best_mu, best_sigma)
                    except:
                        print("No best fit mu and sigma could be found for this combination of size mu and sigma.")
                elif len(sigmas) == 1:
                    best_mu, best_mean = mu_linear_regression(fittable_results, real_Average, directory)
                    best_sigma, best_SD = np.nan, np.nan
                    df_best_fit = pd.DataFrame({'mu' : [best_mu], 'sigma' : [best_sigma], 'length' : ["best fit"], 'average' : [best_mean], 'largest' : ["N/A"], 'smallest' : ["N/A"], 'stdDev' : [best_SD] })
                    include = check_in_range_number(fittable_results, df_best_fit, best_mu, best_sigma, sim)
                elif len(mus) == 1:
                    best_sigma, best_SD = sigma_linear_regression(fittable_results, real_stdDev, directory)
                    best_mu, best_mean = np.nan, np.nan
                    df_best_fit = pd.DataFrame({'mu' : [best_mu], 'sigma' : [best_sigma], 'length' : ["best fit"], 'average' : [best_mean], 'largest' : ["N/A"], 'smallest' : ["N/A"], 'stdDev' : [best_SD] })
                    include = check_in_range_number(fittable_results, df_best_fit, best_mu, best_sigma, sim)
            else: include = False
            if include:
                df_best_fit_all = pd.concat([df_best_fit_all, df_best_fit], ignore_index=True)

    print("Best fit parameters:")
    print(df_best_fit_all)
    print("Full results:")
    print(multi_results)
    
    multi_results_full = pd.concat([df_best_fit_all, multi_results], ignore_index=True)
    with open(os.path.join(directory, 'sim_body_number_statistics.csv'), 'w') as f:  
        multi_results_full.to_csv(f, index = False) 
    print("Saved simulated body number statistics to 'sim_body_number_statistics.csv' in the same directory as your original real body data")

def check_in_range_number(fittable_results, best_mu, best_sigma):
    # Checks if the best fit number mu and sigma are within the range of simulated values; if not, sets "include" to false so we can throw this data out.  
    if best_mu < fittable_results["mu"].min():
        print(f"The best fit mu, {best_mu}, is smaller than the smallest mu you tested.")
        print("Therefore, the best fit estimation is not valid, and will not be saved to the output file.")
        print("Consider generating more simulated data with smaller mu values.")
        include = False
    if best_sigma < fittable_results["sigma"].min():
        print(f"The best fit sigma, {best_sigma}, is smaller than the smallest sigma you tested.")
        print("Therefore, the best fit estimation is not valid, and will not be saved to the output file.")
        print("Consider generating more simulated data with smaller sigma values.")
        include = False
    if best_mu > fittable_results["mu"].max():
        print(f"The best fit mu, {best_mu}, is larger than the largest mu you tested.")
        print("Therefore, the best fit estimation is not valid, and will not be saved to the output file.")
        print("Consider generating more simulated data with larger mu values.")
        include = False
    if best_sigma > fittable_results["sigma"].max():
        print(f"The best fit sigma, {best_sigma}, is larger than the largest sigma you tested.")
        print("Therefore, the best fit estimation is not valid, and will not be saved to the output file.")
        print("Consider generating more simulated data with larger sigma values.")
        include = False
    return include

def check_in_range(multi_results, df_best_fit, best_mu, best_sigma, sim):
    # Check to make sure the mu and sigma values are within the range of the simulated data; otherwise give a message and throw out.
    print("\nHere are the statistics for your simulated data.")
    print(multi_results)
    if best_mu < sim['size_mu'].min(): 
        print(f"The best fit mu, {best_mu}, is smaller than the smallest mu you tested.")
        print("Therefore, the best fit estimation is not valid, and will not be saved to the output file.")
        print("Consider generating more simulated data with smaller mu values.")
        multi_results_full = multi_results
    elif best_sigma < sim['size_sigma'].min():
        print(f"The best fit sigma, {best_sigma}, is smaller than the smallest sigma you tested.")
        print("Therefore, the best fit estimation is not valid, and will not be saved to the output file.")
        print("Consider generating more simulated data with smaller sigma values.")
        multi_results_full = multi_results
    elif best_sigma > sim['size_sigma'].max():
        print(f"The best fit sigma, {best_sigma}, is larger than the largest sigma you tested.")
        print("Therefore, the best fit estimation is not valid, and will not be saved to the output file.")
        print("Consider generating more simulated data with larger sigma values.")
        multi_results_full = multi_results
    elif best_mu > sim['size_mu'].max():
        print(f"The best fit mu, {best_mu}, is larger than the largest mu you tested.")
        print("Therefore, the best fit estimation is not valid, and will not be saved to the output file.")
        print("Consider generating more simulated data with larger mu values.")
        multi_results_full = multi_results
    else:  
        multi_results_full = pd.concat([df_best_fit, multi_results], ignore_index=True)   
        print(f"The mu and sigma that provide the best fit to the average and standard deviation of your real data are: mu={best_mu}, sigma={best_sigma}")
    return multi_results_full

def compare_distribs_area(real, sim):
    sim = sim['area_scaled']
    ks = stats.ks_2samp(real, sim)
    es = stats.epps_singleton_2samp(real, sim)
    #print (f"The Kolmogorov-Smirnov statistic for your two data sets is {ks.statistic:.3f}, and the p-value is {ks.pvalue:.2E}. \n")
    return ks, es

def compare_distribs_number(real, sim):
    sim = sim['number']
    ks = stats.ks_2samp(real, sim)
    es = stats.epps_singleton_2samp(real, sim)
    return ks, es

def multi_compare_area(real, sim, directory):
    print(sim.head())
    mus = sim['size_mu'].value_counts().index.tolist()      # Extracts all of the different values of mu
    mu_list = sorted(mus) 
    print(mu_list)
    sigmas = sim['size_sigma'].value_counts().index.tolist()      # Extracts all of the different values of sigma
    sigma_list = sorted(sigmas)
    print(sigma_list)
    multi_ks_results = pd.DataFrame(columns = ['mu', 'sigma', 'ks', 'pvalue', 'statistic_location'])
    multi_es_results = pd.DataFrame(columns = ['mu', 'sigma', 'es_statistic', 'es_pvalue'])
    for mu in mu_list:
        split_data = sim.loc[sim['size_mu'] == mu]
        for sigma in sigma_list:
            splitter_data = split_data.loc[split_data['size_sigma'] == sigma]
            ks, es = compare_distribs_area(real, splitter_data)
            print (f"The Kolmogorov-Smirnov statistic for your real data vs the simulated data for mu = {mu} and sigma = {sigma}is {ks.statistic:.9f}, and the p-value is {ks.pvalue:.9E}, and the statistic location is {ks.statistic_location}. \n")
            ks_results = pd.DataFrame([[mu, sigma, float(ks.statistic), float(ks.pvalue)]], columns = ['mu', 'sigma', 'ks', 'pvalue'])
            multi_ks_results = pd.concat([multi_ks_results, ks_results], ignore_index= True)
            print (f"The Epps-Singleton statistic for your real data vs the simulated data for mu = {mu} and sigma = {sigma} is {es.statistic:.9f}, and the p-value is {es.pvalue:.9E}. \n")
            es_results = pd.DataFrame([[mu, sigma, float(es.statistic), float(es.pvalue)]], columns = ['mu', 'sigma', 'es_statistic', 'es_pvalue'])
            multi_es_results = pd.concat([multi_es_results, es_results], ignore_index= True)
    
    # Sorting the KS results by ks statistic and making a heatmap
    sorted_ks_results = multi_ks_results.sort_values(by = 'ks')
    print (sorted_ks_results)
    KS_heatmap(sorted_ks_results,  directory = directory)
    
    # Estimating the mu and sigma that would give the lowest possible ks statistic based on fitting a 2nd order polynomial to the data
    sorted_results = sorted_ks_results[["mu", "sigma", "ks"]]
    sorted_results.columns = ["mu", "sigma", "statistic"]
    interpolated_min_statistic, interpolated_mu_min, interpolated_sigma_min = best_fit(sorted_results, directory, stat_name = "KS")
    min_df = pd.DataFrame({'mu': [interpolated_mu_min], 'sigma' : [interpolated_sigma_min], 'ks': [interpolated_min_statistic]})
    full_ks_results = pd.concat([min_df, sorted_ks_results], ignore_index = True)
    
    #Saving the KS results 
    with open(os.path.join(directory, 'ks_results_area.csv'), 'w') as f: 
        full_ks_results.to_csv(f, index = False) 
    print("Kolmogorov-Smirnov results and heatmap saved to the same directory as your original real body data")

    # Sorting the ES results by ES statistic and making a heatmap
    sorted_es_results = multi_es_results.sort_values(by = 'es_statistic') 
    print (sorted_es_results)
    ES_heatmap(sorted_es_results, directory = directory)    
    
    # Estimating the mu and sigma that would give the lowest possible ES statistic based on fitting a 2nd order polynomial to the data
    sorted_results = sorted_es_results[["mu", "sigma", "es_statistic"]]
    sorted_results.columns = ["mu", "sigma", "statistic"]
    interpolated_min_statistic, interpolated_mu_min, interpolated_sigma_min = best_fit(sorted_results, directory, stat_name = "ES")
    min_df = pd.DataFrame({'mu': [interpolated_mu_min], 'sigma' : [interpolated_sigma_min], 'es_statistic': [interpolated_min_statistic]})
    full_es_results = pd.concat([min_df, sorted_es_results], ignore_index = True)

    #Saving the ES results
    with open(os.path.join(directory, 'es_results_area.csv'), 'w') as f:    
        full_es_results.to_csv(f, index = False)
    print("Epps-Singleton results and heatmap saved to the same directory as your original real body data")

    return sorted_ks_results, sorted_es_results 

def multi_compare_number(real, sim, directory):
    size_mu = input("Input the size mu you want to use (run the program in size mode to find this): ")
    size_sigma = input("Input the size sigma you want to use (run the program in size mode to find this): ")
    print(sim.head())
    # print(sim.dtypes)
    sim = sim[(sim['size_mu']) == float(size_mu)]  #filters the data to only include the specified size mu
    print(sim.head()) 
    sim = sim[(sim['size_sigma']) == float(size_sigma)]  #filters the data to only include the specified size sigma
    print(sim.head())
    if len(sim) == 0:
        print("No data found for that size mu and sigma combination - either choose different values or run the program in size mode to find valid values")
        return
    mus = sim['number_mu'].value_counts().index.tolist()      # Extracts all of the different values of number mu
    mu_list = sorted(mus) 
    print(mu_list)
    #print(type(mu_list))
    sigmas = sim['number_sigma'].value_counts().index.tolist()      # Extracts all of the different values of number sigma
    sigma_list = sorted(sigmas)
    print(sigma_list)
    #print(type(sigma_list))
    multi_ks_results = pd.DataFrame(columns = ['mu', 'sigma', 'ks', 'pvalue', 'statistic_location'])
    multi_es_results = pd.DataFrame(columns = ['mu', 'sigma', 'es_statistic', 'es_pvalue'])
    for mu in mu_list:
        split_data = sim.loc[sim['number_mu'] == mu]
        for sigma in sigma_list:
            splitter_data = split_data.loc[split_data['number_sigma'] == sigma]
            #splitter_data.to_csv(os.path.join(directory, f"splitter_data_number_mu{mu}_sigma{sigma}.csv"), index = False)  #Saving the data for each mu and sigma combination to a csv file for verification
            ks, es = compare_distribs_number(real, splitter_data)
            ks_results = pd.DataFrame([[mu, sigma, float(ks.statistic), float(ks.pvalue), float(ks.statistic_location)]], columns = ['mu', 'sigma', 'ks', 'pvalue', 'statistic_location'])
            es_results = pd.DataFrame([[mu, sigma, float(es.statistic), float(es.pvalue)]], columns = ['mu', 'sigma', 'es_statistic', 'es_pvalue'])
            multi_ks_results = pd.concat([multi_ks_results, ks_results], ignore_index= True)
            multi_es_results = pd.concat([multi_es_results, es_results], ignore_index= True)

    #Sorting the KS results by KS statistic and making and saving a heatmap
    sorted_ks_results = multi_ks_results.sort_values(by = 'ks') 
    print (sorted_ks_results)
    KS_heatmap(sorted_ks_results,  directory = directory)

    #Estimating the mu and sigma that would minimize the KS statistic
    sorted_results = sorted_ks_results[["mu", "sigma", "ks"]]
    sorted_results.columns = ["mu", "sigma", "statistic"]
    interpolated_min_statistic, interpolated_mu_min, interpolated_sigma_min = best_fit(sorted_results, directory, stat_name = "KS")
    min_df = pd.DataFrame({'mu': [interpolated_mu_min], 'sigma' : [interpolated_sigma_min], 'ks': [interpolated_min_statistic]})
    
    #Saving the KS results
    full_ks_results = pd.concat([min_df, sorted_ks_results], ignore_index = True)
    with open(os.path.join(directory, 'ks_results_number.csv'), 'w') as f: 
        full_ks_results.to_csv(f, index = False)
    print("Kolmogorov-Smirnov results and heatmap saved to the same directory as your original real body data")

    # Sorting the ES results by ES statistic and making a heatmap
    sorted_es_results = multi_es_results.sort_values(by = 'es_statistic') 
    print (sorted_es_results)
    ES_heatmap(sorted_es_results, directory = directory)  

    # Estimating the mu and sigma that would give the lowest possible ES statistic based on fitting a 2nd order polynomial to the data
    sorted_results = sorted_es_results[["mu", "sigma", "es_statistic"]]
    sorted_results.columns = ["mu", "sigma", "statistic"]
    interpolated_min_statistic, interpolated_mu_min, interpolated_sigma_min = best_fit(sorted_results, directory, stat_name = "ES")
    min_df = pd.DataFrame({'mu': [interpolated_mu_min], 'sigma' : [interpolated_sigma_min], 'es_statistic': [interpolated_min_statistic]})
    
    #Saving the ES results
    full_es_results = pd.concat([min_df, sorted_es_results], ignore_index = True)
    with open(os.path.join(directory, 'es_results_number.csv'), 'w') as f:    
            full_es_results.to_csv(f, index = False)
    print("Epps-Singleton results and heatmap saved to the same directory as your original real body data")

    return sorted_ks_results, sorted_es_results

def KS_heatmap(ks_results, directory):
    ks_pivot = ks_results.pivot(index = 'sigma', columns = 'mu', values = 'ks')
    ks_pivot_nums = ks_pivot.astype(float)  
    plt.figure(figsize = (10,8))
    sns.heatmap(ks_pivot_nums, annot = True, cmap = 'viridis')
    plt.title('KS statistic for different mu and sigma values')
    plt.xlabel('mu values')
    plt.ylabel('sigma values')
    plt.savefig(os.path.join(directory, f"KS_heatmap.png"))
    plt.show()

def ES_heatmap(es_results, directory): 
    es_pivot = es_results.pivot(index = 'sigma', columns = 'mu', values = 'es_statistic')
    es_pivot_nums = es_pivot.astype(float)
    plt.figure(figsize = (10,8))
    sns.heatmap(es_pivot_nums, annot = True, cmap = 'viridis')
    plt.title('Epps-Singleton statistic for different mu and sigma values')
    plt.xlabel('mu values')
    plt.ylabel('sigma values')
    plt.savefig(os.path.join(directory, f"ES_heatmap.png"))
    plt.show()

#Find best fit values of mu and sigma that minimize the statistic
def best_fit(results, directory, stat_name):
    y_predicted, poly_reg_model = predict(results)

    # This code from Gemini 2.5 (via google Colab)
    # Get the coefficients from the trained linear regression model
    coef = poly_reg_model.coef_
    intercept = poly_reg_model.intercept_

    A = np.array([
        [2 * coef[2], coef[3]],
        [coef[3], 2 * coef[4]]
    ])

    B = np.array([
        -coef[0],
        -coef[1]
    ])

    # Solve the system for mu and sigma
    try:
        optimal_mu_sigma = solve(A, B)
        interpolated_mu_min = optimal_mu_sigma[0]
        interpolated_sigma_min = optimal_mu_sigma[1]

        # Now, calculate the predicted statistic at these interpolated mu and sigma values
        # Need to transform these values into polynomial features for prediction
        pr = PolynomialFeatures(degree = 2, include_bias = False)
        # Create a 2D array for a single sample: [[interpolated_mu_min, interpolated_sigma_min]]
        optimal_features = pr.fit_transform([[interpolated_mu_min, interpolated_sigma_min]])

        interpolated_min_statistic = poly_reg_model.predict(optimal_features)[0]

        print(f"Interpolated mu for minimum statistic: {interpolated_mu_min:.4f}")
        print(f"Interpolated sigma for minimum statistic: {interpolated_sigma_min:.4f}")
        print(f"Interpolated minimum statistic: {interpolated_min_statistic:.6f}")

    except np.linalg.LinAlgError:
        print("Could not solve the system of equations. The matrix might be singular.")
    
    print("coeficients of model", coef)
    print ("intercept", intercept)
    #print("predicted y values", y_predicted)
    graph_best_fit(results, poly_reg_model, y_predicted, interpolated_min_statistic, interpolated_mu_min, interpolated_sigma_min, directory, stat_name)
    return interpolated_min_statistic, interpolated_mu_min, interpolated_sigma_min
    
#Order 2 polynomial model, predicting statistic from mu and sigma
def predict(results):
    pr = PolynomialFeatures(degree = 2, include_bias = False)
    s_poly = pr.fit_transform(results[['mu', 'sigma']])
    poly_reg_model = LinearRegression()
    poly_reg_model.fit(s_poly, results['statistic'])
    y_predicted = poly_reg_model.predict(s_poly)
    return y_predicted, poly_reg_model

def linear_regression(multi_results, real_Average, real_stdDev, directory):  # Coded with help from CoPilot
    # Creating the multi-fit linear model
    X = multi_results[['mu', 'sigma']]
    Y = multi_results[['average', 'stdDev']]
    model = LinearRegression()
    model.fit(X, Y)
    
    # Extract coefficients
    intercepts = model.intercept_
    coefs = model.coef_

    # Target values
    average_target = real_Average
    sd_target = real_stdDev

   # Build system of two equations 
    A = coefs

    b = np.array([
        average_target - intercepts[0],
        sd_target - intercepts[1]
    ])

    # Solve for [mu, sigma]
    try:
        mu, sigma = np.linalg.solve(A, b)

    except np.linalg.LinAlgError:
        print("Could not solve the system of equations. The matrix might be singular.")

    # Find the estimated mu and sigma this corresponds to (to check)
    calc_average = intercepts[0] + coefs[0,0]*mu + coefs[0,1]*sigma
    calc_stdDev = intercepts[1] + coefs[1,0]*mu + coefs [1,1]*sigma 

    graph_linear_regression (multi_results, intercepts, coefs, mu, sigma, calc_average, calc_stdDev, directory)

    return mu, sigma, calc_average, calc_stdDev

def mu_linear_regression(multi_results, real_Average, directory): 
    # A simple linear model using only mu and fitting only to the mean values
    X = multi_results[['mu']]
    Y = multi_results[['average']]
    model = LinearRegression()
    model.fit(X, Y)
    
    # Extract coefficients
    intercept = model.intercept_
    coef = model.coef_

    #Solve for best fit mu based on the real mean
    mu = (real_Average - intercept) / coef[0] 
    calc_average = intercept[0] + coef[0] * mu

    graph_mu_linear_regression(multi_results, intercept, coef, mu, calc_average, directory)

    return mu, calc_average

def sigma_linear_regression(multi_results, real_stdDev, directory):
    # A simple linear model using only sigma and fitting only to the standard deviation values
    X = multi_results[['sigma']]
    Y = multi_results[['stdDev']]
    model = LinearRegression()
    model.fit(X, Y)

    # Extract coefficients
    intercept = model.intercept_
    coef = model.coef_

    # Solve for best fit sigma based on the real standard deviation
    sigma = (real_stdDev - intercept) / coef[0]
    calc_stdDev = intercept[0] + coef[0] * sigma

    graph_sigma_linear_regression(multi_results, intercept, coef, sigma, calc_stdDev, directory)

    return sigma, calc_stdDev

def graph_mu_linear_regression(multi_results, intercept, coef, mu, calc_average, directory):
    df = multi_results.astype(float)

# Scatterplot of the simulated data
    plt.figure(figsize=(8, 6))
    plt.scatter(df["mu"], df["average"], color="blue", alpha=0.7)

    # Plot the fitted line
    x_line = np.linspace(df["mu"].min(), df["mu"].max(), 100)
    y_line = intercept[0] + coef[0] * x_line

    plt.plot(x_line, y_line, color="blue", linewidth=2)

    # Add equation text to plot
    equation = f"average = {intercept[0]:.3f} + {coef[0][0]:.3f} x mu"
    plt.text(
        0.05,
        0.95,
        equation,
        transform=plt.gca().transAxes,
        fontsize=11,
        verticalalignment="top",
        bbox=dict(facecolor="white", alpha=0.8)
    )

    # Single point to highlight
    calc_average = intercept + coef[0] * mu

    plt.scatter(
        mu,
        calc_average,
        color="red",
        s=100,
        zorder=5,
        label="Prediction"
    )

    plt.xlabel("mu")
    plt.ylabel("average")
    plt.title("Average vs Mu")
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.savefig(os.path.join(directory, f"Linear Regression plot - Mu.png"))
    plt.show()

def graph_sigma_linear_regression(multi_results, intercept, coef, sigma, calc_stdDev, directory):
    df = multi_results.astype(float)

    # Scatterplot of the simulated data
    plt.figure(figsize=(8, 6))
    plt.scatter(df["sigma"], df["stdDev"], color="blue", alpha=0.7)

    # Plot the fitted line
    x_line = np.linspace(df["sigma"].min(), df["sigma"].max(), 100)
    y_line = intercept[0] + coef[0] * x_line

    plt.plot(x_line, y_line, color="blue", linewidth=2)

    # Add equation text to plot
    equation = f"stdDev = {intercept[0]:.3f} + {coef[0][0]:.3f} x sigma"
    plt.text(
        0.05,
        0.95,
        equation,
        transform=plt.gca().transAxes,
        fontsize=11,
        verticalalignment="top",
        bbox=dict(facecolor="white", alpha=0.8)
    )

    # Single point to highlight
    calc_stdDev = intercept[0] + coef[0][0] * sigma

    plt.scatter(
        sigma,
        calc_stdDev,
        color="red",
        s=100,
        zorder=5,
        label="Prediction"
    )

    plt.xlabel("sigma")
    plt.ylabel("stdDev")
    plt.title("Standard Deviation vs Sigma")
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.savefig(os.path.join(directory, f"Linear Regression plot - Sigma.png"))
    plt.show()

def graph_linear_regression (multi_results, intercepts, coefs, mu, sigma, calc_average, calc_stdDev, directory):
    df = multi_results.astype(float)
     
    # Create grid for regression planes
    mu_min, mu_max = df["mu"].min(), df["mu"].max()
    sigma_min, sigma_max = df["sigma"].min(), df["sigma"].max()

    mu_grid, sigma_grid = np.meshgrid(
        np.linspace(mu_min, mu_max, 30),
        np.linspace(sigma_min, sigma_max, 30)
    )

    # Plane equation:
    # z = b0 + b1*mu + b2*sigma
    mean_plane = (
       intercepts[0]
        + coefs[0,0] * mu_grid
        + coefs[0,1] * sigma_grid
    )

    stdDev_plane = (
       intercepts[1]
        + coefs[1,0] * mu_grid
        + coefs[1,1] * sigma_grid
    )


    # Plot
    fig = plt.figure(figsize=(16, 7))

    # --- Plot 1: mean ---
    ax1 = fig.add_subplot(1, 2, 1, projection="3d")

    ax1.scatter(
        df["mu"],
        df["sigma"],
        df["average"],
        alpha=0.7,
        label="Simulated data"
    )

    ax1.plot_surface(
        mu_grid,
        sigma_grid,
        mean_plane,
        alpha=0.35,
        edgecolor="none"
    )


    # Highlighted point (best fit to real data)
    ax1.scatter(
        mu,
        sigma,
        calc_average,
        color="red",
        s=100,
        marker="o",
        edgecolor="black",
        label="Best Fit"
    )


    ax1.set_title("Regression Plane for mean")
    ax1.set_xlabel("mu")
    ax1.set_ylabel("sigma")
    ax1.set_zlabel("mean")

    mean_eq = (
        f"mean = {intercepts[0]:.3f} "
        f"+ {coefs[0,0]:.3f}·mu "
        f"+ {coefs[0,1]:.3f}·sigma"
    )

    best_fit_mean_desc = (f"best fit mu = {mu:.3f}; best fit sigma = {sigma:.4f}; best fit mean = {calc_average:.0f}")

    ax1.text2D(
        0.05,
        0.95,
        mean_eq,
        transform=ax1.transAxes
    )

    ax1.text2D(
        0.05,
        0.90,
        best_fit_mean_desc,
        transform=ax1.transAxes
    )

    # --- Plot 2: stdDev ---
    ax2 = fig.add_subplot(1, 2, 2, projection="3d")

    ax2.scatter(
        df["mu"],
        df["sigma"],
        df["stdDev"],
        alpha=0.7,
        label="Simulated data"
    )

    ax2.plot_surface(
        mu_grid,
        sigma_grid,
        stdDev_plane,
        alpha=0.35,
        edgecolor="none"
    )


    # Highlighted point
    ax2.scatter(
        mu,
        sigma,
        calc_stdDev,
        color="red",
        s=100,
        marker="o",
        edgecolor="black",
        label="Best Fit"
    )


    ax2.set_title("Regression Plane for stdDev")
    ax2.set_xlabel("mu")
    ax2.set_ylabel("sigma")
    ax2.set_zlabel("stdDev")

    stdDev_eq = (
        f"stdDev = {intercepts[1]:.3f} "
        f"+ {coefs[1,0]:.3f}·mu "
        f"+ {coefs[1,1]:.3f}·sigma"
    )

    best_fit_stdDev_desc = (f"best fit mu = {mu:.3f}; best fit sigma = {sigma:.4f}; best fit stdDev = {calc_stdDev:.0f}")

    ax2.text2D(
        0.05,
        0.95,
        stdDev_eq,
        transform=ax2.transAxes
    )

    ax2.text2D(
        0.05,
        0.90,
        best_fit_stdDev_desc,
        transform=ax2.transAxes
    )

    # Match viewing angles
    ax1.view_init(elev=25, azim=135)
    ax2.view_init(elev=25, azim=135)

    plt.tight_layout()
    plt.savefig(os.path.join(directory, f"Linear Regression plot.png"))
    plt.show()

def graph_best_fit (results, poly_reg_model, y_predicted, interpolated_min_statistic, interpolated_mu_min, interpolated_sigma_min, directory, stat_name):
    df = results.astype(float)

    # Create grid for regression curves
    mu_min, mu_max = df["mu"].min(), df["mu"].max()
    sigma_min, sigma_max = df["sigma"].min(), df["sigma"].max()

    mu_grid, sigma_grid = np.meshgrid(
        np.linspace(mu_min, mu_max, 30),
        np.linspace(sigma_min, sigma_max, 30)
    )

    # Curve equation    
    #z = b0 + b1*mu +b2*sigma + b3*mu**2 + b4*mu*sigma + b5*sigma^2
    intercept = poly_reg_model.intercept_
    coef = poly_reg_model.coef_
    statistic_curve = (
        intercept 
        + coef[0]*mu_grid
        + coef[1]*sigma_grid
        + coef[2]*mu_grid**2
        + coef[3]*sigma_grid*mu_grid   
        + coef[4]*sigma_grid**2
    )

    # Plot
    fig = plt.figure(figsize=(16, 16))
    ax1 = fig.add_subplot(1, 1, 1, projection="3d")

    ax1.scatter(
    df["mu"],
    df["sigma"],
    df["statistic"],
    alpha=0.7,
    label="Simulated data"
    )
    
    ax1.plot_surface(
    mu_grid,
    sigma_grid,
    statistic_curve,
    alpha=0.35,
    edgecolor="none"
    )

    # Highlighted point (best fit to real data)
    ax1.scatter(
    interpolated_mu_min,
    interpolated_sigma_min,
    interpolated_min_statistic,
    color="red",
    s=100,
    marker="o",
    edgecolor="black",
    label=f"Best Fit for {stat_name} statistic"
    )

    curve_eq = (f"statistic = {intercept:.3f} + {coef[0]:.4f}*mu + {coef[1]:.4f}*sigma + {coef[2]:.4f}*mu**2 + {coef[3]:.4f}*mu*sigma + {coef[4]:.4f}sigma**2")  

    best_fit_desc = (f"best fit mu = {interpolated_mu_min:.3f}; best fit sigma = {interpolated_sigma_min:.4f}; statistic of best fit = {interpolated_min_statistic:.5f}")

    ax1.text2D(
        0.05,
        0.95,
        curve_eq,
        transform=ax1.transAxes
    )

    ax1.text2D(
        0.05,
        0.90,
        best_fit_desc,
        transform=ax1.transAxes
    )

    plt.tight_layout()
    plt.savefig(os.path.join(directory, f"polynomial regression plot for {stat_name}.png"))
    plt.show()



def make_graph (real, sim, directory, programMode):
    print(">>Please select an option: ")
    print("[1]: Generate a Q-Q (quantile-quantile) plot)")
    print("[2]: Generate a Violin Plot")
    print("[3]: Generate a Ridgeline Plot")  
    print("[4]: Generate a CDF (cumulative distribution function) plot")    
    if programMode == "1":
        which_graph = input("Which graph would you like to generate to visualize the differences between your real and simulated body size data?")
        if which_graph == "1":
            qqPlot_area(real, sim, directory) 
        elif which_graph == "2":
            violinPlot_area(real, sim, directory, size_mu_list = None, size_sigma_list = None) 
        elif which_graph == "3":
            ridgelinePlot_area(real, sim, directory)  
        elif which_graph == "4":
            cdfPlot_area(real, sim, directory)
        else:
            print("Please choose an option 1 through 4 by typing that number")
    elif programMode == "2":
        which_graph = input("Which graph would you like to generate to visualize the differences between your real and simulated body number data?)")
        if which_graph == "1":
            qqPlot_number(real, sim, directory, size_mu = None, size_sigma = None, number_mu = None, number_sigma = None) 
        elif which_graph == "2":
            violinPlot_number(real, sim, directory, number_mu_list = None, number_sigma_list = None)  
        elif which_graph == "3":
            size_mu, size_sigma, number_mu, number_sigma = get_number_input(sim)
            number_mu_list = [number_mu]
            number_sigma_list = [number_sigma]
            ridgelinePlot_number(real, sim, directory, size_mu = size_mu, size_sigma = size_sigma, number_mu_list = number_mu_list, number_sigma_list = number_sigma_list)   
        elif which_graph == "4":
            cdfPlot_number(real, sim, directory, number_mu_list = None, number_sigma_list = None)
        else:
            print("Please choose an option 1 through 3 by typing that number")

def make_graph_multi(real, sim, directory, programMode):
    print(">>Which graph would you like to generate to visualize the differences between your real and simulated body size data?")
    print("[1]: Generate a set of Q-Q (quantile-quantile) plots of the simulated data vs the real data for different mu and sigma combinations)")
    print("[2]: Generate a set of violin plots of the simulated data vs the real data for  different mu and sigma combinations")
    print("[3]: Generate a Ridgeline Plot with all of the different mu and sigma combinations plotted together")
    print("[4]: Generate a CDF (cumulative distribution function) plot with all of the different mu and sigma combinations plotted together)")   
    which_graph = input(">>Please select an option: ") 
    if programMode == "1":
        mus = sim['size_mu'].value_counts().index.tolist()      # Extracts all of the different values of mu
        sigmas = sim['size_sigma'].value_counts().index.tolist()  # Extracts all of the different values of sigma
        mu_list = sorted(mus) 
        sigma_list = sorted(sigmas)
        if which_graph == "1":
            print("The Q-Q plots will be generated one at a time and automatically saved to the same directory as your original real body data.")
            print("You must close each Q-Q plot to see the next one.")
            for mu in mu_list:
                for sigma in sigma_list:
                    qqPlot_area(real, sim, directory, size_mu = mu, size_sigma = sigma) 
        elif which_graph == "2":
            print("The violin plots will be generated one at a time and automatically saved to the same directory as your original real body data.")
            print("You must close each violin plot to see the next one.")
            violinPlot_area(real, sim, directory, size_mu_list = mu_list, size_sigma_list = sigma_list)
        elif which_graph == "3":
            ridgelinePlot_area(real, sim, directory, size_mu_list = mu_list, size_sigma_list = sigma_list)   
        elif which_graph == "4":
            cdfPlot_area(real, sim, directory, size_mu_list = mu_list, size_sigma_list = sigma_list)
        else:
            print("Please choose an option 1 through 4 by typing that number")
    elif programMode == "2":
        mus = sim['number_mu'].value_counts().index.tolist()      # Extracts all of the different values of mu
        sigmas = sim['number_sigma'].value_counts().index.tolist()  # Extracts all of the different values of sigma
        mu_list = sorted(mus) 
        sigma_list = sorted(sigmas)
        size_mu = None
        size_sigma = None
        while size_mu is None or size_sigma is None:
            size_mu, size_sigma = get_size_input(sim)
        if which_graph == "1":
            print("The Q-Q plots will be generated one at a time and automatically saved to the same directory as your original real body data.")
            print("You must close each Q-Q plot to see the next one.")
            for mu in mu_list:
                for sigma in sigma_list:
                    qqPlot_number(real, sim, directory, size_mu = size_mu, size_sigma = size_sigma, number_mu = mu, number_sigma = sigma) 
        elif which_graph == "2":
            print("The violin plots will be generated one at a time and automatically saved to the same directory as your original real body data.")
            print("You must close each violin plot to see the next one.")
            violinPlot_number(real, sim, directory, size_mu = size_mu, size_sigma = size_sigma, number_mu_list = mu_list, number_sigma_list = sigma_list)
        elif which_graph == "3":
            ridgelinePlot_number(real, sim, directory, size_mu = size_mu, size_sigma = size_sigma, number_mu_list = mu_list, number_sigma_list = sigma_list)   
        elif which_graph == "4":
            cdfPlot_number(real, sim, directory, size_mu = size_mu, size_sigma = size_sigma, number_mu_list = mu_list, number_sigma_list = sigma_list)
        else:
            print("Please choose an option 1 through 4 by typing that number")

def get_size_input(sim):
    print ("available size mus:", list(sim['size_mu'].unique()))
    size_mu = float(input("Input the size_mu you want to use: "))
    if size_mu not in sim['size_mu'].values:
        print("No data found for that size mu - choose a different value ")
        return None, None 
    print ("available size sigmas:", list(sim['size_sigma'].unique()))
    size_sigma = float(input("Input the size_sigma you want to use: "))
    if size_sigma not in sim['size_sigma'].values:
        print("No data found for that size sigma - choose a different value")
        return None, None
    return size_mu, size_sigma

def get_number_input(sim):
    size_mus = list(sim['size_mu'].unique())
    print (f"available size mus: {size_mus}")
    size_mu = float(input("Input the size_mu you want to use: "))
    if size_mu not in sim['size_mu'].values:
        print("No data found for that size mu - choose a different value ")
        return None, None, None, None 
    
    size_sigmas = list(sim['size_sigma'].unique())
    print (f"available size sigmas: {size_sigmas}")
    size_sigma = float(input("Input the size_sigma you want to use: "))
    if size_sigma not in sim['size_sigma'].values:
        print("No data found for that size sigma - choose a different value")
        return None, None, None, None
    
    number_mus = list(sim['number_mu'].unique())
    print (f"available number mus: {number_mus}")
    number_mu = float(input("Input the number_mu you want to use: "))
    if number_mu not in sim['number_mu'].values:
        print("No data found for that number mu - choose a different value")
        return None, None, None, None  
    
    number_sigmas = list(sim['number_sigma'].unique())
    print (f"available number sigmas: {number_sigmas}")
    number_sigma = float(input("Input the number_sigma you want to use: "))
    if number_sigma not in sim['number_sigma'].values:
        print("No data found for that number sigma - choose a different value")
        return None, None, None, None
    
    return size_mu, size_sigma, number_mu, number_sigma

def get_number_only_input(sim):
    print ("available number mus:", list(sim['number_mu'].unique()))
    number_mu = float(input("Input the number_mu you want to use: "))
    if number_mu not in sim['number_mu'].values:
        print("No data found for that number mu - choose a different value")
        return None, None, None, None  
    print ("available number sigmas:", list(sim['number_sigma'].unique()))
    number_sigma = float(input("Input the number_sigma you want to use: "))
    if number_sigma not in sim['number_sigma'].values:
        print("No data found for that number sigma - choose a different value")
        return None, None, None, None
    return number_mu, number_sigma

def qqPlot_area(real, sim, directory, size_mu = None, size_sigma = None):
    if size_mu is None or size_sigma is None:
        print('Choose the values of size_mu and size_sigma you want to use for the Q-Q plot - for example, the values that gave the lowest KS statistic.')
        size_mu, size_sigma = get_size_input(sim)
    else:
        size_mu = size_mu
        size_sigma = size_sigma
    if size_mu is None or size_sigma is None:
        return
    sim = sim[(sim['size_mu']) == float(size_mu)]  #filters the data to only include the specified size mu
    sim = sim[(sim['size_sigma']) == float(size_sigma)]  #filters the data to only include the specified size sigma
    sim = sim['area_scaled']

    if len(sim) > 0:   # To skip any mu and sigma combinations that don't have any simulated data
        plotA = sm.ProbPlot(real)
        plotB = sm.ProbPlot(sim)
        qqplot_2samples(plotA,plotB, line='r', xlabel = 'Quantiles of Experimental Data', ylabel =f'Quantiles for mu = {size_mu}, sigma = {size_sigma}')  
        plt.savefig(os.path.join(directory, f"QQ_plot_area_mu{size_mu}_sigma{size_sigma}.png"))
        plt.show()

def qqPlot_number(real, sim, directory, size_mu = None, size_sigma = None, number_mu = None, number_sigma = None):
    if (number_mu is None or number_sigma is None) and (size_mu is None or size_sigma is None):
        print('Choose the values of size_mu, size_sigma, number_mu, and number_sigma you want to use for the Q-Q plot)')
        print('- for example, the values that gave the lowest KS statistic.')
        size_mu, size_sigma, number_mu, number_sigma = get_number_input(sim)
    elif (number_mu is None or number_sigma is None) and (size_mu is not None or size_sigma is not None):
        size_mu = size_mu
        size_sigma = size_sigma
        number_mu, number_sigma = get_number_only_input(sim)
    else:
        size_mu = size_mu
        size_sigma = size_sigma
        number_mu = number_mu
        number_sigma = number_sigma
    if size_mu is None or size_sigma is None or number_mu is None or number_sigma is None:
        return
    sim = sim[(sim['size_mu']) == float(size_mu)]  #filters the data to only include the specified size mu
    sim = sim[(sim['size_sigma']) == float(size_sigma)]  #filters the data to only include the specified size sigma
    sim = sim[(sim['number_mu']) == float(number_mu)]  #filters the data to only include the specified number mu
    sim = sim[(sim['number_sigma']) == float(number_sigma)]  #filters the data to only include the specified number sigma
    sim = sim['number']

    if len(sim) > 0:   # To skip any mu and sigma combinations that don't have any simulated data
        plotA = sm.ProbPlot(real)
        plotB = sm.ProbPlot(sim)
        qqplot_2samples(plotA,plotB, line='r', xlabel = 'Quantiles of Experimental Data', ylabel =f'Quantiles for number_mu = {number_mu}, sigma = {number_sigma}')  
        plt.savefig(os.path.join(directory, f"QQ_plot_numberMu{number_mu}_numberSigma{number_sigma}.png"))
        plt.show()


def violinPlot_area(real, sim, directory, size_mu_list = None, size_sigma_list = None):
    if size_mu_list is None or size_sigma_list is None:
        print('Choose the values of size_mu and size_sigma you want to use for the violin plot - for example, the values that gave the lowest KS statistic.')
        size_mu, size_sigma = get_size_input(sim)
        size_mu_list = [size_mu]
        size_sigma_list = [size_sigma]
    else:
        size_mu_list = size_mu_list
        size_sigma_list = size_sigma_list     
    if size_mu_list is None or size_sigma_list is None:       
        return
    for size_mu in size_mu_list:
        for size_sigma in size_sigma_list:
            sim_f = sim[(sim['size_mu']) == float(size_mu)]  #filters the data to only include the specified size mu
            sim_f = sim_f[(sim_f['size_sigma']) == float(size_sigma)]  #filters the data to only include the specified size sigma
            sim_f = sim_f['area_scaled']
            data = [real, sim_f]

            if len(sim_f) > 0:   # To skip any mu and sigma combinations that don't have any simulated data
                fig=plt.figure()
                ax = fig.add_subplot(111)   
                sm.graphics.violinplot(data, ax=ax, labels=["Experimental Data", f"mu = {size_mu}, sigma = {size_sigma}"])
                ax.set_xlabel("Data Sets")
                ax.set_ylabel("Body Crossectional Area (square nm)")   
                plt.savefig(os.path.join(directory, f"Violin_plot_size_mu_{size_mu}_size_sigma_{size_sigma}.png"))
                plt.show()

def violinPlot_number(real, sim, directory, size_mu = None, size_sigma = None, number_mu_list = None, number_sigma_list = None):
    if (number_mu_list is None or number_sigma_list is None) and (size_mu is None or size_sigma is None):
        print('Choose the values of size_mu, size_sigma, number_mu, and number_sigma you want to use for the violin plot)')
        print('- for example, the values that gave the lowest KS statistic.')
        size_mu, size_sigma, number_mu, number_sigma = get_number_input(sim)
        number_mu_list = [number_mu]
        number_sigma_list = [number_sigma]
    elif (number_mu_list is not None and number_sigma_list is not None) and (size_mu is None or size_sigma is None):
        print('Choose the values of size_mu and size_sigma you want to use for the violin plot - for example, the values that gave the lowest KS statistic.')
        number_mu_list = number_mu_list
        number_sigma_list = number_sigma_list    
        size_mu, size_sigma = get_size_input(sim)
    else:
        number_mu_list = number_mu_list
        number_sigma_list = number_sigma_list
        size_mu = size_mu
        size_sigma = size_sigma 
    if size_mu is None or size_sigma is None or number_mu_list is None or number_sigma_list is None:
        return
    
    for number_mu in number_mu_list:
        for number_sigma in number_sigma_list:
            sim_f = sim[(sim['size_mu']) == float(size_mu)]  #filters the data to only include the specified size mu
            sim_f = sim_f[(sim_f['size_sigma']) == float(size_sigma)]  #filters the data to only include the specified size sigma
            sim_f = sim_f[(sim_f['number_mu']) == float(number_mu)]  #filters the data to only include the specified number mu
            sim_f = sim_f[(sim_f['number_sigma']) == float(number_sigma)]  #filters the data to only include the specified number sigma
            sim_f = sim_f['number']
            data = [real, sim_f]

            if len(sim_f) > 0:   # To skip any mu and sigma combinations that don't have any simulated data
                fig=plt.figure()
                ax = fig.add_subplot(111)   
                sm.graphics.violinplot(data, ax=ax, labels=["Experimental Data", "Simulated Data"])
                ax.set_title(f"Violin_plot_number_sizeMu{size_mu}_sizeSigma{size_sigma}_numberMu{number_mu}_numberSigma{number_sigma}")
                ax.set_xlabel("Data Sets")
                ax.set_ylabel("Body Number per Slice")   
                plt.savefig(os.path.join(directory, f"Violin_plot_number_sizeMu{size_mu}_sizeSigma{size_sigma}_numberMu{number_mu}_numberSigma{number_sigma}.png"))
                plt.show()

def cdfPlot_area(real, sim, directory, size_mu_list = None, size_sigma_list = None):
    if size_mu_list is None or size_sigma_list is None:
        print('Choose the values of size_mu and size_sigma you want to use for the Q-Q plot - for example, the values that gave the lowest KS statistic.')
        size_mu, size_sigma = get_size_input(sim)
        organization = "0"
    else:
        size_mu, size_sigma = None, None
        print("Would you like each plot to show all values of size_mu, or all values of size_sigma?")
        print("Multiple plots will be generated to cover all of your data in either case")
        print("[1]: Each plot should show all values of mu for a single value of sigma")
        print("[2]: Each plot should show all values of sigma for a single value of mu")
        organization = input()
    if (size_mu is None or size_sigma is None) and (size_mu_list is None or size_sigma_list is None):
        return
    
    if organization == "0":  # Making only a single graph
        sim_f = sim[(sim['size_mu']) == float(size_mu)]  #filters the data to only include the specified size mu
        sim_f = sim_f[(sim_f['size_sigma']) == float(size_sigma)]  #filters the data to only include the specified size sigma
        sim_f = sim_f['area_scaled']

        plt.figure()
        plt.title(f"CDF Plot for mu = {size_mu} and sigma = {size_sigma}")
        plt.xlabel("Body Crossectional Area (square nm)")
        plt.ylabel("Cumulative Probability")
        plt.grid()
        plt.ecdf(real, label = 'Real Data')
        plt.ecdf(sim_f, label = 'Simulated Data')
        plt.legend()
        plt.savefig(os.path.join(directory, f"CDF_plot_area_mu{size_mu}_sigma{size_sigma}.png"))
        plt.show()

    elif organization == "1":  # One plot per sigma, with all mus plotted on each graph
        for size_sigma in size_sigma_list:
            sim_f = sim[(sim['size_sigma']) == float(size_sigma)]  #filters the data to only include the specified size sigma
            plt.figure()
            plt.title(f"CDF Plot for Size Sigma = {size_sigma}")
            plt.xlabel("Body Crossectional Area (square nm)")
            plt.ylabel("Cumulative Probability")
            plt.grid()
            plt.ecdf(real, label = 'Real Data')
            for size_mu in size_mu_list:
                sim_f2 = sim_f[(sim_f['size_mu']) == float(size_mu)]  #filters the data to only include the specified size mu    
                sim_f2 = sim_f2['area_scaled']
                if len(sim_f2) > 0:   # To skip any mu and sigma combinations that don't have any simulated data    
                    plt.ecdf(sim_f2, label = f"Mu = {size_mu}")
            plt.legend()
            plt.savefig(os.path.join(directory, f"CDF_plot_area_sizeSigma{size_sigma}.png"))
            plt.show()

    elif organization == "2":  # One plot per mu, with all sigmas plotted on each graph
        for size_mu in size_mu_list:
            sim_f = sim[(sim['size_mu']) == float(size_mu)]  #filters the data to only include the specified size mu
            plt.figure()
            plt.title(f"CDF Plot for Size Mu = {size_mu}")
            plt.xlabel("Body Crossectional Area (square nm)")
            plt.ylabel("Cumulative Probability")
            plt.grid()
            plt.ecdf(real, label = 'Real Data')
            for size_sigma in size_sigma_list:
                sim_f2 = sim_f[(sim_f['size_sigma']) == float(size_sigma)]  #filters the data to only include the specified size sigma    
                sim_f2 = sim_f2['area_scaled']
                if len(sim_f2) > 0:   # To skip any mu and sigma combinations that don't have any simulated data    
                    plt.ecdf(sim_f2, label = f"Sigma = {size_sigma}")
            plt.legend()
            plt.savefig(os.path.join(directory, f"CDF_plot_area_sizeSigma{size_sigma}.png"))
            plt.show()

def cdfPlot_number(real, sim, directory, size_mu = None, size_sigma = None, number_mu_list = None, number_sigma_list = None):
    organization = None
    if (number_mu_list is None or number_sigma_list is None) and (size_mu is None or size_sigma is None):
        print('Choose the values of size_mu, size_sigma, number_mu, and number_sigma you want to use for the CDF plot - for example, the values that gave the lowest KS statistic.')
        size_mu, size_sigma, number_mu, number_sigma = get_number_input(sim)
        organization = "0"
    elif (number_mu_list is not None and number_sigma_list is not None) and (size_mu is None or size_sigma is None):
        number_mu, number_sigma = None, None
        print('Choose the values of size_mu and size_sigma you want to use for the CDF plot - for example, the values that gave the lowest KS statistic.')
        number_mu_list = number_mu_list
        number_sigma_list = number_sigma_list    
        size_mu, size_sigma = get_size_input(sim)
    else: number_mu, number_sigma = None, None
    if organization == None:
        print("Would you like each plot to show all values of number_mu, or all values of number_sigma?")
        print("Multiple plots will be generated to cover all of your data in either case")
        print("[1]: Each plot should show all values of number_mu for a single value of number_sigma")
        print("[2]: Each plot should show all values of number_sigma for a single value of number_mu")
        organization = input()
    if (size_mu is None or size_sigma is None) or ((number_mu is None or number_sigma is None) and (number_mu_list is None or number_sigma_list is None)):
        return
    
    if organization == "0":  # Making only a single graph
        sim = sim[(sim['size_mu']) == float(size_mu)]  #filters the data to only include the specified size mu
        sim = sim[(sim['size_sigma']) == float(size_sigma)]  #filters the data to only include the specified size sigma
        sim = sim[(sim['number_mu']) == float(number_mu)]  #filters the data to only include the specified number mu
        sim = sim[(sim['number_sigma']) == float(number_sigma)]  #filters the data to only include the specified number sigma
        sim = sim['number']
        plt.figure()
        plt.title(f"CDF Plot for size mu = {size_mu}, size sigma = {size_sigma}, number mu = {number_mu}, and number sigma = {number_sigma}")
        plt.xlabel("Body Number per Slice")
        plt.ylabel("Cumulative Probability")
        plt.grid()
        plt.ecdf(real, label = 'Real Data')
        plt.ecdf(sim, label = 'Simulated Data')
        plt.legend()
        plt.savefig(os.path.join(directory, f"CDF_plot_number_sizeMu{size_mu}_sizeSigma{size_sigma}_numberMu{number_mu}_numberSigma{number_sigma}.png"))
        plt.show()
    
    elif organization == "1":  # One plot per number_sigma, with all number_mus plotted on each graph
        for number_sigma in number_sigma_list:
            sim_f = sim[(sim['number_sigma']) == float(number_sigma)]  #filters the data to only include the specified number sigma
            plt.figure()
            plt.title(f"CDF Plot for Number Sigma = {number_sigma}")
            plt.xlabel("Body Number per Slice")
            plt.ylabel("Cumulative Probability")
            plt.grid()
            plt.ecdf(real, label = 'Real Data')
            for number_mu in number_mu_list:
                sim_f2 = sim_f[(sim_f['number_mu']) == float(number_mu)]  #filters the data to only include the specified number mu    
                sim_f2 = sim_f2['number']
                if len(sim_f2) > 0:   # To skip any mu and sigma combinations that don't have any simulated data    
                    plt.ecdf(sim_f2, label = f"Number Mu = {number_mu}")
            plt.legend()
            plt.savefig(os.path.join(directory, f"CDF_plot_number_numberSigma{number_sigma}.png"))
            plt.show()
        
    elif organization == "2":  # One plot per number_mu, with all number_sigmas plotted on each graph
        for number_mu in number_mu_list:
            sim_f = sim[(sim['number_mu']) == float(number_mu)]  #filters the data to only include the specified number mu
            plt.figure()
            plt.title(f"CDF Plot for Number Mu = {number_mu}")
            plt.xlabel("Body Number per Slice")
            plt.ylabel("Cumulative Probability")
            plt.grid()
            plt.ecdf(real, label = 'Real Data')
            for number_sigma in number_sigma_list:
                sim_f2 = sim_f[(sim_f['number_sigma']) == float(number_sigma)]  #filters the data to only include the specified number sigma    
                sim_f2 = sim_f2['number']
                if len(sim_f2) > 0:   # To skip any mu and sigma combinations that don't have any simulated data    
                    plt.ecdf(sim_f2, label = f"Number Sigma = {number_sigma}")
            plt.legend()
            plt.savefig(os.path.join(directory, f"CDF_plot_number_numberMu{number_mu}.png"))
            plt.show()

def ridgelinePlot_area(real, sim, directory, size_mu_list=None, size_sigma_list=None):
    print("running ridgelinePlot_area")
    if size_mu_list is None or size_sigma_list is None:
        print('Choose the values of size_mu and size_sigma you want to use for the Ridgeline Plot - for example, the values that gave the lowest KS statistic.')
        size_mu, size_sigma = get_size_input(sim)
        size_mu_list = [size_mu]
        size_sigma_list = [size_sigma]
    else:
        size_mu_list = size_mu_list
        size_sigma_list = size_sigma_list
    if size_mu_list is None or size_sigma_list is None:
        return
    headers = []
    real_data_length = len(real)
    max_graph = np.percentile(real, 99)  #Setting the max value for the x-axis to the 99th percentile of the real data.  May increase later if simulated data needs it.  
    min_sim_length = float('inf')
    for size_mu in size_mu_list:
        for size_sigma in size_sigma_list:
            sim_filtered = sim[(sim['size_mu'] == float(size_mu)) & (sim['size_sigma'] == float(size_sigma))]
            sim_length = len(sim_filtered)
            if sim_length < min_sim_length and sim_length > 0:  #Finding the length of the smallest simulated data set that is being plotted to use for resizing
                min_sim_length = sim_length
    if real_data_length > 2*min_sim_length:
        print ("realy data length:", real_data_length)
        print ("min sim length:", min_sim_length)
        print("Your real data has more than twice as many data points as your simulated data - generate more simulated data to use this graphing method")
        return
    elif min_sim_length > 2*real_data_length:
        data_length = min_sim_length
        data = pd.DataFrame(index=range(data_length))  #Creating an empty dataframe with the number of rows equal to the length of the smallest simulated data set
        data['Experimental Data'] = np.resize(real, data_length)  #Resizing the real data to fit the length of the dataframe - this will repeat values if there are fewer real data points than the length of the dataframe. 
    elif real_data_length > min_sim_length:
        data_length = min_sim_length
        print (f"Warning: your real data has more data points than your simulated data - only the first {data_length} data points of your real data will be used for the ridgeline plot")
        print (f"Consider generating more simulated data to use all of your real data in the ridgeline plot") 
        data = pd.DataFrame(index=range(data_length))  #Creating an empty dataframe with the number of rows equal to the length of the smallest simulated data set
        data['Experimental Data'] = real[:data_length]  
    else:
        data_length = real_data_length
        data = pd.DataFrame(index=range(data_length))  #Creating an empty dataframe with the number of rows equal to the length of the smallest simulated data set
        data['Experimental Data'] = real[:data_length] 

    headers.append("Experimental Data")
    for size_mu in size_mu_list:
        for size_sigma in size_sigma_list:
            print(f"size_mu: {size_mu}, size_sigma: {size_sigma}")
            sim_filtered = sim[(sim['size_mu'] == float(size_mu)) & (sim['size_sigma'] == float(size_sigma))]
            sim1 = sim_filtered['area_scaled']
            if len(sim1) > 0:
                sim1np = sim1.to_numpy()
                max_99 = np.percentile(sim1np, 99)  #Setting the max value for the x-axis to the 99th percentile of the largest dataset to avoid outliers dominating the graph.
                if max_99 > max_graph:
                    max_graph = max_99
                data[f"mu = {size_mu}, sigma = {size_sigma}"] = sim1.to_numpy()[:data_length]   #The "to_numpy" is so that it doesn't try to line them up by index, which leads to a lot of NaN's   
                headers.append(f"mu = {size_mu}, sigma = {size_sigma}")
    samples=data.to_numpy().T

    fig = ridgeplot(
        samples=samples,
        bandwidth=40,
        kde_points=np.linspace(0, max_graph, 20),
        colorscale="viridis",
        colormode="row-index",
        opacity=0.6,
        labels=headers,
        spacing=0.5,
        )

    fig.update_layout(
        height=800,
        width=800,
        font_size=12,
        plot_bgcolor="white",
        xaxis_gridcolor="rgba(0, 0, 0, 0.1)",
        yaxis_gridcolor="rgba(0, 0, 0, 0.1)",
        showlegend=False,
    )

    fig.show()
    fig.write_image(os.path.join(directory, f"Ridgeline_plot_area.png"))

def ridgelinePlot_number(real, sim, directory, size_mu = None, size_sigma = None, number_mu_list=None, number_sigma_list=None):
    print("running ridgelinePlot_number")
    print("number_mu_list:", number_mu_list)
    print("number_sigma_list:", number_sigma_list)
    headers = []
    real_data_length = len(real)
    max_graph = np.percentile(real, 99)  #Setting the max value for the x-axis to the 99th percentile of the real data.  May increase later if simulated data needs it.  
    min_sim_length = float('inf')
    sim = sim[(sim['size_mu'] == float(size_mu)) & (sim['size_sigma'] == float(size_sigma))]  #Filtering the data to only include the specified size mu and sigma, since those are not being varied in this graph.
    for number_mu in number_mu_list:
        for number_sigma in number_sigma_list:
            sim_filtered = sim[(sim['number_mu'] == float(number_mu)) & (sim['number_sigma'] == float(number_sigma))]
            sim_length = len(sim_filtered)
            if (sim_length < min_sim_length) and (sim_length > 0):  #Finding the length of the smallest simulated data set that is being plotted to use for resizing the real data if necessary.
                min_sim_length = sim_length
    print ("realy data length:", real_data_length)
    print ("min sim length:", min_sim_length)
    if real_data_length > 2*min_sim_length:
        print("Your real data has more than twice as many data points as your simulated data - generate more simulated data to use this graphing method")
        return
    elif min_sim_length > 2*real_data_length:
        data_length = min_sim_length
        data = pd.DataFrame(index=range(data_length))  #Creating an empty dataframe with the number of rows equal to the length of the smallest simulated data set
        data['Experimental Data'] = np.resize(real, data_length)  #Resizing the real data to fit the length of the dataframe - this will repeat values if there are fewer real data points than the length of the dataframe. 
    elif real_data_length > min_sim_length:
        data_length = min_sim_length
        print (f"Warning: your real data has more data points than your simulated data - only the first {data_length} data points of your real data will be used for the ridgeline plot")
        print (f"Consider generating more simulated data to use all of your real data in the ridgeline plot") 
        data = pd.DataFrame(index=range(data_length))  #Creating an empty dataframe with the number of rows equal to the length of the smallest simulated data set
        data['Experimental Data'] = real[:data_length]  
    else:
        data_length = real_data_length
        data = pd.DataFrame(index=range(data_length))  #Creating an empty dataframe with the number of rows equal to the length of the smallest simulated data set
        data['Experimental Data'] = real[:data_length] 

    headers.append("Experimental Data")
    for number_mu in number_mu_list:
        for number_sigma in number_sigma_list:
            print(f"number_mu: {number_mu}, number_sigma: {number_sigma}")
            sim_filtered = sim[(sim['number_mu'] == float(number_mu)) & (sim['number_sigma'] == float(number_sigma))]
            sim1 = sim_filtered['number']
            if len(sim1) > 0:    # To avoid getting an error from mu and sigma combinations not included in the simulate dataset
                sim1np = sim1.to_numpy()
                max = sim1np.max()  #Setting the max value for the x-axis to the maximum value of the dataset.
                if max > max_graph:
                    max_graph = max
                data[f"mu = {number_mu}, sigma = {number_sigma}"] = sim1.to_numpy()[:data_length]   #The "to_numpy" is so that it doesn't try to line them up by index, which leads to a lot of NaN's   
                headers.append(f"mu = {number_mu}, sigma = {number_sigma}")
    print("headers:", headers)
    print("data:", data)
    print ("max value for x-axis:", max_graph)
    samples=data.to_numpy().T
    print("samples:", samples)

    fig = ridgeplot(
        samples=samples,
        bandwidth=1,
        kde_points=np.linspace(0, max_graph, 100),
        colorscale="viridis",
        colormode="row-index",
        opacity=0.6,
        labels=headers,
        spacing=0.5,
        )

    fig.update_layout(
        height=800,
        width=800,
        font_size=12,
        plot_bgcolor="white",
        xaxis_gridcolor="rgba(0, 0, 0, 0.1)",
        yaxis_gridcolor="rgba(0, 0, 0, 0.1)",
        showlegend=False,
    )

    fig.show()
    fig.write_image(os.path.join(directory, f"Ridgeline_plot_number.png"))

main()
