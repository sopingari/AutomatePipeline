# -*- coding: utf-8 -*-
import subprocess
import time
from datetime import datetime
import os
import sys
import stat
import numpy as np
import csv
import pandas as pd
import logging

    
def main():
    print("Welcome to the Autophagic Vacuole Simulation (AVS) Project")

    while True:
        print(">> Please select from the following options by entering the corresponding number:")
        print("\t[1]: Run the full pipeline with nested iterations over Body Number and Body Size")
        print("\t[2]: Run vacuole_gen alone, with nested iterations but no PIFF file generation, for testing")
        print("\t[3]: Run vacuole_gen alone, with nested iterations but no CC3D, for making PIFF files")
        print("\t[4]: Run CC3D alone in headless mode")
        print("\t[5]: Run Slice Stats alone")
        print("\t[6]: Run AVS Stats alone")
        print("\t[7]: Run Slice Stats on the vacuole_gen output, without CC3D, for testing only")
        print("\t[9]: Read the ReadMe file")
        print("\t[0]: Exit AVS")

        scriptChoice = input()

        if scriptChoice == "0":
            print("--- Now exiting AVS ---")
            sys.exit()
        elif scriptChoice == "1":
            run_pipeline(cc3d = True, PIFF = 1)   # PIFF files will be made for cc3d, but overwritten each run (not saved)
        elif scriptChoice == "2":
            run_pipeline(cc3d = False, PIFF = 0)  
        elif scriptChoice == "3":
            run_pipeline(cc3d = False, PIFF = 2)
        elif scriptChoice == "4":
            run_cc3d_script()
        elif scriptChoice == "5":
            run_SliceStats()
        elif scriptChoice == "6":
            run_AVSStats()
        elif scriptChoice == "7":
            run_pipeline(cc3d = False, PIFF = 1, SliceTest = True) # For testing SliceStats
        elif scriptChoice == "9":
            read_readme()
        else:
            print("--- Invalid Input, please select one of the above options ---")


def load_model_parameters(file_path="Model_Parameters.txt"):
    params = {}
    try:
        with open(file_path, "r") as f:
            for line in f:
                line = line.split("#")[0].strip()
                if "=" in line:
                    key, value = line.split("=", 1)
                    params[key.strip()] = value.strip()
    except FileNotFoundError:
        raise FileNotFoundError(f"Error: Parameters file '{file_path}' not found.")
    except Exception as e:
        raise ValueError(f"Error reading parameters file: {e}")
    
    return params

def get_list_to_iterate_over(start, end, step,endpoint=True):
  n = (end - start) / step
  if endpoint:
    n += 1
  if abs(n-round(n)) > 0.01:
    print(f"Strange step size: start={start} end={end} step={step} n={n}")
  return np.linspace(start, end, num=int(round(n)), endpoint=endpoint, dtype=float).tolist()


def run_pipeline(cc3d = True, PIFF = 1, SliceTest = False):
    """
    Run pipeline iterating over Body Number and Body Size using nested loops.
    """
    
    #Create a directory to store the results in
    runs_dir = 'runs'
    if not os.path.exists(runs_dir):
        os.makedirs(runs_dir)
    
    # Generate a unique overall ID based on the current date and time
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Create a new folder for this run
    run_folder = os.path.join(runs_dir, run_id)
    os.makedirs(run_folder)

    #Intializing headers for vacuole data output csv file
    vac_output_file = os.path.join(run_folder, 'Vacuole_Data_combined.csv')
    
    try:  
        with open(vac_output_file, 'a', newline='') as f:
            writer = csv.writer(f)
            
            # Columns of info to output
            writer.writerow([
                'Run_ID', 
                'timestamp', 
                'Seed', 
                'Grid_Resolution_(dx)', 
                'Body_Radius_Mu', 
                'Body_Radius_Sigma',
                'p-norm value', 
                'Average_Body_Radius', 
                'Standard_Deviation_Body_Radius',
                'Largest_Body_Radius',
                'Average_Distance_from_Origin',
                'Body_Number_Mu',
                'Body_Number_Sigma',  
                'Number_of_Bodies_Requested', 
                'Number_of_Bodies_Placed', 
                'Success_Rate_(%)', 
                'Iterations', 
                'Optimization_Max_Iterations', 
                'Starting_Objective_Function',
                'Ending_Objective_Function',
                'Optimization_Factor',
                'Compactness',
                'Wall_Radius_Mu', 
                'Wall_Radius_Sigma', 
                'Vacuole_Inner_Radius', 
                'Vacuole_Volume', 
                'Total_Body_Volume', 
                'Packing_Density_(%)', 
                'Maximum_Tries',
                'Actual_Tries',
                'Vacuole_x',
                'Vacuole_y',
                'Vacuole_z'
                ]) 
            
    except Exception as e:
        print(f"Error writing completely combined vacuole data CSV file")
        raise
                
    params = load_model_parameters("./attributes/Model_Parameters.txt")
    if not params:
        raise ValueError("Failed to load parameters from Model_Parameters.txt.")


    # Validate loop parameters to prevent infinite loops 
    assert float(params["Body_number_mu_step"]) > 0, "Body_number_mu_step must be > 0"
    assert float(params["Body_number_ending_mu"]) >= float(params["Body_number_starting_mu"]), "Body_number_ending_mu must be >= than Body_number_starting_mu."

    assert float(params["Body_number_sigma_step"]) > 0, "Body_number_sigma_step must be > 0"
    assert float(params["Body_number_ending_sigma"]) >= float(params["Body_number_starting_sigma"]), "Body_number_ending_sigma must be >= than Body_number_starting_sigma."

    assert float(params["Body_radius_mu_step"]) > 0, "Body_radius_mu_step must be > 0"
    assert float(params["Body_radius_ending_mu"]) >= float(params["Body_radius_starting_mu"]), "Body_radius_ending_mu must be >= than Body_radius_starting_mu."

    assert float(params["Body_radius_sigma_step"]) > 0, "Body_radius_sigma_step must be > 0"
    assert float(params["Body_radius_ending_sigma"]) >= float(params["Body_radius_starting_sigma"]), "Body_radius_ending_sigma must be >= than Body_radius_starting_sigma."

    sample_size = int(params["Sample_Size"])
    
    # Parameters for Wall mu and sigma 
    wall_radius_mu = float(params["Wall_Radius_mu"])
    wall_radius_sigma = float(params["Wall_Radius_sigma"])

    # Parameters for Body Number iteration
    mu_body_number_start = float(params["Body_number_starting_mu"])
    mu_body_number_end = float(params["Body_number_ending_mu"])
    mu_body_number_step = float(params["Body_number_mu_step"])
    sigma_body_number_start = float(params["Body_number_starting_sigma"])
    sigma_body_number_end = float(params["Body_number_ending_sigma"])
    sigma_body_number_step = float(params["Body_number_sigma_step"])

    # Parameters for Body Size iteration
    mu_body_size_start = float(params["Body_radius_starting_mu"])
    mu_body_size_end = float(params["Body_radius_ending_mu"])
    mu_body_size_step = float(params["Body_radius_mu_step"])
    sigma_body_size_start = float(params["Body_radius_starting_sigma"])
    sigma_body_size_end = float(params["Body_radius_ending_sigma"])
    sigma_body_size_step = float(params["Body_radius_sigma_step"])
    pvals = float(params["pvals"])
    
    #Parameters for PIFF file generation
    dx = float(params["Scale_Factor"])  #nm per voxel
    show_wall = str(params["Show_Wall"])
    
    #Parameter for number of maximum iterations for optimization
    optimmaxiter = int(params["optimmaxiter"])
    
    # Generate lists of values to iterate over
    mu_body_number_list = get_list_to_iterate_over(mu_body_number_start, mu_body_number_end, mu_body_number_step)
    sigma_body_number_list = get_list_to_iterate_over(sigma_body_number_start, sigma_body_number_end, sigma_body_number_step)
    mu_body_size_list = get_list_to_iterate_over(mu_body_size_start, mu_body_size_end, mu_body_size_step)
    sigma_body_size_list = get_list_to_iterate_over(sigma_body_size_start, sigma_body_size_end, sigma_body_size_step)

    for mu_body_number in mu_body_number_list:
        for sigma_body_number in sigma_body_number_list:
            for mu_body_size in mu_body_size_list:
                for sigma_body_size in sigma_body_size_list:
                    for run_idx in range(sample_size):
                        N_spheroids = int(np.random.lognormal(mean=mu_body_number, sigma=sigma_body_number))
                        # Ross comment: while I understand not wanting long lines (over 80 characters),
                        # it feels more systematic to have it all on one line of output
                        print(f"\nRunning pipeline with N={N_spheroids}, mu_body_number={mu_body_number}, sigma_body_number={sigma_body_number}, "
                          f"mu_body_size={mu_body_size}, sigma_body_size={sigma_body_size}, scale_factor={dx}")
                        timenow = time.strftime("%c")
                        print(f"datetime: {timenow}") # it's nice to say when we started
                        print(f"Sample {run_idx + 1}/{sample_size}")
                        vacuolegenmain(
                            run_folder = run_folder,
                            N_spheroids=N_spheroids,
                            mu_body_number=mu_body_number,
                            sigma_body_number=sigma_body_number,
                            mu_body_size=mu_body_size,
                            sigma_body_size=sigma_body_size,
                            pvals=pvals,
                            wall_radius_mu=wall_radius_mu,
                            wall_radius_sigma=wall_radius_sigma,
                            dx=dx,
                            show_wall=show_wall,
                            optimmaxiter=optimmaxiter,
                            PIFF=PIFF
                            
                        )

                        # Run CompuCell3D simulation (only if called for)
                        if cc3d == True:
                            run_cc3d_script()
                            run_SliceStats()

                        
                        #Test SliceStats (if called for)
                        if SliceTest == True:
                            run_SliceStatsTest()
                        


    print("--- Pipeline execution complete ---")


def vacuolegenmain(run_folder, N_spheroids, mu_body_number, sigma_body_number, mu_body_size, sigma_body_size, pvals, 
                    wall_radius_mu, wall_radius_sigma, dx, show_wall, optimmaxiter, PIFF):
    """
    Run vacuole_gen.py with specified parameters.
    """
    print(f"Running vacuole_gen with N={N_spheroids}, mu_body_size={mu_body_size}, sigma_body_size={sigma_body_size}")

    command = [
        sys.executable,
        "vacuole_gen.py",
        "--run_folder", str(run_folder),
        "--N", str(N_spheroids),
        "--mu", str(mu_body_size),
        "--sigma", str(sigma_body_size),
        "--pvals", str(pvals),
        "--wall_radius_mu", str(wall_radius_mu),
        "--wall_radius_sigma", str(wall_radius_sigma),
        "--mu_body_number", str(mu_body_number),
        "--sigma_body_number", str(sigma_body_number),
        "--dx", str(dx),
        "--show_wall", str(show_wall),
        "--optimmaxiter", str(optimmaxiter),
        "--PIFF", str(PIFF)
        
    ]

    try:
        subprocess.run(command, check=True)
        print(f"vacuole_gen.py executed successfully with N={N_spheroids}, mu={mu_body_size}, sigma={sigma_body_size}")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running vacuole_gen.py: {e}")

def run_vacuole_gen():
    """
    Run vacuole_gen, iterating over Body Number and Body Size using nested loops to create PIFF files, but not squashing them with CC3D or sliceing with slicestats.
    """    
    run_pipeline(cc3d = False)
    
def run_cc3d_script():
    print("Running CC3D simulation using runScript.sh...")
    cc3d_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CompuCell3D')
    run_script = os.path.join(cc3d_folder, 'runScript.sh')

    if not os.path.exists(run_script):
        print("Error: runScript.sh not found in the CompuCell3D folder.")
        return

    st = os.stat(run_script)
    os.chmod(run_script, st.st_mode | stat.S_IEXEC)

    command = [run_script, '-i', './cc3dSimulation/CC3D_sim_for_AVS.cc3d', '-f', '100', '-o', './Output/', '-c', 'infoPrinter']

    try:
        start_time = time.time()
        subprocess.run(command, check=True)
        end_time = time.time()
        print(f"CC3D simulation completed in {end_time - start_time:.2f} seconds.")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running runScript.sh: {e}")


def run_SliceStats():
    print("Running SliceStats...")
    command = [sys.executable, 'SliceStats.py']

    try:
        subprocess.run(command, check=True)
        print("SliceStats executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running SliceStats.py: {e}")

def run_SliceStatsTest():
    print("Running SliceStatsTest...")
    command = [sys.executable, 'SliceStatsTest.py']

    try:
        subprocess.run(command, check=True)
        print("SliceStatsTest executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running SliceStatsTest.py: {e}")


def run_AVSStats():
    print("Running AVSStats...")
    command = [sys.executable, 'AVSStats.py']

    try:
        subprocess.run(command, check=True)
        print("AVSStats executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running AVSStats.py: {e}")


def read_readme():
    print("Reading the ReadMe file...")
    readme_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'README.txt')
    if os.path.exists(readme_file):
        with open(readme_file, 'r') as f:
            print(f.read())
    else:
        print("README file not found.")


if __name__ == "__main__":
    main()
