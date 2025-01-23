# -*- coding: utf-8 -*-
import subprocess
import time
import os
import sys
import stat
import numpy as np

############################################################################################################
#   Autophagic Vacuole Simulation (AVS) Project
#   Updated Pipeline_Interface.py script
#   Last Date Modified: [Current Date]
#
#   This script serves as the master control program for the AVS project, providing a menu-driven interface
#   to run various components of the simulation pipeline.
#
#   Note: The interface (Pipeline_Interface.py), vacuole_gen.py, and ccRunScript.sh are all located within the same src folder.
#         All outputs will remain within the src folder as well.
#
#   Important Note: AVS software was developed and tested primarily on Windows OS. AVS has not (yet) been
#       tested on Mac or Linux Operating Systems.
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

def main():
    print("Welcome to the Autophagic Vacuole Simulation (AVS) Project")

    while True:
        print(">> Please select from the following options by entering the corresponding number:")
        print("\t[1]: Run pipeline based on N input (1 -> N)")
        print("\t[2]: Run vacuole_gen alone (checking outputs)")
        print("\t[3]: Run CC3D alone in headless mode")
        print("\t[4]: Run Slice Stats alone")
        print("\t[5]: Run AVS Stats alone")
        print("\t[6]: Read the ReadMe file")
        print("\t[0]: Exit AVS")

        scriptChoice = input()

        if scriptChoice == "0":
            print("--- Now exiting AVS ---")
            sys.exit()
        elif scriptChoice == "1":
            option_one()
        elif scriptChoice == "2":
            option_two()
        elif scriptChoice == "3":
            option_three()
        elif scriptChoice == "4":
            option_four()
        elif scriptChoice == "5":
            option_five()
        elif scriptChoice == "6":
            option_six()
        else:
            print("--- Invalid Input, please select from options 0 to 6 ---")

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
    
def vacuolegenmain(mu_body_number, sigma_body_number, mu_body_size, sigma_body_size, wall_radius_mu, wall_radius_sigma, fixed_N=None):
    # Generate number of spheroids
    if fixed_N is not None:
        # Use fixed number of spheroids when iterating over Body Size
        N_spheroids = fixed_N
    else:
        # Generate number of spheroids dynamically when iterating over Body Number
        N_spheroids = int(np.random.lognormal(mean=mu_body_number, sigma=sigma_body_number))

    print(f"Running vacuole_gen with N={N_spheroids}, mu_body_size={mu_body_size}, sigma_body_size={sigma_body_size}")

    # Build the command to execute vacuole_gen.py
    command = [
        sys.executable,
        "vacuole_gen.py",
        "--N", str(N_spheroids),               # Number of spheroids
        "--mu", str(mu_body_size),            # Mean for body size
        "--sigma", str(sigma_body_size),       # Sigma for body size
        "--wall_radius_mu", str(wall_radius_mu),  # Wall radius mean
        "--wall_radius_sigma", str(wall_radius_sigma)  # Wall radius sigma        
    ]

    # Execute the command
    src_folder = os.path.dirname(os.path.abspath(__file__))
    try:
        subprocess.run(command, check=True, cwd=src_folder)
        print(f"vacuole_gen.py executed successfully with N={N_spheroids}, mu={mu_body_size}, sigma={sigma_body_size}")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running vacuole_gen.py: {e}")

'''
def run_pipeline(sample_size, mu, sigma):
    """
    Runs the pipeline sample_size times for the given mu and sigma.
    """
    for run_idx in range(1, sample_size + 1):
        print(f"\n--- Pipeline execution {run_idx}/{sample_size} for mu={mu}, sigma={sigma} ---")
        
        # Step 1: Run vacuole_gen with the current mu and sigma
        vacuolegenmain(mu, sigma)
        
        # Step 2: Run the CC3D simulation
        run_cc3d_script()
        
        # Step 3: Run SliceStats (if applicable)
        run_SliceStats()
        
    print(f"--- Completed pipeline for mu={mu}, sigma={sigma} ---")
'''

def option_one():
    """
    Main pipeline logic with the choice to iterate over either Body Number or Body Size.
    """
    params = load_model_parameters("./attributes/Model_Parameters.txt")
    if not params:
        raise ValueError("Failed to load parameters from Model_Parameters.txt.")

    sample_size = int(params["Sample_Size"])
    wall_radius_mu = float(params["Wall_Radius_mu"])
    wall_radius_sigma = float(params["Wall_Radius_sigma"])    

    print("What do you want to iterate over?")
    print("[1] Body Number")
    print("[2] Body Size")
    choice = input("Enter your choice (1/2): ").strip()

    if choice == "1":
        # Iterate over Body Number
        mu_start = float(params["Body_number_starting_mu"])
        mu_end = float(params["Body_number_ending_mu"])
        mu_step = float(params["Body_number_mu_step"])
        sigma_start = float(params["Body_number_starting_sigma"])
        sigma_end = float(params["Body_number_ending_sigma"])
        sigma_step = float(params["Body_number_sigma_step"])

        # Fixed body size parameters
        mu_body_size = float(params["Body_radius_starting_mu"])
        sigma_body_size = float(params["Body_radius_starting_sigma"])

        mu = mu_start
        while mu <= mu_end:
            sigma = sigma_start
            while sigma <= sigma_end:
                print(f"Running pipeline with Body Number: mu={mu}, sigma={sigma}")
                for run_idx in range(sample_size):
                    print(f"Sample {run_idx + 1}/{sample_size} for mu={mu}, sigma={sigma}")
                    vacuolegenmain(
                        mu_body_number=mu,
                        sigma_body_number=sigma,
                        mu_body_size=mu_body_size,
                        sigma_body_size=sigma_body_size,
                        wall_radius_mu=wall_radius_mu,
                        wall_radius_sigma=wall_radius_sigma)
                    run_cc3d_script()
                    #run_SliceStats()
                sigma += sigma_step
            mu += mu_step

    elif choice == "2":
        # Iterate over Body Size
        mu_start = float(params["Body_radius_starting_mu"])
        mu_end = float(params["Body_radius_ending_mu"])
        mu_step = float(params["Body_radius_mu_step"])
        sigma_start = float(params["Body_radius_starting_sigma"])
        sigma_end = float(params["Body_radius_ending_sigma"])
        sigma_step = float(params["Body_radius_sigma_step"])

        # Fixed body number parameters
        mu_body_number = float(params["Body_number_starting_mu"])
        sigma_body_number = float(params["Body_number_starting_sigma"])
        fixed_N = int(np.random.lognormal(mean=mu_body_number, sigma=sigma_body_number))

        mu = mu_start
        while mu <= mu_end:
            sigma = sigma_start
            while sigma <= sigma_end:
                print(f"Running pipeline with Body Size: mu={mu}, sigma={sigma}")
                for run_idx in range(sample_size):
                    print(f"Sample {run_idx + 1}/{sample_size} for mu={mu}, sigma={sigma}")
                    vacuolegenmain(
                        mu_body_number=mu_body_number,
                        sigma_body_number=sigma_body_number,
                        mu_body_size=mu,
                        sigma_body_size=sigma,
                        wall_radius_mu=wall_radius_mu,
                        wall_radius_sigma=wall_radius_sigma,
                        fixed_N=fixed_N)
                    run_cc3d_script()
                    #run_SliceStats()
                sigma += sigma_step
            mu += mu_step

    else:
        print("Invalid choice. Please restart and select either 1 or 2.")

    print("--- Pipeline execution complete ---")

def option_two():
    print("--- Option Two Selected ---")
    # Prompt for N (number of spheroids) required by vacuole_gen.py
    print("Enter the number of spheroids (N) to generate:")
    try:
        N_spheroids = int(input())
        if N_spheroids < 1:
            raise ValueError
    except ValueError:
        print("Invalid input. Please enter a positive integer for N.")
        return

    vacuolegenmain()
    print("--- Option Two Complete ---")

def option_three():
    print("--- Option Three Selected ---")
    run_cc3d_script()
    print("--- Option Three Complete ---")

def option_four():
    print("--- Option Four Selected ---")
    run_SliceStats()
    print("--- Option Four Complete ---")

def option_five():
    print("--- Option Five Selected ---")
    run_AVSStats()
    print("--- Option Five Complete ---")

def option_six():
    print("--- Option Six Selected ---")
    read_readme()
    print("--- Option Six Complete ---")

def vacuole_gen_main(N_spheroids):
    print("Running vacuole_gen to generate initial conditions...")

    # Build the command to execute vacuole_gen.py with the required --N argument
    command = [
        sys.executable,  # Path to the Python interpreter
        'vacuole_gen.py',  # Path to the vacuole_gen.py script
        '--N', str(N_spheroids)
    ]

    # Execute the command in the src folder
    src_folder = os.path.dirname(os.path.abspath(__file__))  # Get the directory of AVS.py (src folder)
    try:
        subprocess.run(command, check=True, cwd=src_folder)
        print("vacuole_gen.py has been executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running vacuole_gen.py: {e}")

def run_cc3d_script():
    print("Running CC3D simulation using runScript.sh...")

    # Get the path to the CompuCell3D folder
    cc3d_folder = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'CompuCell3D')
    run_script = os.path.join(cc3d_folder, 'runScript.sh')

    # Check if runScript.sh exists in the CompuCell3D folder
    if not os.path.exists(run_script):
        print("Error: runScript.sh not found in the CompuCell3D folder.")
        return

    # Make sure runScript.sh is executable
    st = os.stat(run_script)
    os.chmod(run_script, st.st_mode | stat.S_IEXEC)

    # Define the arguments to pass to runScript.sh
    input_file = './cc3dSimulation/CC3D_sim_for_AVS.cc3d'
    frames = '100'
    output_dir = './Output/'
    console_printer = 'infoPrinter'

    # Command to run with arguments
    command = [run_script, '-i', input_file, '-f', frames, '-o', output_dir, '-c', console_printer]

    # Run runScript.sh with the necessary arguments in the CompuCell3D folder
    try:
        start_time = time.time()
        subprocess.run(command, check=True, cwd=cc3d_folder)
        end_time = time.time()
        print(f"CC3D simulation completed in {end_time - start_time:.2f} seconds.")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running runScript.sh: {e}")
        
def run_SliceStats():
    print("--- Option Four Selected ---")
    print("Running SliceStats.py...")
    
    # Build the command to execute SliceStats.py
    command = [
        sys.executable,  # Path to the Python interpreter
        'SliceStats.py',  # Path to the SliceStats.py script
    ]
    
    # Define the folder where SliceStats.py is located
    src_folder = os.path.dirname(os.path.abspath(__file__))
    
    try:
        # Run SliceStats.py in the src folder
        subprocess.run(command, check=True, cwd=src_folder)
        print("SliceStats.py has been executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running SliceStats.py: {e}")
    except FileNotFoundError:
        print("SliceStats.py not found. Please ensure it exists in the src folder.")
        
    print("--- Option Four Complete ---")   
    

def run_AVSStats():
    print("Running AVSStats...")
    command = [
        sys.executable,
        'AVSStats.py'
    ]

    src_folder = os.path.dirname(os.path.abspath(__file__))

    try:
        # Run AVSStats.py in the src folder
        subprocess.run(command, check=True, cwd=src_folder)
        print("AVSStats has been executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running AVSStats: {e}")
    except FileNotFoundError:
        print("AVSStats.py not found. Please ensure it exists in the src folder.")


def read_readme():
    print("Reading the ReadMe file...")
    src_folder = os.path.dirname(os.path.abspath(__file__))
    readme_file = os.path.join(src_folder, 'README.txt')  # Adjust the path if necessary
    if os.path.exists(readme_file):
        with open(readme_file, 'r') as f:
            content = f.read()
            print(content)
    else:
        print("README file not found.")

if __name__ == "__main__":
    main()
