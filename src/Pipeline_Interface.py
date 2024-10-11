# -*- coding: utf-8 -*-
import subprocess
import time
import os
import sys

from tkinter import Tk
from tkinter.filedialog import askopenfilename, askdirectory

import vacuole_gen2  # Updated import statement

############################################################################################################
#   Autophagic Vacuole Simulation (AVS) Project
#   Re-engineered AVS.py script
#   Author: [Your Name]
#   Last Date Modified: [Current Date]
#
#   This script serves as the master control program for the AVS project, providing a menu-driven interface
#   to run various components of the simulation pipeline.
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

# paramsFile is used to keep track of several variables used by multiple scripts.
paramsFile = "attributes/Model_Parameters.txt"

def main():
    print("Welcome to the Autophagic Vacuole Simulation (AVS) Project")
    print("Would you like to run this program with File Explorer selection enabled? [y/n]")
    explorerEnabled = False

    optionSelection = input()
    if optionSelection.lower() == 'y':
        explorerEnabled = True

    while True:
        print(">> Please select from the following options by entering the corresponding number:")
        print("\t[1]: Run pipeline based on N input (1 -> N)")
        print("\t[2]: Run vacuole_gen alone (checking outputs)")
        print("\t[3]: Run CC3D alone in headless mode")
        print("\t[4]: Run Slice Stats alone")
        print("\t[5]: Run AVS Stats alone")
        print("\t[6]: Update Model_Parameters.txt")
        print("\t[7]: Read the ReadMe file")
        print("\t[0]: Exit AVS")

        scriptChoice = input()

        if scriptChoice == "0":
            print("--- Now exiting AVS ---")
            sys.exit()
        elif scriptChoice == "1":
            option_one(explorerEnabled)
        elif scriptChoice == "2":
            option_two(explorerEnabled)
        elif scriptChoice == "3":
            option_three(explorerEnabled)
        elif scriptChoice == "4":
            option_four()
        elif scriptChoice == "5":
            option_five()
        elif scriptChoice == "6":
            option_six()
        elif scriptChoice == "7":
            option_seven()
        else:
            print("--- Invalid Input, please select from options 0 to 7 ---")

def option_one(fileSelectOpt):
    print("--- Option One Selected ---")
    print("Enter the number of runs (N) you want to perform (1 -> N):")
    try:
        num_runs = int(input())
        if num_runs < 1:
            raise ValueError
    except ValueError:
        print("Invalid input. Please enter a positive integer.")
        return

    # Ask for a base output directory
    if fileSelectOpt:
        print("Select base output directory for all runs:")
        Tk().withdraw()
        base_output_dir = askdirectory()
    else:
        print("Enter base output directory for all runs:")
        base_output_dir = input()

    # Ensure base output directory exists
    if not os.path.exists(base_output_dir):
        os.makedirs(base_output_dir)

    # Collect CC3D executable path once
    if fileSelectOpt:
        print("Select your CC3D executable (e.g., 'runScript.bat' or 'cc3d.sh'):")
        Tk().withdraw()
        cc3d_executable = askopenfilename()
    else:
        print("Enter the full path to your CC3D executable (e.g., 'runScript.bat' or 'cc3d.sh'):")
        cc3d_executable = input()

    for i in range(1, num_runs+1):
        print(f"\n--- Running pipeline iteration {i} ---")
        iteration_dir = os.path.join(base_output_dir, f"Run_{i}")
        if not os.path.exists(iteration_dir):
            os.makedirs(iteration_dir)

        # Step 1: Run vacuole_gen to generate the initial conditions, prompting for N
        vacuole_gen_main(fileSelectOpt, output_dir=iteration_dir, prompt_for_N=True)

        # Step 2: Run CC3D simulation using the PIF file only
        run_cc3d_pif_only(cc3d_executable, iteration_dir)

    print("--- Option One Complete ---")

def option_two(fileSelectOpt):
    print("--- Option Two Selected ---")
    if fileSelectOpt:
        print("Select output directory for vacuole_gen outputs:")
        Tk().withdraw()
        output_dir = askdirectory()
    else:
        print("Enter output directory for vacuole_gen outputs:")
        output_dir = input()

    vacuole_gen_main(fileSelectOpt, output_dir=output_dir, prompt_for_N=True)
    print("--- Option Two Complete ---")

def option_three(fileSelectOpt):
    print("--- Option Three Selected ---")
    run_cc3d_model(fileSelectOpt)
    print("--- Option Three Complete ---")

def option_four():
    print("--- Option Four Selected ---")
    print("Slice Stats functionality is not implemented yet.")
    print("--- Option Four Complete ---")

def option_five():
    print("--- Option Five Selected ---")
    print("AVS Stats functionality is not implemented yet.")
    print("--- Option Five Complete ---")

def option_six():
    print("--- Option Six Selected ---")
    update_model_parameters()
    print("--- Option Six Complete ---")

def option_seven():
    print("--- Option Seven Selected ---")
    print("ReadMe functionality is not implemented yet.")
    print("--- Option Seven Complete ---")

def vacuole_gen_main(fileSelectOpt, output_dir, prompt_for_N=False):
    print("Running vacuole_gen to generate initial conditions...")
    # Collect parameters from the user or Model_Parameters.txt
    params = read_model_parameters()

    if prompt_for_N:
        print(f"Current number of spheroids (N_SPHEROIDS) is {int(params.get('N_SPHEROIDS', 50))}.")
        print("Do you want to use this value? [y/n]")
        choice = input().lower()
        if choice == 'n':
            print("Enter the number of spheroids (N_SPHEROIDS) to generate:")
            try:
                params['N_SPHEROIDS'] = int(input())
            except ValueError:
                print("Invalid input. Using existing value.")
        else:
            print("Using existing value of N_SPHEROIDS.")

    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Build the command to execute vacuole_gen2.py
    command = [
        sys.executable,  # Path to the Python interpreter
        'vacuole_gen2.py',  # Path to the vacuole_gen2.py script
        '--N', str(int(params.get('N_SPHEROIDS', 50))),
        '--x_max', str(float(params.get('X_MAX', 120.0))),
        '--y_max', str(float(params.get('Y_MAX', 120.0))),
        '--z_max', str(float(params.get('Z_MAX', 50.0))),
        '--min_radius', str(float(params.get('MIN_RADIUS', 3.0))),
        '--max_radius', str(float(params.get('MAX_RADIUS', 8.0))),
        '--wall_outer_radius', str(float(params.get('WALL_OUTER_RADIUS', 40.0))),
        '--wall_thickness', str(float(params.get('WALL_THICKNESS', 2.0))),
        '--dx', str(float(params.get('DX', 1.0))),
        '--max_tries', str(int(params.get('MAX_TRIES', 1000))),
        '--output', os.path.join(output_dir, 'initial_conditions.piff')
    ]

    # Execute the command
    try:
        subprocess.run(command, check=True)
        print(f"vacuole_gen outputs have been saved to {output_dir}")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running vacuole_gen2.py: {e}")

def run_cc3d_pif_only(cc3d_executable, working_dir):
    print("Running CC3D simulation in headless mode using PIF file...")

    # Ensure that CC3D simulation outputs are stored in the specified output directory
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)

    # Run CC3D with the PIF file
    pif_file = os.path.join(working_dir, 'initial_conditions.piff')

    # Create a minimal CC3D project directory
    cc3d_project_dir = os.path.join(working_dir, 'CC3D_Project')
    if not os.path.exists(cc3d_project_dir):
        os.makedirs(cc3d_project_dir)

    # Copy the PIF file into the CC3D project directory
    import shutil
    shutil.copy(pif_file, os.path.join(cc3d_project_dir, 'initial.piff'))

    # Create a minimal model file within the project directory
    minimal_model_file = os.path.join(cc3d_project_dir, 'minimal_model.cc3d')
    with open(minimal_model_file, 'w') as f:
        f.write('''<?xml version="1.0" encoding="UTF-8"?>
<CompuCell3D>
    <Metadata>
        <NumberOfSteps>1000</NumberOfSteps>
    </Metadata>
    <Potts>
        <Dimensions x="200" y="200" z="200"/>
        <Steps>1000</Steps>
        <Temperature>10.0</Temperature>
        <NeighborOrder>1</NeighborOrder>
    </Potts>
    <Plugin Name="CellType">
        <CellType TypeName="Medium" TypeId="0"/>
        <CellType TypeName="Body" TypeId="1"/>
        <CellType TypeName="Wall" TypeId="2"/>
    </Plugin>
    <Steppable Type="PIFInitializer">
        <PIFName>initial.piff</PIFName>
    </Steppable>
</CompuCell3D>
''')

    # Run CC3D simulation
    start_time = time.time()
    subprocess.run([cc3d_executable, "--exitWhenDone", "-i", minimal_model_file], cwd=cc3d_project_dir)
    end_time = time.time()
    print(f"CC3D simulation completed in {end_time - start_time:.2f} seconds.")

    # Move CC3D outputs back to the working directory
    cc3d_output_dir = os.path.join(cc3d_project_dir, 'Simulation')
    if os.path.exists(cc3d_output_dir):
        dest_output_dir = os.path.join(working_dir, 'CC3D_Outputs')
        shutil.move(cc3d_output_dir, dest_output_dir)
        print(f"CC3D outputs have been saved to {dest_output_dir}")
    else:
        print("CC3D output directory not found.")

def run_cc3d_model(fileSelectOpt):
    print("Running CC3D simulation in headless mode using model.cc3d...")
    if fileSelectOpt:
        print("Select your CC3D executable (e.g., 'runScript.bat' or 'cc3d.sh'):")
        Tk().withdraw()
        cc3d_executable = askopenfilename()
        print("Select your 'model.cc3d' file:")
        Tk().withdraw()
        cc3d_model_file = askopenfilename()
    else:
        print("Enter the full path to your CC3D executable (e.g., 'runScript.bat' or 'cc3d.sh'):")
        cc3d_executable = input()
        print("Enter the full path to your 'model.cc3d' file:")
        cc3d_model_file = input()

    # Ask for output directory
    if fileSelectOpt:
        print("Select output directory for CC3D outputs:")
        Tk().withdraw()
        output_dir = askdirectory()
    else:
        print("Enter output directory for CC3D outputs:")
        output_dir = input()

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Run CC3D simulation
    start_time = time.time()
    subprocess.run([cc3d_executable, "--exitWhenDone", "-i", cc3d_model_file], cwd=output_dir)
    end_time = time.time()
    print(f"CC3D simulation completed in {end_time - start_time:.2f} seconds.")
    print(f"CC3D outputs have been saved to {output_dir}")

def update_model_parameters():
    print("Updating Model_Parameters.txt with new parameters.")
    params = {}
    print("\n>> Enter the number of spheroids (N_SPHEROIDS):")
    params['N_SPHEROIDS'] = float(input())
    print(">> Enter X_MAX:")
    params['X_MAX'] = float(input())
    print(">> Enter Y_MAX:")
    params['Y_MAX'] = float(input())
    print(">> Enter Z_MAX:")
    params['Z_MAX'] = float(input())
    print(">> Enter MIN_RADIUS:")
    params['MIN_RADIUS'] = float(input())
    print(">> Enter MAX_RADIUS:")
    params['MAX_RADIUS'] = float(input())
    print(">> Enter WALL_OUTER_RADIUS:")
    params['WALL_OUTER_RADIUS'] = float(input())
    print(">> Enter WALL_THICKNESS:")
    params['WALL_THICKNESS'] = float(input())
    print(">> Enter DX (grid resolution):")
    params['DX'] = float(input())
    print(">> Enter MAX_TRIES:")
    params['MAX_TRIES'] = int(input())

    # Write parameters to Model_Parameters.txt
    with open(paramsFile, 'w') as f:
        for key, value in params.items():
            f.write(f"{key}:{value}\n")
    print("Model_Parameters.txt has been updated.")

def read_model_parameters():
    params = {}
    try:
        with open(paramsFile, 'r') as f:
            for line in f:
                if ':' in line:
                    key, value = line.strip().split(':', 1)
                    params[key.strip()] = float(value.strip())
    except FileNotFoundError:
        print(f"{paramsFile} not found. Using default parameters.")
    return params

if __name__ == "__main__":
    main()
