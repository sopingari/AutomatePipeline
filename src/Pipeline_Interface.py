# -*- coding: utf-8 -*-
import subprocess
import time
import os
import sys
import stat

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

def option_one():
    print("--- Option One Selected ---")
    print("Enter the number of runs (N) you want to perform (1 -> N):")
    try:
        num_runs = int(input())
        if num_runs < 1:
            raise ValueError
    except ValueError:
        print("Invalid input. Please enter a positive integer.")
        return

    # Prompt for N (number of spheroids) required by vacuole_gen.py
    print("Enter the number of spheroids (N) to generate:")
    try:
        N_spheroids = int(input())
        if N_spheroids < 1:
            raise ValueError
    except ValueError:
        print("Invalid input. Please enter a positive integer for N.")
        return
    use_model_parameters = ask_use_model_parameters()

    for i in range(1, num_runs+1):
        print(f"\n--- Running pipeline iteration {i} ---")

        # Step 1: Run vacuole_gen to generate the initial conditions
        vacuole_gen_main(N_spheroids, use_model_parameters)

        # Step 2: Run CC3D simulation
        run_cc3d_script()

    print("--- Option One Complete ---")

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
    # Ask if the user wants to use model_parameters.txt
    use_model_parameters = ask_use_model_parameters()

    vacuole_gen_main(N_spheroids, use_model_parameters)
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

def ask_use_model_parameters():
    while True:
        print("Do you want to use parameters from model_parameters.txt? (y/n):")
        choice = input().strip().lower()
        if choice in ['y', 'yes']:
            return True
        elif choice in ['n', 'no']:
            return False
        else:
            print("Invalid input. Please enter 'y' or 'n'.")
            
def read_parameters(filename):
    parameters = {}
    try:
        with open(filename, 'r') as f:
            for line in f:
                # Remove comments and whitespace
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                if '#' in line:
                # Remove inline comments
                    line = line.split('#', 1)[0].strip()
                # Split parameter name and value
                if '=' in line:
                    name, value = line.split('=', 1)
                    name = name.strip()
                    value = value.strip()
                    # Store as string; we can handle types later if needed
                    parameters[name] = value
    except FileNotFoundError:
        print(f"Parameter file {filename} not found. Proceeding with default parameters.")
    return parameters

def map_parameters(parameters):
    # Mapping from text file parameter names to vacuole_gen.py's expected command-line argument names
    parameter_mapping = {
        'Wall_Radius_mu': 'wall_Radius_Mu',
        'Wall_Radius_sigma': 'wall_Radius_Sigma',
        'Body_radius_mu': 'mu',
        'Body_radius_sigma': 'sigma',
        'Body_number': 'bodies',
        # Add more mappings if necessary
    }

    mapped_parameters = {}
    for param_name, value in parameters.items():
        if param_name in parameter_mapping:
            method_param = parameter_mapping[param_name]
            mapped_parameters[method_param] = value
    return mapped_parameters


def vacuole_gen_main(N_spheroids, use_model_parameters):
    print("\nRunning vacuole_gen to generate initial conditions...")

    # Build the base command to execute vacuole_gen.py
    command = [
        sys.executable,          # Path to the Python interpreter
        'vacuole_gen.py'         # Path to the vacuole_gen.py script
    ]

    if use_model_parameters:
       # Path to model_parameters.txt in attributes folder
        model_params_path = os.path.join('attributes', 'Model_Parameters.txt')

        # Read parameters from model_parameters.txt
        parameters = read_parameters(model_params_path)
        method_parameters = map_parameters(parameters)

        if method_parameters:
            # Add mapped parameters as command-line arguments
            for key, value in method_parameters.items():
                command.extend([f'--{key}', str(value)])
            print("Using parameters from model_parameters.txt:")
            for key, value in method_parameters.items():
                print(f"  {key} = {value}")
        else:
            print("No matching parameters found in model_parameters.txt. Running with default parameters.")
    else:
        # Add the '--N' argument to the command
        command.extend(['--N', str(N_spheroids)])
        print(f"Using --N argument with value: {N_spheroids}")

    # Execute the command in the src folder
    src_folder = os.path.dirname(os.path.abspath(__file__))  # Get the directory of pipeline_interface.py
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
