# -*- coding: utf-8 -*-
import subprocess
import time
import os
import sys
import stat
import numpy as np

def main():
    print("Welcome to the Autophagic Vacuole Simulation (AVS) Project")

    while True:
        print(">> Please select from the following options by entering the corresponding number:")
        print("\t[1]: Run pipeline with nested iterations over Body Number and Body Size")
        print("\t[2]: Run CC3D alone in headless mode")
        print("\t[3]: Run Slice Stats alone")
        print("\t[4]: Run AVS Stats alone")
        print("\t[5]: Read the ReadMe file")
        print("\t[0]: Exit AVS")

        scriptChoice = input()

        if scriptChoice == "0":
            print("--- Now exiting AVS ---")
            sys.exit()
        elif scriptChoice == "1":
            run_pipeline()
        elif scriptChoice == "2":
            option_three()
        elif scriptChoice == "3":
            option_four()
        elif scriptChoice == "4":
            option_five()
        elif scriptChoice == "5":
            option_six()
        else:
            print("--- Invalid Input, please select from options 0 to 5 ---")


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


def run_pipeline():
    """
    Run pipeline iterating over Body Number and Body Size using nested loops.
    """
    params = load_model_parameters("./attributes/Model_Parameters.txt")
    if not params:
        raise ValueError("Failed to load parameters from Model_Parameters.txt.")

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

    mu_body_number = mu_body_number_start
    while mu_body_number <= mu_body_number_end:
        sigma_body_number = sigma_body_number_start
        while sigma_body_number <= sigma_body_number_end:
            mu_body_size = mu_body_size_start
            while mu_body_size <= mu_body_size_end:
                sigma_body_size = sigma_body_size_start
                while sigma_body_size <= sigma_body_size_end:
                    N_spheroids = int(np.random.lognormal(mean=mu_body_number, sigma=sigma_body_number))
                    print(f"\nRunning pipeline with N={N_spheroids}, mu_body_number={mu_body_number}, sigma_body_number={sigma_body_number}")
                    print(f"mu_body_size={mu_body_size}, sigma_body_size={sigma_body_size}")

                    for run_idx in range(sample_size):
                        print(f"Sample {run_idx + 1}/{sample_size}")
                        vacuolegenmain(
                            N_spheroids=N_spheroids,
                            mu_body_number=mu_body_number,
                            sigma_body_number=sigma_body_number,
                            mu_body_size=mu_body_size,
                            sigma_body_size=sigma_body_size,
                            wall_radius_mu=wall_radius_mu,
                            wall_radius_sigma=wall_radius_sigma
                        )

                        # Run CompuCell3D simulation
                        run_cc3d_script()

                    sigma_body_size += sigma_body_size_step
                mu_body_size += mu_body_size_step
            sigma_body_number += sigma_body_number_step
        mu_body_number += mu_body_number_step

    print("--- Pipeline execution complete ---")


def vacuolegenmain(N_spheroids, mu_body_number, sigma_body_number, mu_body_size, sigma_body_size, wall_radius_mu, wall_radius_sigma):
    """
    Run vacuole_gen.py with specified parameters.
    """
    print(f"Running vacuole_gen with N={N_spheroids}, mu_body_size={mu_body_size}, sigma_body_size={sigma_body_size}")

    command = [
        sys.executable,
        "vacuole_gen.py",
        "--N", str(N_spheroids),
        "--mu", str(mu_body_size),
        "--sigma", str(sigma_body_size),
        "--wall_radius_mu", str(wall_radius_mu),
        "--wall_radius_sigma", str(wall_radius_sigma)
    ]

    try:
        subprocess.run(command, check=True)
        print(f"vacuole_gen.py executed successfully with N={N_spheroids}, mu={mu_body_size}, sigma={sigma_body_size}")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running vacuole_gen.py: {e}")


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


def run_SliceStats():
    print("Running SliceStats...")
    command = [sys.executable, 'SliceStats.py']

    try:
        subprocess.run(command, check=True)
        print("SliceStats executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while running SliceStats.py: {e}")


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
