import numpy as np
import math
import argparse
import numpy as np
import math
import argparse
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import logging
from datetime import datetime
import os
import shutil


#Set up logging
def setup_logging(run_folder):
    log_file = os.path.join(run_folder, f'statistics.log')
    logging.basicConfig(filename=log_file, level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')

def generate_random_direction():
    """
    Generates a random unit vector direction within the positive octant.
    
    Returns:
        tuple: (x, y, z) components of the unit vector, all non-negative.
    """
    phi = random.uniform(0, math.pi / 2)  # Azimuthal angle between 0 and π/2
    costheta = random.uniform(0, 1)       # Cos(theta) between 0 and 1
    theta = math.acos(costheta)           # Polar angle between 0 and π/2
    
    x = math.sin(theta) * math.cos(phi)
    y = math.sin(theta) * math.sin(phi)
    z = math.cos(theta)
    
    return (x, y, z)


def generate_random_point_adjacent(existing_spheroid, new_radius, buffer=0.01):
    """
    Generates a random point adjacent to an existing spheroid.
    
    Parameters:
        existing_spheroid (dict): The spheroid to place the new spheroid next to.
        new_radius (float): Radius of the new spheroid.
        buffer (float): Additional space to ensure no overlap.
    
    Returns:
        tuple: (x, y, z) coordinates for the new spheroid.
    """
    direction = generate_random_direction()
    distance = existing_spheroid['radius'] + new_radius + buffer
    x = existing_spheroid['x'] + direction[0] * distance
    y = existing_spheroid['y'] + direction[1] * distance
    z = existing_spheroid['z'] + direction[2] * distance
    return (x, y, z)

def check_overlap(new_spheroid, spheroids):
    """
    Checks if the new spheroid overlaps with any existing spheroids.
    
    Parameters:
        new_spheroid (dict): The spheroid to be placed.
        spheroids (list): List of existing spheroids.
    
    Returns:
        bool: True if overlaps with any existing spheroid, False otherwise.
    """
    for spheroid in spheroids:
        dx = new_spheroid['x'] - spheroid['x']
        dy = new_spheroid['y'] - spheroid['y']
        dz = new_spheroid['z'] - spheroid['z']
        distance = math.sqrt(dx**2 + dy**2 + dz**2)
        if distance < (new_spheroid['radius'] + spheroid['radius']):
            return True
    return False

def is_within_wall(new_spheroid, wall_center, inner_radius):
    """
    Checks if the new spheroid is entirely within the inner radius of the wall.
    
    Parameters:
        new_spheroid (dict): The spheroid to be placed.
        wall_center (tuple): The (x, y, z) coordinates of the wall's center.
        inner_radius (float): The inner radius of the wall.
    
    Returns:
        bool: True if entirely within the wall, False otherwise.
    """
    dx = new_spheroid['x'] - wall_center[0]
    dy = new_spheroid['y'] - wall_center[1]
    dz = new_spheroid['z'] - wall_center[2]
    distance = math.sqrt(dx**2 + dy**2 + dz**2)
    return (distance + new_spheroid['radius']) <= inner_radius

def generate_spheroids(N, center, inner_radius, min_radius, max_radius, max_tries=1000):
    """
    Generates N spheroids randomly placed within a central cluster inside the wall.
    
    Parameters:
        N (int): Number of spheroids to generate.
        center (tuple): The (x, y, z) coordinates of the central cluster.
        inner_radius (float): The inner radius of the wall.
        min_radius (float): Minimum radius of spheroids.
        max_radius (float): Maximum radius of spheroids.
        max_tries (int): Maximum attempts to place a spheroid without overlap.
    
    Returns:
        list: List of spheroid dictionaries with 'x', 'y', 'z', and 'radius'.
    """
    spheroids = []
    
    # Generate list of spheroid radii and sort them descendingly for better packing
    spheroid_radii = [random.uniform(min_radius, max_radius) for _ in range(N)]
    spheroid_radii.sort(reverse=True)
    
    for i, radius in enumerate(spheroid_radii):
        placed = False
        for attempt in range(max_tries):
            if not spheroids:
                # Place the first spheroid at the center
                x, y, z = center
            else:
                # Select a random existing spheroid to place adjacent to
                existing = random.choice(spheroids)
                x, y, z = generate_random_point_adjacent(existing, radius)
            
            new_spheroid = {'x': x, 'y': y, 'z': z, 'radius': radius}
            
            if not check_overlap(new_spheroid, spheroids):
                if is_within_wall(new_spheroid, center, inner_radius):
                    spheroids.append(new_spheroid)
                    placed = True
                    break
        if not placed:
            print(f"Warning: Could not place spheroid {i+1} after {max_tries} attempts.")

    logging.info(f"Generated {len(spheroids)} spheroids out of {N} requested.")
    return spheroids

def generate_piff_file(spheroids, wall_center, wall_outer_radius, wall_thickness, dx=1.0, filename='output.piff'):
    """
    Generates a PIFF file with spheroids and surrounding wall.
    
    Parameters:
        spheroids (list): List of spheroid dictionaries.
        wall_center (tuple): The (x, y, z) coordinates of the wall's center.
        wall_outer_radius (float): The outer radius of the wall.
        wall_thickness (float): The thickness of the wall.
        dx (float): The resolution for grid boxes.
        filename (str): The name of the output PIFF file.
    """
    piff_lines = []
    cell_id = 1  # Start CellID from 1
    
    # Generate boxes for spheroids (Body cells)
    for spheroid in spheroids:
        x0, y0, z0 = spheroid['x'], spheroid['y'], spheroid['z']
        R = spheroid['radius']
        # Define bounding box
        x_min = x0 - R
        x_max = x0 + R
        y_min = y0 - R
        y_max = y0 + R
        z_min = z0 - R
        z_max = z0 + R
        
        # Generate grid within bounding box
        x_vals = np.arange(x_min, x_max, dx)
        y_vals = np.arange(y_min, y_max, dx)
        z_vals = np.arange(z_min, z_max, dx)
        
        for x in x_vals:
            for y in y_vals:
                for z in z_vals:
                    # Calculate center of the voxel
                    voxel_center = (x + dx / 2, y + dx / 2, z + dx / 2)
                    distance_sq = ((voxel_center[0] - x0) ** 2 +
                                  (voxel_center[1] - y0) ** 2 +
                                  (voxel_center[2] - z0) ** 2)
                    if distance_sq <= R ** 2:
                        # line = f"{cell_id} Body {x:.2f} {x+dx:.2f} {y:.2f} {y+dx:.2f} {z:.2f} {z+dx:.2f}"
                        line = f"{cell_id} Body {int(x)} {int(x+dx)} {int(y)} {int(y+dx)} {int(z)} {int(z+dx)}"
                        piff_lines.append(line)
        cell_id += 1  # Increment CellID for the next spheroid
    
    # Generate boxes for the hollow sphere (Wall)
    wall_cell_id = cell_id  # Assign a unique CellID for the wall
    wall_inner_radius = wall_outer_radius - wall_thickness
    
    # Define bounding box for the wall
    x_min = wall_center[0] - wall_outer_radius
    x_max = wall_center[0] + wall_outer_radius
    y_min = wall_center[1] - wall_outer_radius
    y_max = wall_center[1] + wall_outer_radius
    z_min = wall_center[2] - wall_outer_radius
    z_max = wall_center[2] + wall_outer_radius
    
    # Generate grid within bounding box
    x_vals = np.arange(x_min, x_max, dx)
    y_vals = np.arange(y_min, y_max, dx)
    z_vals = np.arange(z_min, z_max, dx)
    
    for x in x_vals:
        for y in y_vals:
            for z in z_vals:
                # Calculate center of the voxel
                voxel_center = (x + dx / 2, y + dx / 2, z + dx / 2)
                distance_sq = ((voxel_center[0] - wall_center[0]) ** 2 +
                              (voxel_center[1] - wall_center[1]) ** 2 +
                              (voxel_center[2] - wall_center[2]) ** 2)
                if wall_inner_radius ** 2 <= distance_sq <= wall_outer_radius ** 2:
                    # line = f"{wall_cell_id} Wall {x:.2f} {x+dx:.2f} {y:.2f} {y+dx:.2f} {z:.2f} {z+dx:.2f}"
                    line = f"{wall_cell_id} Wall {int(x)} {int(x+dx)} {int(y)} {int(y+dx)} {int(z)} {int(z+dx)}"
                    piff_lines.append(line)
    
    # Write the PIFF content to a file
    with open(filename, 'w') as f:
        for line in piff_lines:
            f.write(line + '\n')
    
    print(f"PIFF file '{filename}' generated with {len(spheroids)} spheroids and surrounding wall.")
    logging.info(f"PIFF file '{filename}' generated with {len(spheroids)} spheroids and surrounding wall.")

def visualize_spheroids_and_wall(spheroids, wall_center, wall_outer_radius, wall_thickness, x_max, y_max, z_max):
    """
    Visualizes the spheroids and the surrounding wall in 3D.
    
    Parameters:
        spheroids (list): List of spheroid dictionaries.
        wall_center (tuple): The (x, y, z) coordinates of the wall's center.
        wall_outer_radius (float): The outer radius of the wall.
        wall_thickness (float): The thickness of the wall.
        x_max (float): Maximum x-axis limit for visualization.
        y_max (float): Maximum y-axis limit for visualization.
        z_max (float): Maximum z-axis limit for visualization.
    """
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot the hollow wall as two wireframes (outer and inner spheres)
    plot_hollow_sphere(ax, wall_center, wall_outer_radius, wall_thickness, color='lightblue', alpha=0.2)
    
    # Plot the spheroids
    for spheroid in spheroids:
        plot_sphere(ax, spheroid['x'], spheroid['y'], spheroid['z'],
                    spheroid['radius'], color='red', alpha=0.6)
    
    # Set plot limits
    ax.set_xlim(0, x_max)
    ax.set_ylim(0, y_max)
    ax.set_zlim(0, z_max)
    
    # Labels
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
    # Equal aspect ratio
    ax.set_box_aspect([x_max, y_max, z_max])
    
    plt.title('Randomly Placed and Clustered Spheroids within Surrounding Wall')
    plt.show()

def plot_sphere(ax, x_center, y_center, z_center, radius, color='r', alpha=0.6):
    """
    Plots a single sphere.
    
    Parameters:
        ax (Axes3D): The 3D axes to plot on.
        x_center (float): X-coordinate of the sphere's center.
        y_center (float): Y-coordinate of the sphere's center.
        z_center (float): Z-coordinate of the sphere's center.
        radius (float): Radius of the sphere.
        color (str): Color of the sphere.
        alpha (float): Transparency of the sphere.
    """
    u = np.linspace(0, 2 * np.pi, 30)
    v = np.linspace(0, np.pi, 30)
    x = x_center + radius * np.outer(np.cos(u), np.sin(v))
    y = y_center + radius * np.outer(np.sin(u), np.sin(v))
    z = z_center + radius * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, color=color, alpha=alpha, linewidth=0, shade=True)

def plot_hollow_sphere(ax, center, outer_radius, wall_thickness, color='lightblue', alpha=0.2):
    """
    Plots a hollow sphere (wall) using wireframes for outer and inner boundaries.
    
    Parameters:
        ax (Axes3D): The 3D axes to plot on.
        center (tuple): The (x, y, z) coordinates of the wall's center.
        outer_radius (float): The outer radius of the wall.
        wall_thickness (float): The thickness of the wall.
        color (str): Color of the wall wireframes.
        alpha (float): Transparency of the wireframes.
    """
    inner_radius = outer_radius - wall_thickness
    u = np.linspace(0, 2 * np.pi, 30)
    v = np.linspace(0, np.pi, 30)
    
    # Outer sphere
    x_outer = center[0] + outer_radius * np.outer(np.cos(u), np.sin(v))
    y_outer = center[1] + outer_radius * np.outer(np.sin(u), np.sin(v))
    z_outer = center[2] + outer_radius * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_wireframe(x_outer, y_outer, z_outer, color=color, alpha=alpha, linewidth=0.5)
    
    # Inner sphere
    x_inner = center[0] + inner_radius * np.outer(np.cos(u), np.sin(v))
    y_inner = center[1] + inner_radius * np.outer(np.sin(u), np.sin(v))
    z_inner = center[2] + inner_radius * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_wireframe(x_inner, y_inner, z_inner, color=color, alpha=alpha, linewidth=0.5)


def main(args):
    # Generate a unique run ID based on the current date and time
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Create a new folder for this run
    runs_dir = 'runs'
    if not os.path.exists(runs_dir):
        os.makedirs(runs_dir)
    run_folder = os.path.join(runs_dir, run_id)
    os.makedirs(run_folder)

    setup_logging(run_folder)

    logging.info(f"Starting new run with ID: {run_id}")
    logging.info(f"Arguments: {args}")

    # Extract grid dimensions from arguments
    x_max = args.x_max
    y_max = args.y_max
    z_max = args.z_max
    N_SPHEROIDS = args.N
    MIN_RADIUS = args.min_radius
    MAX_RADIUS = args.max_radius
    WALL_OUTER_RADIUS = args.wall_outer_radius
    WALL_THICKNESS = args.wall_thickness
    DX = args.dx
    MAX_TRIES = args.max_tries
    filename = args.output

    # Ensure z_max is sufficient to prevent negative z coordinates
    required_z_max = 2 * (WALL_OUTER_RADIUS - WALL_THICKNESS + MAX_RADIUS)
    if z_max < required_z_max:
        print(f"Warning: z_max ({z_max}) is less than required ({required_z_max}). Setting z_max to {required_z_max}.")
        z_max = required_z_max
    
    # Define wall center based on grid dimensions (ensure it's within positive octant)
    wall_center = (x_max / 2, y_max / 2, z_max / 2)
    
    # Calculate inner radius of the wall
    wall_inner_radius = WALL_OUTER_RADIUS - WALL_THICKNESS
    
    # Generate spheroids
    spheroids = generate_spheroids(N_SPHEROIDS, wall_center, wall_inner_radius, MIN_RADIUS, MAX_RADIUS, max_tries=MAX_TRIES)
    
    if len(spheroids) < N_SPHEROIDS:
        print(f"Only placed {len(spheroids)} out of {N_SPHEROIDS} spheroids.")
    
    
    #we will do a function call here that writes the spheroids to a csv file 
    #spheroids, wall_center, WALL_OUTER_RADIUS, WALL_THICKNESS, 
    #compute and write the apb volume,apb st dev,apb mean, apb radius 
    #every spheroid gets its own row including the vacuole  



    # Generate PIFF file
    generate_piff_file(spheroids, wall_center, WALL_OUTER_RADIUS, WALL_THICKNESS, dx=DX, filename=filename)

    # Save a copy of the PIFF file in the run folder
    shutil.copy(filename, os.path.join(run_folder, filename))
    logging.info(f"Saved copy of PIFF file in run folder: {run_folder}")

    logging.info(f"Run {run_id} completed successfully.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Generate PIFF file with N randomly clustered spheroids inside a surrounding hollow Wall.')
    parser.add_argument('--N', type=int, required=True, help='Number of internal spheroids to generate')
    parser.add_argument('--x_max', type=float, default=120.0, help='Maximum x dimension of the space (default: 120)')
    parser.add_argument('--y_max', type=float, default=120.0, help='Maximum y dimension of the space (default: 120)')
    parser.add_argument('--z_max', type=float, default=50.0, help='Maximum z dimension of the space (default: 50)')
    parser.add_argument('--min_radius', type=float, default=3.0, help='Minimum radius for spheroids (default: 3)')
    parser.add_argument('--max_radius', type=float, default=8.0, help='Maximum radius for spheroids (default: 8)')
    parser.add_argument('--wall_outer_radius', type=float, default=40.0, help='Outer radius of the wall (default: 40)')
    parser.add_argument('--wall_thickness', type=float, default=2.0, help='Thickness of the wall (default: 2)')
    parser.add_argument('--dx', type=float, default=1.0, help='Resolution for grid boxes (default: 1.0)')
    parser.add_argument('--max_tries', type=int, default=1000, help='Maximum attempts to place each spheroid (default: 1000)')
    parser.add_argument('--output', type=str, default='output.piff', help='Output PIFF file name (default: output.piff)')
    
    args = parser.parse_args()
    main(args)
