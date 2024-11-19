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
import json
import pandas as pd
from scipy.stats import ortho_group
from scipy.optimize import minimize
from scipy.spatial import distance_matrix

def setup_logging(run_folder, seed):
    # Create log file path
    log_file = os.path.join(run_folder, f'run.log')
    
    # Configure logging
    logging.basicConfig(
        filename=log_file,
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    # Create statistics file path
    stats_file = os.path.join(run_folder, f'statistics.json')
    
    # Log the seed
    logging.info(f"Random seed: {seed}")
    
    return stats_file

def cartesian_product(a,b,c=[0]):
    return np.array(np.meshgrid(a,b,c)).T.reshape(-1,3)

def total_dist_to_pt(pos_array_as_vec, apb_df, ctr_to_use:np.ndarray=None):
    pos_array = np.reshape(pos_array_as_vec,(len(apb_df),-1))
    if ctr_to_use is None:
        ctr_to_use=np.zeros_like(pos_array[1,:])
    dists = np.sqrt(np.sum((pos_array-ctr_to_use)**2,axis=1))
    return np.sum(dists)

def move_along_dir(df, pos_array, ii, direction, stepsize, placeVacuole=False):
    dist_from_origin = (df.loc[0,"r"]+df.loc[ii,"r"])*(1.0+1e-8)
    pvals = df['p']
    pvals = pvals[0:ii]
    pvals_e = np.expand_dims(pvals,axis=1)
    apvals = (pvals[0:ii] + df.loc[ii,"p"])/2
    apvals_e = np.expand_dims(apvals,axis=1)
    
    stopflag = False
    distmins = df.head(ii)["r"]+df.loc[ii,"r"]
    
    while not stopflag:
        pos = direction*dist_from_origin
        dists = (np.sum(np.abs(pos_array[0:ii,:]-pos)**apvals_e,axis=1))**(1.0/apvals)
        if placeVacuole:
            close_enough = (dists + df.head(ii)["r"]) <= vacRad
            stopflag = not all(close_enough)
        else:
            far_enough = dists >= distmins
            stopflag = all(far_enough)       
        
        dist_from_origin += stepsize
    
    if placeVacuole:
        dist_from_origin -= 2*stepsize
        pos = direction*dist_from_origin
    else:
        dist_from_origin -= stepsize
    return (pos,dist_from_origin)

def inter_APB_boundary_distances(pos_array_as_vec, apb_df):
    pos_array = np.reshape(pos_array_as_vec,(len(apb_df),-1))
    cur_dist_mat = distance_matrix(pos_array,pos_array)
    rvec = apb_df['r'].to_numpy()
    min_dist_mat = rvec[:,None]+rvec
    d = cur_dist_mat - min_dist_mat
    dists_as_vec = d[np.triu_indices(d.shape[0],k=1)]
    return dists_as_vec

def get_directions_hcp(bodies, ndim=3, rng=None, unitvectors=True, exclude_origin=True, random_rotate=True):
    side_length = int(np.ceil(bodies**(1.0/ndim)))
    hsl = int(np.ceil(side_length/2))
    hsl = hsl + (hsl+1)%2
    
    ii_list = range(-hsl,hsl+1)
    jj_list = ii_list
    if ndim==2:
        kk_list = [0]
    else:
        kk_list = ii_list
    
    iijjkk = cartesian_product(ii_list,jj_list,kk_list)
    if exclude_origin:
        keepthese = np.sum(iijjkk**2,axis=1) > 0
        iijjkk = iijjkk[keepthese,:]
    coords = np.zeros_like(iijjkk,dtype=float)
    
    coords[:,0] = 2*iijjkk[:,0] + (iijjkk[:,1]+iijjkk[:,2])%2
    coords[:,1] = np.sqrt(3)*(iijjkk[:,1]+(1/3.0)*(iijjkk[:,2] % 2))
    if ndim==2:
        coords = np.delete(coords,2,axis=1)
    else:
        coords[:,2] = 2*np.sqrt(6)/3*iijjkk[:,2]
    
    if random_rotate:
        scipy_randomGen = ortho_group
        scipy_randomGen.random_state=rng
        dirmat = ortho_group.rvs(ndim)
        coords_rotated = dirmat @ coords.transpose()
        coords = coords_rotated.transpose()
    
    dists_true = np.sqrt(np.sum(coords**2,axis=1))
    dists_perturbed = dists_true + rng.uniform(low=0,high=0.001)
    dists_order = np.argsort(dists_perturbed)
    dir_array = coords[dists_order,:]
    orig_lengths = np.expand_dims(np.sqrt(np.sum(dir_array**2,axis=1)),axis=1)
    
    if unitvectors:
        dir_array = dir_array / orig_lengths
    
    return dir_array, orig_lengths

def get_directions_ortho_random(bodies, ndim=3, rng=None, unitvectors=True):
    dir_array = np.zeros((0,ndim))
    for ii in range(bodies):
        scipy_randomGen = ortho_group
        scipy_randomGen.random_state=rng
        dirmat = scipy_randomGen.rvs(ndim)
        dir_array = np.append(dir_array,dirmat,axis=0)
    
    orig_lengths = np.expand_dims(np.sqrt(np.sum(dir_array**2,axis=1)),axis=1)
    return dir_array, orig_lengths

def genBalls3(bodies=20, wall_Radius_Mu=6.8, wall_Radius_Sigma=0.34, mu=5, sigma=0.2, iterations=4,
              ndim=3, rng=None, method='hcp', plist:list=None, pcterrcap=10.0, 
              posOctant=True, optimmaxiter=100, maxVacuoleIterations=100):
    
    r_normals = rng.standard_normal(bodies)*sigma+mu
    r = np.exp(r_normals)
    
    if plist is None:
        plist = [2.0]*bodies
    elif len(plist)==1:
        plist = plist*bodies
    
    stepsize = (pcterrcap/100.0)*np.min(r)
    
    if method=='hcp':
        dirmat,orig_lengths = get_directions_hcp(bodies,ndim,rng,unitvectors=True)
    elif method=='ortho_random':
        dirmat,orig_lengths = get_directions_ortho_random(bodies,ndim,rng)
    else:
        print(f'unknown method:{method}')
        return (None,None,None)
    
    dirmat_safe = dirmat.copy()
    
    ii = 0
    d = {'bodynum': ii, 'r': r[ii],'p':plist[ii]}
    df = pd.DataFrame(data=[d])
    
    pos_array = np.zeros(shape=(bodies,3))
    
    for ii0 in range(bodies-1):
        ii = ii0+1
        d = {'bodynum': ii, 'r': r[ii],'p':plist[ii]}
        df = pd.concat([df,pd.DataFrame(data=[d])],ignore_index=True)
        
        ndir_to_try = min(iterations,dirmat.shape[0])
        pos_list = []
        dist_hist = np.inf*np.ones((ndir_to_try,))
        
        for jj in range(ndir_to_try):
            dir=np.zeros((3,),dtype=float)
            dir[0:ndim]=dirmat[jj,:]
            pos,dist_from_origin = move_along_dir(df,pos_array,ii,dir,stepsize)
            pos_list.append(pos)
            dist_hist[jj] = dist_from_origin
        
        jjbest = np.argmin(dist_hist)
        pos = pos_list[jjbest]
        pos_array[ii,:] = pos
        dirmat = np.delete(dirmat,jjbest,axis=0)
        orig_lengths = np.delete(orig_lengths,jjbest,axis=0)
    
    if optimmaxiter > 0:
        optim_method='trust-constr'
        myoptions={'maxiter':optimmaxiter}
        optim_df = df
        x0 = np.ravel(pos_array)
        cons = {'type': 'ineq', 'args': (optim_df,), 'fun': inter_APB_boundary_distances}
        
        ofv_original = total_dist_to_pt(x0,optim_df)
        print("ofv_original:",ofv_original)
        
        res = minimize(total_dist_to_pt, x0, args=optim_df,method=optim_method,constraints=cons,options=myoptions)
        
        pos_array_safe = pos_array.copy()
        pos_array_as_vec = res['x']
        pos_array = np.reshape(pos_array_as_vec,(len(optim_df),-1))
        ofv_final = total_dist_to_pt(pos_array_as_vec,optim_df)
        print("ofv_final:",ofv_final)
    
    df['bodyType'] = "APB"
    
    # Calculate means before shifting
    pos_array_means = np.mean(pos_array, axis=0)
    
    if posOctant:
        pos_array_shifted_to_origin = pos_array - pos_array_means
        
        vacp = 2.0
        dists_to_origin = (np.sum(np.abs(pos_array_shifted_to_origin)**vacp,axis=1))**(1.0/vacp)
        dists_plus_rad = dists_to_origin + r
        vacRad = 0
        iterCount = 0
        
        while any(dists_plus_rad > vacRad) and (iterCount < maxVacuoleIterations):
            r_normals = rng.standard_normal(1)*wall_Radius_Sigma+wall_Radius_Mu
            vacRad = np.exp(r_normals)
            iterCount += 1
            
        if(iterCount >= maxVacuoleIterations):
            print("Warning: Maxed Out on Vacuole Size Iterations.")
            print(f"Iterations:{iterCount}")
            print(f"{dists_plus_rad}")
            print(f"{wall_Radius_Sigma}")
            print(f"{wall_Radius_Mu}")
            vacRad = -9999
        
        pos_array_shifted_to_1st_octant = pos_array_shifted_to_origin + vacRad
        r_and_pos_array = np.hstack((np.expand_dims(r,axis=1),pos_array_shifted_to_1st_octant))
        
        vac_pos_array = np.ones_like(pos_array[0,:])*vacRad
        vac_r_and_pos = np.hstack((vacRad,vac_pos_array))
        r_and_pos_array_w_vac = np.vstack((vac_r_and_pos,r_and_pos_array))
        pos_array = np.delete(r_and_pos_array_w_vac,0,axis=1)
        
        d = {'bodynum': -100, 'bodyType':'Vacuole', 'r': vacRad,'p':vacp}
        df_just_vac = pd.DataFrame(data=[d])
        df = pd.concat([df_just_vac,df],ignore_index=True)
    
    dfpos = pd.DataFrame.from_records(pos_array)
    dfpos.rename(columns={0: 'x', 1: 'y', 2: 'z'}, inplace=True)
    df = pd.concat([df.reset_index(drop=True), dfpos], axis=1)
    
    return(df,pos_array,r_and_pos_array,dirmat_safe)

def calculate_spheroid_metrics(spheroid, spheroids_df, vacuole, max_radius):
    """
    Calculate additional metrics for a single spheroid.
    """
    volume = (4/3) * np.pi * (spheroid.r**3)
    
    distance_from_center = np.sqrt(
        (spheroid.x - vacuole['x'])**2 + 
        (spheroid.y - vacuole['y'])**2 + 
        (spheroid.z - vacuole['z'])**2
    )
    
    neighbors = []
    for other in spheroids_df.itertuples():
        if other.Index == spheroid.Index:
            continue
        
        dx = spheroid.x - other.x
        dy = spheroid.y - other.y
        dz = spheroid.z - other.z
        center_distance = np.sqrt(dx**2 + dy**2 + dz**2)
        surface_distance = center_distance - (spheroid.r + other.r)
        
        neighbors.append({
            'center_distance': center_distance,
            'surface_distance': surface_distance
        })
    
    neighbors.sort(key=lambda x: x['center_distance'])
    
    local_radius = 2.5 * max_radius
    local_neighbors = [n for n in neighbors if n['center_distance'] <= local_radius]
    local_density = len(local_neighbors) / ((4/3) * np.pi * local_radius**3)
    
    distance_to_wall = abs(distance_from_center - vacuole['r'])
    
    return {
        'volume': volume,
        'distance_from_center': distance_from_center,
        'nearest_neighbor_center_distance': neighbors[0]['center_distance'] if neighbors else None,
        'nearest_neighbor_surface_distance': neighbors[0]['surface_distance'] if neighbors else None,
        'num_neighbors': len(local_neighbors),
        'local_density': local_density,
        'distance_to_wall': distance_to_wall
    }


def write_combined_csv(run_folder, run_id, args, df):
    """
    Write a single CSV file containing both summary and detailed information.
    """
    import csv
    from datetime import datetime
    
    # Separate vacuole and spheroids
    vacuole = df[df['bodyType'] == 'Vacuole'].iloc[0]
    spheroids = df[df['bodyType'] == 'APB']
    
    # Calculate global statistics
    total_spheroid_volume = sum((4/3) * np.pi * (s['r']**3) for _, s in spheroids.iterrows())
    vacuole_volume = (4/3) * np.pi * vacuole['r']**3
    success_rate = len(spheroids) / args.N * 100
    
    # Create combined CSV file
    output_file = os.path.join(run_folder, f'{run_id}_combined.csv')
    
    try:
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            
            # Write run summary section
            writer.writerow(['=== Run Summary ==='])
            writer.writerow(['Parameter', 'Value'])
            summary_data = [
                ('Run ID', run_id),
                ('Timestamp', datetime.now().isoformat()),
                ('Seed', args.seed),
                ('Number of Spheroids Requested', args.N),
                ('Number of Spheroids Placed', len(spheroids)),
                ('Success Rate (%)', success_rate),
                ('Minimum Radius', args.min_radius),
                ('Maximum Radius', args.max_radius),
                ('Vacuole Radius', float(vacuole['r'])),
                ('Wall Thickness', args.wall_thickness),
                ('Grid Resolution (dx)', args.dx),
                ('Maximum Tries', args.max_tries),
                ('Total Spheroid Volume', float(total_spheroid_volume)),
                ('Vacuole Volume', float(vacuole_volume)),
                ('Total Volume', float(total_spheroid_volume + vacuole_volume)),
                ('Distribution Mean (mu)', np.log((args.min_radius + args.max_radius) / 2)),
                ('Distribution Sigma', args.sigma),
                ('Wall Radius Mean (mu)', np.log(args.wall_outer_radius)),
                ('Wall Radius Sigma', 0.1),
                ('Iterations', args.iterations),
                ('Optimization Max Iterations', args.optimmaxiter)
            ]
            writer.writerows(summary_data)
            
            # Add spacing between sections
            writer.writerow([])
            writer.writerow(['=== Detailed Cell Information ==='])
            
            # Write headers for detailed section
            detailed_headers = [
                'cell_id', 'type', 'x', 'y', 'z', 'radius', 'volume', 
                'distance_from_center', 'nearest_neighbor_center_distance',
                'nearest_neighbor_surface_distance', 'num_neighbors', 
                'local_density', 'distance_to_wall', 'p_value'
            ]
            writer.writerow(detailed_headers)
            
            # Write spheroid data
            for i, spheroid in enumerate(spheroids.itertuples(), 1):
                metrics = calculate_spheroid_metrics(
                    spheroid=spheroid,
                    spheroids_df=spheroids,
                    vacuole=vacuole,
                    max_radius=args.max_radius
                )
                
                row_data = [
                    i,  # cell_id
                    'Body',  # type
                    float(spheroid.x),
                    float(spheroid.y),
                    float(spheroid.z),
                    float(spheroid.r),
                    float(metrics['volume']),
                    float(metrics['distance_from_center']),
                    float(metrics['nearest_neighbor_center_distance']) if metrics['nearest_neighbor_center_distance'] is not None else None,
                    float(metrics['nearest_neighbor_surface_distance']) if metrics['nearest_neighbor_surface_distance'] is not None else None,
                    metrics['num_neighbors'],
                    float(metrics['local_density']),
                    float(metrics['distance_to_wall']),
                    float(spheroid.p)
                ]
                writer.writerow(row_data)
            
            # Write vacuole data
            vacuole_row = [
                len(spheroids) + 1,  # cell_id
                'Vacuole',  # type
                float(vacuole['x']),
                float(vacuole['y']),
                float(vacuole['z']),
                float(vacuole['r']),
                float(vacuole_volume),
                0.0,  # distance_from_center
                None,  # nearest_neighbor_center_distance
                None,  # nearest_neighbor_surface_distance
                None,  # num_neighbors
                None,  # local_density
                0.0,  # distance_to_wall
                float(vacuole['p'])
            ]
            writer.writerow(vacuole_row)
            
            # Add statistical summary section
            writer.writerow([])
            writer.writerow(['=== Statistical Summary ==='])
            
            # Calculate additional statistics
            avg_radius = spheroids['r'].mean()
            std_radius = spheroids['r'].std()
            avg_distance = spheroids.apply(
                lambda row: np.sqrt(row['x']**2 + row['y']**2 + row['z']**2), 
                axis=1
            ).mean()
            
            stats_data = [
                ('Average Spheroid Radius', float(avg_radius)),
                ('Standard Deviation of Radius', float(std_radius)),
                ('Average Distance from Origin', float(avg_distance)),
                ('Total Number of Spheroids', len(spheroids)),
                ('Packing Density (%)', float(total_spheroid_volume / vacuole_volume * 100))
            ]
            writer.writerows(stats_data)
        
        logging.info(f"Generated combined CSV file: {output_file}")
        return output_file
        
    except Exception as e:
        logging.error(f"Error writing combined CSV file: {str(e)}")
        raise

def generate_piff_file(df, dx=1.0, filename='output.piff'):
    """
    Generates a PIFF file with spheroids and surrounding wall.
    """
    piff_lines = []
    cell_id = 1

    # Separate vacuole and spheroids
    vacuole = df[df['bodyType'] == 'Vacuole'].iloc[0]
    spheroids = df[df['bodyType'] == 'APB']

    # Generate boxes for spheroids (Body cells)
    for index, spheroid in spheroids.iterrows():
        x0, y0, z0 = spheroid['x'], spheroid['y'], spheroid['z']
        R = spheroid['r']
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
                    voxel_center = (x + dx/2, y + dx/2, z + dx/2)
                    distance_sq = ((voxel_center[0] - x0)**2 +
                                (voxel_center[1] - y0)**2 +
                                (voxel_center[2] - z0)**2)
                    if distance_sq <= R**2:
                        line = f"{cell_id} Body {int(x)} {int(x+dx)} {int(y)} {int(y+dx)} {int(z)} {int(z+dx)}"
                        piff_lines.append(line)
        cell_id += 1

    # Generate boxes for the vacuole (Wall)
    wall_cell_id = cell_id
    x0, y0, z0 = vacuole['x'], vacuole['y'], vacuole['z']
    R_outer = vacuole['r']
    R_inner = R_outer - vacuole['r'] * 0.05  # Wall thickness

    # Define bounding box for the wall
    x_min = x0 - R_outer
    x_max = x0 + R_outer
    y_min = y0 - R_outer
    y_max = y0 + R_outer
    z_min = z0 - R_outer
    z_max = z0 + R_outer

    # Generate grid within bounding box
    x_vals = np.arange(x_min, x_max, dx)
    y_vals = np.arange(y_min, y_max, dx)
    z_vals = np.arange(z_min, z_max, dx)

    for x in x_vals:
        for y in y_vals:
            for z in z_vals:
                voxel_center = (x + dx/2, y + dx/2, z + dx/2)
                distance_sq = ((voxel_center[0] - x0)**2 +
                            (voxel_center[1] - y0)**2 +
                            (voxel_center[2] - z0)**2)
                if R_inner**2 <= distance_sq <= R_outer**2:
                    line = f"{wall_cell_id} Wall {int(x)} {int(x+dx)} {int(y)} {int(y+dx)} {int(z)} {int(z+dx)}"
                    piff_lines.append(line)

    with open(filename, 'w') as f:
        for line in piff_lines:
            f.write(line + '\n')

    print(f"PIFF file '{filename}' generated with {len(spheroids)} spheroids and surrounding wall.")
    logging.info(f"PIFF file '{filename}' generated with {len(spheroids)} spheroids and surrounding wall.")

def main(args):
    # Generate a unique run ID
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")

    # Create run directory
    runs_dir = 'runs'
    if not os.path.exists(runs_dir):
        os.makedirs(runs_dir)
    run_folder = os.path.join(runs_dir, run_id)
    os.makedirs(run_folder)

    # Setup logging
    stats_file = setup_logging(run_folder, args.seed)

    # Set random seed
    if args.seed is None:
        args.seed = random.randint(1, 1000000)
    random.seed(args.seed)
    np.random.seed(args.seed)
    rng = np.random.default_rng(args.seed)

    logging.info(f"Starting new run with ID: {run_id}")
    logging.info(f"Arguments: {args}")

    # Extract parameters
    args.wall_outer_radius = 60.0  # Fixed value

    # Generate spheroids
    df, pos_array, r_and_pos_array, dirmat_safe = genBalls3(
        bodies=args.N,
        wall_Radius_Mu=np.log(args.wall_outer_radius),
        wall_Radius_Sigma=0.1,
        mu=np.log((args.min_radius + args.max_radius) / 2),
        sigma=args.sigma,
        iterations=args.iterations,
        ndim=3,
        rng=rng,
        method='hcp',
        plist=None,
        pcterrcap=10.0,
        posOctant=True,
        optimmaxiter=args.optimmaxiter,
        maxVacuoleIterations=100
    )

    # Generate output files
    csv_file = write_combined_csv(run_folder, run_id, args, df)
    csv_filename = args.csv_output if args.csv_output else 'output.csv'
    df.to_csv(csv_filename, index=False)
    generate_piff_file(df, dx=args.dx, filename=args.output)

    # Save copies of files
    shutil.copy(args.output, os.path.join(run_folder, args.output))
    shutil.copy(csv_file, os.path.join(run_folder, os.path.basename(csv_file)))
    logging.info(f"Saved copy of PIFF file and CSV file in run folder: {run_folder}")

    # Copy to CompuCell3D directory
    cc3d = './CompuCell3D/cc3dSimulation/Simulation'
    shutil.copy(args.output, os.path.join(cc3d, args.output))

    logging.info(f"Run {run_id} completed successfully.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Generate PIFF file with N randomly clustered spheroids inside a surrounding hollow Wall.')
    parser.add_argument('--N', type=int, required=True, help='Number of internal spheroids to generate')
    parser.add_argument('--min_radius', type=float, default=3.0, help='Minimum radius for spheroids (default: 3)')
    parser.add_argument('--max_radius', type=float, default=8.0, help='Maximum radius for spheroids (default: 8)')
    parser.add_argument('--wall_outer_radius', type=float, default=40.0, help='Outer radius of the wall (default: 40)')
    parser.add_argument('--wall_thickness', type=float, default=2.0, help='Thickness of the wall (default: 2)')
    parser.add_argument('--dx', type=float, default=1.0, help='Resolution for grid boxes (default: 1.0)')
    parser.add_argument('--max_tries', type=int, default=1000, help='Maximum attempts to place each spheroid (default: 1000)')
    parser.add_argument('--output', type=str, default='output.piff', help='Output PIFF file name (default: output.piff)')
    parser.add_argument('--seed', type=int, help='Random seed for reproducibility (default: random)')
    parser.add_argument('--sigma', type=float, default=0.2, help='Standard deviation for log-normal distribution of spheroid radii')
    parser.add_argument('--iterations', type=int, default=4, help='Number of iterations for direction selection')
    parser.add_argument('--optimmaxiter', type=int, default=100, help='Maximum iterations for optimization')
    parser.add_argument('--csv_output', type=str, help='Output CSV file name')

    args = parser.parse_args()
    main(args)
