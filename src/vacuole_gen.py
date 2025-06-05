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
import xml.etree.ElementTree as ET
import csv
import pandas as pd
from scipy.stats import ortho_group
from scipy.optimize import minimize
from scipy.spatial import distance_matrix

#TAKEN FROM GENBALLS07

def cartesian_product(a,b,c=[0]):
  return np.array(np.meshgrid(a,b,c)).T.reshape(-1,3)

def total_dist_to_pt(pos_array_as_vec,apb_df,ctr_to_use:np.ndarray=None):
  # all we use apb_df for is getting the number of points,
  # so we can reshape the input to be a matrix.
  pos_array = np.reshape(pos_array_as_vec,(len(apb_df),-1))
  if ctr_to_use is None:
    ctr_to_use=np.zeros_like(pos_array[1,:])
  dists = np.sqrt(np.sum((pos_array-ctr_to_use)**2,axis=1))
  #print("obj func: pos_array",pos_array)
  #print("obj func: sum(dists to origin)", np.sum(dists))
  return np.sum(dists)
  #  orig_lengths = np.expand_dims(np.sqrt(np.sum(dir_array**2,axis=1)),axis=1)

def move_along_dir(df,pos_array,ii,direction,stepsize,placeVacuole=False):
  # the direction (dir) passed in should usually be a unit vector,
  # but it's not actually required, so we don't check for it.
  dist_from_origin = (df.loc[0,"r"]+df.loc[ii,"r"])*(1.0+1e-8)
  pvals = df['p']
  # In summer 2023, Connor W. and Ross played with this, using this Desmos,
  # https://www.desmos.com/calculator/7kvz3hdynt
  # and found:
    # If you have two APBs, one with p-value p1 and one with p2,
    # taking min(p1,p2) is "not safe": it's easy to get a "false negative" for a collision--we might end up placing APBs that overlap.
    # taking max(p1,p2) is "safe": we can't have a "false negative" for a collision--we'll never accidentally place APBs that actually overlap.
    # In between those two is the arithmetic average, (p1+p2)/2, which we didn't prove is safe but seemed to work most of the time: "mostly safe"
    # To get the tightest clustering, we want to use the smallest possible combined pval, so (p1+p2)/2 is better/more compact than max(p1,p2)
    # We could get slightly better clustering by using the geometric mean, which is always <= the arithmetic mean; see
    # https://en.wikipedia.org/wiki/Inequality_of_arithmetic_and_geometric_means.
  #
  pvals = pvals[0:ii]
  pvals_e = np.expand_dims(pvals,axis=1)
  # now take the average of each already-placed APB's p-value with the p-value of the APB we're trying to place now.
  apvals = (pvals[0:ii] + df.loc[ii,"p"])/2 # arithmetic mean
  apvals_e = np.expand_dims(apvals,axis=1)

  # without the *(1.0+1e-8), sometimes it would still seem that the first
  # placement we tried overlapped the two balls by 1e-16 or so, which is just
  # roundoff error.
  stopflag = False
  # We only need distmins when usings an APB and not a Vacuole but not worth an if statement here
  distmins =  df.head(ii)["r"]+df.loc[ii,"r"] # each individual placed APB's radius, plus radius of this APB
  #print("distmins:")
  #print(distmins)
  while not stopflag:
    pos = direction*dist_from_origin
    # Compute distance from this APB's current trial position (pos) to all others that we've already placed:
    # dists = np.sqrt(np.sum((pos_array[0:ii,:]-pos)**2,axis=1))
    # This next line uses the p-value of each body to determine the distance between that body and
    # the current position estimate. While that's okay at this point, it gets weird
    # when we use the distance later and compare it to distmins: adding radii only
    # really works for p=2, circles/spheres.
    # For other values of p, it can (I think) over- or under-estimate the gap between
    # the APB and the candidate position.
    # But doing it "right" is extremely difficult mathematically--it can involve solving
    # high-order polynomials, essentially (technically that's for ellipses with p=2, but
    # this problem is just as hard)

    # This old code uses the p-values of the already-placed APBs and ignores the p-value of the one we're trying to place:
    #dists = (np.sum(np.abs(pos_array[0:ii,:]-pos)**pvals_e,axis=1))**(1.0/pvals)
    # Now let's also incorporate the p-value of the APB we're trying to place:
    dists = (np.sum(np.abs(pos_array[0:ii,:]-pos)**apvals_e,axis=1))**(1.0/apvals)
    if placeVacuole:
        close_enough = (dists + df.head(ii)["r"]) <= vacRadInner # aross15 updating from vacRad to vacRadInner
        # note that we do still want to use df's "r" column which is the APB's _outer_ radius,
        # and compare those with the vacuole's _inner_ radius.
        stopflag = not all(close_enough)
    else:
        # Make a set of logicals (True/False), to say whether this APB
        # is now far enough from each other APB:
        far_enough = dists >= distmins
        stopflag = all(far_enough)       

    dist_from_origin += stepsize
  # Done with while loop now  
  if placeVacuole:
    dist_from_origin -= 2*stepsize #we moved the vacuole just far enough to be unfeasible so move it back one step to feasibility
    pos = direction*dist_from_origin
  else:
    dist_from_origin -= 1*stepsize  # subtracting stepsize because we added it at the end of the last iteration of the loop.
    # no need to recompute position
  return (pos,dist_from_origin) 


def inter_APB_boundary_distances(pos_array_as_vec,apb_df):
  pos_array = np.reshape(pos_array_as_vec,(len(apb_df),-1))
  cur_dist_mat = distance_matrix(pos_array,pos_array) # between centers, not boundaries.
  #print(cur_dist_mat)
  #print("type:",type(apb_df))
  #print(apb_df)
  rvec = apb_df['r'].to_numpy()
  #print(rvec)
  # adding two vectors (or, the same vector twice) to get all pairs of sums:
  # https://stackoverflow.com/questions/47829946/how-to-sum-two-vectors-to-obtain-a-matrix-of-sums-of-all-pairs
  min_dist_mat =  rvec[:,None]+rvec
  #print(min_dist_mat)
  d = cur_dist_mat - min_dist_mat
  # But, we don't want the diagonal--who cares about distances from an APB to itself?
  # And, the matrix is symmetric, so we only care about the upper triangle, or the lower triangle, not both.
  # Let's extract the upper triangle part:
  # https://stackoverflow.com/questions/8905501/extract-upper-or-lower-triangular-part-of-a-numpy-matrix/44395030
  # a[np.triu_indices(3, k = 1)]
  dists_as_vec = d[np.triu_indices(d.shape[0],k=1)]
  #print("constr: pos_array",pos_array)
  #print(f"constr: amin(dists) {np.amin(dists_as_vec):.10e}")
  return dists_as_vec

def get_directions_hcp(bodies,ndim=3,rng=None,unitvectors=True,exclude_origin=True,random_rotate=True):

  side_length= int( np.ceil( bodies**(1.0/ndim) ) )
  hsl = int( np.ceil(side_length/2) ) # half-side-length
  hsl = hsl + (hsl+1)%2 # if hsl is even, hsl+1 is odd, so (mod 2) gives us +1
  ii_list = range(-hsl,hsl+1) # the +1 is so we include hsl, not just stop at hsl-1.
  jj_list = ii_list # don't need to .copy it since we won't be changing ii_list or jj_list or kk_list.
  if ndim==2:
    kk_list = [0]
  else:
    kk_list = ii_list

  iijjkk = cartesian_product(ii_list,jj_list,kk_list)
  if exclude_origin:
    keepthese = np.sum(iijjkk**2,axis=1) > 0
    iijjkk = iijjkk[keepthese,:]
  coords = np.zeros_like(iijjkk,dtype=float) #without dtype=float,
  # it set the array up with integers, which induced severe rounding trouble!

  coords[:,0] = 2*iijjkk[:,0] + (iijjkk[:,1]+iijjkk[:,2])%2
  coords[:,1] = np.sqrt(3)*(iijjkk[:,1]+(1/3.0)*(iijjkk[:,2] % 2))
  if ndim==2: # if 2-d, return only 2 columns; delete the "3rd" (index=2) column.
    coords = np.delete(coords,2,axis=1)
  else: # create the info that should go there:
    coords[:,2] = 2*np.sqrt(6)/3*iijjkk[:,2]

  if random_rotate:

    scipy_randomGen = ortho_group

    scipy_randomGen.random_state=rng

    dirmat = ortho_group.rvs(ndim)

    coords_rotated = dirmat @ coords.transpose() # not sure this works mathematically or in python!
    coords = coords_rotated.transpose()

  dists_true = np.sqrt(np.sum(coords**2,axis=1))
  dists_perturbed = dists_true + rng.uniform(low=0,high=0.001)
  dists_order = np.argsort(dists_perturbed)
  dir_array = coords[dists_order,:]
  #print(dir_array.shape) # debugging
  orig_lengths = np.expand_dims(np.sqrt(np.sum(dir_array**2,axis=1)),axis=1)
  if unitvectors:
    dir_array = dir_array / orig_lengths # turn each into a unit vector.

  return dir_array, orig_lengths


def get_directions_ortho_random(bodies,ndim=3,rng=None,unitvectors=True):

  dir_array = np.zeros((0,ndim))
  for ii in range(bodies):
    scipy_randomGen = ortho_group
    scipy_randomGen.random_state=rng
    dirmat = scipy_randomGen.rvs(ndim)
    dir_array = np.append(dir_array,dirmat,axis=0)

  orig_lengths = np.expand_dims(np.sqrt(np.sum(dir_array**2,axis=1)),axis=1)
  return dir_array, orig_lengths # all orig_lengths should be 1, hopefully!


def genBalls3(bodies=20, wall_Radius_Mu=6.8, wall_Radius_Sigma=0.34, mu=5, sigma=0.2, iterations=4,ndim=3,rng=None,method='hcp',plist:list=None,pcterrcap=10.0,posOctant=True,optimmaxiter=100,maxVacuoleIterations=100):

  # generate body radii:
  r_normals = rng.standard_normal(bodies)*sigma+mu # one for each APB
  r = np.exp(r_normals) # turn the normals into log-normals in nanometers

  if plist is None:
    plist = [args.pvals]*bodies # repeats pvals (the p-norm value) as many times as is specified.
  elif len(plist)==1:
    plist = plist*bodies

  stepsize = (pcterrcap/100.0)*np.min(r)

  if method=='hcp':
    dirmat,orig_lengths=get_directions_hcp(bodies,ndim,rng,unitvectors=True)
  elif method=='ortho_random':
    dirmat,orig_lengths=get_directions_ortho_random(bodies,ndim,rng)
  else:
    print(f'unknown method:{method}')
    return (None,None,None)
  dirmat_safe = dirmat.copy() # not sure we'll need it, but can't hurt.

  ii = 0
  d = {'bodynum': ii, 'r': r[ii],'p':plist[ii]}
  df = pd.DataFrame(data=[d])

  pos_array = np.zeros(shape=(bodies,3))  # always 3 dim: (x,y,z) but if ndim=2, then one coordinate will just be 0.

  for ii0 in range(bodies-1):
    ii = ii0+1
    #print(ii)
    d = {'bodynum': ii, 'r': r[ii],'p':plist[ii]}

    #MV  pd <- df
    df = pd.concat([df,pd.DataFrame(data=[d])],ignore_index=True)
    #print(df)

    ndir_to_try = min(iterations,dirmat.shape[0])
    pos_list = []
    dist_hist = np.inf*np.ones((ndir_to_try,))
    #print(dirmat.shape) # debugging
    for jj in range(ndir_to_try):
      #print(jj,ndir_to_try)
      dir=np.zeros((3,),dtype=float)
      dir[0:ndim]=dirmat[jj,:]
      #print(dir,orig_lengths[jj]) # debugging
      pos,dist_from_origin = move_along_dir(df,pos_array,ii,dir,stepsize)
      pos_list.append(pos)
      dist_hist[jj] = dist_from_origin
    # done with the for loop. Which one was the best?
    jjbest = np.argmin(dist_hist)
    pos = pos_list[jjbest]

    # MV - Desired ABP Placement, Recored This Somewhere
    pos_array[ii,:] = pos
    # and delete this direction from the dirmat:
    dirmat = np.delete(dirmat,jjbest,axis=0)
    orig_lengths = np.delete(orig_lengths,jjbest,axis=0)
    #print(pos_array)

  if optimmaxiter > 0: # since if maxiter==0, then don't do the optimization!
    optim_method='trust-constr' # interior-point method
    myoptions={'maxiter':optimmaxiter}
    
    # calculate compactness (distance from every body to every other body) of the starting bodies
    #print(pos_array)
    body_starting_distances = 0
    for body1 in pos_array:
      for body2 in pos_array:
        body_distance = math.sqrt((body1[0] - body2[0])**2 + (body1[1] - body2[1])**2 + (body1[2] - body2[2])**2) 
        body_starting_distances += body_distance
    #print ("body_starting_distances = ", body_starting_distances)
    
    optim_df = df
    x0 = np.ravel(pos_array) # turn into a 1-dimensional array, since that's what minimize works with.
    cons= ({'type': 'ineq', 'args': (optim_df,), 'fun': inter_APB_boundary_distances})
    # "inequality means that it is to be non-negative"

    # print the objective function value (ofv) for the current positions:
    ofv_original = total_dist_to_pt(x0,optim_df)
    print("ofv_original:",ofv_original)
    #scipy.optimize.minimize(fun, x0, args=(), method='trust-constr', hess=None, hessp=None, bounds=None, constraints=(), tol=None, callback=None, options={'grad': None, 'xtol': 1e-08, 'gtol': 1e-08, 'barrier_tol': 1e-08, 'sparse_jacobian': None, 'maxiter': 1000, 'verbose': 0, 'finite_diff_rel_step': None, 'initial_constr_penalty': 1.0, 'initial_tr_radius': 1.0, 'initial_barrier_parameter': 0.1, 'initial_barrier_tolerance': 0.1, 'factorization_method': None, 'disp': False})
    res = minimize(total_dist_to_pt, x0, args=optim_df,method=optim_method,constraints=cons,options=myoptions) 
    # res stands for results, put into pos_array
    pos_array_safe = pos_array.copy()
    pos_array_as_vec = res['x']
    pos_array = np.reshape(pos_array_as_vec,(len(optim_df),-1))
    
    # calculate compactness (distance from every body to every other body) of the ending bodies
    #print(pos_array)
    body_ending_distances = 0
    for body1 in pos_array:
      for body2 in pos_array:
        body_distance = math.sqrt((body1[0] - body2[0])**2 + (body1[1] - body2[1])**2 + (body1[2] - body2[2])**2) 
        body_ending_distances += body_distance
    #print ("body_ending_distances = ", body_ending_distances)
    compactness = (body_starting_distances / body_ending_distances)-1   #0 if no compaction, positive if compacted
    #print ("compactness = ", compactness)
    
    #
    ofv_final = total_dist_to_pt(pos_array_as_vec,optim_df)
    print("ofv_final:",ofv_final)
    #Done with optimization step
    
  else:
    ofv_original = 'No optimization'
    ofv_final = 'No optimization'
    
       
  df['bodyType'] = "APB"  #MV - Dataframe Label
    
  # We'll want this if we don't do the posOctant shift:
  r_and_pos_array = np.hstack((np.expand_dims(r,axis=1),pos_array))

  if posOctant: 
    pos_array_means = np.mean(pos_array,axis=0) # mean of each column
    pos_array_shifted_to_origin = pos_array - pos_array_means
    # Generate a vacuole INNER radius that will include all of the Autophagic Bodies
    # We'll suppose that the vacuole has p=2
    vacp = 2.0
    dists_to_origin = (np.sum(np.abs(pos_array_shifted_to_origin)**vacp,axis=1))**(1.0/vacp)
    dists_plus_rad = dists_to_origin + r
    vacRadInner = 0 # aross15 changing this from just vacRad to vacRadInner to be more clear
    iterCount = 0
    while any(dists_plus_rad > vacRadInner)and(iterCount<maxVacuoleIterations):
        #generate new vacRadInner according to a lognormal
        r_normals = rng.standard_normal(1)[0]*wall_Radius_Sigma+wall_Radius_Mu #aross15 adding the [0] to get just a scalar, not an np.array
        vacRadInner = np.exp(r_normals) # turn the normals into log-normals
        iterCount += 1
    if(iterCount>=maxVacuoleIterations):
        print("Warning: Maxed Out on Vacuole Size Iterations.")
        print(f"Iterations:{iterCount}")
        print(f"{dists_plus_rad}")
        print(f"{wall_Radius_Sigma}")
        print(f"{wall_Radius_Mu}")
        vacRad = -9999

    vacRadOuter = 1.05*vacRadInner

    # Done with checking, now shift everything:
    pos_array_shifted_to_1st_octant = pos_array_shifted_to_origin + vacRadOuter
    r_and_pos_array = np.hstack((np.expand_dims(r,axis=1),pos_array_shifted_to_1st_octant))

    # Then add the vacuole wall as the first item, and radius in the first column.
    # from R:    cellSize=2*vacRad ; output_spheregen <- rbind(cellSize/2, output_spheregen)
    vac_pos_array = np.ones_like(pos_array[0,:])*vacRadOuter #  x, y, and z if needed.

    vac_r_and_pos = np.hstack((vacRadOuter,vac_pos_array)) # we need to record vacRadOuter not vacRadInner to make sure that it all fits in the positive octant
    r_and_pos_array_w_vac = np.vstack((vac_r_and_pos,r_and_pos_array))
    pos_array = np.delete(r_and_pos_array_w_vac,0,axis=1) # copy everything except 1st col, which is r.
    d = {'bodynum': -100, 'bodyType':'Vacuole', 'rOuter': vacRadOuter,'p':vacp, 'rInner':vacRadInner} #sbackues updated 'r' to 'rOuter'
    df_just_vac = pd.DataFrame(data=[d])
    #MV Pandas Concatenation of Vac Parameter and Dataframe
    df = pd.concat([df_just_vac,df],ignore_index=True)
  dfpos = pd.DataFrame.from_records(pos_array)

  # rename the columns
  dfpos.rename(columns={0: 'x', 1: 'y', 2: 'z'}, inplace=True)
  # add the resulting columns to the main dataframe:
  df = pd.concat([df.reset_index(drop=True), dfpos], axis=1)

  return(df,pos_array,r_and_pos_array,dirmat_safe,iterCount, ofv_original, ofv_final, compactness)


def log_statistics(args, df):
    """
    Logs statistics about the generated spheroids and wall.

    Parameters:
        args: The argparse arguments.
        df (DataFrame): DataFrame containing spheroid and vacuole data.
    """
    # Extract spheroid and vacuole data from df
    spheroids = df[df['bodyType'] == 'APB']
    vacuole = df[df['bodyType'] == 'Vacuole'].iloc[0]

    stats = {
        "seed": args.seed,
        "num_spheroids_requested": args.N,
        "num_spheroids_placed": len(spheroids),
        "vacuole_outer_radius": vacuole['rOuter'], # aross15 making it more clear
        "vacuole_inner_radius": vacuole['rInner'],    #sbackues
        "dx": args.dx,
        "optimmaxiter": args.optimmaxiter
    }

    # Calculate additional statistics
    spheroid_radii = spheroids['r'].values
    stats["mean_spheroid_radius"] = float(spheroid_radii.mean())
    stats["min_spheroid_radius"] = float(spheroid_radii.min())
    stats["max_spheroid_radius"] = float(spheroid_radii.max())

    # Calculate total volume of spheroids
    total_volume = ((4/3) * np.pi * (spheroid_radii ** 3)).sum()
    stats["total_spheroid_volume"] = float(total_volume)
    
    # Vacuole volume
    vacuole_volume = (4/3) * np.pi * vacuole['rInner'] ** 3
    stats["vacuole_volume"] = float(vacuole_volume.item()) if isinstance(vacuole_volume, np.ndarray) else float(vacuole_volume)

    # Log the statistics
    logging.info(f"Statistics: {stats}")

#Set up logging
def setup_logging(runs_dir, seed):
    log_file = os.path.join(runs_dir, f'run.log')
    logging.basicConfig(filename=log_file, level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s')

    # Create a separate statistics file
    stats_file = os.path.join(runs_dir, f'statistics.json')
    
    # Log the seed
    logging.info(f"Random seed: {seed}")

    return stats_file

#genballs
def generate_piff_file(df, dx, show_wall, filename='output.piff'):
    """
    Generates a PIFF file with spheroids and surrounding wall.

    Parameters:
        df (DataFrame): DataFrame containing spheroid and vacuole data.
        dx (float): The resolution for grid boxes.
        filename (str): The name of the output PIFF file.
    """
    
    df.to_csv('dataframe_df.csv')
    piff_lines = []
    cell_id = 1  # Start CellID from 1
    max_voxel_value = -float('inf')  # Initialize max voxel value as negative infinity

    # Separate vacuole and spheroids
    vacuole = df[df['bodyType'] == 'Vacuole'].iloc[0]
    spheroids = df[df['bodyType'] == 'APB']

    # Generate boxes for spheroids (Body cells)
    # dx is the scale factor, in nm per voxel.  Dividing by dx converts nm to voxels.  
    for index, spheroid in spheroids.iterrows():
        x0, y0, z0 = spheroid['x'] / dx, spheroid['y'] / dx, spheroid['z'] / dx   #Center of spheriod in voxels 
        R = spheroid['r'] / dx  #R = radius of spheroid in voxels 
        pvals = spheroid['p']    #pvals = p-value (p-norm); 2=spherical
        #print("pvals =", pvals)
        # Define bounding box (in voxels) 
        x_min = x0 - R
        x_max = x0 + R
        y_min = y0 - R
        y_max = y0 + R
        z_min = z0 - R
        z_max = z0 + R

        # Ensure all inputs to np.arange are scalars
        x_min = float(x_min.item()) if isinstance(x_min, np.ndarray) else float(x_min)
        x_max = float(x_max.item()) if isinstance(x_max, np.ndarray) else float(x_max)
        y_min = float(y_min.item()) if isinstance(y_min, np.ndarray) else float(y_min)
        y_max = float(y_max.item()) if isinstance(y_max, np.ndarray) else float(y_max)
        z_min = float(z_min.item()) if isinstance(z_min, np.ndarray) else float(z_min)
        z_max = float(z_max.item()) if isinstance(z_max, np.ndarray) else float(z_max)
        dx = float(dx.item()) if isinstance(dx, np.ndarray) else float(dx)

        # Generate values using np.arange.  All values in voxels.  
        x_vals = np.arange(x_min, x_max, 1)
        y_vals = np.arange(y_min, y_max, 1)
        z_vals = np.arange(z_min, z_max, 1)      

        for x in x_vals:
            for y in y_vals:
                for z in z_vals:
                    # Calculate center of the voxel
                    voxel_center = (x + 0.5, y + 0.5, z + 0.5)
                    distance_sq = ((voxel_center[0] - x0) ** pvals +
                                   (voxel_center[1] - y0) ** pvals +
                                   (voxel_center[2] - z0) ** pvals)
                    if distance_sq <= R ** pvals:
                        line = f"{cell_id} Body {int(x)} {int(x)} {int(y)} {int(y)} {int(z)} {int(z)}"
                        piff_lines.append(line)
                        
                        # Track the maximum voxel value
                        max_voxel_value = max(max_voxel_value, x, y, z)
                        
        cell_id += 1  # Increment CellID for the next spheroid

    #Generate the wall, if desired
    if show_wall.lower() == "true":

      # Generate boxes for the vacuole (Wall)
      wall_cell_id = cell_id  # Assign a unique CellID for the wall
      
      # Scale physical coordinates (nm) and radii down to grid units (voxels).
      x0, y0, z0 = vacuole['x'] / dx, vacuole['y'] / dx, vacuole['z'] / dx   #Center of the vacuol in voxels
      R_outer = vacuole['rOuter'] / dx
      R_inner = vacuole['rInner'] / dx #There is both an outer and inner because the vacuole is a hollow sphere.  In voxels.  
  
      # Define bounding box for the wall
      x_min = x0 - R_outer
      x_max = x0 + R_outer
      y_min = y0 - R_outer
      y_max = y0 + R_outer
      z_min = z0 - R_outer
      z_max = z0 + R_outer
  
      # Ensure all inputs to np.arange are scalars
      x_min = float(x_min.item()) if isinstance(x_min, np.ndarray) else float(x_min)
      x_max = float(x_max.item()) if isinstance(x_max, np.ndarray) else float(x_max)
      y_min = float(y_min.item()) if isinstance(y_min, np.ndarray) else float(y_min)
      y_max = float(y_max.item()) if isinstance(y_max, np.ndarray) else float(y_max)
      z_min = float(z_min.item()) if isinstance(z_min, np.ndarray) else float(z_min)
      z_max = float(z_max.item()) if isinstance(z_max, np.ndarray) else float(z_max)
      dx = float(dx.item()) if isinstance(dx, np.ndarray) else float(dx)
  
      # Generate the grid points for voxel placement, the dx step size ensures that the grid points are spaced correctly
      x_vals = np.arange(x_min, x_max, 1)
      y_vals = np.arange(y_min, y_max, 1)
      z_vals = np.arange(z_min, z_max, 1)
  
      for x in x_vals:
          for y in y_vals:
              for z in z_vals:
                  # Calculate center of the voxel
                  voxel_center = (x + 0.5, y + 0.5, z + 0.5)
                  distance_sq = ((voxel_center[0] - x0) ** 2 +
                                 (voxel_center[1] - y0) ** 2 +
                                 (voxel_center[2] - z0) ** 2)
                  if R_inner ** 2 <= distance_sq <= R_outer ** 2:
                      line = f"{wall_cell_id} Wall {int(x)} {int(x)} {int(y)} {int(y)} {int(z)} {int(z)}"
                      piff_lines.append(line)
  
                      # Track the maximum voxel value
                      max_voxel_value = max(max_voxel_value, x, y, z)
                    
    # Write the PIFF content to a file
    with open(filename, 'w') as f:
        for line in piff_lines:
            f.write(line + '\n')

    if show_wall.lower() == "true": 
      print(f"PIFF file '{filename}' generated with {len(spheroids)} spheroids and surrounding wall.")
      logging.info(f"PIFF file '{filename}' generated with {len(spheroids)} spheroids and surrounding wall.")
    elif show_wall.lower() == "false":
      print (f"PIFF file '{filename}' generated with {len(spheroids)} spheroids.")
      logging.info(f"PIFF file '{filename}' generated with {len(spheroids)} spheroids.")
    print(f"Greatest voxel value found: {max_voxel_value}")
    logging.info(f"Greatest voxel value found: {max_voxel_value}")
    

    xml_file_path = './CompuCell3D/cc3dSimulation/Simulation/clustertest.xml'
    update_dimensions_in_xml(xml_file_path, max_voxel_value + 3)
    
    #return piff_lines

def update_dimensions_in_xml(xml_file_path, greatest_voxel_value):
    # Parse the XML file
    tree = ET.parse(xml_file_path)
    root = tree.getroot()

    # Find the <Dimensions> element
    for elem in root.iter('Dimensions'):
        # Update the x, y, and z attributes with the greatest voxel value
        elem.set('x', str(greatest_voxel_value))
        elem.set('y', str(greatest_voxel_value))
        elem.set('z', str(greatest_voxel_value))
        print(f"Updated dimensions to: x={greatest_voxel_value}, y={greatest_voxel_value}, z={greatest_voxel_value}")

    # Save the changes back to the XML file
    tree.write(xml_file_path)

def calculate_spheroid_metrics(spheroid, spheroids_df, vacuole, max_radius=8.0):
    """
    Calculate additional metrics for a single spheroid.
    """
    volume = (4/3) * np.pi * (spheroid.r**3)
    
    distance_from_center = np.sqrt(
        (spheroid.x - vacuole['x'])**2 + 
        (spheroid.y - vacuole['y'])**2 + 
        (spheroid.z - vacuole['z'])**2
    )
       
    distance_to_wall = abs(distance_from_center - vacuole['rInner'])
    
    
    
    return {
        'volume': volume,
        'distance_from_center': distance_from_center,
        'distance_to_wall': distance_to_wall
    }

def write_combined_csv(run_folder, run_id, args, df):
    """
    Write a single CSV file containing both summary and detailed information.
    """
    
    # Separate vacuole and spheroids
    vacuole = df[df['bodyType'] == 'Vacuole'].iloc[0]
    spheroids = df[df['bodyType'] == 'APB']
    
    # Calculate global statistics
    # Vacuole volume based on inner radius; the hollow volume that the spheres can be in is what matters #sbackues
    total_spheroid_volume = sum((4/3) * np.pi * (s['r']**3) for _, s in spheroids.iterrows())
    vacuole_volume = (4/3) * np.pi * (vacuole['rInner'].item() if isinstance(vacuole['rInner'], np.ndarray) else vacuole['rInner'])**3
    success_rate = len(spheroids) / args.N * 100
    
    
    # Create output file path
    output_file = os.path.join(run_folder, f'{run_id}_combined.csv')
    
    try:
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            
            # === Run Summary Section ===
            writer.writerow(['=== Run Summary ==='])
            writer.writerow(['Parameter', 'Value'])
            summary_data = [
                ('Run ID', run_id),
                ('Timestamp', datetime.now().isoformat()),
                ('Seed', args.seed),
                ('Number of Spheroids Requested', args.N),
                ('Number of Spheroids Placed', len(spheroids)),
                ('Success Rate (%)', success_rate),
                ('Body Radius Mu', args.mu),
                ('Body Radius Sigma', args.sigma),
                ('Wall Radius Mu', args.wall_radius_mu),
                ('Wall Radius Sigma', args.wall_radius_sigma),
                ('Body Number Mu', args.mu_body_number),
                ('Body Number Sigma', args.sigma_body_number), 
                ('Vacuole Inner Radius', float(vacuole['rInner'].item() if isinstance(vacuole['rInner'], np.ndarray) else vacuole['rInner'])), #sbackues
                ('Grid Resolution (dx)', args.dx),
                ('Maximum Tries', args.max_tries),
                ('Total Spheroid Volume', float(total_spheroid_volume)),
                ('Vacuole Volume', float(vacuole_volume)),
                ('Total Volume', float(total_spheroid_volume + vacuole_volume)),
                ('Iterations', args.iterations),
                ('Optimization Max Iterations', args.optimmaxiter)
            ]
            writer.writerows(summary_data)
            
            # === Detailed Cell Information Section ===
            writer.writerow([])
            writer.writerow(['=== Detailed Cell Information ==='])
            writer.writerow([
                'cell_id', 'type', 'x', 'y', 'z', 'radius', 'volume',
                'distance_from_center', 'nearest_neighbor_center_distance',
                'nearest_neighbor_surface_distance', 'num_neighbors',
                'local_density', 'distance_to_wall', 'p_value'
            ])
            
            # Process and write spheroid data first
            for idx, spheroid in enumerate(spheroids.itertuples(), 1):
                metrics = calculate_spheroid_metrics(spheroid, spheroids, vacuole)
                
                writer.writerow([
                    idx,  # cell_id
                    'Body',  # type
                    spheroid.x,
                    spheroid.y,
                    spheroid.z,
                    spheroid.r,
                    metrics['volume'],
                    metrics['distance_from_center'],
                    metrics['distance_to_wall'],
                    spheroid.p
                ])
            
            # Write vacuole data last
            writer.writerow([
                len(spheroids) + 1,  # cell_id
                'Vacuole',  # type
                vacuole['x'],
                vacuole['y'],
                vacuole['z'],
                vacuole['rInner'],
                vacuole_volume,
                0.0,  # distance_from_center
                '',  # nearest_neighbor_center_distance
                '',  # nearest_neighbor_surface_distance
                '',  # num_neighbors
                '',  # local_density
                0.0,  # distance_to_wall
                vacuole['p']
            ])
            
            # === Statistical Summary Section ===
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

def write_body_size_combined_csv(runs_dir, run_id, args, df):
    """
    Write a single CSV file for all of the runs containing body size information.    
    """
    
    # Separate vacuole and spheroids
    vacuole = df[df['bodyType'] == 'Vacuole'].iloc[0]
    spheroids = df[df['bodyType'] == 'APB']
    
    # Create combined output file in the runs directory for all bodies, to check body size distributions
    body_output_file = os.path.join(runs_dir, 'Body_Size_combined.csv')
        
    try:
        with open(body_output_file, 'a', newline='') as f:
            writer = csv.writer(f)
            
            # === Detailed Cell Information Section ===
            writer.writerow([
                'Run ID', 
                'Cell ID', 
                'type', 
                'x', 
                'y', 
                'z', 
                'radius', 
                'volume',
                'distance_from_center',
                'nearest_neighbor_center_distance',
                'nearest_neighbor_surface_distance',
                'num_neighbors',
                'local_density',
                'distance_to_wall', 
                'p_value', 
                'Body_Radius_Mu', 
                'Body_Radius_Sigma'
            ])
            
            # Process and write spheroid data first
            for idx, spheroid in enumerate(spheroids.itertuples(), 1):
                metrics = calculate_spheroid_metrics(spheroid, spheroids, vacuole)
                
                writer.writerow([
                    run_id,
                    idx,  # cell_id
                    'Body',  # type
                    spheroid.x,
                    spheroid.y,
                    spheroid.z,
                    spheroid.r,
                    metrics['volume'],
                    metrics['distance_from_center'],
                    metrics['distance_to_wall'],
                    spheroid.p,
                    args.mu,
                    args.sigma
                ])
        logging.info(f"Updated combined body size CSV file: {body_output_file}")
        return body_output_file
              
    except Exception as e:
        logging.error(f"Error writing completely combined body size CSV file")
        raise

def write_vacuole_data_csv(runs_dir, run_id, args, df, iterCount, ofv_original, ofv_final, compactness):
    """
    Write a single CSV file for all of the runs containing information on the vacuole size and body number.  Puts it in the "runs" folder.  
    """
    
    # Separate vacuole and spheroids
    vacuole = df[df['bodyType'] == 'Vacuole'].iloc[0]
    spheroids = df[df['bodyType'] == 'APB']
    
    # Calculate global statistics
    # Vacuole volume based on inner radius; the hollow volume that the spheres can be in is what matters #sbackues
    total_spheroid_volume = sum((4/3) * np.pi * (s['r']**3) for _, s in spheroids.iterrows())
    vacuole_volume = (4/3) * np.pi * (vacuole['rInner'].item() if isinstance(vacuole['rInner'], np.ndarray) else vacuole['rInner'])**3
    success_rate = len(spheroids) / args.N * 100
    
    # Calculate additional statistics
    avg_radius = spheroids['r'].mean()
    std_radius = spheroids['r'].std()
    largest_radius = spheroids['r'].max()
    avg_distance = spheroids.apply(
        lambda row: np.sqrt(row['x']**2 + row['y']**2 + row['z']**2), 
        axis=1
    ).mean()
    if (ofv_original) == 'No optimization':
      optim_factor = 'No optimization'  
    else:
      optim_factor = (ofv_original/ofv_final) - 1
    
            
    # Create combined output file in the runs directory for all runs
    vac_output_file = os.path.join(runs_dir, 'Vacuole_Data_combined.csv')
    
    try:  
        with open(vac_output_file, 'a', newline='') as f:
            writer = csv.writer(f)

            # Column headers were created in Pipeline_Interface
            # outputting a line of actual info
            writer.writerow([
                run_id,
                datetime.now().isoformat(),
                args.seed,
                args.dx,
                args.mu,    #body radius mu
                args.sigma, #body radius sigma
                args.pvals, #p-norm value
                float(avg_radius),
                float(std_radius),
                float(largest_radius),
                float(avg_distance),
                args.mu_body_number,
                args.sigma_body_number,
                args.N,
                len(spheroids),
                success_rate,
                args.iterations,     #Body placement attempts
                args.optimmaxiter,   #Max Number of iterations of the clustering code
                ofv_original,        #ofv before clustering  
                ofv_final,           #ofv after clustering (should be less)  
                optim_factor,
                compactness,
                args.wall_radius_mu,
                args.wall_radius_sigma,
                float(vacuole['rInner'].item() if isinstance(vacuole['rInner'], np.ndarray) else vacuole['rInner']),
                float(vacuole_volume),
                float(total_spheroid_volume),
                float(total_spheroid_volume / vacuole_volume * 100),
                args.max_tries,
                iterCount,
                float(vacuole['x'].item() if isinstance(vacuole['x'], np.ndarray) else vacuole['x']),  # Add x-coordinate
                float(vacuole['y'].item() if isinstance(vacuole['y'], np.ndarray) else vacuole['y']),  # Add y-coordinate
                float(vacuole['z'].item() if isinstance(vacuole['z'], np.ndarray) else vacuole['z'])   # Add z-coordinate
                ])                                           
                
        logging.info(f"Updated combined vacuole data CSV file: {vac_output_file}")
        return vac_output_file
              
    except Exception as e:
        logging.error(f"Error writing completely combined vacuole data CSV file")
        raise


def main(args):

    # Generate a unique run ID based on the current date and time
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    run_folder = args.run_folder

    # Setup logging and get stats file path
    stats_file = setup_logging(run_folder, args.seed)

    if args.seed is None:
        args.seed = random.randint(1, 1000000)  # Generate a random seed if not provided
    random.seed(args.seed)
    np.random.seed(args.seed)
    rng = np.random.default_rng(args.seed)

    logging.info(f"Starting new run with ID: {run_id}")
    logging.info(f"Arguments: {args}")

    # Extract parameters from arguments
    N_SPHEROIDS = args.N
    BODY_MU = args.mu
    BODY_SIGMA = args.sigma
    WALL_RADIUS_MU = args.wall_radius_mu
    WALL_RADIUS_SIGMA = args.wall_radius_sigma
    DX = args.dx
    MAX_TRIES = args.max_tries
    filename = args.output
    PIFF = args.PIFF

    # Generate spheroids using genBalls3
    df, pos_array, r_and_pos_array, dirmat_safe, iterCount, ofv_original, ofv_final, compactness = genBalls3(
        bodies=N_SPHEROIDS,
        wall_Radius_Mu=WALL_RADIUS_MU,
        wall_Radius_Sigma=WALL_RADIUS_SIGMA,  
        mu = BODY_MU,
        sigma = BODY_SIGMA,
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
  
    # Log statistics
    log_statistics(args, df)
    
    # Write the combined CSVs directly to the run folder
    write_body_size_combined_csv(run_folder, run_id, args, df)
    write_vacuole_data_csv(run_folder, run_id, args, df, iterCount, ofv_original, ofv_final, compactness)
    
    # If desired, generate PIFF file and save it to the simulation folder (for cc3d use)
    if PIFF == 1 or PIFF == 2:
        generate_piff_file(df, dx=args.dx, show_wall = args.show_wall, filename=filename)

        cc3d = './CompuCell3D/cc3dSimulation/Simulation'
        if os.path.exists(cc3d):
            shutil.copy(filename, os.path.join(cc3d, filename))
            logging.info(f"Copied PIFF file to CC3D simulation folder")
        
   
    #If desired Save a copy of the PIFF file in a subfolder within the run folder for later inspection
    # Also add statistics for that run as a csv to that folder.  
    if PIFF == 2:
        run_subfolder = os.path.join(run_folder, run_id)
        os.makedirs(run_subfolder)
        piff_copy_path = os.path.join(run_subfolder, filename)
        shutil.copy(filename, piff_copy_path)
        logging.info(f"Saved copy of PIFF file in run folder: {piff_copy_path}")
        write_combined_csv(run_subfolder, run_id, args, df)

    logging.info(f"Run {run_id} completed successfully.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Generate PIFF file with N randomly clustered spheroids inside a surrounding hollow Wall.')
    parser.add_argument('--run_folder', type=str, required=True, help='Folder to save the results into')
    parser.add_argument('--N', type=int, required=True, help='Number of internal spheroids to generate')
    parser.add_argument('--mu', type=float, required=True, help='Log-normal mean for spheroid radii')
    parser.add_argument('--sigma', type=float, required=True, help='Log-normal sigma for spheroid radii') 
    parser.add_argument('--pvals', type = float, required = False, help = 'p-value (p-norm) for controlling body sphericity.  2 = spherical')
    parser.add_argument("--wall_radius_mu", type=float, required=True, help="Log-normal mean for wall radius")
    parser.add_argument("--wall_radius_sigma", type=float, required=True, help="Log-normal sigma for wall radius")    
    parser.add_argument('--mu_body_number', type=float, required=True, help='Log-normal mean for body number')   
    parser.add_argument('--sigma_body_number', type=float, required=True, help='Log-normal sigma for body number')     
    parser.add_argument('--dx', type=float, required=True, help='Resolution for grid boxes, in nm per pixel')
    parser.add_argument('--show_wall', type=str, required=True, help='Whether or not to show the vacuole wall in the PIFF file')
    parser.add_argument('--optimmaxiter', type=int, required=True, help='Maximum iterations for optimization') 
    parser.add_argument('--wall_outer_radius', type=float, default=40.0, help='Outer radius of the wall (default: 40)')
    parser.add_argument('--wall_thickness', type=float, default=2.0, help='Thickness of the wall for visualization (default: 2) - not used in the PIFF file')
    parser.add_argument('--max_tries', type=int, default=1000, help='Maximum attempts to place each spheroid (default: 1000)')
    parser.add_argument('--output', type=str, default='output.piff', help='Output PIFF file name (default: output.piff)')
    parser.add_argument('--seed', type=int, help='Random seed for reproducibility (default: random)')
    parser.add_argument('--iterations', type=int, default=4, help='Number of iterations for direction selection')
    parser.add_argument('--PIFF', type=int, default=1, help='0=no PIFF, 1=PIFF overwritten, 2=PIFF saved')
    
    

    args = parser.parse_args()
    main(args)

