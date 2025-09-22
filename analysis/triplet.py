from __future__ import print_function
import matplotlib.pyplot as plt
import itertools
import mdtraj as md
import os
import numpy as np
import pandas as pd

def calculate_angle(o1_coords, o2_coords, o3_coords):
    """Calculate the angle between three oxygen atoms."""
    # Calculate vectors
    v1 = o2_coords - o1_coords
    v2 = o3_coords - o1_coords

    # Normalize vectors
    v1_normalized = v1 / np.linalg.norm(v1)
    v2_normalized = v2 / np.linalg.norm(v2)

    # Calculate the cosine of the angle using the dot product
    cos_theta = np.dot(v1_normalized, v2_normalized)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)  # Clip values to avoid NaNs due to floating-point errors

    # Calculate the angle in radians
    angle = np.arccos(cos_theta)

    # Convert the angle to degrees
    angle_deg = np.degrees(angle)

    return angle_deg

def triplet(DIRECTORY_PATH, XTC_FILE_NAME, GRO_FILE_NAME, gibbs_location, OUTPUT_FILE_PATH)
# Construct the full path to the XTC and GRO files
    xtc_file = os.path.join(DIRECTORY_PATH, folder, XTC_FILE_NAME)
    gro_file = os.path.join(DIRECTORY_PATH, folder, GRO_FILE_NAME)

    # Load XTC file
    t = md.load(xtc_file, top=gro_file)

    
    ## THIS IS WHEN GIBB'S DIVIDNG SURFACE IS DEFINED WITH OTHER FUNCTION
    ## THIS IS JUST TO REDUCE THE COMPUTATIONAL TIME
    
    z_value = gibbs_location  
  
    angles_data = []
    iter_max = t.n_frames
    for i in range(0, iter_max):
        START = i*5
        STOP = i*5+1
        STRIDE = 1

        chunk = t[START:STOP:STRIDE]

        # Precompute atom indices
        o_indices = np.array([a.index for a in chunk.topology.atoms if a.element.symbol == 'O' and a.residue.name == 'HOH' and 1.5 <= chunk.xyz[0, a.index, 0] <= 3.5 
                                and 1.5 <= chunk.xyz[0, a.index, 1] <= 3.5 and chunk.xyz[0, a.index, 2] <= z_value[0]+0.2])
        if len(o_indices) == 0:
            continue  # Skip frame if no valid O atoms found

        # Get coordinates of O atoms
        o_coords = chunk.xyz[0, o_indices, :]

        # Find triplet angles for each O atom
        for idx, o1_coords in enumerate(o_coords):
            # Find nearby O atoms within a distance of 0.33 using compute_neighbors with periodic boundary conditions
            nearby_o_indices = md.compute_neighbors(chunk, 0.33, np.array([o_indices[idx]]), haystack_indices=o_indices, periodic=True)[0]

            # Ensure there are at least four nearby O atoms
            if len(nearby_o_indices) < 2:
                continue   

            # Get the coordinates of the four closest O atoms
            closest_o_indices = nearby_o_indices[:2]
            closest_o_coords = chunk.xyz[0, closest_o_indices, :]

            # Ensure nearby O atoms are from HOH residue
            if all(chunk.topology.atom(oi).residue.name == 'HOH' for oi in closest_o_indices):
                # Calculate the angles for all unique pairs of the four nearest O atoms
                for j in range(2):
                    for k in range(j + 1, 2):
                        o2_coords = closest_o_coords[j]
                        o3_coords = closest_o_coords[k]
                        angle = calculate_angle(o1_coords, o2_coords, o3_coords)
                        angles_data.append(angle)
                
            
        print('   folder:', folder, '   frame number:', i + 1, '   Angles:', len(angles_data))


    file_path = OUTPUT_FILE_PATH
    
    # Save results to CSV file
    df = pd.DataFrame(angles_data, columns=['Angle'])
    output_filename = f'FILE_NAME.csv'
    output_file = os.path.join(file_path, output_filename)
    df.to_csv(output_file, index=False)
