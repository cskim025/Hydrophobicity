from __future__ import print_function
import matplotlib.pyplot as plt
import itertools
import mdtraj as md
import os
import numpy as np
import pandas as pd

def calculate_angle_with_z_axis(traj, n_index, h_index):
    """Calculate the angle between the N-H vector and the z-axis."""
    # Get coordinates of the atoms
    n_coords = traj.xyz[0, n_index, :]
    h_coords = traj.xyz[0, h_index, :]

    # Compute vector from N to H
    nh_vector = h_coords - n_coords

    # Define the z-axis direction vector
    z_axis = np.array([0, 0, 1])

    # Normalize the NH vector
    nh_vector_normalized = nh_vector / np.linalg.norm(nh_vector)

    # Calculate the cosine of the angle using the dot product
    cos_theta = np.dot(nh_vector_normalized, z_axis)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)  # Clip values to avoid NaNs due to floating-point errors

    # Calculate the angle in radians
    angle = np.arccos(cos_theta)

    # Convert the angle to degrees
    angle_deg = np.degrees(angle)

    return angle_deg

#folders = ['FILE_1', 'FILE_2]

def angle_calc(FILE_PATH, FILE_NAME, folders, z_value, START_TIME, ATOM_NUM_AFT_N, OUTFILE_NAME):
  
  for folder in folders:
      # Construct the full path to the XTC and GRO files
      xtc_file = os.path.join(r'FILE_PATH', folder, FILE_NAME + '.xtc')
      gro_file = os.path.join(r'FILE_PATH', folder, FILE_NAME + '.gro')
  
      # Load XTC file
      t = md.load(xtc_file, top=gro_file)
  
      angles_data = []
  
      for i in range(START_TIME, t.n_frames):
          START = i
          STOP = i + 1
          STRIDE = 100
  
          chunk = t[START:STOP:STRIDE]
  
          len_var = len([a.index for a in chunk.topology.atoms if a.element.symbol == 'N'])
  
          angle_tempo = []
  
          for j in range(0, len_var):
              # get the indices of N and its first H
              n_index = [a.index for a in chunk.topology.atoms if a.element.symbol == 'N'][j]
              target_index = [a.index for a in chunk.topology.atoms][n_index+ATOM_NUM_AFT_N]
  
              # Compute the angle
              angle = calculate_angle_with_z_axis(chunk, n_index, target_index)
  
              angle_tempo.append(angle)
  
          angles_data.append(angle_tempo)
  
          print('   folder:', folder, '   frame number:', i+1, '   The angle:', angle_tempo)
  
      # Save results to CSV file
      df = pd.DataFrame(angles_data)
      output_filename = OUTFILE_NAME+'.csv'
      df.to_csv(output_filename, index=False)

