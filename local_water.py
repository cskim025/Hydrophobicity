import os
import mdtraj as md
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import os
import mdtraj as md
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


# List of folders
folders = ['FILE_1', 'FILE_2']

def local_water(DIRECTORY_PATH, XTC_FILE_NAME, GRO_FILE_NAME, gibbs_location, OUTPUT_FILE_PATH)

  for folder in folders:
      # Construct the full path to the XTC and GRO files
      xtc_file = os.path.join(DIRECTORY_PATH, folder, XTC_FILE_NAME+'.xtc')
      gro_file = os.path.join(DIRECTORY_PATH, folder, GRO_FILE_NAME+'.gro')
  
      # Load XTC file
      t = md.load(xtc_file, top=gro_file)
  
      whole_water = []
  
      # Determine slab height
      table, bonds = t.topology.to_dataframe()
      z_location = gibbs_location

      iter_max = t.n_frames
      for i in range(0, iter_max):
          START = i
          STOP = i + 1
          STRIDE = 100
  
          chunk = t[START:STOP:STRIDE]
  
          # Extract Oxygen coordinates below the z_location threshold
          o_indices = np.array([a.index for a in chunk.topology.atoms if a.element.symbol == 'O' and a.residue.name == 'HOH' and chunk.xyz[0, a.index, 2] <= z_location + 0.3])
          o_coords = chunk.xyz[0, o_indices, :]
  
          # Create slab grid coordinates
          x = np.array(range(0, 51)) * 0.1
          y = np.array(range(0, 51)) * 0.1
          X, Y = np.meshgrid(x, y, indexing='xy')
          Z = np.full(X.ravel().shape, z_location+0.3)
          slab_coords = np.column_stack((X.ravel(), Y.ravel(), Z))
  
          # Use cKDTree for efficient nearest-neighbor search
          tree = cKDTree(o_coords)
          neighbors_count = tree.query_ball_point(slab_coords, r=0.1)
  
          count_water = np.array([len(neighbors) for neighbors in neighbors_count])
          whole_water.append(count_water)
  
          print(f'   folder: {folder}, frame number: {i + 1}, length: {len(whole_water)}')
  
      # Compute the average number of water molecules at each slab coordinate
      average_water = np.mean(whole_water, axis=0)
  
      # Create the result DataFrame
      coordinate_water = pd.DataFrame({
          'X': X.ravel(),
          'Y': Y.ravel(),
          'Number': average_water
      })
  
      file_path = OUTPUT_FILE_PATH
      output_filename = f'{folder}_Water.csv'
      output_file = os.path.join(file_path, output_filename)
  
      # Save to CSV
      coordinate_water.to_csv(output_file, index=False)
