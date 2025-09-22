import os
import mdtraj as md
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


folders = ['FILE_1', 'FILE_2']

def angle_calc(FILE_PATH, FILE_NAME, folders, START_TIME, OUTFILE_NAME, local_density):
  
  # Initialize an empty DataFrame for hexagonal factors
  hexagonal_factor_df = pd.DataFrame(columns=folders)
  
  for folder in folders:
      xtc_file = os.path.join(FILE_PATH, folder, FILE_NAME+'.xtc')
      gro_file = os.path.join(FILE_PATH, folder, FILE_NAME+'.gro')

      local_density_data = pd.DataFrame(local_density)
      
      # Load trajectory
      t = md.load(xtc_file, top=gro_file)
      table, bonds = t.topology.to_dataframe()
  
      molecule_list = sorted(list(set(table['resName'])))
      molecule_list.remove("HOH")
  
      polar_name = molecule_list[0]
      polar_table = table[table["resSeq"] == min(table.loc[table["resName"] == molecule_list[0], "resSeq"])]
      polar_carbon_index = polar_table.loc[max(polar_table[polar_table["element"] == "C"].index), "name"]
      polar_carbon_indices = chunk.topology.select(f"(name == {polar_carbon_index}) and (resname == {polar_name})")
  
      nonpolar_name = molecule_list[1]
      nonpolar_table = table[table["resSeq"] == min(table.loc[table["resName"] == molecule_list[1], "resSeq"])]
      nonpolar_carbon_index = nonpolar_table.loc[max(nonpolar_table[nonpolar_table["element"] == "C"].index), "name"]
      nonpolar_carbon_indices = chunk.topology.select(f"(name == {nonpolar_carbon_index}) and (resname == {nonpolar_name})")
  
      selected_indices = np.array([*polar_carbon_indices, *nonpolar_carbon_indices])
  
      # Extract positions of selected atoms
      selected_positions = t.xyz[:, selected_indices, :]  # Shape: (n_frames, n_selected_atoms, 3)
  
      # Compute average position across all frames (mean of all frames)
      avg_positions = selected_positions[START_TIME:] # Shape: (n_selected_atoms, 3)
      
      # Perform Voronoi tessellation for the last frame
      avg_positions_2d = avg_positions[:, :, :2]
      frame_positions_2d = np.mean(avg_positions_2d, axis=0)
      vor = Voronoi(frame_positions_2d)
  
      # Plot Voronoi
      fig, ax = plt.subplots(figsize=(5, 4.7))
      voronoi_plot_2d(vor, ax=ax, show_vertices=False, line_colors='blue', linewidth=1.5)
  
      # Extract Voronoi bounds
      vor_min_x, vor_max_x = 1, 4
      vor_min_y, vor_max_y = 1, 4
  
      # Extract hexbin bounds from data
      data_min_x, data_max_x = 1, 4
      data_min_y, data_max_y = 1, 4
  
      # Set combined axis limits with padding
      padding = 0.1
  
      ax.set_xlim(1, 4)
      ax.set_ylim(1, 4)
  
      # Extract x, y, and density values from the DataFrame
      x = local_density_data['X']
      y = local_density_data['Y']
      density_values = local_density_data['Number'] / 0.12887112887112886 #3.662466876656167
  
      # Plot hexbin for local density data with dynamic bounds
      hb = ax.hexbin(x, y, C=density_values, gridsize=30, cmap='viridis', linewidths=0.5, alpha=0.5,
                      vmin=0.4, vmax=0.6)
  
      # Add a color bar for the density values
      cbar = plt.colorbar(hb, ax=ax)
      cbar.set_label('Local Density Value', fontsize=14)
  
      # Plot Voronoi centers with polar (blue) and nonpolar (red) markers
      ax.scatter(frame_positions_2d[:len(polar_carbon_indices), 0], frame_positions_2d[:len(polar_carbon_indices), 1], marker='o', color='blue', zorder=3, s=20, alpha=1)
      ax.scatter(frame_positions_2d[len(polar_carbon_indices):, 0], frame_positions_2d[len(polar_carbon_indices):, 1], marker='o', color='red', zorder=3, s=20, alpha=1)
  
      # Customize axis labels
      ax.set_xlabel('X Position (nm)', fontsize=20)
      ax.set_ylabel('Y Position (nm)', fontsize=20)
  
      # Final display
      plt.tight_layout()
  
      # Specify output location
      output_folder = FILE_PATH
      output_folder = os.path.normpath(output_folder)
      output_file = os.path.join(output_folder, OUTFILE_NAME+".png")
  
      plt.savefig(output_file, dpi=500, bbox_inches='tight', format='png')  # High resolution PNG
  
      print(f"Figure saved as '{output_file}'")
  
      plt.show()
