## This file is to use in Linux and High Performance Computing environment using Slurm 


from __future__ import print_function
import matplotlib.pyplot as plt
import itertools
import mdtraj as md
import os
import numpy as np
import pandas as pd
from optparse import OptionParser

class number_counter_z:

    def __init__(self, in_file = None, in_path = None, out_path = None, cal_time = 4):

        self.in_file = in_file
        self.in_path = in_path
        self.out_path = out_path
        self.cal_time = cal_time

    def compute( self ):

        traj = self.load_md_trajactory( self.in_path, self.in_file )
        cal_time = self.cal_time

        SOL_num_here, MOH_num_here = self.number_counting(traj, cal_time)

        return SOL_num_here, MOH_num_here

    def check_path(path):

        import os
        
        if not os.path.exists(path):
            raise RuntimeError('\n %s does not exist. Check your input' %(path))

    def load_md_trajactory(self, path, filenames):
        
        if path[-1] == '/':
            xtc_file = path + filenames + '.gro'
            gro_file = path + filenames + '.xtc'
        else:
            xtc_file = path + '/' + filenames + '.gro'
            gro_file = path + '/' + filenames + '.xtc'
        
        self.check_path(xtc_file)
        self.check_path(gro_file)

        traj = md.load(xtc_file, top = gro_file)

        return traj

    def number_counting(traj, cal_time):
        
        # Get z-dimension of the box
        box_vectors = traj.unitcell_lengths
        box_height = np.mean(box_vectors[:, 2])

        # Bin parameters
        bin_width = 0.02
        n_bins = int(np.ceil(box_height / bin_width))
        z_centers = np.linspace(bin_width / 2, box_height - bin_width / 2, n_bins)

        # Time slicing
        START = cal_time * 1000
        STOP = traj.n_frames
        STRIDE = 100
        chunk = traj[START:STOP:STRIDE]

        # Atom selection
        table, bonds = traj.topology.to_dataframe()
        SOL_index = table[table['resName'] == 'HOH'].index

        # Get coordinates and compute COMs
        coords = chunk.xyz[:, SOL_index, :]
        n_frames, n_atoms, _ = coords.shape
        assert n_atoms % 3 == 0

        coords_grouped = coords.reshape(n_frames, -1, 3, 3)
        masses = np.array([15.999, 1.008, 1.008])
        mass_sum = np.sum(masses)
        SOL_com = np.sum(coords_grouped * masses[None, None, :, None], axis=2) / mass_sum  # shape: (n_frames, n_molecules, 3)

        # Bin COMs along z-axis
        hist_all = np.zeros((n_frames, n_bins), dtype=int)

        for i in range(n_frames):
            z_vals = SOL_com[i, :, 2]  # z-coordinates of all COMs in frame i
            hist, _ = np.histogram(z_vals, bins=n_bins, range=(0, box_height))
            hist_all[i] = hist  # store per-frame histogram

        avg_density = np.mean(hist_all, axis=0)

        return avg_density

            
if __name__ == "__main__":
    
    use = 'Usage: %prog [options]'
    parser = OptionParser(usage=use)
    parser.add_option('-f', '--in', dest='infile', action='store', type='string', help='input_file', default='.')
    parser.add_option('-w', '--wd', dest='inpath', action='store', type='string', help='working directory', default='.')
    parser.add_option('-p', '--op', dest='outpath', action='store', type='string', help='output_directory', default='.')
    parser.add_option('-t', '--time', dest='caltime', action='store', type='string', help='number', default='4')

    (options, args) = parser.parse_args()
    args = { 'in_file': options.infile,
             'in_path': options.inpath,
             'out_path': options.outpath,
             'cal_time': options.caltime}
    
    out_file = args['in_file'].split('_')[0]
    counting = number_counter( **args )
    SOL_number = counting.compute()
    
    SOL_number = pd.DataFrame(SOL_number)

    SOL_number.to_csv("SOL_num.csv", header = False, index= False)
