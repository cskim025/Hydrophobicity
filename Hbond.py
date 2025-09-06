from __future__ import print_function
import matplotlib.pyplot as plt
import itertools
import mdtraj as md
import os
import numpy as np
import pandas as pd

# Modified version of code written by Dr. Dallin
def _compute_bounded_geometry( traj,
                                triplets,
                                distance_indices = [ 0, 2 ],
                                angle_indices = [ 1, 0, 2 ] ):
                                  
        # Calculate the requested distances
        distances = md.compute_distances( traj,
                                            triplets[ :, distance_indices ],
                                            periodic = True )

        # Calculate angles using the law of cosines
        abc_pairs = zip( angle_indices, angle_indices[1:] + angle_indices[:1] )
        abc_distances = []

        # calculate distances (if necessary)
        for abc_pair in abc_pairs:
            if set( abc_pair ) == set( distance_indices ):
                abc_distances.append( distances )
            else:
                abc_distances.append( md.compute_distances( traj, triplets[ :, abc_pair ], ) )

        a, b, c = abc_distances
        cosines = ( a ** 2 + b ** 2 - c ** 2 ) / ( 2 * a * b )
        np.clip(cosines, -1, 1, out=cosines) # avoid NaN error
        angles = np.arccos(cosines)

        return distances, angles

def _get_bond_triplets( traj, z_value, polar_name, acc_or_don):
    def get_donors(e0, e1):
        # Find all matching bonds
        elems = set((e0, e1))

        if acc_or_don == 'N' or acc_or_don == 'O':
            atoms = [ (one, two) for one, two in traj.topology.bonds
                        if set((one.element.symbol, two.element.symbol)) == elems and one.residue.name == 'HOH' and
                    (traj.xyz[0, one.index, 2] <= z_value)]
        else:
            atoms = [ (one, two) for one, two in traj.topology.bonds
                    if set((one.element.symbol, two.element.symbol)) == elems]

        # Get indices for the remaining atoms
        indices = []
        for a0, a1 in atoms:
            pair = (a0.index, a1.index)
            # make sure to get the pair in the right order, so that the index
            # for e0 comes before e1
            if a0.element.symbol == e1:
                pair = pair[::-1]
            indices.append(pair)

        return indices

    #nh_donors = get_donors('N', 'H')
    if acc_or_don == 'N' or acc_or_don == 'O':
        donors = get_donors('O', 'H')
    else:
        char = list(acc_or_don)
        char_1 = char[0]
        char_2 = char[1]

        donors = get_donors(char_1, char_2)

    xh_donors = np.array(donors)

    if len(xh_donors) == 0:
        # if there are no hydrogens or protein in the trajectory, we get
        # no possible pairs and return nothing
        return np.zeros((0, 3), dtype=int)

    if acc_or_don == 'N' or acc_or_don == 'O':

        acceptor_elements = frozenset((acc_or_don)) #'C' 
        acceptors = [ a.index for a in traj.topology.atoms
                        if a.element.symbol in acceptor_elements and a.residue.name == polar_name]
    else:
        acceptor_elements = frozenset(('O')) #'C' 
        acceptors = [ a.index for a in traj.topology.atoms
                    if a.element.symbol in acceptor_elements and one.residue.name == 'HOH' and 
                    (traj.xyz[0, one.index, 2] <= z_value)]

    # Make acceptors a 2-D numpy array
    acceptors = np.array(acceptors)[:, np.newaxis]

    # Generate the cartesian product of the donors and acceptors
    xh_donors_repeated = np.repeat(xh_donors, acceptors.shape[0], axis=0)
    acceptors_tiled = np.tile(acceptors, (xh_donors.shape[0], 1))
    bond_triplets = np.hstack((xh_donors_repeated, acceptors_tiled))

    # Filter out self-bonds
    self_bond_mask = (bond_triplets[:, 0] == bond_triplets[:, 2])
    return bond_triplets[np.logical_not(self_bond_mask), :]
    

