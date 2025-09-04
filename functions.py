from __future__ import print_function
import matplotlib.pyplot as plt
import itertools
import mdtraj as md
import os
import numpy as np
import pandas as pd

class HydroAnaly:
  
  def _compute_bounded_geometry( traj,
                                  triplets,
                                  distance_indices = [ 0, 2 ],
                                  angle_indices = [ 1, 0, 2 ] ):
          '''this function computes the distances between the atoms involved in
          the hydrogen bonds and the H-donor...acceptor angle using the law of
          cosines.
  
          Inputs
          ------
          traj : md.traj
          triplets : np.array, shape[n_possible_hbonds, 3], dtype=int
              An array containing the indices of all possible hydrogen bonding triplets
          distance_indices : [LIST], [ donor_index, acceptor_index ], default = [ 0, 2 ]
              A list containing the position indices of the donor and acceptor atoms
          angle_indices : [LIST], [ h_index, donor_index, acceptor_index ], default = [ 1, 0, 2 ]
              A list containing the position indices of the H, donor, and acceptor
              atoms. Default is H-donor...acceptor angle
  
          Outputs
          -------
          distances : np.array, shape[n_possible_hbonds, 1], dtype=float
              An array containing the distance between the donor and acceptor atoms
          angles : np.array, shape[n_possible_hbonds, 1], dtype=float
              An array containing the triplet angle between H-donor...acceptor atoms
          '''
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
  
          # Law of cosines calculation to find the H-Donor...Acceptor angle
          #            c**2 = a**2 + b**2 - 2*a*b*cos(C)
          #                        acceptor
          #                          /\
          #                         /  \
          #                      c /    \ b
          #                       /      \
          #                      /______(_\
          #                     H    a     donor
          a, b, c = abc_distances
          cosines = ( a ** 2 + b ** 2 - c ** 2 ) / ( 2 * a * b )
          np.clip(cosines, -1, 1, out=cosines) # avoid NaN error
          angles = np.arccos(cosines)
  
          return distances, angles
