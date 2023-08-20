# =============================================
# Author: Benjamin R. Gilbert
# Email: brg4@illinois.edu
# =============================================

import numpy as np

import time
import os

import importlib
import sys

#sys.path.insert(1,'./LAMMPS_helpers')
#sys.path.append('./LAMMPS_helpers/')
import LAMMPS_helpers.read_LAMMPS_DNA as r_L_D
importlib.reload(r_L_D)

in_dir = '/home/ben/Data/btree_chromo/disentanglement_testing/volume_expansion/'
in_label = 'volume_expansion'
#in_dir = '/home/ben/Data/btree_chromo/loop_compaction_testing_new/loops_25/'
#in_label = 'loops_25'
#in_dir = '/home/ben/Data/btree_chromo/small_system_Jan/'
#in_label = 'woloops_wotopo'
# in_dir = '/home/ben/Data/btree_chromo/diffusion_testing_v4/'
# in_label = 'diffusion_bond_bending_twisting_v4'
#in_dir = '/home/ben/Data/btree_chromo/full_chromo_movie/'
#in_label = 'full_chromo_run0'

out_dir = '/home/ben/Data/btree_chromo/disentanglement_testing/volume_expansion/'
out_label = 'volume_expansion'
# out_dir = '/home/ben/Data/btree_chromo/loop_compaction_testing_new/loops_25/'
# out_label = 'loops_25'
# out_dir = '/home/been/Data/btree_chromo/small_system_Jan/'
# out_label = 'woloops_wotopo'
# out_dir = '/home/ben/Data/btree_chromo/diffusion_testing_v4/'
# out_label = 'diffusion_bond_bending_twisting_v4'
#out_dir = '/home/ben/Data/btree_chromo/full_chromo_movie/'
#out_label = 'full_chromo_run0'

rep_specification = 'range'

if rep_specification == 'range':
    min_rep = 1
    max_rep = 5
    reps = np.arange(min_rep,max_rep+1,dtype=np.int32)
elif rep_specification == 'manual':
    reps = [43]
    reps = np.array(reps,dtype=np.int32)

N_reps = reps.shape[0]
    

timestep_options = ['Last','FirstLast','Slices','All']
timestep_select = timestep_options[3]
slices = [0,50000,100000,150000,200000,250000]


if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

for i_rep in range(N_reps):

    rep = reps[i_rep]
    
    rep_label = '_rep' + str(rep).zfill(5)

    in_file = in_dir + in_label + rep_label  + '.lammpstrj'

    print('LAMMPS trajectory: ' + in_file)

    traj = r_L_D.read_traj(in_file,in_timesteps=timestep_select,slices=slices)

    out_file = out_dir + out_label + rep_label + '.pkl_traj'

    print('Pickle trajectory: ' + out_file)
    
    r_L_D.write_pickle_traj(out_file,traj)
