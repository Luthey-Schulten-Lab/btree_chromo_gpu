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

# in_dir = '/home/ben/Data/btree_chromo/disentanglement_testing/wextraloops_wotopo/'
# in_label = 'wextraloops_wotopo'
in_dir = '/mnt/HDD_linux/Data/btree_chromo/loop_compaction_testing/loops_0/'
in_label = 'loops_0'
min_rep = 1
max_rep = 1

# out_dir = '/home/ben/Data/btree_chromo/disentanglement_testing/wextraloops_wotopo/'
# out_label = 'wextraloops_wotopo'
out_dir = '/mnt/HDD_linux/Data/btree_chromo/loop_compaction_testing/loops_0/'
out_label = 'loops_0'

timestep_options = ['Last','FirstLast','Slices','All']
timestep_select = timestep_options[3]
slices = [0,50000,100000,150000,200000,250000]


if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

rep_count = 0

for rep in range(min_rep,max_rep+1):

    rep_label = '_rep' + str(rep).zfill(5)

    in_file = in_dir + in_label + rep_label  + '.lammpstrj'

    print('LAMMPS trajectory: ' + in_file)

    traj = r_L_D.read_traj(in_file,in_timesteps=timestep_select,slices=slices)

    out_file = out_dir + out_label + rep_label + '.pkl_traj'

    print('Pickle trajectory: ' + out_file)
    
    r_L_D.write_pickle_traj(out_file,traj)
