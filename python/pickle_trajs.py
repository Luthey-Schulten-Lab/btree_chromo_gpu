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

in_dir = '/home/ben/Data/btree_chromo/diffusion_testing/'
in_label = 'diffusion_no_interactions'
min_rep = 45
max_rep = 50

out_dir = '/home/ben/Data/btree_chromo/diffusion_testing/'
out_label = 'diffusion_no_interactions'

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
