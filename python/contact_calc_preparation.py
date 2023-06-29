# =============================================
# Author: Benjamin R. Gilbert
# Email: brg4@illinois.edu
# =============================================

import numpy as np

import time
import os

import importlib
import sys

import LAMMPS_helpers.read_LAMMPS_DNA as r_L_D
importlib.reload(r_L_D)

in_dir = '/home/ben/Data/btree_chromo/disentanglement_testing/wextraloops_wtopo/'
in_label = 'wextraloops_wtopo'
min_rep = 1
max_rep = 5

N_stages = 10
i_stage = 9
dt = 50000
T_total = 20000000
T_stage = T_total//N_stages
print(T_stage)

out_dir = '/home/ben/Data/contact_maps/btree_chromo/disentanglement_stages/stage_' + str(i_stage) + '/'
out_label = 'stage_' + str(i_stage)

if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

rep_count = 1

#ts_contacts = np.arange(15000000,20000000+1,250000,dtype=np.int32)
#ts_contacts = np.arange(15000000,20000000+1,50000,dtype=np.int32)
ts_contacts = np.arange(i_stage*T_stage+dt,(i_stage+1)*T_stage+1,dt,dtype=np.int32)
print(ts_contacts)
print(ts_contacts.shape)

write_coords_flag = True

for rep in range(min_rep,max_rep+1):

    rep_label = '_rep' + str(rep).zfill(5)

    in_file = in_dir + in_label + rep_label + '.pkl_traj'

    print('LAMMPS trajectory: '+ in_file)

    if write_coords_flag == True:
        
        traj = r_L_D.read_pickle_traj(in_file)

        for i_t in range(ts_contacts.shape[0]):

            t_target = str(ts_contacts[i_t])

            print(t_target)
    
            N = traj[t_target]['DNA']['N']
            x = traj[t_target]['DNA']['x']

            x = np.reshape(x,(N,3),order='F')

            out_file = out_dir + 'x_' + out_label +\
                '_rep' + str(rep_count).zfill(5) + '.bin'  # output file
            print(out_file)

            with open(out_file,'wb') as f:

                x.tofile(f)

            rep_count += 1
