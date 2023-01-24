# =============================================
# Author: Benjamin R. Gilbert
# Email: brg4@illinois.edu
# =============================================

import numpy as np

import time
import os

import importlib
import sys
import glob

#sys.path.insert(1,'./LAMMPS_helpers')
#sys.path.append('./LAMMPS_helpers/')
import LAMMPS_helpers.read_LAMMPS_DNA as r_L_D
importlib.reload(r_L_D)
import LAMMPS_helpers.fork_partitions as f_p
importlib.reload(f_p)
import segregation_calculations.disentanglement_plotting as de_p
importlib.reload(de_p)

in_dir = '/home/ben/Data/btree_chromo/partition_testing/'
in_label = 'partition_test'
f_p_label = 'test_f_p'
min_rep = 1
max_rep = 3

t_d = np.arange(500000,2500000+1,50000,dtype=np.int32)
R = 4.0*17.0
t_d = np.sort(t_d)
print(t_d)

out_dir = '/home/ben/Data/btree_chromo/partition_testing/'
out_label = 'partition_test'

if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

rep_count = 0

for rep in range(min_rep,max_rep+1):

    # read the trajectory

    rep_label = '_rep' + str(rep).zfill(5)

    in_file = in_dir + in_label + rep_label + '.pkl_traj'

    print('LAMMPS trajectory: ' + in_file)

    traj = r_L_D.read_pickle_traj(in_file)

    # get the set of fork partition files

    f_p_files = in_dir + f_p_label + rep_label + '_*.dat'

    f_p_files = glob.glob(f_p_files)

    ts_f_ps = dict()

    ts_f_ps['ts'] = np.zeros((len(f_p_files)),dtype=np.int32)
    ts_f_ps['fork_partitions'] = []

    for i in range(len(f_p_files)):

        t = f_p_files[i].split('_')
        t = t[-1]
        t = t.split('.')
        t = t[0]
        t = t.strip('t')

        ts_f_ps['ts'][i] = int(t)

        f_ps = f_p.fork_partitions()
        ts_f_ps['fork_partitions'].append(f_ps)

    ts_arg_sorted = np.argsort(ts_f_ps['ts'])

    for i_t in range(ts_arg_sorted.shape[0]):

        ts_f_ps['fork_partitions'][i_t].parse_f_p_file(f_p_files[ts_arg_sorted[i_t]])

    ts_f_ps['ts'] = np.sort(ts_f_ps['ts'])

    print(ts_f_ps['ts'])
    print(ts_f_ps['fork_partitions'])
    
    for i_t in range(t_d.shape[0]):

        i_t_m = t_d.shape[0] - 1 - i_t
        
        ts_f_p_selected = np.argwhere(ts_f_ps['ts']<t_d[i_t_m])[-1][0]

        f_ps = ts_f_ps['fork_partitions'][ts_f_p_selected]

        x = traj[str(t_d[i_t_m])]['DNA']['x']
        print(str(t_d[i_t_m]))
        print(x.shape[0])
        
        forks, d_o_d = de_p.disentanglement(f_ps,x,R)
        
        if i_t == 0:
            final_forks = forks
            d_o_d_forks_ts = np.zeros((len(final_forks),t_d.shape[0]),dtype=np.double)

        for i_fork in range(len(forks)):

            for j_fork in range(len(final_forks)):

                if forks[i_fork] == final_forks[j_fork]:

                    d_o_d_forks_ts[j_fork,i_t_m] = d_o_d[i_fork]

    print(d_o_d_forks_ts)

        

        
            

        
        
        
