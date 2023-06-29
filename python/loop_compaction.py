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
import chromosome_structure.chromosome_structure as c_s
importlib.reload(c_s)

conditions = np.array([0,25,50,100],dtype=np.int32)
N_conditions = conditions.shape[0]
in_base_dir = '/home/ben/Data/btree_chromo/loop_compaction_testing_new/'
in_base_label = 'loops_'
min_rep = 1
max_rep = 1

out_dir = '/home/ben/Documents/svn/Minimal_Cell_Chromosome_Organization_2022/tex_figures/looping_summary/raw_plots/'
out_label = 'loop_compaction_v2'

write_RoG = True
read_RoG = True

if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

RoG_file = out_dir + out_label + '_RoG.pickle'

ts_RoG = np.arange(0,4000000+1,50000,dtype=np.int32)
print(ts_RoG)
freq = 1
window = 100

if write_RoG:

    RoG = c_s.new_RoG(N_conditions,conditions,
                      max_rep-min_rep+1,
                      ts_RoG.shape[0],ts_RoG,
                      freq,window)

    for i_condition in range(N_conditions):

        in_file_temp = in_base_dir + in_base_label + str(conditions[i_condition]) + '/'
        in_file_temp += (in_base_label + str(conditions[i_condition]))

        for rep in range(min_rep,max_rep+1):

            # read the trajectory

            rep_label = '_rep' + str(rep).zfill(5)

            in_file = in_file_temp + rep_label + '.pkl_traj'

            print('LAMMPS trajectory: ' + in_file)

            traj = r_L_D.read_pickle_traj(in_file)

            RoG = c_s.fill_conditions_reps_RoG(RoG,traj,i_condition,rep-min_rep)
        
    c_s.write_RoG(RoG_file,RoG)

if read_RoG == True:

    RoG = c_s.read_RoG(RoG_file)
            
RoG_file = out_dir + out_label + '_RoG.png'

c_s.plot_RoG(RoG_file,RoG)
        
        
        
