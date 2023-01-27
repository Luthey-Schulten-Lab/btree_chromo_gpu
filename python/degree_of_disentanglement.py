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
import segregation_calculations.disentanglement_plotting as de_plot
importlib.reload(de_plot)

in_dir = '/home/ben/Data/btree_chromo/partition_testing/'
in_label = 'partition_noloops_test'
f_p_label = 'test_noloops_f_p'
min_rep = 1
max_rep = 5

out_dir = '/home/ben/Documents/svn/Minimal_Cell_Chromosome_Organization_2022/tex_figures/disentanglement/raw_plots/'
out_label = 'noloops'

write_DoD = True
read_DoD = True

if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

DoD_file = out_dir + out_label + '_DoD.pickle'

ts_DoD = np.arange(500000,2500000+1,50000,dtype=np.int32)
R = 4.0*17.0
ts_DoD = np.sort(ts_DoD)

if write_DoD:

    rep_count = 0

    for rep in range(min_rep,max_rep+1):

        # read the trajectory

        rep_label = '_rep' + str(rep).zfill(5)

        in_file = in_dir + in_label + rep_label + '.pkl_traj'

        print('LAMMPS trajectory: ' + in_file)

        traj = r_L_D.read_pickle_traj(in_file)

        # get the set of fork partition files

        f_p_files = in_dir + f_p_label + rep_label + '_*.dat'
      
        if rep == min_rep:

            DoD = de_plot.new_DoD(R,
                                  max_rep-min_rep+1,
                                  ts_DoD.shape[0],
                                  ts_DoD)

        DoD = de_plot.fill_reps_DoD(DoD,
                                    traj,
                                    f_p_files,
                                    rep-min_rep)
        
    DoD['d_reps_forks_ts'] = np.where(DoD['d_reps_forks_ts']>0.5,
                                      2.0*(DoD['d_reps_forks_ts']-0.5),
                                      0.0)
        
    de_plot.write_DoD(DoD_file,DoD)

if read_DoD == True:

    DoD = de_plot.read_DoD(DoD_file)
            
DoD_file = out_dir + out_label + '_DoD.pdf'

de_plot.plot_DoD(DoD_file,DoD)
        
        
        
