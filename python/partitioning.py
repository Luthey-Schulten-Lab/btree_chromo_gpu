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
import segregation_calculations.partition_plotting as par_plot
importlib.reload(par_plot)

in_dir = '/home/ben/Data/btree_chromo/disentanglement_testing/volume_expansion/'
in_label = 'volume_expansion'
f_p_label = 'volume_expansion_f_p'
min_rep = 1
max_rep = 5

out_dir = '/home/ben/Data/btree_chromo/disentanglement_testing/volume_expansion/raw_plots/'
out_label = 'volume_expansion'

write_par_CoM = True
read_par_CoM = True

if not os.path.isdir(out_dir):
    os.makedirs(out_dir)

par_CoM_file = out_dir + out_label + '_par_CoM.pickle'

ts_par_CoM = np.arange(500000,20000000+1,250000,dtype=np.int32)
print(ts_par_CoM)
R = 4.0*17.0
ts_par_CoM = np.sort(ts_par_CoM)

R_bdry = 900.0

if write_par_CoM:

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

            par_CoM = par_plot.new_par_CoM(R,
                                  max_rep-min_rep+1,
                                  ts_par_CoM.shape[0],
                                  ts_par_CoM)

        par_CoM = par_plot.fill_reps_par_CoM(par_CoM,
                                    traj,
                                    f_p_files,
                                    rep-min_rep)
        
    # par_CoM['d_reps_forks_ts'] = np.where(par_CoM['d_reps_forks_ts']>0.5,
    #                                   2.0*(par_CoM['d_reps_forks_ts']-0.5),
    #                                   0.0)
        
    par_plot.write_par_CoM(par_CoM_file,par_CoM)

if read_par_CoM == True:

    par_CoM = par_plot.read_par_CoM(par_CoM_file)
            
#par_CoM_file = out_dir + out_label + '_par_CoM.pdf'
par_CoM_file = out_dir + out_label + '_par_CoM.png'

par_plot.plot_par_CoM_scaled(par_CoM_file,par_CoM,R_bdry)
        
        
        
