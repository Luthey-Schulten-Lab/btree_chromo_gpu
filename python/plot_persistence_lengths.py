# =============================================
# Author: Benjamin R. Gilbert
# Email: brg4@illinois.edu
# =============================================

import numpy as np

import time
import os

import imp
import sys

sys.path.insert(1,'./LAMMPS_helpers')
import read_LAMMPS_DNA as r_L_D
imp.reload(r_L_D)
sys.path.insert(1,'./mechanical_calculations')
import persistence_length_plotting as pl_plot
imp.reload(pl_plot)

in_dir = '/home/ben/Data/btree_chromo/diffusion_testing/'
in_label = 'diffusion_bond_bending_twisting'
min_rep = 1
max_rep = 5
N_reps = max_rep - min_rep + 1

out_dir = '/home/ben/Documents/svn/Minimal_Cell_Chromosome_Organization_2022/tex_figures/DNA_mechanical_properties/raw_plots/'
out_label = 'persistence_lengths'

s_min = 1
s_max = 10

if not os.path.isdir(out_dir):
    os.makedirs(out_dir)


for rep in range(min_rep,max_rep+1):

    in_file = in_dir + in_label + '_rep'+str(rep).zfill(5) + '.pkl_traj'

    print('LAMMPS trajectory: '+str(in_file))

    traj = r_L_D.read_pickle_traj(in_file)
    
    if (rep == min_rep):
        
        disp_corr = pl_plot.new_polymer_correlations(s_min,
                                                     s_max,
                                                     N_reps,
                                                     traj['timesteps'].shape[0])

    pl_plot.fill_displacement_correlations_circular(disp_corr,
                                                    traj,
                                                    rep-min_rep)

    print(disp_corr['corr'][:,0,0])
