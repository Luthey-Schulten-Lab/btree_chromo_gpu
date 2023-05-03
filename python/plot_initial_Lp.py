# =============================================
# Author: Benjamin R. Gilbert
# Email: brg4@illinois.edu
# =============================================

import numpy as np

import time
import os

import importlib
import sys

import mechanical_calculations.initial_persistence_length_plotting as i_Lp
importlib.reload(i_Lp)

in_dir = '/home/ben/Data/sc_chains/lp_testing/'
in_label = 'x_chain_lp_test_alt'

MC_steps = [0,1000,10000,50000,100000,250000,500000,1000000]
MC_steps = np.array(MC_steps,dtype=np.int32)
N_MC = MC_steps.shape[0]

min_rep = 1
max_rep = 20
N_reps = max_rep - min_rep + 1

s_min = 0
s_max = 25
Lp_true = 450.0
l_0 = 34.0

out_dir = '/home/ben/Documents/svn/Minimal_Cell_Chromosome_Organization_2022/tex_figures/sc_chain_Lp/raw_plots/'
out_label = 'Lp_MC_normal_alt'

if not os.path.isdir(out_dir):
    os.makedirs(out_dir)


corr_file = out_dir + 'poly_corr_' + out_label + '.pkl'

write_corr_flag = True
read_corr_flag = True

if write_corr_flag == True:

    poly_corr = i_Lp.new_initial_polymer_correlations(s_min,s_max,
                                                      MC_steps,
                                                      N_reps,
                                                      N_MC,
                                                      Lp_true,l_0)
    for i_MC in range(N_MC):

        MC_label = '_' + str(MC_steps[i_MC]) + 'MC'

        rep_count = 0

        for rep in range(min_rep,max_rep+1):

            rep_label = '_rep' + str(rep).zfill(5)

            in_file = in_dir + in_label + MC_label + rep_label + '.bin'

            print(in_file)

            with open(in_file,'rb') as f:

                x = np.fromfile(f,dtype=np.double,count=-1)

            x = x.reshape((3,x.shape[0]//3),order='F').T

            # print(x.shape)

            poly_corr = i_Lp.fill_bond_correlations_circular(poly_corr,
                                                             x,
                                                             i_MC,
                                                             rep_count)

            rep_count += 1



    i_Lp.write_poly_corr(corr_file,poly_corr)

if read_corr_flag == True:

    poly_corr = i_Lp.read_poly_corr(corr_file)


fig_file = out_dir + out_label + '_Lp_plot.png'

i_Lp.plot_initial_Lp(poly_corr,fig_file)
