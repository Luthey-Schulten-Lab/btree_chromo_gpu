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
sys.path.insert(1,'./diffusion_calculations')
import diffusion_plotting as d_plot
imp.reload(d_plot)

in_dir = '/home/ben/Data/btree_chromo/diffusion_testing_v4/'
in_label = 'diffusion_bond_bending_twisting_v4'

#out_dir = '/home/ben/Documents/svn/Minimal_Cell_Chromosome_Organization_2022/tex_figures/diffusion/raw_plots/'
out_dir = '/home/ben/Data/btree_chromo/diffusion_testing_v4/'
out_label = 'diffusion_bond_bending_twisting_v4'

rep_specification = 'range'

if rep_specification == 'range':
    min_rep = 1
    max_rep = 50
    reps = np.arange(min_rep,max_rep+1,dtype=np.int32)
elif rep_specification == 'manual':
    reps = [1,2,3,4,5,6,8,9,10,11,12,13,15,16,17,18,19,20,22,23,24,25,26,27,29,30,31,32,33,34,36,37,38,39,40,41,43,44,45,46,47,48,49]
    reps = np.array(reps,dtype=np.int32)

N_reps = reps.shape[0]

fig_filetype = ".png"

write_shell_diff = True
read_shell_diff = True

if not os.path.isdir(out_dir) and write_shell_diff:
    os.makedirs(out_dir)

shell_diff_file = out_dir + out_label + '_shell_diff.pickle'

dt = 0.1 # timestep - [ns]

# radii of shells - [A]
radii = np.array([1500,2000],dtype=np.int32)

if write_shell_diff == True:

    rep_count = 0

    for i_rep in range(N_reps):

        rep = reps[i_rep]

        in_file = in_dir + in_label + '_rep'+str(rep).zfill(5) + '.pkl_traj'

        print('LAMMPS trajectory: '+str(in_file))

        traj = r_L_D.read_pickle_traj(in_file)

        if i_rep == 0:

            shell_MSD = d_plot.new_shell_MSD(radii,
                                             N_reps,
                                             traj['timesteps'].shape[0],
                                             traj['timesteps'])

        shell_MSD = d_plot.fill_reps_shell_MSDs(shell_MSD,traj,i_rep)

    d_plot.write_shell_diff(shell_diff_file,shell_MSD)

if read_shell_diff == True:

    shell_MSD = d_plot.read_shell_diff(shell_diff_file)

    
shell_MSD = d_plot.shell_MSD_calc_diffusion(shell_MSD,dt,[100,200])
    

shell_MSD_fig_file = out_dir + out_label + '_shellMSD'

d_plot.plot_shell_MSD_reps(shell_MSD_fig_file,
                           fig_filetype,
                           shell_MSD,
                           dt,
                           [50,100],
                           'DNA')

d_plot.plot_shell_MSD_reps(shell_MSD_fig_file,
                           fig_filetype,
                           shell_MSD,
                           dt,
                           [50,100],
                           'ribo')

shell_law_fig_file = out_dir + out_label + '_law'

d_plot.plot_shell_law_reps(shell_law_fig_file,
                           fig_filetype,
                           shell_MSD,
                           dt,
                           [50,100],
                           'DNA')

d_plot.plot_shell_law_reps(shell_law_fig_file,
                           fig_filetype,
                           shell_MSD,
                           dt,
                           [50,100],
                           'ribo')

shell_diff_vs_density_fig_file = out_dir + out_label + '_DiffvDensity'

d_plot.plot_shell_diff_vs_density_reps(shell_diff_vs_density_fig_file,
                                       fig_filetype,
                                       shell_MSD,
                                       'DNA')

# d_plot.plot_shell_diff_vs_density_reps(shell_diff_vs_density_fig_file,
#                                        fig_filetype,
#                                        shell_MSD,
#                                        'ribo')

shell_law_vs_density_fig_file = out_dir + out_label + '_LawvDensity'

d_plot.plot_shell_law_vs_density_reps(shell_law_vs_density_fig_file,
                                      fig_filetype,
                                      shell_MSD,
                                      'DNA')

# d_plot.plot_shell_law_vs_density_reps(shell_law_vs_density_fig_file,
#                                       fig_filetype,
#                                       shell_MSD,
#                                       'ribo')

shell_diff_vs_diff_fig_file = out_dir + out_label + '_DiffvDiff'

d_plot.plot_shell_diff_vs_diff_reps(shell_diff_vs_diff_fig_file,
                                    fig_filetype,
                                    shell_MSD)

shell_law_vs_law_fig_file = out_dir + out_label + '_LawvLaw'

d_plot.plot_shell_law_vs_law_reps(shell_law_vs_law_fig_file,
                                  fig_filetype,
                                  shell_MSD)
