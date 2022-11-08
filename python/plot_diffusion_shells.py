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

# in_dir = '/home/ben/Data/LAMMPS/5000mono_testcase/5000mono_init/'
# in_label = '5000mono_init'
in_dir = '/home/ben/Data/LAMMPS/s1c15_503ribos/s1c15_503ribos_init/'
in_label = 's1c15_503ribos_init'
min_rep = 1
max_rep = 50

# out_dir = '/home/ben/Data/LAMMPS/5000mono_testcase/5000mono_init/'
# out_label = '5000mono_init'
out_dir = '/home/ben/Data/LAMMPS/s1c15_503ribos/s1c15_503ribos_init/'
out_label = 's1c15_503ribos_init'

timestep_options = ['Last','FirstLast','Slices','All']
timestep_select = timestep_options[3]
slices = [0,50000,100000,150000,200000,250000]

write_shell_diff = True
read_shell_diff = False

if not os.path.isdir(out_dir) and write_flag:
    os.makedirs(out_dir)

shell_diff_file = out_dir + out_label + '_shell_diff.pickle'

radii = np.array([500,1000,1500,2000],dtype=np.int32)

if write_shell_diff == True:

    rep_count = 0

    for rep in range(min_rep,max_rep+1):

        in_file = in_dir + in_label + '_rep'+str(rep).zfill(5) + '.lammpstrj'

        print('LAMMPS trajectory: '+str(in_file))

        traj = r_L_D.read_traj(in_file,in_timesteps=timestep_select,slices=slices)

        if rep == min_rep:

            shell_MSD = d_plot.new_shell_MSD(radii,
                                             max_rep-min_rep+1,
                                             traj['timesteps'].shape[0],
                                             traj['timesteps'])

        shell_MSD = d_plot.fill_reps_shell_MSDs(shell_MSD,traj,rep-min_rep)

    #print(shell_MSD)

    d_plot.write_shell_diff(shell_diff_file,shell_MSD)

if read_shell_diff == True:

    shell_MSD = d_plot.read_shell_diff(shell_diff_file)

shell_MSD = d_plot.shell_MSD_calc_diffusion(shell_MSD,0.1,[25000,51000])
    

shell_MSD_fig_file = out_dir + out_label + '_shellMSD.pdf'

d_plot.plot_shell_MSD_reps(shell_MSD_fig_file,
                           shell_MSD,
                           0.1,
                           [500,1000],
                           'DNA')

d_plot.plot_shell_MSD_reps(shell_MSD_fig_file,
                           shell_MSD,
                           0.1,
                           [500,1000],
                           'ribo')

shell_law_fig_file = out_dir + out_label + '_law.pdf'

d_plot.plot_shell_law_reps(shell_law_fig_file,
                           shell_MSD,
                           0.1,
                           [500,1000],
                           'DNA')

d_plot.plot_shell_law_reps(shell_law_fig_file,
                           shell_MSD,
                           0.1,
                           [500,1000],
                           'ribo')

shell_diff_vs_density_fig_file = out_dir + out_label + '_DiffvDensity.pdf'

d_plot.plot_shell_diff_vs_density_reps(shell_diff_vs_density_fig_file,
                           shell_MSD,
                           'DNA')

d_plot.plot_shell_diff_vs_density_reps(shell_diff_vs_density_fig_file,
                                       shell_MSD,
                                       'ribo')

shell_law_vs_density_fig_file = out_dir + out_label + '_LawvDensity.pdf'

d_plot.plot_shell_law_vs_density_reps(shell_law_vs_density_fig_file,
                           shell_MSD,
                           'DNA')

d_plot.plot_shell_law_vs_density_reps(shell_law_vs_density_fig_file,
                                       shell_MSD,
                                       'ribo')

shell_diff_vs_diff_fig_file = out_dir + out_label + '_DiffvDiff.pdf'

d_plot.plot_shell_diff_vs_diff_reps(shell_diff_vs_diff_fig_file,
                                    shell_MSD)

shell_law_vs_law_fig_file = out_dir + out_label + '_LawvLaw.pdf'

d_plot.plot_shell_law_vs_law_reps(shell_law_vs_law_fig_file,
                                  shell_MSD)
