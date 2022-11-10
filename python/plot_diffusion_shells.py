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

in_dir = '/home/ben/Data/btree_chromo/diffusion_testing/'
in_label = 'diffusion_no_interactions'
min_rep = 1
max_rep = 10

out_dir = '/home/ben/Documents/svn/Minimal_Cell_Chromosome_Organization_2022/tex_figures/diffusion/raw_plots/'
out_label = 'diffusion_no_interactions'

timestep_options = ['Last','FirstLast','Slices','All']
timestep_select = timestep_options[3]
slices = [0,50000,100000,150000,200000,250000]

write_shell_diff = False
read_shell_diff = True

if not os.path.isdir(out_dir) and write_shell_diff:
    os.makedirs(out_dir)

shell_diff_file = out_dir + out_label + '_shell_diff.pickle'

dt = 0.1 # timestep - [ns]

# radii of shells - [A]
radii = np.array([1000,2000],dtype=np.int32)

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

else:

    quit()

shell_MSD = d_plot.shell_MSD_calc_diffusion(shell_MSD,dt,[50000,100000])
    

shell_MSD_fig_file = out_dir + out_label + '_shellMSD.pdf'

d_plot.plot_shell_MSD_reps(shell_MSD_fig_file,
                           shell_MSD,
                           dt,
                           [50,100],
                           'DNA')

d_plot.plot_shell_MSD_reps(shell_MSD_fig_file,
                           shell_MSD,
                           dt,
                           [50,100],
                           'ribo')

shell_law_fig_file = out_dir + out_label + '_law.pdf'

d_plot.plot_shell_law_reps(shell_law_fig_file,
                           shell_MSD,
                           dt,
                           [50,100],
                           'DNA')

d_plot.plot_shell_law_reps(shell_law_fig_file,
                           shell_MSD,
                           dt,
                           [50,100],
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
