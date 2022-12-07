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
sys.path.insert(1,'./RDF_calculations')
import RDF_plotting as rdf_plot
imp.reload(rdf_plot)

in_dir = '/home/ben/Data/btree_chromo/diffusion_testing/'
in_label = 'diffusion_bond_bending_twisting'
min_rep = 1
max_rep = 4

out_dir = '/home/ben/Documents/svn/Minimal_Cell_Chromosome_Organization_2022/tex_figures/RDF/raw_plots/'
out_label = 'RDF_bond_bending_twisting'

write_shell_RDF = True
read_shell_RDF = False

if not os.path.isdir(out_dir) and write_shell_diff:
    os.makedirs(out_dir)

shell_RDF_file = out_dir + out_label + '_shell_RDF.pickle'

dt = 0.1 # timestep - [ns]

# number of modes
N_modes = 20
Rc = 5000.0

# radii of shells - [A]
radii = np.array([1500,2000],dtype=np.int32)

# timesteps for RDF
dts = 10000
ts_RDF = np.arange(1010000,2000000,dts,dtype=np.int32)

if write_shell_RDF == True:

    rep_count = 0

    for rep in range(min_rep,max_rep+1):

        in_file = in_dir + in_label + '_rep'+str(rep).zfill(5) + '.pkl_traj'

        print('LAMMPS trajectory: '+str(in_file))

        traj = r_L_D.read_pickle_traj(in_file)

        if rep == min_rep:

            shell_RDF = rdf_plot.new_shell_RDF(N_modes,
                                               Rc,
                                               radii,
                                               max_rep-min_rep+1,
                                               ts_RDF.shape[0],
                                               ts_RDF)

        shell_RDF = rdf_plot.fill_reps_shell_RDFs(shell_MSD,traj,rep-min_rep)

    rdf_plot.write_shell_diff(shell_diff_file,shell_RDF)

if read_shell_diff == True:

    shell_RDF = rdf_plot.read_shell_RDF(shell_RDF_file)

    
