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

# in_dir = '/home/ben/Data/btree_chromo/diffusion_testing/'
# in_label = 'diffusion_bond_bending_twisting'
in_dir = '/home/ben/Data/btree_chromo/diffusion_testing_v4/'
in_label = 'diffusion_bond_bending_twisting_v4'

# out_dir = '/home/ben/Documents/svn/Minimal_Cell_Chromosome_Organization_2022/tex_figures/RDF/raw_plots/'
# out_label = 'no_interactions'
out_dir = '/home/ben/Data/btree_chromo/diffusion_testing_v4/'
out_label = 'bond_bending_twisting_v4'

rep_specification = 'range'

if rep_specification == 'range':
    min_rep = 1
    max_rep = 50
    reps = np.arange(min_rep,max_rep+1,dtype=np.int32)
elif rep_specification == 'manual':
    #reps = [6,13,20,27,34,41,49]
    reps = [1,2,3,4,5,6,8,9,10,11,12,13,15,16,17,18,19,20,22,23,24,25,26,27,29,30,31,32,33,34,36,37,38,39,40,41,44,45,46,47,48,49]
    reps = np.array(reps,dtype=np.int32)

N_reps = reps.shape[0]

write_single_shell_RDFs = False
merge_shell_RDFs = True
read_shell_RDF = True

if not os.path.isdir(out_dir) and write_single_shell_RDFs:
    os.makedirs(out_dir)

merged_shell_RDF_file = out_dir + out_label + '_merged_shell_RDF.pickle'

dt = 0.1 # timestep - [ns]

# number of modes
N_modes = 100
Rc = 500.0

# radii of shells - [A]
radii = np.array([1500,2000],dtype=np.int32)

# timesteps for RDF
#dts = 10000
#ts_RDF = np.arange(1990000,2000000+dts,dts,dtype=np.int32)
dts = 100000
ts_RDF = np.arange(1000000,2000000+dts,dts,dtype=np.int32)

if write_single_shell_RDFs == True:

    for i_rep in range(N_reps):

        rep = reps[i_rep]

        in_file = in_dir + in_label + '_rep'+str(rep).zfill(5) + '.pkl_traj'

        shell_RDF_file = out_dir + out_label + '_shell_RDF_rep'+str(rep).zfill(5) + '.pickle'

        print('LAMMPS trajectory: '+str(in_file))

        traj = r_L_D.read_pickle_traj(in_file)

        shell_RDF = rdf_plot.new_single_shell_RDF(N_modes,
                                                  Rc,
                                                  radii,
                                                  ts_RDF.shape[0],
                                                  ts_RDF)

        shell_RDF = rdf_plot.fill_shell_RDF(shell_RDF,traj)

        print('writing RDF file:')
        rdf_plot.write_shell_RDF(shell_RDF_file,shell_RDF)
        print(shell_RDF_file)

if merge_shell_RDFs == True:

    merged_shell_RDF = rdf_plot.new_merged_shell_RDF(N_modes,
                                                     Rc,
                                                     radii,
                                                     N_reps,
                                                     ts_RDF.shape[0],
                                                     ts_RDF)

    for i_rep in range(N_reps):

        rep = reps[i_rep]

        single_shell_RDF_file = out_dir + out_label + '_shell_RDF_rep'+str(rep).zfill(5) + '.pickle'

        print(single_shell_RDF_file)

        single_shell_RDF = rdf_plot.read_shell_RDF(single_shell_RDF_file)

        merged_shell_RDF = rdf_plot.fill_merged_shell_RDF(merged_shell_RDF,
                                                          single_shell_RDF,
                                                          i_rep)

    rdf_plot.write_shell_RDF(merged_shell_RDF_file,merged_shell_RDF)

if read_shell_RDF == True:

    shell_RDF = rdf_plot.read_shell_RDF(merged_shell_RDF_file)

    shell_RDFs_file = out_dir + out_label + '_shell_RDFs.png'

    rdf_plot.plot_shell_RDFs(shell_RDFs_file,shell_RDF,1000)

    
