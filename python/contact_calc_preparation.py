# =============================================
# Author: Benjamin R. Gilbert
# Email: brg4@illinois.edu
# =============================================

import numpy as np

import time
import os

import implib
import sys

import LAMMPS_helpers.read_LAMMPS_DNA as r_L_D
importlib.reload(r_L_D)

in_dir = '/home/ben/Data/LAMMPS/s1c15_503ribos/s1c15_503ribos_m25000_ml5000_mr5000/'
in_label = 's1c15_503ribos_m25000_ml5000_mr5000'
min_rep = 1
max_rep = 20

out_dir = '/home/ben/Data/contact_maps/s1c15_503ribos_LAMMPS/s1c15_503ribos_m25000_ml5000_mr5000/'
out_label = 's1c15_503ribos_m25000_ml5000_mr5000'

write_flag = True
write_coords_flag = True
timestep_options = ['Last','FirstLast','Slices','All']
timestep_select = timestep_options[3]

if not os.path.isdir(out_dir) and write_flag:
    os.makedirs(out_dir)

# size of chromosome
genome_size = 54338

genome_features = {}
unif_gene_size = 1000

gene_starts = np.arange(1,genome_size,unif_gene_size)
gene_ends = gene_starts + (unif_gene_size - 1)
gene_ends[-1] = min(gene_ends[-1],genome_size-1)

N_genes = gene_starts.shape[0]

for i_gene in range(N_genes):
    
    genome_features['gene {:d}'.format(i_gene)] = [gene_starts[i_gene],gene_ends[i_gene]]

# create mother chromosome
mother = t_c.chromosome(genome_size,0,genome_features)

mother.read_replication_state(in_dir+in_label+'_chromo_state.txt')

mother.partition_units()

CG_file = out_dir + out_label + '_CGinfo.txt'

if write_flag:
    
    mother.write_replication_state(out_dir+out_label+'_chromo_state.txt')
    mother.write_CG_file(CG_file,50)

rep_count = 0

for rep in range(min_rep,max_rep+1):

    in_file = in_dir + in_label + '_rep'+str(rep).zfill(5) + '.lammpstrj'

    print('LAMMPS trajectory: '+str(in_file))

    if write_coords_flag == True:

        print(timestep_select)
        
        traj = r_L_D.read_traj(in_file,in_timesteps=timestep_select)

        # t_target = str(traj['timesteps'][-1])

        # out_file = out_dir + 'x_' + out_label +\
        #     '_t' + str(t_target.zfill(len(str(traj['timesteps'][-1])))) +\
        #     '_rep' + str(rep).zfill(5) + '.bin'  # output file
        # print(out_file)

        if rep == min_rep:
            N_t = traj['timesteps'].shape[0]
            rep_count += (min_rep-1)*N_t

        for i_t in range(N_t):

            t_target = str(traj['timesteps'][i_t])

            print(t_target)
    
            N = traj[t_target]['DNA']['N']
            x = traj[t_target]['DNA']['x']

            x = np.reshape(x,(N,3),order='F')

            out_file = out_dir + 'x_' + out_label +\
                '_rep' + str(rep_count).zfill(5) + '.bin'  # output file
            print(out_file)

            with open(out_file,'wb') as f:

                x.tofile(f)

            rep_count+=1
