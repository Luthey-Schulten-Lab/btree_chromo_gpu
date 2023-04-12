# =============================================
# Author: Benjamin R. Gilbert
# Email: brg4@illinois.edu
# =============================================

import numpy as np

import time
import os

import importlib
import sys

import contact_calculations.contact_plotting as c_plot
importlib.reload(c_plot)

in_dir_p = '/mnt/HDD_linux/Data/contact_maps/btree_chromo/disentanglement_testing/wextraloops_wtopo/'
in_label_p = 'wextraloops_wtopo'

in_dir_m = '/mnt/HDD_linux/Data/contact_maps/btree_chromo/disentanglement_testing/wextraloops_wotopo/'
in_label_m = 'wextraloops_wotopo'

in_label_alt = in_label_p
CGinfo_file = in_dir_p + in_label_alt + '_CGmap.dat'

out_dir = '/mnt/HDD_linux/Documents/svn/Minimal_Cell_Chromosome_Organization_2022/tex_figures/replication_contact_maps/raw_plots/'
out_label = 'wextraloops_wtopo_vs_wextraloops_wotopo'

mat = c_plot.read_matrix(in_dir,in_label,True)

CGinfo = c_plot.read_CGinfo_file(CGinfo_file)

mat_p = c_plot.read_matrix(in_dir_p,'mapped_'+in_label_p,True)
dm_p = np.mean(np.diagonal(mat_p))
print('dm_p={:.3f}\n'.format(dm_p))
mat_m = c_plot.read_matrix(in_dir_m,'mapped_'+in_label_m,True)
dm_m = np.mean(np.diagonal(mat_m))
print('dm_m={:.3f}\n'.format(dm_m))

mat = mat_m - mat_p

thresh = [1.0E-6,1.0E-2]

#norm_flag = 'manual'
norm_flag = 'scaled'

c_plot.plot_mat_diff(out_dir,out_label,mat,CGinfo,[80,80],True,False)
