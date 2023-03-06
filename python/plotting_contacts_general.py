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

in_dir = '/mnt/HDD_linux/Data/contact_maps/btree_chromo/disentanglement_testing/woloops_wotopo/'
in_label = 'woloops_wotopo'

in_label_alt = in_label
CGinfo_file = in_dir + in_label_alt + '_CGmap.dat'

out_dir = '/mnt/HDD_linux/Documents/svn/Minimal_Cell_Chromosome_Organization_2022/tex_figures/replication_contact_maps/raw_plots/'
out_label = 'woloops_wotopo'

mat = c_plot.read_matrix(in_dir,in_label,True)

CGinfo = c_plot.read_CGinfo_file(CGinfo_file)

mapped_mat = c_plot.read_matrix(in_dir,'mapped_'+in_label,True)

thresh = [1.0E-6,1.0E-1]

#norm_flag = 'manual'
norm_flag = 'scaled'

c_plot.plot_mat(out_dir,out_label,mat,CGinfo,[87,87],norm_flag,thresh,True,True)

c_plot.plot_mat_mapped(out_dir,out_label,mapped_mat,CGinfo,[87,87],norm_flag,thresh,True,True)

#c_plot.plot_mat_mapped_dev(out_dir,out_label,mapped_mat,CGinfo,[87,87],True,True)

#c_plot.plot_mapping(out_dir,out_label,CGinfo,[87,87],True,True)
