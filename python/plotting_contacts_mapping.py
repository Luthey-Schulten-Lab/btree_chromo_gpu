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

in_dir = '/mnt/HDD_linux/Data/btree_chromo/theta_example/'
in_label = 'theta_example_scaled'

in_label_alt = in_label
CGinfo_file = in_dir + in_label_alt + '_CGmap_t3.dat'

out_dir = '/mnt/HDD_linux/Documents/svn/Minimal_Cell_Chromosome_Organization_2022/tex_figures/contact_calc_summary/raw_plots/'
out_label = 'theta_example_scaled'

CGinfo = c_plot.read_CGinfo_file(CGinfo_file)

checkerboard_length = 20

c_plot.plot_mapping_simple(out_dir,out_label,CGinfo,checkerboard_length,[80,80],True,True)
