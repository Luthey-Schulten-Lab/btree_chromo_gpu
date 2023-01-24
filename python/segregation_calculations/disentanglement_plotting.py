import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as colormaps
import matplotlib.patches as patches
from matplotlib.lines import Line2D

from mpl_toolkits.axes_grid1 import make_axes_locatable

import os

import pickle


# import matplotlib.font_manager as font_manager
# font_manager._rebuild()

#print(matplotlib.get_cachedir())

# use LaTeX fonts in the plot
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage[helvet]{sfmath}'

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'

mm = 1/25.4
l_w_ratio = 1/1.618

def disentanglement(f_ps,x,R):

    R2 = np.power(R,2.0)

    N_forks = f_ps.get_N_forks()
    forks = f_ps.get_fork_list()
    print(N_forks)
    print(forks)

    d_o_d = np.zeros((N_forks),dtype=np.double)

    for i_fork in range(N_forks):

        f_p = f_ps.get_fork_partition(i_fork)

        N_l = f_p.get_left_size()
        N_r = f_p.get_right_size()
        a_l = f_p.get_full_left_mono_array()
        a_r = f_p.get_full_right_mono_array()
        # print(forks[i_fork])
        # print(a_l.shape)
        # print(a_r.shape)

        phi_bar_l, phi_bar_r = step_proximity_fraction(x,
                                                       N_l,
                                                       N_r,
                                                       a_l,
                                                       a_r,
                                                       R2)

        #d_o_d[i_fork] = (phi_bar_l + phi_bar_r)/2.0
        d_o_d[i_fork] = 2.0*phi_bar_l*phi_bar_r/(phi_bar_l + phi_bar_r)
        
    return forks, d_o_d

def step_proximity_fraction(x,N_l,N_r,a_l,a_r,R2):

    phi_bar_l = 0.0
    phi_bar_r = 0.0

    n_s_l = step_proximity_same(x,N_l,a_l,R2)
    n_s_r = step_proximity_same(x,N_r,a_r,R2)

    n_o_l, n_o_r = step_proximity_opposite(x,N_l,a_l,N_r,a_r,R2)

    alpha = 1.0 #N_l/N_r
    phi_l = np.divide(n_s_l,(n_s_l + alpha*n_o_l))
    phi_r = np.divide(n_s_r,(n_s_r + n_o_r/alpha))

    phi_bar_l = np.mean(phi_l)
    phi_bar_r = np.mean(phi_r)

    # print("phi_bar_l = " + str(phi_bar_l))
    # print("phi_bar_r = " + str(phi_bar_r))
    
    return phi_bar_l, phi_bar_r

def step_proximity_same(x,N,a,R2):

    xt = x[a,:]

    n_s = np.zeros((N),dtype=np.int32)

    for i in range(N-1):
        
        xc = xt[i,:]
        xnc = xt[i+1:,:]

        dx = xnc - np.tile(xc,(xnc.shape[0],1))

        d2 = np.sum(np.power(dx,2.0),axis=1)

        p = np.argwhere(d2<R2) + (i + 1)

        n_s[p] += 1
        n_s[i] += p.shape[0]

    return n_s

def step_proximity_opposite(x,N_l,a_l,N_r,a_r,R2):

    n_o_l = np.zeros((N_l),dtype=np.int32)
    n_o_r = np.zeros((N_r),dtype=np.int32)

    x_l = x[a_l,:]
    x_r = x[a_r,:]

    for i in range(N_l):

        xc = x_l[i,:]
        xnc = x_r

        dx = xnc - np.tile(xc,(xnc.shape[0],1))

        d2 = np.sum(np.power(dx,2.0),axis=1)

        p = np.argwhere(d2<R2)

        n_o_r[p] += 1
        n_o_l[i] += p.shape[0]
    
    return n_o_l, n_o_r
