import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as colormaps
import matplotlib.patches as patches
from matplotlib.lines import Line2D

from mpl_toolkits.axes_grid1 import make_axes_locatable


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

def new_polymer_correlations(s_min,s_max,N_reps,N_t):

    polymer_correlations = dict()

    polymer_correlations['N_reps'] = N_reps
    polymer_correlations['N_t'] = N_t

    polymer_correlations['s'] = np.arange(s_min,s_max+1,dtype=np.int32)
    N_s = s_max - s_min + 1
    polymer_correlations['N_s'] = N_s

    polymer_correlations['corr'] = np.zeros((N_s,N_reps,N_t),dtype=np.double)

    return polymer_correlations

def fill_displacement_correlations_circular(poly_corr,traj,i_rep):

    for i_t in range(poly_corr['N_t']):

        t = traj['t'][i_t]
        x = traj[t]['DNA']['x']

        dx = x - np.roll(x,1,axis=0)

        udx = dx/np.sqrt(np.sum(np.power(dx,2.0),axis=1))

        for i_s in range(poly_corr['N_s']):

            s = poly_corr['s'][i_s]

            poly_corr['corr'][i_s,i_rep,i_t] = np.dot(udx.T,np.roll(udx,s,axis=0))
    
    return poly_corr

def plot_diffusion(fig_file,t,MSD,lims):
    
    fig_size = [87,87]

    fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

    ax = plt.gca()

    ax.set_xlabel(r'$t$ - Time [ns]', fontsize=8)
    ax.set_ylabel(r'$\langle|\mathbf{x}(t)-\mathbf{x}(0)|^2\rangle$ - MSD [\AA$^2$]', fontsize=8)
    #ax.set_title(r'MSD of 10~bp DNA monomers',fontsize=10)

    temp_marker_style = dict(marker='.', markersize=8,fillstyle='none')

    ax.plot(t,MSD,c='k',lw=1.0,alpha=0.35,**temp_marker_style)

    lims = np.array(lims,dtype=np.float32)

    multi_fits_linear(ax,lims,t,MSD)

    ax.tick_params(axis='x',labelsize=8)
    ax.tick_params(axis='y',labelsize=8)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    plt.tight_layout()

    ax.legend(fontsize=6)

    #plt.grid()

    fig.savefig(fig_file,dpi=300)

    return
    
