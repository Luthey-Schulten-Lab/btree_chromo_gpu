import numpy as np

import pickle

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

def new_initial_polymer_correlations(s_min,s_max,
                                     MC_steps,
                                     N_reps,N_MC,
                                     Lp_true,l_0):

    polymer_correlations = dict()

    polymer_correlations['Lp_true'] = Lp_true
    polymer_correlations['l_0'] = l_0

    polymer_correlations['MC_steps'] = MC_steps
    polymer_correlations['N_reps'] = N_reps
    polymer_correlations['N_MC'] = N_MC

    polymer_correlations['s'] = np.arange(s_min,s_max+1,dtype=np.int32)
    N_s = s_max - s_min + 1
    polymer_correlations['N_s'] = N_s

    polymer_correlations['corr'] = np.zeros((N_MC,N_reps,N_s),dtype=np.double)

    return polymer_correlations

def fill_bond_correlations_circular(poly_corr,x,i_MC,i_rep):

    dx = x - np.roll(x,1,axis=0)

    udx = dx*np.tile(np.reciprocal(np.sqrt(np.sum(np.power(dx,2.0),axis=1))),(3,1)).T

    for i_s in range(poly_corr['N_s']):

        s = poly_corr['s'][i_s]

        s_udx = np.roll(udx,s,axis=0)

        y = np.zeros(udx.shape[0],dtype=np.double)

        for i in range(udx.shape[0]):

            y[i] = np.dot(udx[i],s_udx[i])

        poly_corr['corr'][i_MC,i_rep,i_s] = np.mean(y)
    
    return poly_corr

def write_poly_corr(poly_corr_file,poly_corr):

    with open(poly_corr_file,'wb') as f:
        pickle.dump(poly_corr,f)

    return

def read_poly_corr(poly_corr_file):

    with open(poly_corr_file,'rb') as f:
        poly_corr = pickle.load(f)

    return poly_corr


def plot_initial_Lp(poly_corr,fig_file):

    fig_size = [174,87]

    cmap = colormaps.get_cmap('plasma')
    c_space_lower_lim = 0.0
    c_space_upper_lim = 0.7
    c_space = np.linspace(c_space_lower_lim,
                          c_space_upper_lim,
                          poly_corr['N_MC'])

    fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

    ax = plt.gca()

    ax.set_xlabel(r'$s$ - monomer displacement', fontsize=8)
    ax.set_ylabel(r'$\langle\hat{\mathbf{t}}_{i}\cdot\hat{\mathbf{t}}_{i+s}\rangle$ - bond correlation', fontsize=8)

    tick_length = 4.0
    tick_width = 2.0
    ax.tick_params(labelsize=9,
                   length=tick_length,
                   width=tick_width,
                   direction='out',
                   left=True,
                   right=False,
                   bottom=True,
                   top=False,
                   which='major')

    ax.tick_params(labelsize=9,
                   length=tick_length/1.5,
                   width=tick_width/1.5,
                   direction='out',
                   left=True,
                   right=False,
                   bottom=True,
                   top=False,
                   which='minor')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_linewidth(2.0)
    ax.spines['bottom'].set_linewidth(2.0)

    x = poly_corr['s']
    y = np.exp(-x*poly_corr['l_0']/poly_corr['Lp_true'])
    ax.plot(x,y,
            linewidth=2.5,
            color='black',
            alpha=0.5,
            linestyle='dashed',
            zorder=-2)

    for i_MC in range(poly_corr['N_MC']):

        for i_rep in range(poly_corr['N_reps']):

            y = poly_corr['corr'][i_MC,i_rep,:]

            ax.plot(x,y,
                    color=cmap(c_space[i_MC]),
                    linewidth=0.3,
                    alpha=0.4,zorder=-1)


        y = np.mean(poly_corr['corr'][i_MC,:,:],axis=0)
        
        ax.plot(x,y,
                color='white',
                linewidth=2.0,
                alpha=1.0,
                zorder=0)

        
        temp_label = r''+'{:d}'.format(poly_corr['MC_steps'][i_MC])
        
        ax.plot(x,y,
                color=cmap(c_space[i_MC]),
                linewidth=1.5,
                alpha=1.0,
                zorder=1,
                label=temp_label)


    ax.legend(title=r'\# MC Steps',
              title_fontsize=10,
              fontsize=8)

    #plt.grid()

    fig.savefig(fig_file,dpi=300)
