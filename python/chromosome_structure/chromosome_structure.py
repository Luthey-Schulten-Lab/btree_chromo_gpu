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

def new_RoG(N_conditions,conditions,N_reps,N_t,t,freq,window):

    RoG = dict()

    RoG['N_conditions'] = N_conditions
    RoG['conditions'] = conditions
    RoG['N_reps'] = N_reps
    RoG['N_t'] = N_t
    RoG['t'] = t
    RoG['freq'] = freq
    RoG['window'] = window
    RoG['R'] = np.zeros((N_conditions,N_reps,N_t),dtype=np.double)

    return RoG

def fill_conditions_reps_RoG(RoG,traj,i_condition,i_rep):

    # loop over the timesteps for the disentangle calculations
    for i_t in range(RoG['N_t']):

        t = str(RoG['t'][i_t])

        # get the coordinates of the DNA at the timestep
        x = traj[t]['DNA']['x']

        R = calc_windowed_RoG(x,RoG['freq'],RoG['window'])

        RoG['R'][i_condition,i_rep,i_t] = R
    
    return RoG

def write_RoG(RoG_file,RoG):

    with open(RoG_file,'wb') as f:
        pickle.dump(RoG,f)

    return

def read_RoG(RoG):

    with open(RoG,'rb') as f:
        RoG = pickle.load(f)

    return RoG

def calc_windowed_RoG(x,f,w):

    N = x.shape[0]

    i = 0
    j = 0
    M = N//w
    if N%M > 0:
        M += 1
    R = np.zeros((M),dtype=np.double)

    ds = np.arange(0,w,1,dtype=np.int32)

    while i < x.shape[0]:

        s = i + ds
        s = np.remainder(s,N)

        xs = x[s,:]

        xm = np.mean(xs,axis=0)

        dxm = xs - xm

        R[j] = np.sqrt(np.mean(np.sum(np.power(dxm,2.0),axis=1)))
        
        i += w
        j += 1

    R = np.mean(R)

    return R

def plot_RoG(fig_file,RoG):
    
    fig_size = [1.4*80,85]

    cmap = colormaps.get_cmap('plasma')
    c_space_lower_lim = 0.0
    c_space_upper_lim = 0.7
    c_space = np.linspace(c_space_lower_lim,
                          c_space_upper_lim,
                          RoG['N_conditions'])

    fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

    ax = plt.gca()

    ax.set_xlabel(r'$t$ - simulation time', fontsize=10)
    ax.set_ylabel(r'radius of gyration ' +
                  '({:d} monomers)'.format(RoG['window']) +
                  ' - [\AA]', fontsize=10)

    #temp_marker_style = dict(marker='.', markersize=8,fillstyle='none')

    t_trans = RoG['t']/RoG['t'][-1]

    for i_condition in range(RoG['N_conditions']):

        temp_color = cmap(c_space[i_condition])

        y = np.mean(RoG['R'][i_condition,:,:],axis=0)
        
        ax.plot(t_trans,y,
                lw=2.0,
                alpha=1.0,
                c=temp_color)

        temp_label = '{:d} loops'.format(RoG['conditions'][i_condition])
        ax.annotate(r''+temp_label,
                    xy=(t_trans[-1]+0.02,y[-1]),
                    xycoords='data',
                    va='center',
                    ha='left',
                    color=temp_color,
                    fontsize=8,
                    clip_on=False)

    ax.set_xticks(ticks=[0.0,0.25,0.5,0.75,1.0],minor=False)
    ax.set_xticks(ticks=[0.125,0.375,0.625,0.875],minor=True)

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

    ax.set_xticklabels(labels=[r'$0$',r'$T_f/4$',r'$T_f/2$',r'$3T_f/4$',r'$T_f$'])

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_linewidth(2.0)
    ax.spines['bottom'].set_linewidth(2.0)

    ax.set_xlim(xmin=0.0,xmax=1.1*t_trans[-1])

    #ax.legend(fontsize=6)

    plt.tight_layout()

    print(fig_file)

    fig.savefig(fig_file,dpi=300)

    plt.close()

    return
