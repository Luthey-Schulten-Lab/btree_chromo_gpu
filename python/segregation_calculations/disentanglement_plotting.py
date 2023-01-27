import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as colormaps
import matplotlib.patches as patches
from matplotlib.lines import Line2D

from mpl_toolkits.axes_grid1 import make_axes_locatable

import os

import pickle
import glob

from . import fork_partitions as f_p


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

def new_DoD(R,N_reps,N_t,t):

    DoD = dict()

    DoD['R'] = R
    DoD['N_reps'] = N_reps
    DoD['N_t'] = N_t
    DoD['t'] = t

    DoD['forks'] = []
    DoD['N_forks'] = 0

    return DoD

def fill_reps_DoD(DoD,traj,f_p_files,i_rep):

    # use glob to get all of the fork_partition files
    f_p_files = glob.glob(f_p_files)

    # create a dictionary of the timesteps with new fork_partitions
    ts_f_ps = dict()
    ts_f_ps['ts'] = np.zeros((len(f_p_files)),dtype=np.int32)
    ts_f_ps['fork_partitions'] = []

    for i in range(len(f_p_files)):

        t = f_p_files[i].split('_')
        t = t[-1]
        t = t.split('.')
        t = t[0]
        t = t.strip('t')

        ts_f_ps['ts'][i] = int(t)

        f_ps = f_p.fork_partitions()
        #f_ps = fork_partitions()
        # append an empty fork_partition
        ts_f_ps['fork_partitions'].append(f_ps)

    # find args to sort the timesteps of the fork_partitions
    ts_arg_sorted = np.argsort(ts_f_ps['ts'])

    # loop over the sorted args and parse the fork_partition files
    for i_t in range(ts_arg_sorted.shape[0]):
        ts_f_ps['fork_partitions'][i_t].parse_f_p_file(f_p_files[ts_arg_sorted[i_t]])

    # now sort the timetseps
    ts_f_ps['ts'] = np.sort(ts_f_ps['ts'])

    # loop over the timesteps for the disentangle calculations
    for i_t in range(DoD['N_t']):

        # increment backwards to begin with the most replicated state
        i_t_m = DoD['N_t'] - 1 - i_t

        # select fork_partition from most recent update
        ts_f_p_selected = np.argwhere(ts_f_ps['ts']<DoD['t'][i_t_m])[-1][0]
        f_ps = ts_f_ps['fork_partitions'][ts_f_p_selected]

        # get the coordinates of the DNA at the timestep
        x = traj[str(DoD['t'][i_t_m])]['DNA']['x']

        # get the forks and calculate the DoD
        forks, d = disentanglement(f_ps,x,DoD['R'])

        # perform actions for first iteration
        if i_t == 0:
            # store the final set of forks
            final_forks = forks
            # allocate array with results
            d_forks_ts = np.zeros((len(final_forks),DoD['N_t']),
                                  dtype=np.double)

        for i_fork in range(len(forks)):

            for j_fork in range(len(final_forks)):

                if forks[i_fork] == final_forks[j_fork]:

                    d_forks_ts[j_fork,i_t_m] = d[i_fork]

    if i_rep == 0:

        DoD['N_forks'] = len(final_forks)
        DoD['forks'] = final_forks

        DoD['d_reps_forks_ts'] = np.zeros((DoD['N_reps'],
                                           DoD['N_forks'],
                                           DoD['N_t']),
                                          dtype=np.double)

        DoD['d_reps_forks_ts'][0,:,:] = d_forks_ts

    else:

        for i_fork in range(len(final_forks)):

            k = -1

            for j_fork in range(DoD['N_forks']):

                if final_forks[i_fork] == DoD['forks'][j_fork]:

                    k = j_fork
                    break

            if k != -1:

                DoD['d_reps_forks_ts'][i_rep,k,:] = d_forks_ts[i_fork,:]

            else:

                DoD['N_forks'] += 1
                DoD['forks'].append(final_forks[i_fork])

                # add new elements to array
                zero_d_forks_ts = np.zeros((DoD['N_reps'],
                                            DoD['N_t']),
                                           dtype=np.double)

                np.append(DoD['d_reps_forks_ts'],
                          zero_d_forks_ts,
                          axis=1)

                DoD['d_reps_forks_ts'][i_rep,DoD['N_forks']-1,:] = d_forks_ts[i_fork,:]
    
    return DoD

def write_DoD(DoD_file,DoD):

    with open(DoD_file,'wb') as f:
        pickle.dump(DoD,f)

    return

def read_DoD(DoD_file):

    with open(DoD_file,'rb') as f:
        DoD = pickle.load(f)

    return DoD

def disentanglement(f_ps,x,R):

    R2 = np.power(R,2.0)

    N_forks = f_ps.get_N_forks()
    forks = f_ps.get_fork_list()

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

        phi_bar_l, phi_bar_r = tanh_proximity_fraction(x,
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

    alpha = N_l/N_r
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

def tanh_proximity_fraction(x,N_l,N_r,a_l,a_r,R2):

    phi_bar_l = 0.0
    phi_bar_r = 0.0

    n_s_l = tanh_proximity_same(x,N_l,a_l,R2)
    n_s_r = tanh_proximity_same(x,N_r,a_r,R2)

    n_o_l, n_o_r = tanh_proximity_opposite(x,N_l,a_l,N_r,a_r,R2)

    alpha = N_l/N_r
    phi_l = np.divide(n_s_l,(n_s_l + alpha*n_o_l))
    phi_r = np.divide(n_s_r,(n_s_r + n_o_r/alpha))

    phi_bar_l = np.mean(phi_l) 
    phi_bar_r = np.mean(phi_r)

    # print("phi_bar_l = " + str(phi_bar_l))
    # print("phi_bar_r = " + str(phi_bar_r))
    
    return phi_bar_l, phi_bar_r

def tanh_proximity_same(x,N,a,R2):

    R2i = 1.0/R2

    xt = x[a,:]

    n_s = np.zeros((N),dtype=np.double)

    for i in range(N-1):
        
        xc = xt[i,:]
        xnc = xt[i+1:,:]

        dx = xnc - np.tile(xc,(xnc.shape[0],1))

        d2 = np.sum(np.power(dx,2.0),axis=1)

        y = np.sqrt(R2i*d2)
        y = 1.0 + np.tanh(-y)
        
        p = np.arange(i+1,N,1,dtype=np.int32)

        n_s[p] += y
        n_s[i] += np.sum(y)

    return n_s

def tanh_proximity_opposite(x,N_l,a_l,N_r,a_r,R2):

    R2i = 1.0/R2

    n_o_l = np.zeros((N_l),dtype=np.double)
    n_o_r = np.zeros((N_r),dtype=np.double)

    x_l = x[a_l,:]
    x_r = x[a_r,:]

    for i in range(N_l):

        xc = x_l[i,:]
        xnc = x_r

        dx = xnc - np.tile(xc,(xnc.shape[0],1))

        d2 = np.sum(np.power(dx,2.0),axis=1)

        y = np.sqrt(R2i*d2)
        y = 1.0 + np.tanh(-y)

        n_o_r += y
        n_o_l[i] += np.sum(y)
    
    return n_o_l, n_o_r


def plot_DoD(fig_file,DoD):

    cmap = colormaps.get_cmap('winter')
    c_space = np.linspace(0.0,1.0,DoD['N_forks'])

    fig_size = [87,87]

    fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

    ax = plt.gca()

    ax.set_xlabel(r'timestep', fontsize=12)
    ax.set_ylabel(r'degree of disentanglement', fontsize=12)


    tick_length = 4.0
    tick_width = 2.0
    ax.tick_params(labelsize=9,
                   length=tick_length,
                   width=tick_width,
                   direction='out',
                   left=True,
                   right=False,
                   bottom=True,
                   top=False)

    #ax.spines['right'].set_visible(False)
    #ax.spines['top'].set_visible(False)
    ax.spines['left'].set_linewidth(2.0)
    ax.spines['bottom'].set_linewidth(2.0)
    ax.spines['right'].set_linewidth(2.0)
    ax.spines['top'].set_linewidth(2.0)
    

    ax.set_xlim(xmin=0,xmax=DoD['t'][-1])
    ax.set_ylim(ymin=0.0,ymax=1.0)

    ax.grid(which='major',axis='both',zorder=-4)

    d_means = np.mean(DoD['d_reps_forks_ts'],axis=0)

    for i_fork in range(DoD['N_forks']):

        temp_color = cmap(c_space[i_fork])

        for i_rep in range(DoD['N_reps']):

            ax.plot(DoD['t'],
                    DoD['d_reps_forks_ts'][i_rep,i_fork,:],
                    lw=1.0,
                    alpha=0.4,
                    ls='-',
                    c=temp_color,
                    zorder=-2)

        ax.plot(DoD['t'],
                d_means[i_fork,:],
                lw=2.0,
                alpha=1.0,
                ls='-',
                c='w',
                zorder=1)
            
        ax.plot(DoD['t'],
                d_means[i_fork,:],
                lw=1.5,
                alpha=1.0,
                ls='-',
                c=temp_color,
                zorder=2)

        temp_str = DoD['forks'][i_fork]

        x_fork_label = 1.01
        y_fork_label = d_means[i_fork,-1]
        ax.annotate(text=r''+temp_str,
                    xy=(x_fork_label,y_fork_label),
                    xycoords='axes fraction',
                    color=temp_color,
                    va='center',
                    ha='left')

        

    plt.tight_layout()

    fig.savefig(fig_file,dpi=300)

    plt.close()

    return
