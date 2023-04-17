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

def new_par_CoM(R,N_reps,N_t,t):

    par_CoM = dict()

    par_CoM['N_reps'] = N_reps
    par_CoM['N_t'] = N_t
    par_CoM['t'] = t

    par_CoM['forks'] = []
    par_CoM['N_forks'] = 0

    return par_CoM

def fill_reps_par_CoM(par_CoM,traj,f_p_files,i_rep):

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
    for i_t in range(par_CoM['N_t']):

        # increment backwards to begin with the most replicated state
        i_t_m = par_CoM['N_t'] - 1 - i_t
        print('t = ' + str(par_CoM['t'][i_t_m]))

        # select fork_partition from most recent update
        ts_f_p_selected = np.argwhere(ts_f_ps['ts']<par_CoM['t'][i_t_m])[-1][0]
        print(ts_f_ps['ts'])
        print(ts_f_ps['ts'][ts_f_p_selected])
        f_ps = ts_f_ps['fork_partitions'][ts_f_p_selected]

        # get the coordinates of the DNA at the timestep
        x = traj[str(par_CoM['t'][i_t_m])]['DNA']['x']

        # get the forks and calculate the par_CoM
        forks, d, a = distance_CoM(f_ps,x)

        # perform actions for first iteration
        if i_t == 0:
            # store the final set of forks
            final_forks = forks
            # allocate array with results
            d_forks_ts = np.zeros((len(final_forks),par_CoM['N_t']),
                                  dtype=np.double)

            a_forks_ts = 0.5*np.ones((len(final_forks),par_CoM['N_t']),
                                     dtype=np.double)

        for i_fork in range(len(forks)):

            for j_fork in range(len(final_forks)):

                if forks[i_fork] == final_forks[j_fork]:

                    d_forks_ts[j_fork,i_t_m] = d[i_fork]
                    a_forks_ts[j_fork,i_t_m] = a[i_fork]

    if i_rep == 0:

        par_CoM['N_forks'] = len(final_forks)
        par_CoM['forks'] = final_forks

        par_CoM['d_reps_forks_ts'] = np.zeros((par_CoM['N_reps'],
                                           par_CoM['N_forks'],
                                           par_CoM['N_t']),
                                          dtype=np.double)

        par_CoM['a_reps_forks_ts'] = 0.5*np.ones((par_CoM['N_reps'],
                                                  par_CoM['N_forks'],
                                                  par_CoM['N_t']),
                                                 dtype=np.double)

        par_CoM['d_reps_forks_ts'][0,:,:] = d_forks_ts
        par_CoM['a_reps_forks_ts'][0,:,:] = a_forks_ts

    else:

        for i_fork in range(len(final_forks)):

            k = -1

            for j_fork in range(par_CoM['N_forks']):

                if final_forks[i_fork] == par_CoM['forks'][j_fork]:

                    k = j_fork
                    break

            if k != -1:

                par_CoM['d_reps_forks_ts'][i_rep,k,:] = d_forks_ts[i_fork,:]
                par_CoM['a_reps_forks_ts'][i_rep,k,:] = a_forks_ts[i_fork,:]

            else:

                par_CoM['N_forks'] += 1
                par_CoM['forks'].append(final_forks[i_fork])

                # add new elements to array
                zero_d_forks_ts = np.zeros((par_CoM['N_reps'],
                                            par_CoM['N_t']),
                                           dtype=np.double)

                zero_a_forks_ts = 0.5*np.ones((par_CoM['N_reps'],
                                               par_CoM['N_t']),
                                              dtype=np.double)

                np.append(par_CoM['d_reps_forks_ts'],
                          zero_d_forks_ts,
                          axis=1)

                np.append(par_CoM['a_reps_forks_ts'],
                          zero_a_forks_ts,
                          axis=1)

                par_CoM['d_reps_forks_ts'][i_rep,par_CoM['N_forks']-1,:] = d_forks_ts[i_fork,:]
                par_CoM['a_reps_forks_ts'][i_rep,par_CoM['N_forks']-1,:] = a_forks_ts[i_fork,:]
    
    return par_CoM

def write_par_CoM(par_CoM_file,par_CoM):

    with open(par_CoM_file,'wb') as f:
        pickle.dump(par_CoM,f)

    return

def read_par_CoM(par_CoM_file):

    with open(par_CoM_file,'rb') as f:
        par_CoM = pickle.load(f)

    return par_CoM

def distance_CoM(f_ps,x):

    N_forks = f_ps.get_N_forks()
    forks = f_ps.get_fork_list()

    d_CoM = np.zeros((N_forks),dtype=np.double)
    a = 0.5*np.ones((N_forks),dtype=np.double)

    for i_fork in range(N_forks):

        f_p = f_ps.get_fork_partition(i_fork)

        N_l = f_p.get_left_size()
        N_r = f_p.get_right_size()
        a_l = f_p.get_full_left_mono_array()
        a_r = f_p.get_full_right_mono_array()

        xm_l = np.mean(x[a_l,:],axis=0)
        xm_r = np.mean(x[a_r,:],axis=0)

        d_CoM[i_fork] = np.linalg.norm(xm_l-xm_r)

        q = float(N_l)/(N_l+N_r)

        p = np.array([-q,0.0,3.0,-2.0],dtype=np.double)
        p = np.polynomial.Polynomial(p)
        roots = p.roots()

        # print(roots)

        a[i_fork] = roots[1]
        
    return forks, d_CoM, a


def plot_par_CoM(fig_file,par_CoM,R_bdry):

    L_part = (3.0/4.0)*R_bdry

    cmap = colormaps.get_cmap('RdPu')
    c_space_lower_lim = 0.5
    c_space_upper_lim = 1.0
    c_space = np.linspace(c_space_lower_lim,
                          c_space_upper_lim,
                          par_CoM['N_forks'])

    #fig_size = [87,87]
    fig_size = [70,35]

    fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

    ax = plt.gca()

    #ax.set_xlabel(r'timestep', fontsize=12)
    #ax.set_ylabel(r'degree of disentanglement', fontsize=12)

    t_trans = par_CoM['t']/par_CoM['t'][-1]
    
    ax.set_xticks(ticks=[0.0,0.5,1.0],minor=False)
    ax.set_xticks(ticks=[0.25,0.75],minor=True)

    tick_length = 4.0
    tick_width = 2.0
    ax.tick_params(labelsize=9,
                   length=tick_length,
                   width=tick_width,
                   direction='out',
                   left=True,
                   right=True,
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

    ax.set_xticklabels(labels=[r'$0$',r'$T_f/2$',r'$T_f$'])

    #ax.spines['right'].set_visible(False)
    #ax.spines['top'].set_visible(False)
    ax.spines['left'].set_linewidth(2.0)
    ax.spines['bottom'].set_linewidth(2.0)
    ax.spines['right'].set_linewidth(2.0)
    ax.spines['top'].set_linewidth(2.0)
    

    ax.set_xlim(xmin=0.0,xmax=t_trans[-1])
    ax.set_ylim(ymin=0.0,ymax=0.6)


    #ax.grid(which='both',axis='both',zorder=-4)

    d_means = np.mean(par_CoM['d_reps_forks_ts'],axis=0)

    for i_fork in range(par_CoM['N_forks']):

        temp_color = cmap(c_space[i_fork])

        for i_rep in range(par_CoM['N_reps']):

            ax.plot(t_trans,
                    par_CoM['d_reps_forks_ts'][i_rep,i_fork,:]/L_part,
                    lw=1.0,
                    alpha=0.4,
                    ls='-',
                    c=temp_color,
                    zorder=-2)

        ax.plot(t_trans,
                d_means[i_fork,:]/L_part,
                lw=2.0,
                alpha=1.0,
                ls='-',
                c='w',
                zorder=1)
            
        ax.plot(t_trans,
                d_means[i_fork,:]/L_part,
                lw=1.5,
                alpha=1.0,
                ls='-',
                c=temp_color,
                zorder=2)

        temp_str = par_CoM['forks'][i_fork]

        # label at start
        
        i_label = 0
        for i in range(t_trans.shape[0]):
            if d_means[i_fork,i] > 0.0:
                break
            else:
                i_label = i

        x_fork_label = t_trans[i_label] - 0.02
        y_fork_label = 0.08

        # label at end
        # x_fork_label = 1.01
        # y_fork_label = d_means[i_fork,-1]
        
        ax.annotate(text=r''+temp_str,
                    xy=(x_fork_label,y_fork_label),
                    xycoords='axes fraction',
                    color=temp_color,
                    va='center',
                    ha='right')

        

    plt.tight_layout()

    fig.savefig(fig_file,dpi=300)

    plt.close()

    return

def plot_par_CoM_scaled(fig_file,par_CoM,R_bdry):

    #L_part = (3.0/4.0)*R_bdry

    cmap = colormaps.get_cmap('RdPu')
    c_space_lower_lim = 0.5
    c_space_upper_lim = 1.0
    c_space = np.linspace(c_space_lower_lim,
                          c_space_upper_lim,
                          par_CoM['N_forks'])

    #fig_size = [87,87]
    fig_size = [70,35]

    fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

    ax = plt.gca()

    #ax.set_xlabel(r'timestep', fontsize=12)
    #ax.set_ylabel(r'degree of disentanglement', fontsize=12)

    t_trans = par_CoM['t']/par_CoM['t'][-1]
    
    ax.set_xticks(ticks=[0.0,0.5,1.0],minor=False)
    ax.set_xticks(ticks=[0.25,0.75],minor=True)

    tick_length = 4.0
    tick_width = 2.0
    ax.tick_params(labelsize=9,
                   length=tick_length,
                   width=tick_width,
                   direction='out',
                   left=True,
                   right=True,
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

    ax.set_xticklabels(labels=[r'$0$',r'$T_f/2$',r'$T_f$'])

    #ax.spines['right'].set_visible(False)
    #ax.spines['top'].set_visible(False)
    ax.spines['left'].set_linewidth(2.0)
    ax.spines['bottom'].set_linewidth(2.0)
    ax.spines['right'].set_linewidth(2.0)
    ax.spines['top'].set_linewidth(2.0)
    

    ax.set_xlim(xmin=0.0,xmax=t_trans[-1])
    ax.set_ylim(ymin=0.0,ymax=0.6)

    a = np.copy(par_CoM['a_reps_forks_ts'])

    ar = np.divide(np.power(a,2.0),2.0*a+1.0)
    al = np.divide(np.power(a-1.0,2.0),3.0-2.0*a)

    L_part = 3.0*R_bdry*(ar+al)

    print(np.mean(L_part[:,0,:],axis=0))

    d_trans = np.divide(par_CoM['d_reps_forks_ts'],L_part)


    #ax.grid(which='both',axis='both',zorder=-4)

    d_means = np.mean(d_trans,axis=0)

    for i_fork in range(par_CoM['N_forks']):

        temp_color = cmap(c_space[i_fork])

        for i_rep in range(par_CoM['N_reps']):

            ax.plot(t_trans,
                    d_trans[i_rep,i_fork,:],
                    lw=1.0,
                    alpha=0.4,
                    ls='-',
                    c=temp_color,
                    zorder=-2)

        ax.plot(t_trans,
                d_means[i_fork,:],
                lw=2.0,
                alpha=1.0,
                ls='-',
                c='w',
                zorder=1)
            
        ax.plot(t_trans,
                d_means[i_fork,:],
                lw=1.5,
                alpha=1.0,
                ls='-',
                c=temp_color,
                zorder=2)

        temp_str = par_CoM['forks'][i_fork]

        # label at start
        
        i_label = 0
        for i in range(t_trans.shape[0]):
            if d_means[i_fork,i] > 0.0:
                break
            else:
                i_label = i

        x_fork_label = t_trans[i_label] - 0.02
        y_fork_label = 0.08

        # label at end
        # x_fork_label = 1.01
        # y_fork_label = d_means[i_fork,-1]
        
        ax.annotate(text=r''+temp_str,
                    xy=(x_fork_label,y_fork_label),
                    xycoords='axes fraction',
                    color=temp_color,
                    va='center',
                    ha='right')

        

    plt.tight_layout()

    fig.savefig(fig_file,dpi=300)

    plt.close()

    return
