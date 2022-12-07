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


def monic_cheb_coeffs(N_modes,x):

    c = np.zeros((N_modes),dtype=np.double)

    half = 1.0/2.0
    quarter = 1.0/4.0

    c[0] = 1.0
    c[1] = x
    c[2] = x*c[1] - half*c[0]

    for i in range(3,N_modes):
        c[i] = x*c[i-1] - quarter*c[i-2]

    return c


def new_shell_RDF(N_modes,Rc,radii,N_reps,N_t,t):

    shell_RDF = dict()

    shell_RDF['N_modes'] = N_modes
    shell_RDF['Rc'] = Rc

    shell_RDF['N_reps'] = N_reps
    shell_RDF['N_t'] = N_t

    N_shells = radii.shape[0]

    shell_RDF['N_shells'] = N_shells

    shell_lims = np.zeros((N_shells,2),dtype=np.int32)

    for i_rad in range(N_shells):

        if i_rad == 0:

            shell_lims[i_rad,0] = 0
            shell_lims[i_rad,1] = radii[i_rad]

        else:

            shell_lims[i_rad,0] = radii[i_rad-1]
            shell_lims[i_rad,1] = radii[i_rad]

    shell_RDF['t'] = t

    shell_RDF['shell_lims'] = shell_lims

    shell_volumes = (4.0/3.0)*np.pi*(np.power(shell_lims[:,1],3.0)-np.power(shell_lims[:,0],3.0))

    shell_RDF['shell_volumes'] = shell_volumes
    shell_RDF['total_volume'] = np.sum(shell_volumes)

    shell_RDF['ribo_counts'] = np.zeros((N_shells,N_reps),dtype=np.int32)

    shell_RDF['coeffs'] = np.zeros((N_shells,N_reps,N_t,N_modes),dtype=np.double)

    return shell_RDF

def fill_reps_shell_RDFs(shell_RDF,traj,i_rep):

    t0 = str(shell_RDF['t'][0])
    
    x0_ribo = traj[t0]['ribo']['x']

    x0_R2 = np.sum(np.power(x0_ribo,2),axis=1)

    print(x0_R2.shape)

    shell_assignments = []

    for i_shell in range(shell_RDF['N_shells']):

        shell_assignments.append(np.argwhere((x0_R2 >= np.power(shell_RDF['shell_lims'][i_shell,0],2.0)) &
                                             (x0_R2 < np.power(shell_RDF['shell_lims'][i_shell,1],2.0))).flatten())

        shell_RDF['ribo_counts'][i_shell,i_rep] = shell_assignments[i_shell].shape[0]

    for i_t in range(1,shell_RDF['N_t']):

        t_temp = str(shell_RDF['t'][i_t])
        x_DNA = traj[t_temp]['DNA']['x']
        x_ribo = traj[t_temp]['ribo']['x']

        rho = traj[t_temp]['DNA']['N']/shell_RDF['total_volume']

        for i_shell in range(shell_RDF['N_shells']):

            if shell_MSD['ribo_counts'][i_shell,i_rep] > 0:

                x_ribo_shell = x_ribo[shell_assignments[i_shell],:]

                for i_ribo in range(x_ribo_shell.shape[0]):

                    xr = x_ribo_shell[i,:]

                    dx = x_DNA - np.tile((x_DNA.shape[0],1))

                    r = np.sum(np.power(dx,2.0),axis=1)

                    r = r[np.argwhere(r < shell_RDF['Rc'])]

                    z = (2*r - shell_RDF['Rc'])/(shell_RDF['Rc'])

                    


                

    return shell_RDF


def write_shell_RDF(shell_RDF_file,shell_RDF):

    with open(shell_RDF_file,'wb') as f:
        pickle.dump(shell_RDF,f)

    return

def read_shell_RDF(shell_RDF_file):

    with open(shell_RDF_file,'rb') as f:
        shell_MSD = pickle.load(f)

    return shell_RDF
