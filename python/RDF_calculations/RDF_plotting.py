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


def cheb_modes_eval(N_modes,x):

    T = np.zeros((N_modes),dtype=np.double)

    T[0] = 1.0
    T[1] = x

    for i in range(2,N_modes):
        T[i] = 2*x*T[i-1] - T[i-2]

    return T

def cheb_ortho_W(N_modes):

    W = np.ones((N_modes),dtype=np.double)

    W[0] = np.pi

    W[1:N_modes] = np.pi/2.0

    return W


def new_single_shell_RDF(N_modes,Rc,radii,N_t,t):

    shell_RDF = dict()

    shell_RDF['N_modes'] = N_modes
    shell_RDF['Rc'] = Rc
    shell_RDF['local_volume'] = (4.0/3.0)*np.pi*np.power(Rc,3.0)

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

    shell_RDF['ribo_counts'] = np.zeros((N_shells),dtype=np.int32)

    shell_RDF['coeffs'] = np.zeros((N_shells,N_t,N_modes),dtype=np.double)
    shell_RDF['total_densities'] = np.zeros((N_t),dtype=np.double)
    shell_RDF['local_densities'] = np.zeros((N_shells,N_t),dtype=np.double)

    return shell_RDF


def fill_shell_RDF(shell_RDF,traj):

    t0 = str(shell_RDF['t'][0])
    
    x0_ribo = traj[t0]['ribo']['x']

    x0_R2 = np.sum(np.power(x0_ribo,2),axis=1)

    print(x0_R2.shape)

    W = cheb_ortho_W(shell_RDF['N_modes'])

    shell_assignments = []

    for i_shell in range(shell_RDF['N_shells']):

        shell_assignments.append(np.argwhere((x0_R2 >= np.power(shell_RDF['shell_lims'][i_shell,0],2.0)) &
                                             (x0_R2 < np.power(shell_RDF['shell_lims'][i_shell,1],2.0))).flatten())

        shell_RDF['ribo_counts'][i_shell] = shell_assignments[i_shell].shape[0]

    for i_t in range(1,shell_RDF['N_t']):

        t_temp = str(shell_RDF['t'][i_t])
        x_DNA = traj[t_temp]['DNA']['x']
        x_ribo = traj[t_temp]['ribo']['x']

        rho = traj[t_temp]['DNA']['N']/shell_RDF['total_volume']

        shell_RDF['total_densities'][i_t] = rho

        print('rho = '+str(rho))

        for i_shell in range(shell_RDF['N_shells']):

            print("i_shell = "+str(i_shell))
            
            if shell_RDF['ribo_counts'][i_shell] > 0:

                x_ribo_shell = x_ribo[shell_assignments[i_shell],:]

                N_DNA_per_ribo = 0
                N_ribo_per_shell = x_ribo_shell.shape[0]

                temp_coeffs_ribo = np.zeros((shell_RDF['N_modes']),dtype=np.double)
                
                for i_ribo in range(N_ribo_per_shell):

                    xr = x_ribo_shell[i_ribo,:]

                    dx = x_DNA - xr

                    rsqrd = np.sum(np.power(dx,2.0),axis=1)

                    r = np.sqrt(rsqrd)

                    Rc_filter = np.argwhere(r < shell_RDF['Rc'])

                    r = r[Rc_filter]

                    z = (2*r - shell_RDF['Rc'])/(shell_RDF['Rc'])
                    z_denom = np.power(1.0+z,-2.0)

                    w = np.reciprocal(np.sqrt(1.0-np.power(z,2.0)))

                    i_N_DNA_per_ribo = r.shape[0]

                    temp_coeffs_DNA = np.zeros((shell_RDF['N_modes']),dtype=np.double)
                    
                    for i_DNA in range(i_N_DNA_per_ribo):

                        c = cheb_modes_eval(shell_RDF['N_modes'],z[i_DNA])

                        c = w[i_DNA]*c

                        c = z_denom[i_DNA]*c

                        temp_coeffs_DNA += c

                    temp_coeffs_DNA = temp_coeffs_DNA/i_N_DNA_per_ribo

                    temp_coeffs_ribo += temp_coeffs_DNA

                    N_DNA_per_ribo += i_N_DNA_per_ribo

                N_DNA_per_ribo = N_DNA_per_ribo/N_ribo_per_shell
                temp_coeffs_ribo = temp_coeffs_ribo/N_ribo_per_shell

                local_rho = N_DNA_per_ribo/shell_RDF['local_volume']
                print('local_rho = '+str(local_rho))
                shell_RDF['local_densities'][i_shell,i_t] = local_rho

                a = 8.0/3.0
                a = a*(local_rho/rho)
                temp_coeffs_ribo = a*temp_coeffs_ribo
                temp_coeffs_ribo = np.divide(temp_coeffs_ribo,W)

                shell_RDF['coeffs'][i_shell,i_t,:] = temp_coeffs_ribo

    return shell_RDF


def new_merged_shell_RDF(N_modes,Rc,radii,N_reps,N_t,t):

    shell_RDF = dict()

    shell_RDF['N_modes'] = N_modes
    shell_RDF['Rc'] = Rc
    shell_RDF['local_volume'] = (4.0/3.0)*np.pi*np.power(Rc,3.0)

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
    shell_RDF['total_densities'] = np.zeros((N_reps,N_t),dtype=np.double)
    shell_RDF['local_densities'] = np.zeros((N_shells,N_reps,N_t),dtype=np.double)

    return shell_RDF


def fill_merged_shell_RDF(merged_shell_RDF,single_shell_RDF,i_rep):

    merged_shell_RDF['ribo_counts'][:,i_rep] = single_shell_RDF['ribo_counts']
    merged_shell_RDF['coeffs'][:,i_rep,:,:] = single_shell_RDF['coeffs']
    merged_shell_RDF['total_densities'][i_rep,:] = single_shell_RDF['total_densities']
    merged_shell_RDF['local_densities'][:,i_rep,:] = single_shell_RDF['local_densities']

    return merged_shell_RDF


def new_shell_RDF(N_modes,Rc,radii,N_reps,N_t,t):

    shell_RDF = dict()

    shell_RDF['N_modes'] = N_modes
    shell_RDF['Rc'] = Rc
    shell_RDF['local_volume'] = (4.0/3.0)*np.pi*np.power(Rc,3.0)

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
    shell_RDF['total_densities'] = np.zeros((N_reps,N_t),dtype=np.double)
    shell_RDF['local_densities'] = np.zeros((N_shells,N_reps,N_t),dtype=np.double)

    return shell_RDF

def fill_reps_shell_RDFs(shell_RDF,traj,i_rep):

    t0 = str(shell_RDF['t'][0])
    
    x0_ribo = traj[t0]['ribo']['x']

    x0_R2 = np.sum(np.power(x0_ribo,2),axis=1)

    print(x0_R2.shape)

    W = cheb_ortho_W(shell_RDF['N_modes'])

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

        shell_RDF['total_densities'][i_rep,i_t] = rho

        print('rho = '+str(rho))

        for i_shell in range(shell_RDF['N_shells']):

            print("i_shell = "+str(i_shell))
            
            if shell_RDF['ribo_counts'][i_shell,i_rep] > 0:

                x_ribo_shell = x_ribo[shell_assignments[i_shell],:]

                N_DNA_per_ribo = 0
                N_ribo_per_shell = x_ribo_shell.shape[0]

                temp_coeffs_ribo = np.zeros((shell_RDF['N_modes']),dtype=np.double)
                
                for i_ribo in range(N_ribo_per_shell):

                    xr = x_ribo_shell[i_ribo,:]

                    dx = x_DNA - xr

                    rsqrd = np.sum(np.power(dx,2.0),axis=1)

                    r = np.sqrt(rsqrd)

                    Rc_filter = np.argwhere(r < shell_RDF['Rc'])

                    r = r[Rc_filter]

                    z = (2*r - shell_RDF['Rc'])/(shell_RDF['Rc'])
                    z_denom = np.power(1.0+z,-2.0)

                    w = np.reciprocal(np.sqrt(1.0-np.power(z,2.0)))

                    i_N_DNA_per_ribo = r.shape[0]

                    temp_coeffs_DNA = np.zeros((shell_RDF['N_modes']),dtype=np.double)
                    
                    for i_DNA in range(i_N_DNA_per_ribo):

                        c = cheb_modes_eval(shell_RDF['N_modes'],z[i_DNA])

                        c = w[i_DNA]*c

                        c = z_denom[i_DNA]*c

                        temp_coeffs_DNA += c

                    temp_coeffs_DNA = temp_coeffs_DNA/i_N_DNA_per_ribo

                    temp_coeffs_ribo += temp_coeffs_DNA

                    N_DNA_per_ribo += i_N_DNA_per_ribo

                N_DNA_per_ribo = N_DNA_per_ribo/N_ribo_per_shell
                temp_coeffs_ribo = temp_coeffs_ribo/N_ribo_per_shell

                local_rho = N_DNA_per_ribo/shell_RDF['local_volume']
                print('local_rho = '+str(local_rho))
                shell_RDF['local_densities'][i_shell,i_rep,i_t] = local_rho

                a = 8.0/3.0
                a = a*(local_rho/rho)
                temp_coeffs_ribo = a*temp_coeffs_ribo
                temp_coeffs_ribo = np.divide(temp_coeffs_ribo,W)

                shell_RDF['coeffs'][i_shell,i_rep,i_t,:] = temp_coeffs_ribo

    return shell_RDF


def write_shell_RDF(shell_RDF_file,shell_RDF):

    with open(shell_RDF_file,'wb') as f:
        pickle.dump(shell_RDF,f)

    return

def read_shell_RDF(shell_RDF_file):

    with open(shell_RDF_file,'rb') as f:
        shell_RDF = pickle.load(f)

    return shell_RDF


def eval_cheb_RDF(Rc,a,x):

    z = (2*x-Rc)/Rc

    y = np.polynomial.chebyshev.chebval(z,a)

    return y

def plot_shell_RDFs(fig_file,shell_RDF,res):

    cmap = colormaps.get_cmap('tab10')

    print(shell_RDF['total_volume'])
    
    #fig_size = [87,87]
    #fig_size = [174,174]
    fig_size = [174,87]

    fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

    ax = plt.gca()
    sigma = np.power(2.0,1/6.0)*11.7

    ax.set_xlabel(r'$r$ - radial distance [nm]', fontsize=12)
    ax.set_ylabel(r'$g_{ribo-DNA}(r)$ - R.D.F.', fontsize=12)


    tick_length = 4.0
    tick_width = 2.0
    ax.tick_params(axis='x',labelsize=9,
                   length=tick_length,
                   width=tick_width)
    ax.tick_params(axis='y',labelsize=9,
                   length=tick_length,
                   width=tick_width)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_linewidth(2.0)
    ax.spines['bottom'].set_linewidth(2.0)

    # ax.set_xlim(xmin=0.0,xmax=shell_RDF['Rc']//10)
    ax.set_xlim(xmin=0.0,xmax=40)
    ax.set_ylim(ymin=0.0,ymax=1.1) # restore this for Figure 5
    # ax.set_ylim(ymin=0.0,ymax=1.8)
    
    
    for i_shell in range(shell_RDF['N_shells']):

        temp_color = cmap(i_shell%shell_RDF['N_shells'])
        temp_label = r'{:d} - {:d} nm'.format(shell_RDF['shell_lims'][i_shell,0]//10,
                                                                         shell_RDF['shell_lims'][i_shell,1]//10)

        x = np.linspace(0.0,shell_RDF['Rc'],res,dtype=np.double)

        a = np.mean(shell_RDF['coeffs'][i_shell,:,:,:],axis=0)
        a = np.mean(a,axis=0)

        print(a.shape)


        y = eval_cheb_RDF(shell_RDF['Rc'],a,x)
        sigma_val = eval_cheb_RDF(shell_RDF['Rc'],a,10.0*sigma)

        # lw 2.0
        ax.plot(x/10,y,
                lw=2.5,
                alpha=1.0,
                ls='-',
                c=temp_color,
                label=temp_label)

        x_sigma = np.array([sigma,sigma,0.0],dtype=np.double)
        y_sigma = np.array([0.0,sigma_val,sigma_val],dtype=np.double)

        # lw 1.5
        ax.plot(x_sigma,y_sigma,
                lw=2.0,
                alpha=0.75,
                ls='--',
                c=temp_color,
                zorder=-1)


    sigma_text = r'$2^{1/6}\sigma_{DNA-ribo}=$'+'${:.1f}$ nm'.format(sigma)
    ax.annotate(sigma_text,
                xy=(sigma+1.0,0.1),
                fontsize=12)
    ax.legend(fontsize=8,
              loc='lower right')

    plt.tight_layout()

    fig.savefig(fig_file,dpi=300)

    plt.close()
    
    return
