import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as colormaps
import matplotlib.patches as patches
from matplotlib.lines import Line2D

from mpl_toolkits.axes_grid1 import make_axes_locatable

import os

import pandas as pd

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


def read_MSDs(in_dir,in_label,rep_flag,min_rep,max_rep):

    MSDs = dict()

    N_reps = max_rep - min_rep + 1

    if rep_flag == False:

        in_file = in_dir + in_label + '_MSD.json'

        with open(in_file,'r') as f:
    
            df = pd.read_json(f,lines=True)

        ts = df['timestep'].to_numpy(dtype=np.int32)

        for key in df.keys()[1:]:

            MSDs[key] = df[key].to_numpy(dtype=np.float32)

    else:

        for rep in range(min_rep,max_rep+1):

            in_file = in_dir + in_label + '_rep'+str(rep).zfill(5) + '_MSD.json'

            with open(in_file,'r') as f:
    
                df = pd.read_json(f,lines=True)

            if rep == min_rep:

                ts = df['timestep'].to_numpy(dtype=np.int32)

                for key in df.keys()[1:]:

                    MSDs[key] = df[key].to_numpy(dtype=np.float32)

            else:

                for key in df.keys()[1:]:

                    MSDs[key] += df[key].to_numpy(dtype=np.float32)

        for key in MSDs.keys():

            MSDs[key] = MSDs[key]/N_reps
    
    return ts, MSDs

def multi_fits(ax,locus_size,regimes,contact_sum):
    
    for regime in reversed(regimes):

        X=np.ones((regime[1]-regime[0],2),dtype=np.double,order='F')

        z=(np.arange(regime[0],regime[1])*(locus_size/1000.0))
        if regime[0] == 0:
            z += 1

        X[:,1]=np.log(z)

        y=np.log(contact_sum[regime[0]:regime[1]])

        y=np.matmul(X.T,y)

        X=np.matmul(X.T,X)

        beta=np.matmul(np.linalg.inv(X),y)

        print(beta[1])

        y=beta[0]+beta[1]*np.log(z)

        y=np.exp(y)

        if regime[0] == 0:
            z[0]=0
            
        z=np.rint(z).astype(np.int64)
            
        temp_label=r'{:d}-{:d} kbp: '.format(z[0],z[-1])+'$s=-$'+"{:.3f}".format(-beta[1])
        
        ax.plot(z,y,ls="-",c='w',lw=2.0,alpha=0.65)
        ax.plot(z,y,ls="--",lw=1.5,alpha=0.9,label=temp_label)
    
    return

def multi_fits_linear(ax,lims,x,y):
    
    A_to_micro = 1.0E+4 # A/mu
    ns_to_s = 1.0E+9 # ns/s
    
    N_lims = lims.shape[0]
    
    for i_lim in range(N_lims):
        
        i_ll = np.argwhere(x>=lims[i_lim,0])[0][0]
        i_ul = np.argwhere(x<lims[i_lim,1])[-1][0]
        
        print(i_ll)
        print(i_ul)
        
        ll = x[i_ll]
        ul = x[i_ul]
        
        x_temp = x[i_ll:i_ul+1]
        y_temp = y[i_ll:i_ul+1]
        
        
        X = np.ones((x_temp.shape[0],2),dtype=np.float32)
        X[:,1] = x_temp
        
        print(X.shape)
        print(y_temp.shape)
        
        ls_out = np.linalg.lstsq(X,y_temp)
        
        beta = ls_out[0]
        
        print(beta) 
            
        temp_label=r'{:.1E} to {:.1E} [ns]: '.format(int(x_temp[0]),int(x_temp[-1]))+'$D=${:.3E} [$\mu$m$^2$/s]'.format((beta[1]*ns_to_s)/(6.0*(A_to_micro**2.0)))
        
        ax.plot(x_temp,beta[1]*x_temp+beta[0],ls="-",c='w',lw=2.0,alpha=0.65)
        ax.plot(x_temp,beta[1]*x_temp+beta[0],ls="--",lw=1.5,alpha=0.9,label=temp_label)
    
    return

def multi_fits_loglog(ax,lims,x,y):
    
    N_lims = lims.shape[0]
    
    for i_lim in range(N_lims):
        
        i_ll = np.argwhere(x>=lims[i_lim,0])[0][0]
        i_ul = np.argwhere(x<lims[i_lim,1])[-1][0]
        
        print(i_ll)
        print(i_ul)
        
        ll = x[i_ll]
        ul = x[i_ul]
        
        x_temp_nl = x[i_ll:i_ul+1]/x[0]
        x_temp = np.log(x[i_ll:i_ul+1]/x[0])
        y_temp = np.log(y[i_ll:i_ul+1]/y[0])
        
        
        X = np.ones((x_temp.shape[0],2),dtype=np.float32)
        X[:,1] = x_temp
        
        print(X.shape)
        print(y_temp.shape)
        
        ls_out = np.linalg.lstsq(X,y_temp)
        
        beta = ls_out[0]
        
        print('beta=\n')
        print(beta) 
            
        temp_label=r'$\alpha=$ {:.3E}'.format(beta[1])
        
        ax.plot(x_temp_nl,np.exp(beta[0]+beta[1]*x_temp),ls="-",c='w',lw=2.0,alpha=0.65)
        ax.plot(x_temp_nl,np.exp(beta[0]+beta[1]*x_temp),ls="--",lw=1.5,alpha=0.9,label=temp_label)
    
    return

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
    
def plot_diffusion_laws(fig_file,t,MSD,lims):
    
    fig_size = [87,87]

    fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

    ax = plt.gca()

    ax.set_xlabel(r'$\log[t/t_0]$', fontsize=8)
    ax.set_ylabel(r'$\log\big[$MSD$(t)/$MSD$(t_0)\big]$', fontsize=8)
    #ax.set_title(r'MSD of 10~bp DNA monomers',fontsize=10)

    print(t)
    print(MSD)

    temp_marker_style = dict(marker='.', markersize=8,fillstyle='none')

    ax.plot(t/t[0],MSD/MSD[0],c='k',lw=1.0,alpha=0.35,**temp_marker_style)

    multi_fits_loglog(ax,lims,t,MSD)

    ax.tick_params(axis='x',labelsize=8)
    ax.tick_params(axis='y',labelsize=8)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    ax.set_xscale('log')
    ax.set_yscale('log')

    plt.tight_layout()

    ax.legend(title=r'$\frac{MSD(t)}{MSD(t_0)}\propto\big(t/t_0\big)^{\alpha}$',title_fontsize=7,fontsize=6)

    #plt.grid()

    fig.savefig(fig_file,dpi=300)

def new_shell_MSD(radii,N_reps,N_t,t):

    shell_MSD = dict()

    shell_MSD['N_reps'] = N_reps
    shell_MSD['N_t'] = N_t

    N_shells = radii.shape[0]

    shell_MSD['N_shells'] = N_shells

    shell_lims = np.zeros((N_shells,2),dtype=np.int32)

    for i_rad in range(N_shells):

        if i_rad == 0:

            shell_lims[i_rad,0] = 0
            shell_lims[i_rad,1] = radii[i_rad]

        else:

            shell_lims[i_rad,0] = radii[i_rad-1]
            shell_lims[i_rad,1] = radii[i_rad]

    shell_MSD['t'] = t

    shell_MSD['shell_lims'] = shell_lims

    shell_volumes = (4.0/3.0)*np.pi*(np.power(shell_lims[:,1],3.0)-np.power(shell_lims[:,0],3.0))

    shell_MSD['shell_volumes'] = shell_volumes

    shell_MSD['DNA_density_reps'] = np.zeros((N_shells,N_reps),dtype=np.float32)
    shell_MSD['ribo_density_reps'] = np.zeros((N_shells,N_reps),dtype=np.float32)
    

    shell_MSD['DNA_counts'] = np.zeros((N_shells,N_reps),dtype=np.int32)
    shell_MSD['ribo_counts'] = np.zeros((N_shells,N_reps),dtype=np.int32)

    shell_MSD['DNA_MSDs_reps'] = np.zeros((N_shells,N_reps,N_t),dtype=np.float32)
    shell_MSD['ribo_MSDs_reps'] = np.zeros((N_shells,N_reps,N_t),dtype=np.float32)

    shell_MSD['DNA_BrownianDiff'] = np.zeros((N_shells,N_reps,2),dtype=np.float32)
    shell_MSD['ribo_BrownianDiff'] = np.zeros((N_shells,N_reps,2),dtype=np.float32)
    
    shell_MSD['DNA_AnomalousDiff'] = np.zeros((N_shells,N_reps,2),dtype=np.float32)
    shell_MSD['ribo_AnomalousDiff'] = np.zeros((N_shells,N_reps,2),dtype=np.float32)

    return shell_MSD

def fill_reps_shell_MSDs(shell_MSD,traj,i_rep):

    t0 = str(shell_MSD['t'][0])

    x0_DNA = traj[t0]['DNA']['x']

    x0_R2 = np.sum(np.power(x0_DNA,2),axis=1)

    print(x0_R2.shape)

    shell_assignments = []

    for i_shell in range(shell_MSD['N_shells']):

        shell_assignments.append(np.argwhere((x0_R2 >= np.power(shell_MSD['shell_lims'][i_shell,0],2.0)) &
                                             (x0_R2 < np.power(shell_MSD['shell_lims'][i_shell,1],2.0))).flatten())

        shell_MSD['DNA_counts'][i_shell,i_rep] = shell_assignments[i_shell].shape[0]

    for i_t in range(1,shell_MSD['N_t']):

        x_DNA = traj[str(shell_MSD['t'][i_t])]['DNA']['x']

        for i_shell in range(shell_MSD['N_shells']):

            if shell_MSD['DNA_counts'][i_shell,i_rep] > 0:

                dx_DNA = x_DNA[shell_assignments[i_shell],:] - x0_DNA[shell_assignments[i_shell],:]
                # print(dx_DNA.shape)
                
                shell_MSD['DNA_MSDs_reps'][i_shell,i_rep,i_t] = np.mean(np.sum(np.power(dx_DNA,2.0),axis=1))

    
    x0_ribo = traj[t0]['ribo']['x']

    x0_R2 = np.sum(np.power(x0_ribo,2),axis=1)

    print(x0_R2.shape)

    shell_assignments = []

    for i_shell in range(shell_MSD['N_shells']):

        shell_assignments.append(np.argwhere((x0_R2 >= np.power(shell_MSD['shell_lims'][i_shell,0],2.0)) &
                                             (x0_R2 < np.power(shell_MSD['shell_lims'][i_shell,1],2.0))).flatten())

        shell_MSD['ribo_counts'][i_shell,i_rep] = shell_assignments[i_shell].shape[0]

    for i_t in range(1,shell_MSD['N_t']):

        x_ribo = traj[str(shell_MSD['t'][i_t])]['ribo']['x']

        for i_shell in range(shell_MSD['N_shells']):

            if shell_MSD['ribo_counts'][i_shell,i_rep] > 0:

                dx_ribo = x_ribo[shell_assignments[i_shell],:] - x0_ribo[shell_assignments[i_shell],:]
                # print(dx_ribo.shape)

                shell_MSD['ribo_MSDs_reps'][i_shell,i_rep,i_t] = np.mean(np.sum(np.power(dx_ribo,2.0),axis=1))

    shell_MSD['DNA_density_reps'][:,i_rep] = shell_MSD['DNA_counts'][:,i_rep]/shell_MSD['shell_volumes']
    shell_MSD['ribo_density_reps'][:,i_rep] = shell_MSD['ribo_counts'][:,i_rep]/shell_MSD['shell_volumes']

    return shell_MSD

def shell_MSD_calc_diffusion(shell_MSD,dt,lims):

    x = shell_MSD['t']*dt
    
    # i_ll = np.argwhere(x>=lims[0])[0][0]
    # i_ul = np.argwhere(x<lims[1])[-1][0]

    i_ll = lims[0]
    i_ul = lims[1]
        
    ll = x[i_ll]
    ul = x[i_ul]
        
    x_temp = x[i_ll:i_ul+1]

    shell_MSD['fit_range'] = x_temp

    for i_shell in range(shell_MSD['N_shells']):

        for i_rep in range(shell_MSD['N_reps']):

            if shell_MSD['DNA_counts'][i_shell,i_rep] > 0:

                y_temp = shell_MSD['DNA_MSDs_reps'][i_shell,i_rep,i_ll:i_ul+1]

                X = np.ones((x_temp.shape[0],2),dtype=np.float32)
                X[:,1] = x_temp
        
                ls_out = np.linalg.lstsq(X,y_temp)

                shell_MSD['DNA_BrownianDiff'][i_shell,i_rep,:] = np.copy(ls_out[0])

            if shell_MSD['ribo_counts'][i_shell,i_rep] > 0:

                y_temp = shell_MSD['ribo_MSDs_reps'][i_shell,i_rep,i_ll:i_ul+1]

                X = np.ones((x_temp.shape[0],2),dtype=np.float32)
                X[:,1] = x_temp
        
                ls_out = np.linalg.lstsq(X,y_temp)

                shell_MSD['ribo_BrownianDiff'][i_shell,i_rep,:] = np.copy(ls_out[0])

    x = shell_MSD['t'][1:]*dt
    x0 = x[0]
    print(x0)
    
    # i_ll = np.argwhere(x>=lims[0])[0][0]
    # i_ul = np.argwhere(x<lims[1])[-1][0]
        
    # ll = x[i_ll]
    # ul = x[i_ul]
        
    x_temp = x[(i_ll-1):i_ul]/x0
    shell_MSD['fit_range_law'] = x_temp
    x_temp = np.log(x_temp)

    print(x_temp.shape)

    for i_shell in range(shell_MSD['N_shells']):

        for i_rep in range(shell_MSD['N_reps']):

            if shell_MSD['DNA_counts'][i_shell,i_rep] > 0:

                y_temp = shell_MSD['DNA_MSDs_reps'][i_shell,i_rep,i_ll:i_ul+1]/shell_MSD['DNA_MSDs_reps'][i_shell,i_rep,1]
                y_temp = np.log(y_temp)
                print(i_ll)
                print(i_ul)
                print(y_temp.shape)

                X = np.ones((x_temp.shape[0],2),dtype=np.float32)
                X[:,1] = x_temp
        
                ls_out = np.linalg.lstsq(X,y_temp)

                shell_MSD['DNA_AnomalousDiff'][i_shell,i_rep,:] = np.copy(ls_out[0])

            if shell_MSD['ribo_counts'][i_shell,i_rep] > 0:

                y_temp = shell_MSD['ribo_MSDs_reps'][i_shell,i_rep,i_ll:i_ul+1]/shell_MSD['ribo_MSDs_reps'][i_shell,i_rep,1]
                y_temp = np.log(y_temp)

                X = np.ones((x_temp.shape[0],2),dtype=np.float32)
                X[:,1] = x_temp
        
                ls_out = np.linalg.lstsq(X,y_temp)

                shell_MSD['ribo_AnomalousDiff'][i_shell,i_rep,:] = np.copy(ls_out[0])

    return shell_MSD

def plot_shell_MSD_reps(fig_file,shell_MSD,dt,lims,species_select):

    A_to_micro = 1.0E+4 # A/mu
    ns_to_s = 1.0E+9 # ns/s

    cmap = colormaps.get_cmap('tab10')

    fig_size = [87,87]
    #fig_size = [174,174]

    fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

    ax = plt.gca()

    ax.set_xlabel(r'${\tau}$ - Time [ns]', fontsize=12)
    ax.set_ylabel(r'$\langle|\mathbf{x}(t+{\tau})-\mathbf{x}(t)|^2\rangle$ - MSD [\AA$^2$]', fontsize=12)
    #ax.set_title(r'MSD of 10~bp DNA monomers',fontsize=10)

    temp_marker_style = dict(marker='.', markersize=6,fillstyle='none')

    if species_select == 'DNA':
        MSD = shell_MSD['DNA_MSDs_reps']
        counts = shell_MSD['DNA_counts']
        BrownianDiff = shell_MSD['DNA_BrownianDiff']
    elif species_select == 'ribo':
        MSD = shell_MSD['ribo_MSDs_reps']
        counts = shell_MSD['ribo_counts']
        BrownianDiff = shell_MSD['ribo_BrownianDiff']

    t = shell_MSD['t']*dt

    for i_shell in range(shell_MSD['N_shells']):

        temp_color = cmap(i_shell%shell_MSD['N_shells'])

        for i_rep in range(shell_MSD['N_reps']):

            if counts[i_shell,i_rep] > 0:
                #ax.plot(t,MSD[i_shell,i_rep],c=temp_color,lw=1.0,alpha=0.3,**temp_marker_style,zorder=0)
                ax.plot(t,MSD[i_shell,i_rep],c=temp_color,lw=1.0,alpha=0.3,zorder=0)

    lims = np.array(lims,dtype=np.float32)

    # multi_fits_linear(ax,lims,t,MSD)

    ax.tick_params(axis='x',labelsize=9)
    ax.tick_params(axis='y',labelsize=9)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    legend_elements = []

    for i_shell in range(shell_MSD['N_shells']):

        valid_reps = np.argwhere(counts[i_shell] > 0).flatten()

        beta = BrownianDiff[i_shell,valid_reps,:]

        print(beta.shape)

        mean_beta = np.mean(beta,axis=0)
        stddev_beta = np.std(beta,axis=0)

        print(mean_beta.shape)
        print(stddev_beta.shape)

        temp_color = cmap(i_shell%shell_MSD['N_shells'])
        #temp_label = r'Shell \#'+str(i_shell+1) + ', '
        temp_label = r'{:d} - {:d} nm: '.format(shell_MSD['shell_lims'][i_shell,0]//10,
                                                                         shell_MSD['shell_lims'][i_shell,1]//10)
        # temp_label += '$D=${:.3E}\pm${:.3E} [$\mu$m$^2$/s]'.format((mean_beta[1]*ns_to_s)/(6.0*(A_to_micro**2.0)),
        #                                                            (stddev_beta[1]*ns_to_s)/(6.0*(A_to_micro**2.0)))
        #temp_label += '$D=$ {:.3E} [$\mu$m$^2$/s]'.format((mean_beta[1]*ns_to_s)/(6.0*(A_to_micro**2.0)))
        temp_label += '$D=$ {:.3E}'.format((mean_beta[1]*ns_to_s)/(6.0*(A_to_micro**2.0)))

        x_fit = shell_MSD['fit_range']
        y_fit = mean_beta[1]*shell_MSD['fit_range']+mean_beta[0]
        
        ax.plot(x_fit,y_fit,
                lw=3.0,
                alpha=0.65,
                ls='-',
                c='w')
        ax.plot(x_fit,y_fit,
                lw=2.5,
                alpha=0.95,
                ls='--',
                label=temp_label,
                c=temp_color)

        # legend_elements.append(Line2D([0], [0],
        #                               color=temp_color,
        #                               label=temp_label,
        #                               lw=2.0))

    #ax.legend(handles=legend_elements,fontsize=6)
    #ax.legend(fontsize=6)
    ax.legend(title=r'Diffusion Constant, $D$ - [$\mu$m$^2$/s]',title_fontsize=9,fontsize=8)

    plt.tight_layout()

    #plt.grid()

    new_fig_file = fig_file.strip('.pdf')
    new_fig_file = new_fig_file + '_' + species_select + '.pdf'

    fig.savefig(new_fig_file,dpi=300)

    return

def plot_shell_law_reps(fig_file,shell_MSD,dt,lims,species_select):

    A_to_micro = 1.0E+4 # A/mu
    ns_to_s = 1.0E+9 # ns/s

    cmap = colormaps.get_cmap('tab10')

    fig_size = [87,87]
    #fig_size = [174,174]

    fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

    ax = plt.gca()

    ax.set_xlabel(r'$\log[(t+{\tau})/t]$', fontsize=12)
    ax.set_ylabel(r'$\log\big[$MSD$(t+{\tau})/$MSD$(t)\big]$', fontsize=12)
    #ax.set_title(r'MSD of 10~bp DNA monomers',fontsize=10)

    temp_marker_style = dict(marker='.', markersize=6,fillstyle='none')

    t = shell_MSD['t'][1:]*dt
    
    if species_select == 'DNA':
        MSD = shell_MSD['DNA_MSDs_reps'][:,:,1:]
        counts = shell_MSD['DNA_counts']
        AnomalousDiff = shell_MSD['DNA_AnomalousDiff']
    elif species_select == 'ribo':
        MSD = shell_MSD['ribo_MSDs_reps'][:,:,1:]
        counts = shell_MSD['ribo_counts']
        AnomalousDiff = shell_MSD['ribo_AnomalousDiff']

    for i_shell in range(shell_MSD['N_shells']):

        temp_color = cmap(i_shell%shell_MSD['N_shells'])

        for i_rep in range(shell_MSD['N_reps']):

            if counts[i_shell,i_rep] > 0:
                #ax.plot(t,MSD[i_shell,i_rep],c=temp_color,lw=1.0,alpha=0.3,**temp_marker_style,zorder=0)
                ax.plot(t/t[0],MSD[i_shell,i_rep,:]/MSD[i_shell,i_rep,0],c=temp_color,lw=1.0,alpha=0.3,zorder=0)

    lims = np.array(lims,dtype=np.float32)

    # multi_fits_linear(ax,lims,t,MSD)

    ax.tick_params(axis='x',labelsize=9)
    ax.tick_params(axis='y',labelsize=9)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    ax.set_xscale('log')
    ax.set_yscale('log')

    legend_elements = []

    for i_shell in range(shell_MSD['N_shells']):

        valid_reps = np.argwhere(counts[i_shell] > 0).flatten()

        beta = AnomalousDiff[i_shell,valid_reps,:]

        print(beta.shape)

        mean_beta = np.mean(beta,axis=0)
        stddev_beta = np.std(beta,axis=0)

        print(mean_beta.shape)
        print(stddev_beta.shape)

        temp_color = cmap(i_shell%shell_MSD['N_shells'])
        #temp_label = r'Shell \#'+str(i_shell+1) + ', '
        temp_label = r'{:d} - {:d} nm: '.format(shell_MSD['shell_lims'][i_shell,0]//10,
                                                                         shell_MSD['shell_lims'][i_shell,1]//10)
        # temp_label += '$D=${:.3E}\pm${:.3E} [$\mu$m$^2$/s]'.format((mean_beta[1]*ns_to_s)/(6.0*(A_to_micro**2.0)),
        #                                                            (stddev_beta[1]*ns_to_s)/(6.0*(A_to_micro**2.0)))
        #temp_label += '${\alpha}$ {:.3E}'.format(mean_beta[1])
        #temp_label += ' {:.3E}'.format(mean_beta[1])
        temp_label += '$\\alpha=$ {:.3E}'.format(mean_beta[1])

        x_fit = shell_MSD['fit_range_law']
        y_fit = np.exp(mean_beta[1]*np.log(shell_MSD['fit_range_law'])+mean_beta[0])
        
        ax.plot(x_fit,y_fit,
                lw=3.0,
                alpha=0.65,
                ls='-',
                c='w')
        ax.plot(x_fit,y_fit,
                lw=2.5,
                alpha=0.95,
                ls='--',
                label=temp_label,
                c=temp_color)

        # legend_elements.append(Line2D([0], [0],
        #                               color=temp_color,
        #                               label=temp_label,
        #                               lw=2.0))

    #ax.legend(handles=legend_elements,fontsize=6)
    ax.legend(title=r'$\frac{MSD(t+{\tau})}{MSD(t)}\propto\big(\frac{t+{\tau}}{t}\big)^{\alpha}$',title_fontsize=9,fontsize=8)

    plt.tight_layout()

    #plt.grid()

    new_fig_file = fig_file.strip('.pdf')
    new_fig_file = new_fig_file + '_' + species_select + '.pdf'

    fig.savefig(new_fig_file,dpi=300)

    return

def plot_shell_diff_vs_density_reps(fig_file,shell_MSD,species_select):

    A_to_micro = 1.0E+4 # A/mu
    ns_to_s = 1.0E+9 # ns/s

    cmap = colormaps.get_cmap('tab10')

    fig_size = [87,87]
    #fig_size = [174,174]

    fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

    ax = plt.gca()

    ax.set_xlabel(r'$n_{DNA}$ - Density [\#/$\mu$m$^3$]', fontsize=12)
    ax.set_ylabel(r'$D$ - Brownian Diffusion Constant [$\mu$m$^2$/s]', fontsize=12)
    #ax.set_title(r'MSD of 10~bp DNA monomers',fontsize=10)

    temp_marker_style = dict(marker='.', markersize=6,fillstyle='none')

    if species_select == 'DNA':
        density = shell_MSD['DNA_density_reps']
        counts = shell_MSD['DNA_counts']
        BrownianDiff = shell_MSD['DNA_BrownianDiff']
    elif species_select == 'ribo':
        density = shell_MSD['DNA_density_reps']
        counts = shell_MSD['ribo_counts']
        BrownianDiff = shell_MSD['ribo_BrownianDiff']

    ax.tick_params(axis='x',labelsize=9)
    ax.tick_params(axis='y',labelsize=9)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    #temp_marker_style = dict(marker='.', markersize=8,fillstyle='none')

    for i_shell in range(shell_MSD['N_shells']):

        valid_reps = np.argwhere(counts[i_shell] > 0).flatten()

        beta = BrownianDiff[i_shell,valid_reps,:]

        mean_beta = np.mean(beta,axis=0)
        stddev_beta = np.std(beta,axis=0)

        temp_color = cmap(i_shell%shell_MSD['N_shells'])
        #temp_label = r'Shell \#'+str(i_shell+1) + ', '
        temp_label = r'{:d} - {:d} nm: '.format(shell_MSD['shell_lims'][i_shell,0]//10,
                                                                         shell_MSD['shell_lims'][i_shell,1]//10)
        x_temp = density[i_shell,valid_reps]*(A_to_micro**3.0)
        y_temp = (beta[:,1]*ns_to_s)/(6.0*(A_to_micro**2.0))

        cov = np.cov(x_temp,y_temp,dtype=np.float32)

        invcov = np.linalg.inv(cov)

        w, v = np.linalg.eig(cov)

        theta = np.rad2deg(np.arctan2(v[1,0],v[0,0]))
        
        xm = np.mean(x_temp)
        ym = np.mean(y_temp)

        cov_ellipse = patches.Ellipse((xm,ym),2*np.sqrt(w[0]),2*np.sqrt(w[1]),angle=theta,color=temp_color,alpha=0.25,zorder=0)

        ax.add_patch(cov_ellipse)

        temp_label += r'$r_{{x,y}}=$ {:.1E}'.format(cov[0,1]/np.sqrt(cov[0,0]*cov[1,1]))
        
        ax.scatter(x_temp,y_temp,
                   label=temp_label,
                   color=temp_color,
                   alpha=0.75,zorder=1)


    #ax.legend(handles=legend_elements,fontsize=6)
    ax.legend(fontsize=8)

    plt.tight_layout()

    #plt.grid()

    new_fig_file = fig_file.strip('.pdf')
    new_fig_file = new_fig_file + '_' + species_select + '.pdf'

    fig.savefig(new_fig_file,dpi=300)

    return

def plot_shell_law_vs_density_reps(fig_file,shell_MSD,species_select):

    A_to_micro = 1.0E+4 # A/mu
    ns_to_s = 1.0E+9 # ns/s

    cmap = colormaps.get_cmap('tab10')

    fig_size = [87,87]
    #fig_size = [174,174]

    fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

    ax = plt.gca()

    ax.set_xlabel(r'$n_{DNA}$ - Density [\#/$\mu$m$^3$]', fontsize=12)
    ax.set_ylabel(r'$\alpha$ - Anomalous Diffusion Law', fontsize=12)
    #ax.set_title(r'MSD of 10~bp DNA monomers',fontsize=10)

    temp_marker_style = dict(marker='.', markersize=6,fillstyle='none')

    if species_select == 'DNA':
        density = shell_MSD['DNA_density_reps']
        counts = shell_MSD['DNA_counts']
        AnomalousDiff = shell_MSD['DNA_AnomalousDiff']
    elif species_select == 'ribo':
        density = shell_MSD['DNA_density_reps']
        counts = shell_MSD['ribo_counts']
        AnomalousDiff = shell_MSD['ribo_AnomalousDiff']

    ax.tick_params(axis='x',labelsize=9)
    ax.tick_params(axis='y',labelsize=9)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    #temp_marker_style = dict(marker='.', markersize=8,fillstyle='none')

    for i_shell in range(shell_MSD['N_shells']):

        valid_reps = np.argwhere(counts[i_shell] > 0).flatten()

        beta = AnomalousDiff[i_shell,valid_reps,:]

        mean_beta = np.mean(beta,axis=0)
        stddev_beta = np.std(beta,axis=0)

        temp_color = cmap(i_shell%shell_MSD['N_shells'])
        #temp_label = r'Shell \#'+str(i_shell+1) + ', '
        temp_label = r'{:d} - {:d} nm: '.format(shell_MSD['shell_lims'][i_shell,0]//10,
                                                                         shell_MSD['shell_lims'][i_shell,1]//10)
        x_temp = density[i_shell,valid_reps]*(A_to_micro**3.0)
        y_temp = beta[:,1]

        cov = np.cov(x_temp,y_temp,dtype=np.float32)

        w, v = np.linalg.eig(cov)

        theta = np.rad2deg(np.arctan2(v[1,0],v[0,0]))
        
        xm = np.mean(x_temp)
        ym = np.mean(y_temp)

        cov_ellipse = patches.Ellipse((xm,ym),2*np.sqrt(w[0]),2*np.sqrt(w[1]),angle=theta,color=temp_color,alpha=0.25,zorder=0)

        ax.add_patch(cov_ellipse)

        temp_label += r'$r_{{x,y}}=$ {:.1E}'.format(cov[0,1]/np.sqrt(cov[0,0]*cov[1,1]))
        
        ax.scatter(x_temp,y_temp,
                   label=temp_label,
                   color=temp_color,
                   alpha=0.75,zorder=1)

    #ax.legend(handles=legend_elements,fontsize=6)
    ax.legend(fontsize=8)

    plt.tight_layout()

    #plt.grid()

    new_fig_file = fig_file.strip('.pdf')
    new_fig_file = new_fig_file + '_' + species_select + '.pdf'

    fig.savefig(new_fig_file,dpi=300)

    return

def plot_shell_diff_vs_diff_reps(fig_file,shell_MSD):

    A_to_micro = 1.0E+4 # A/mu
    ns_to_s = 1.0E+9 # ns/s

    cmap = colormaps.get_cmap('tab10')

    fig_size = [87,87]
    #fig_size = [174,174]

    fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

    ax = plt.gca()

    ax.set_xlabel(r'$D_{{DNA}}$ - Diffusion Constant [$\mu$m$^2$/s]', fontsize=12)
    ax.set_ylabel(r'$D_{{ribo}}$ - Diffusion Constant [$\mu$m$^2$/s]', fontsize=12)
    #ax.set_title(r'MSD of 10~bp DNA monomers',fontsize=10)

    temp_marker_style = dict(marker='.', markersize=6,fillstyle='none')

    density = shell_MSD['DNA_density_reps']
    DNA_counts = shell_MSD['DNA_counts']
    DNA_BrownianDiff = shell_MSD['DNA_BrownianDiff']
    ribo_counts = shell_MSD['ribo_counts']
    ribo_BrownianDiff = shell_MSD['ribo_BrownianDiff']

    ax.tick_params(axis='x',labelsize=9)
    ax.tick_params(axis='y',labelsize=9)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    #temp_marker_style = dict(marker='.', markersize=8,fillstyle='none')
    cov_ellipses = []

    for i_shell in range(shell_MSD['N_shells']):

        DNA_valid_reps = np.argwhere(DNA_counts[i_shell] > 0).flatten()
        ribo_valid_reps = np.argwhere(ribo_counts[i_shell] > 0).flatten()

        valid_reps = np.intersect1d(DNA_valid_reps,ribo_valid_reps)
        
        DNA_beta = DNA_BrownianDiff[i_shell,valid_reps,:]
        ribo_beta = ribo_BrownianDiff[i_shell,valid_reps,:]

        temp_color = cmap(i_shell%shell_MSD['N_shells'])
        #temp_label = r'Shell \#'+str(i_shell+1) + ', '
        temp_label = r'{:d} - {:d} nm: '.format(shell_MSD['shell_lims'][i_shell,0]//10,
                                                                         shell_MSD['shell_lims'][i_shell,1]//10)
        x_temp = (DNA_beta[:,1]*ns_to_s)/(6.0*(A_to_micro**2.0))
        y_temp = (ribo_beta[:,1]*ns_to_s)/(6.0*(A_to_micro**2.0))

        cov = np.cov(x_temp,y_temp,dtype=np.float32)

        w, v = np.linalg.eig(cov)

        theta = np.rad2deg(np.arctan2(v[1,0],v[0,0]))
        
        xm = np.mean(x_temp)
        ym = np.mean(y_temp)

        cov_ellipse = patches.Ellipse((xm,ym),
                                      2*np.sqrt(w[0]),
                                      2*np.sqrt(w[1]),
                                      angle=theta,
                                      facecolor=None,
                                      edgecolor='white',
                                      fill=False,
                                      alpha=0.95,
                                      linewidth=2.25,
                                      zorder=2)

        cov_ellipses.append(cov_ellipse)
        
        cov_ellipse = patches.Ellipse((xm,ym),
                                      2*np.sqrt(w[0]),
                                      2*np.sqrt(w[1]),
                                      angle=theta,
                                      facecolor=None,
                                      edgecolor=temp_color,
                                      fill=False,
                                      alpha=0.95,
                                      linewidth=1.25,
                                      linestyle='-',
                                      zorder=3)

        cov_ellipses.append(cov_ellipse)

        temp_label += r'$r_{{x,y}}=$ {:.1E}'.format(cov[0,1]/np.sqrt(cov[0,0]*cov[1,1]))
        
        ax.scatter(x_temp,y_temp,
                   label=temp_label,
                   color=temp_color,
                   alpha=0.75,zorder=1)

    for i in range(len(cov_ellipses)):
        ax.add_patch(cov_ellipses[i])

    #ax.legend(handles=legend_elements,fontsize=6)
    ax.legend(fontsize=8)

    plt.tight_layout()

    #plt.grid()

    fig.savefig(fig_file,dpi=300)

    return

def plot_shell_law_vs_law_reps(fig_file,shell_MSD):

    A_to_micro = 1.0E+4 # A/mu
    ns_to_s = 1.0E+9 # ns/s

    cmap = colormaps.get_cmap('tab10')

    fig_size = [87,87]
    #fig_size = [174,174]

    fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

    ax = plt.gca()

    ax.set_xlabel(r'$\alpha_{{DNA}}$ - Anomalous Diffusion Law', fontsize=12)
    ax.set_ylabel(r'$\alpha_{{ribo}}$ - Anomalous Diffusion Law', fontsize=12)
    #ax.set_title(r'MSD of 10~bp DNA monomers',fontsize=10)

    temp_marker_style = dict(marker='.', markersize=6,fillstyle='none')

    density = shell_MSD['DNA_density_reps']
    DNA_counts = shell_MSD['DNA_counts']
    DNA_AnomalousDiff = shell_MSD['DNA_AnomalousDiff']
    ribo_counts = shell_MSD['ribo_counts']
    ribo_AnomalousDiff = shell_MSD['ribo_AnomalousDiff']

    ax.tick_params(axis='x',labelsize=9)
    ax.tick_params(axis='y',labelsize=9)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    #temp_marker_style = dict(marker='.', markersize=8,fillstyle='none')
    cov_ellipses = []

    for i_shell in range(shell_MSD['N_shells']):

        DNA_valid_reps = np.argwhere(DNA_counts[i_shell] > 0).flatten()
        ribo_valid_reps = np.argwhere(ribo_counts[i_shell] > 0).flatten()

        valid_reps = np.intersect1d(DNA_valid_reps,ribo_valid_reps)
        
        DNA_beta = DNA_AnomalousDiff[i_shell,valid_reps,:]
        ribo_beta = ribo_AnomalousDiff[i_shell,valid_reps,:]

        temp_color = cmap(i_shell%shell_MSD['N_shells'])
        #temp_label = r'Shell \#'+str(i_shell+1) + ', '
        temp_label = r'{:d} - {:d} nm: '.format(shell_MSD['shell_lims'][i_shell,0]//10,
                                                                         shell_MSD['shell_lims'][i_shell,1]//10)
        x_temp = DNA_beta[:,1]
        y_temp = ribo_beta[:,1]

        cov = np.cov(x_temp,y_temp,dtype=np.float32)

        w, v = np.linalg.eig(cov)

        theta = np.rad2deg(np.arctan2(v[1,0],v[0,0]))
        
        xm = np.mean(x_temp)
        ym = np.mean(y_temp)

        cov_ellipse = patches.Ellipse((xm,ym),
                                      2*np.sqrt(w[0]),
                                      2*np.sqrt(w[1]),
                                      angle=theta,
                                      facecolor=None,
                                      edgecolor='white',
                                      fill=False,
                                      alpha=0.95,
                                      linewidth=2.25,
                                      zorder=2)

        cov_ellipses.append(cov_ellipse)
        
        cov_ellipse = patches.Ellipse((xm,ym),
                                      2*np.sqrt(w[0]),
                                      2*np.sqrt(w[1]),
                                      angle=theta,
                                      facecolor=None,
                                      edgecolor=temp_color,
                                      fill=False,
                                      alpha=0.95,
                                      linewidth=1.25,
                                      linestyle='-',
                                      zorder=3)

        cov_ellipses.append(cov_ellipse)
        

        temp_label += r'$r_{{x,y}}=$ {:.1E}'.format(cov[0,1]/np.sqrt(cov[0,0]*cov[1,1]))
        
        ax.scatter(x_temp,y_temp,
                   label=temp_label,
                   color=temp_color,
                   alpha=0.75,zorder=1)

    for i in range(len(cov_ellipses)):
        ax.add_patch(cov_ellipses[i])

    #ax.legend(handles=legend_elements,fontsize=6)
    ax.legend(fontsize=8)

    plt.tight_layout()

    #plt.grid()

    fig.savefig(fig_file,dpi=300)

    return

def write_shell_diff(shell_diff_file,shell_MSD):

    with open(shell_diff_file,'wb') as f:
        pickle.dump(shell_MSD,f)

    return

def read_shell_diff(shell_diff_file):

    with open(shell_diff_file,'rb') as f:
        shell_MSD = pickle.load(f)

    return shell_MSD
