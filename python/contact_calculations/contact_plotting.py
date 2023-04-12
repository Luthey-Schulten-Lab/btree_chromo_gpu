import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as colormaps
import matplotlib.patches as patches

from mpl_toolkits.axes_grid1 import make_axes_locatable

import os

import imp


plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r'\usepackage[helvet]{sfmath}'

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'

mm = 1/25.4
l_w_ratio = 1/1.618

def SK_norm(mat,tol,N,maxiter):

    e = np.ones((N),dtype=np.float32)

    i = 0

    error = 100.0

    Dinv = np.diag(np.reciprocal(np.matmul(mat,e)))
    x = np.matmul(Dinv,e)

    error = np.linalg.norm(x-e)

    while (i < maxiter) and (error > tol):

        Dinv = np.diag(np.reciprocal(np.matmul(mat,x)))
        x = np.matmul(Dinv,e)
        error = np.linalg.norm(x-e)
        print(str(i)+' '+str(error))
        i += 1

    mat = np.matmul(np.diag(x),np.matmul(mat,np.diag(x)))

    return mat

def read_matrix(in_dir,in_label,norm_flag):

    # read the pairs binary file
    with open(in_dir+'pairs_'+in_label+'.bin','rb') as f:
    
        pairs = np.fromfile(f,dtype=np.int32,count=-1)

    # reshape the pairs array
    pairs = pairs.reshape((2,pairs.shape[0]//2),order='F').T
    N_pairs = pairs.shape[0]

    # read the matrix binary file depending on the requested normalization
    if norm_flag == True:

        val_file = in_dir+'rval_norm_'+in_label+'.bin'

    else:

        val_file = in_dir+'rval_'+in_label+'.bin'

    print(val_file)

    with open(val_file,'rb') as f:

        vals = np.fromfile(f,dtype=np.float32,count=-1)

    print(vals.shape[0])

    # map the matrix back to the symmetric matrix from the pairs
    N = pairs[-1,0]
    print(N)
    pairs = pairs - 1
    mat = np.zeros((N,N),dtype=np.float32)
    mat[pairs[:,0],pairs[:,1]] = vals
    mat += (mat.T - np.diag(np.diagonal(mat)))

    return mat

def read_CGinfo_file(CGinfo_file):

    CGinfo = dict()

    with open(CGinfo_file,'r') as f:

        temp_line = f.readline()
        CGinfo['N_base'] = int(temp_line.split('=')[1].strip())
        temp_line = f.readline()
        CGinfo['N'] = int(temp_line.split('=')[1].strip())
        temp_line = f.readline()
        CGinfo['f_CG'] = int(temp_line.split('=')[1].strip())
        temp_line = f.readline()
        CGinfo['N_base_CG'] = int(temp_line.split('=')[1].strip())
        temp_line = f.readline()
        CGinfo['N_CG'] = int(temp_line.split('=')[1].strip())

        f.readline()
        f.readline()

        temp_line = f.readline()
        CGinfo['N_ter'] = int(temp_line.split('=')[1].strip())

        CGinfo['ter_chromos'] = []
        CGinfo['ter_ranges'] = []

        for i_ter in range(CGinfo['N_ter']):

            temp_line = f.readline()
            temp_line = temp_line.strip()
            temp_line = temp_line.split(',')

            CGinfo['ter_chromos'].append(temp_line[0])
            CGinfo['ter_ranges'].append([int(temp_line[1]),int(temp_line[2])])

        CGinfo['ter_ranges'] = np.array(CGinfo['ter_ranges'],dtype=np.int32)
    
    with open(CGinfo_file,'r') as f:

        CGinfo['map'] = np.loadtxt(f,\
                                   delimiter=',',\
                                   skiprows=10+CGinfo['N_ter'],\
                                   dtype=np.int32,usecols=(0,1)) - 1
        

    return CGinfo

def map_matrix(mat,CGinfo):

    mapped_mat = np.copy(mat[:CGinfo['N_base_CG'],:CGinfo['N_base_CG']])

    for i in range(CGinfo['N_base_CG'],CGinfo['N_CG']):

        mapped_mat[CGinfo['map'][i,1],CGinfo['map'][i,1]] += mat[i,i]
        
        for j in range(0,i):

            mapped_mat[CGinfo['map'][i,1],CGinfo['map'][j,1]] += mat[i,j]
            mapped_mat[CGinfo['map'][j,1],CGinfo['map'][i,1]] += mat[i,j]

    return mapped_mat

def plot_mat(out_dir,out_label,mat,CGinfo,fig_size,norm_flag,thresh,CGinfo_flag,overlay_flag,dm_flag):

    # fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))
    fig = plt.figure(figsize=(1.15*fig_size[0]*mm,fig_size[1]*mm))

    ax = plt.gca()

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)


    # if norm_flag == 'scaled':
    #     map_norm = matplotlib.colors.Normalize(vmin=0.0, vmax=thresh)
    # elif norm_flag == 'log':
    #     map_norm = matplotlib.colors.LogNorm(vmin=np.min(mat), vmax=thresh)
    # elif norm_flag == 'logmax':
    #     map_norm = matplotlib.colors.LogNorm(vmin=np.min(mat), vmax=np.max(mat))
    # else:
    #     map_norm = matplotlib.colors.Normalize(vmin=0.0, vmax=np.max(mat))

    if norm_flag == 'scaled':
        map_norm = matplotlib.colors.Normalize(vmin=0.0, vmax=thresh[1])
    elif norm_flag == 'log':
        if np.min(mat) <= 0.0:
            temp_min = mat.flatten()
            temp_min.sort()
            temp_min = temp_min[np.argwhere(temp_min>0.0)[0]]
            mat = np.where(mat==0.0,temp_min,mat)
        else:
            temp_min = np.min(mat)
        temp_max = np.max(mat)
        print('max = '+str(temp_max))
        map_norm = matplotlib.colors.LogNorm(vmin=temp_min, vmax=thresh)
    elif norm_flag == 'logmax':
        if np.min(mat) <= 0.0:
            temp_min = mat.flatten()
            temp_min.sort()
            temp_min = temp_min[np.argwhere(temp_min>0.0)[0]]
            mat = np.where(mat==0.0,temp_min,mat)
        else:
            temp_min = np.min(mat)
        temp_min = np.percentile(mat,10)
        temp_max = np.max(mat)
        print('min = '+str(temp_min))
        print('max = '+str(temp_max))
        map_norm = matplotlib.colors.LogNorm(vmin=temp_min, vmax=temp_max)
    elif norm_flag == 'manual':
        map_norm = matplotlib.colors.LogNorm(vmin=thresh[0], vmax=thresh[1])
    else:
        map_norm = matplotlib.colors.Normalize(vmin=0.0, vmax=np.max(mat))

    im = ax.imshow(mat,cmap='binary',norm=map_norm)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)

    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label(label=r'Contact Frequency', fontsize=7)
    cbar.ax.tick_params(labelsize=6)


    if dm_flag == True:
        dm = np.mean(np.diagonal(mat))
        print(dm)
        arrow = dict()
        arrow['width'] = 2
        arrow['headwidth'] = 5
        arrow['headlength'] = 7
        arrow['shrink'] = 0.0
        arrow['ec'] = 'black'
        arrow['fc'] = 'white'
        
        cax.annotate(r'$\frac{\sum_{i}\mathsf{A}_{ii}}{N}=$'+' {:.3f}'.format(dm),
                     xy=(0.0,dm),
                     xytext=(-0.5,dm),
                     xycoords='data',
                     arrowprops=arrow,
                     ha='right',
                     va='center',
                     fontsize=8)


    s_expand = 0.015*float(CGinfo['N_CG'])

    ax.set_xlim(xmin=0.0-s_expand,xmax=float(CGinfo['N_CG'])+s_expand)
    ax.set_ylim(ymin=float(CGinfo['N_CG'])+s_expand,ymax=1.0-s_expand)


    if CGinfo_flag == True:
        ticks = np.arange(0.0,float(CGinfo['N_CG'])/CGinfo['N_base_CG']+0.02,0.5,dtype=np.float32)
        tick_labels = []

        tick_labels.append(r'{:.1f} (ter)'.format(ticks[0]))
        tick_labels.append(r'{:.1f} (ori)'.format(ticks[1]))
        tick_labels.append(r'{:.1f} (ter)'.format(ticks[2]))
        
        if ticks.shape[0] > 3:
            for i in range(3,ticks.shape[0]):
                tick_labels.append(r'{:.1f}'.format(ticks[i]))
        ticks = ticks*CGinfo['N_base_CG']

    if overlay_flag == True:

        overlay_color = 'red'

        for i_ter in range(CGinfo['N_ter']):
        # Create a Rectangle patch
            xy_anchor = (CGinfo['ter_ranges'][i_ter][0]-1,CGinfo['ter_ranges'][i_ter][0]-1)
            w = CGinfo['ter_ranges'][i_ter][1] - CGinfo['ter_ranges'][i_ter][0] + 1
            h = CGinfo['ter_ranges'][i_ter][1] - CGinfo['ter_ranges'][i_ter][0] + 1

            print(xy_anchor)
            print(w)
            print(h)
            rect = patches.Rectangle(xy=xy_anchor,width=w,height=h,linewidth=1.75,edgecolor='white',facecolor='none')
            ax.add_patch(rect)
            rect = patches.Rectangle(xy=xy_anchor,width=w,height=h,linewidth=1.0,edgecolor=overlay_color,facecolor='none')
            ax.add_patch(rect)

            sl = 0.85
            su = 1.0
            x = np.array([float(xy_anchor[0])+sl*float(w)/2.0,float(xy_anchor[0])+su*float(w)/2.0],dtype=np.float32)
            y = np.array([float(xy_anchor[1])+(2.0-sl)*float(h)/2.0,float(xy_anchor[1])+(2.0-su)*float(h)/2.0],dtype=np.float32)
            #ax.plot(x,y,linewidth=0.25,color=overlay_color)

            # x = np.array([float(xy_anchor[0])+float(w)/2.0,float(xy_anchor[0])+float(w)/2.0],dtype=np.float32)
            # y = np.array([float(xy_anchor[1])+0.7*float(h)/2.0,float(xy_anchor[1])+0.85*float(h)/2.0],dtype=np.float32)
            # ax.plot(x,y,linewidth=0.25,color=overlay_color)

            # if i_ter == 0:
            #     xy_anchor = (CGinfo['ter_ranges'][i_ter][0],CGinfo['ter_ranges'][i_ter][1])
            #     ax.annotate(text=r''+CGinfo['ter_chromos'][i_ter],xy=xy_anchor,va='bottom',ha='left',fontsize=6,color=overlay_color)
            # else:
            #     xy_anchor = (CGinfo['ter_ranges'][i_ter][0],CGinfo['ter_ranges'][i_ter][1])
            #     ax.annotate(text=r''+CGinfo['ter_chromos'][i_ter],xy=xy_anchor,va='bottom',ha='right',fontsize=6,color=overlay_color)

            xy_anchor = (CGinfo['ter_ranges'][i_ter][0],CGinfo['ter_ranges'][i_ter][1])
            ax.annotate(text=r''+CGinfo['ter_chromos'][i_ter],xy=xy_anchor,va='bottom',ha='left',fontsize=6,color=overlay_color)

    print(tick_labels)
    print(ticks)
    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels,fontsize=6,ha='center')
    ax.set_yticks(ticks)
    ax.set_yticklabels(tick_labels,fontsize=6,rotation=90,va='center')

    ax.set_xlabel(r'Total DNA content',fontsize=8)
    ax.set_ylabel(r'Total DNA content',fontsize=8)

    # plt.tight_layout()

    #fig.savefig(out_dir+out_label+'_contacts.pdf',dpi=300)
    fig.savefig(out_dir+out_label+'_contacts.png',dpi=300)

    return

def plot_mat_mapped(out_dir,out_label,mat,CGinfo,fig_size,norm_flag,thresh,CGinfo_flag,overlay_flag):

    fig = plt.figure(figsize=(1.15*fig_size[0]*mm,fig_size[1]*mm))

    ax = plt.gca()

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    # if norm_flag == 'scaled':
    #     map_norm = matplotlib.colors.Normalize(vmin=0.0, vmax=thresh)
    # elif norm_flag == 'log':
    #     map_norm = matplotlib.colors.LogNorm(vmin=np.min(mat), vmax=thresh)
    # elif norm_flag == 'logmax':
    #     map_norm = matplotlib.colors.LogNorm(vmin=np.min(mat), vmax=np.max(mat))
    # else:
    #     map_norm = matplotlib.colors.Normalize(vmin=0.0, vmax=np.max(mat))

    if norm_flag == 'scaled':
        map_norm = matplotlib.colors.Normalize(vmin=0.0, vmax=thresh[1])
    elif norm_flag == 'log':
        if np.min(mat) <= 0.0:
            temp_min = mat.flatten()
            temp_min.sort()
            temp_min = temp_min[np.argwhere(temp_min>0.0)[0]]
            mat = np.where(mat==0.0,temp_min,mat)
        else:
            temp_min = np.min(mat)
        temp_max = np.max(mat)
        print('max = '+str(temp_max))
        map_norm = matplotlib.colors.LogNorm(vmin=temp_min, vmax=thresh)
    elif norm_flag == 'logmax':
        if np.min(mat) <= 0.0:
            temp_min = mat.flatten()
            temp_min.sort()
            temp_min = temp_min[np.argwhere(temp_min>0.0)[0]]
            mat = np.where(mat==0.0,temp_min,mat)
        else:
            temp_min = np.min(mat)
        temp_min = np.percentile(mat,10)
        temp_max = np.max(mat)
        print('min = '+str(temp_min))
        print('max = '+str(temp_max))
        map_norm = matplotlib.colors.LogNorm(vmin=temp_min, vmax=temp_max)
    elif norm_flag == 'manual':
        map_norm = matplotlib.colors.LogNorm(vmin=thresh[0], vmax=thresh[1])
    else:
        map_norm = matplotlib.colors.Normalize(vmin=0.0, vmax=np.max(mat))

    im = ax.imshow(mat,cmap='binary',norm=map_norm)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)

    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label(label=r'Contact Frequency', fontsize=7)
    cbar.ax.tick_params(labelsize=6)

    s_expand = 0.015*float(CGinfo['N_base_CG'])

    ax.set_xlim(xmin=0.0-s_expand,xmax=float(CGinfo['N_base_CG'])+s_expand)
    ax.set_ylim(ymin=float(CGinfo['N_base_CG'])+s_expand,ymax=1.0-s_expand)

    
    if CGinfo_flag == True:
        ticks = np.arange(0.0,1.0+0.02,0.5,dtype=np.float32)
        tick_labels = []
        tick_labels.append(r'{:.1f} (ter)'.format(ticks[0]))
        tick_labels.append(r'{:.1f} (ori)'.format(ticks[1]))
        tick_labels.append(r'{:.1f} (ter)'.format(ticks[2]))
        ticks = ticks*CGinfo['N_base_CG']

    if overlay_flag == True:

        overlay_color = 'red'
        s = 0.95

        for i_ter in range(CGinfo['N_ter']):
            
            if i_ter > 0:

                xmap = float(CGinfo['map'][CGinfo['ter_ranges'][i_ter][0]-1,1])
                ymap = float(CGinfo['map'][CGinfo['ter_ranges'][i_ter][1]-1,1])
                w = float(CGinfo['ter_ranges'][i_ter][1] - CGinfo['ter_ranges'][i_ter][0] + 1)
                h = float(CGinfo['ter_ranges'][i_ter][1] - CGinfo['ter_ranges'][i_ter][0] + 1)

                x = np.array([xmap,xmap,xmap+s*w],dtype=np.float32)
                y = np.array([ymap-s*h,ymap,ymap],dtype=np.float32)

                ax.plot(x,y,linewidth=1.5,color='white',alpha=1.0)
                ax.plot(x,y,linewidth=1.0,color=overlay_color,alpha=1.0)

    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels,fontsize=6,ha='left')
    ax.set_yticks(ticks)
    ax.set_yticklabels(tick_labels,fontsize=6,rotation=90,va='center')

    ax.set_xlabel(r'Mapped DNA content',fontsize=8)
    ax.set_ylabel(r'Mapped DNA content',fontsize=8)

    #plt.tight_layout()

    fig.savefig(out_dir+out_label+'_contacts_mapped.pdf',dpi=300)

    return

def plot_mat_abs(out_dir,out_label,mat,CGinfo,fig_size,norm_flag,thresh,CGinfo_flag,overlay_flag):

    fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

    ax = plt.gca()

    if norm_flag == 'scaled':
        map_norm = matplotlib.colors.Normalize(vmin=0.0, vmax=thresh)
    elif norm_flag == 'log':
        if np.min(mat) <= 0.0:
            temp_min = mat.flatten()
            temp_min.sort()
            temp_min = temp_min[np.argwhere(temp_min>0.0)[0]]
            mat = np.where(mat==0.0,temp_min,mat)
        else:
            temp_min = np.min(mat)
        temp_max = np.max(mat)
        print('max = '+str(temp_max))
        map_norm = matplotlib.colors.LogNorm(vmin=temp_min, vmax=thresh)
    elif norm_flag == 'logmax':
        if np.min(mat) <= 0.0:
            temp_min = mat.flatten()
            temp_min.sort()
            temp_min = temp_min[np.argwhere(temp_min>0.0)[0]]
            mat = np.where(mat==0.0,temp_min,mat)
        else:
            temp_min = np.min(mat)
        temp_min = np.percentile(mat,10)
        temp_max = np.max(mat)
        print('min = '+str(temp_min))
        print('max = '+str(temp_max))
        map_norm = matplotlib.colors.LogNorm(vmin=temp_min, vmax=temp_max)
    else:
        map_norm = matplotlib.colors.Normalize(vmin=0.0, vmax=np.max(mat))

    im = ax.imshow(mat,cmap='plasma',norm=map_norm)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)

    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label(label=r'Contact Frequency', fontsize=6)
    cbar.ax.tick_params(labelsize=4)


    if CGinfo_flag == True:
        #ticks = np.arange(0.0,float(CGinfo['N_CG'])/CGinfo['N_base_CG']+0.05,0.5,dtype=np.float32)
        half_tick = (CGinfo['N_CG']//2)*10*CGinfo['f_CG']//1000
        print(half_tick)
        ticks = np.arange(-half_tick,half_tick+1,half_tick,dtype=np.int32)
        tick_labels = []

        tick_labels.append(r'{:3d} (ter)'.format(ticks[0]))
        tick_labels.append(r'{:3d} (ori)'.format(ticks[1]))
        tick_labels.append(r'{:3d} (ter)'.format(ticks[2]))
        
        if ticks.shape[0] > 3:
            for i in range(ticks.shape[0]):
                tick_labels.append(r'{:.1f}'.format(ticks[i]))

    if overlay_flag == True:

        overlay_color = 'limegreen'

        for i_ter in range(CGinfo['N_ter']):
        # Create a Rectangle patch
            xy_anchor = (CGinfo['ter_ranges'][i_ter][0]-1,CGinfo['ter_ranges'][i_ter][0]-1)
            w = CGinfo['ter_ranges'][i_ter][1] - CGinfo['ter_ranges'][i_ter][0] + 1
            h = CGinfo['ter_ranges'][i_ter][1] - CGinfo['ter_ranges'][i_ter][0] + 1

            print(xy_anchor)
            print(w)
            print(h)
            rect = patches.Rectangle(xy=xy_anchor,width=w,height=h,linewidth=0.25,edgecolor=overlay_color,facecolor='none')
            ax.add_patch(rect)

            sl = 0.85
            su = 1.0
            x = np.array([float(xy_anchor[0])+sl*float(w)/2.0,float(xy_anchor[0])+su*float(w)/2.0],dtype=np.float32)
            y = np.array([float(xy_anchor[1])+(2.0-sl)*float(h)/2.0,float(xy_anchor[1])+(2.0-su)*float(h)/2.0],dtype=np.float32)
            #ax.plot(x,y,linewidth=0.25,color=overlay_color)

            # x = np.array([float(xy_anchor[0])+float(w)/2.0,float(xy_anchor[0])+float(w)/2.0],dtype=np.float32)
            # y = np.array([float(xy_anchor[1])+0.7*float(h)/2.0,float(xy_anchor[1])+0.85*float(h)/2.0],dtype=np.float32)
            # ax.plot(x,y,linewidth=0.25,color=overlay_color)

            if i_ter == 0:
                xy_anchor = (CGinfo['ter_ranges'][i_ter][0],CGinfo['ter_ranges'][i_ter][1])
                ax.annotate(text=r''+CGinfo['ter_chromos'][i_ter],xy=xy_anchor,va='bottom',ha='left',fontsize=6,color=overlay_color)
            else:
                xy_anchor = (CGinfo['ter_ranges'][i_ter][0],CGinfo['ter_ranges'][i_ter][1])
                ax.annotate(text=r''+CGinfo['ter_chromos'][i_ter],xy=xy_anchor,va='bottom',ha='right',fontsize=6,color=overlay_color)


    print(ticks)
    ticks = (ticks+half_tick)
    print(ticks)
    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels,fontsize=7,ha='left')
    ax.set_yticks(ticks)
    ax.set_yticklabels(tick_labels,fontsize=7,rotation=90,va='center')

    ax.set_xlabel(r'Genome Position [kbp]',fontsize=8)
    ax.set_ylabel(r'Genome Position [kbp]',fontsize=8)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    plt.tight_layout()

    fig.savefig(out_dir+out_label+'_contacts.pdf',dpi=300)
    fig.savefig(out_dir+out_label+'_contacts.png',dpi=300)

    return

def contact_mat_alignment(N,mat):
    
    aligned_mat = np.zeros((N,N),dtype=mat.dtype)
    
    aligned_mat[0,:] = mat[0,:]
    
    for i in range(1,N):
        aligned_mat[i,:] = np.roll(mat[i,:],-i,axis=0)
    
    return aligned_mat

def multi_fits(ax,locus_size,regimes,contact_sum):
    
    for regime in reversed(regimes):

        X=np.ones((regime[1]-regime[0],2),dtype=np.double,order='F')

        z=(np.arange(regime[0],regime[1])*(locus_size/1000.0))
        if regime[0] == 0:
            z += 1

        X[:,1]=np.log(z)

        y=np.log(contact_sum[regime[0]:regime[1]])
        y_old = np.copy(y)

        y=np.matmul(X.T,y)

        X=np.matmul(X.T,X)

        beta=np.matmul(np.linalg.inv(X),y)

        print(beta[1])

        y=beta[0]+beta[1]*np.log(z)

        y=np.exp(y)

        SS_tot = np.sum(np.power(y_old-np.mean(y_old),2.0))
        SS_res = np.sum(np.power(np.log(y)-y_old,2.0))

        R_sqrd = 1.0 - SS_res/SS_tot

        if regime[0] == 0:
            z[0]=0
            
        z=np.rint(z).astype(np.int64)
            
        temp_label=r'{:d}-{:d} kbp: '.format(z[0],z[-1])+'$s=$ '+"{:.3f}".format(-beta[1]) + \
            ', $R^2=$ {:.2f}'.format(R_sqrd)
            
        
        ax.plot(z,y,ls="-",c='w',lw=2.0,alpha=0.65)
        ax.plot(z,y,ls="--",lw=1.5,alpha=0.9,label=temp_label)
    
    return

def plot_contact_law(out_dir,out_label,mat,CGinfo,fig_size):

    aligned_mat = contact_mat_alignment(CGinfo['N_base_CG'],mat)

    contact_sum = np.mean(aligned_mat,axis=0)

    #print(contact_sum)
    print(np.sum(contact_sum))

    fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

    ax = plt.gca()

    ax.set_xlabel(r'Genomic Distance [kbp] - $d$', fontsize=8)
    ax.set_ylabel(r'Contact Probability - $P(d)$', fontsize=7)

    temp_marker_style = dict(marker='.', markersize=8,fillstyle='none')

    lower_lim = 1
    upper_lim = 151
    locus_size = CGinfo['f_CG']*10

    za=(np.arange(lower_lim,upper_lim)*(locus_size/1000.0))

    ax.plot(za,contact_sum[lower_lim:upper_lim],c='k',lw=1.0,alpha=0.35,**temp_marker_style)

    regimes = [[1,11],[11,51],[51,101]]

    multi_fits(ax,locus_size,regimes,contact_sum)

    ax.set_xlim(left=lower_lim*(locus_size/1000.0),right=upper_lim*(locus_size/1000.0))

    ax.set_yscale('log')
    ax.set_xscale('log')

    ax.legend(fontsize=7,loc='upper right')

    ax.tick_params(axis='x',labelsize=7)
    ax.tick_params(axis='y',labelsize=7)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    plt.tight_layout()

    #plt.grid()

    fig.savefig(out_dir+out_label+'_contact_law.pdf',dpi=300)
    fig.savefig(out_dir+out_label+'_contact_law.png',dpi=300)

    return

def plot_mapping(out_dir,out_label,CGinfo,fig_size,CGinfo_flag,overlay_flag):

    fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

    ax = plt.gca()

    N_plots = CGinfo['N_ter']*(CGinfo['N_ter']-1)//2 + CGinfo['N_ter']-1

    cmap_color_select = colormaps.get_cmap('gist_rainbow')
    color_selections = []

    print('N_plots='+str(N_plots))

    for i_plot in range(N_plots):

        color_selections.append(cmap_color_select(float(i_plot)/N_plots))


    i_plot = 0

    temp_ter_ranges = np.array(CGinfo['ter_ranges'],dtype=np.int32)
    ter_sorted = np.flip(np.argsort(temp_ter_ranges[:,1]-temp_ter_ranges[:,0]))
    temp_ter_ranges = temp_ter_ranges[ter_sorted,:]
    print(temp_ter_ranges)
        
    for i_ter in range(CGinfo['N_ter']):

        if i_ter > 0:
            k_ter = i_ter
        else:
            k_ter = 1
        
        for j_ter in range(k_ter,CGinfo['N_ter']):

            new_img = np.zeros((CGinfo['N_CG'],CGinfo['N_CG']),dtype=np.int32)

            for i in range(temp_ter_ranges[i_ter,0]-1,temp_ter_ranges[i_ter,1]-1):

                for j in range(temp_ter_ranges[j_ter,0]-1,temp_ter_ranges[j_ter,1]-1):

                    new_img[i,j] = 1
                    new_img[CGinfo['map'][i,1],CGinfo['map'][j,1]] = 1
                    new_img[j,i] = 1
                    new_img[CGinfo['map'][j,1],CGinfo['map'][i,1]] = 1

            cmap = matplotlib.colors.ListedColormap(['white',color_selections[i_plot]])
            norm = matplotlib.colors.BoundaryNorm([0,1],cmap.N)

            ax.imshow(new_img,cmap=cmap,norm=norm,alpha=0.75*new_img.astype(dtype=np.float32),zorder=i_plot)

            i_plot += 1


    if CGinfo_flag == True:
        ticks = np.arange(0.0,float(CGinfo['N_CG'])/CGinfo['N_base_CG']+0.02,0.5,dtype=np.float32)
        tick_labels = []

        tick_labels.append(r'{:.1f} (ter)'.format(ticks[0]))
        tick_labels.append(r'{:.1f} (ori)'.format(ticks[1]))
        tick_labels.append(r'{:.1f} (ter)'.format(ticks[2]))
        
        if ticks.shape[0] > 3:
            for i in range(3,ticks.shape[0]):
                tick_labels.append(r'{:.1f}'.format(ticks[i]))
        ticks = ticks*CGinfo['N_base_CG']

    if overlay_flag == True:

        overlay_color = 'black'

        for i_ter in range(CGinfo['N_ter']):
        # Create a Rectangle patch
            xy_anchor = (CGinfo['ter_ranges'][i_ter][0]-1,CGinfo['ter_ranges'][i_ter][0]-1)
            w = CGinfo['ter_ranges'][i_ter][1] - CGinfo['ter_ranges'][i_ter][0] + 1
            h = CGinfo['ter_ranges'][i_ter][1] - CGinfo['ter_ranges'][i_ter][0] + 1

            # print(xy_anchor)
            # print(w)
            # print(h)
            rect = patches.Rectangle(xy=xy_anchor,width=w,height=h,linewidth=1.0,edgecolor='white',facecolor='none',zorder=N_plots)
            ax.add_patch(rect)
            rect = patches.Rectangle(xy=xy_anchor,width=w,height=h,linewidth=0.85,edgecolor=overlay_color,facecolor='none',zorder=N_plots+1)
            ax.add_patch(rect)

            sl = 0.85
            su = 1.0
            x = np.array([float(xy_anchor[0])+sl*float(w)/2.0,float(xy_anchor[0])+su*float(w)/2.0],dtype=np.float32)
            y = np.array([float(xy_anchor[1])+(2.0-sl)*float(h)/2.0,float(xy_anchor[1])+(2.0-su)*float(h)/2.0],dtype=np.float32)
            #ax.plot(x,y,linewidth=0.25,color=overlay_color)

            # x = np.array([float(xy_anchor[0])+float(w)/2.0,float(xy_anchor[0])+float(w)/2.0],dtype=np.float32)
            # y = np.array([float(xy_anchor[1])+0.7*float(h)/2.0,float(xy_anchor[1])+0.85*float(h)/2.0],dtype=np.float32)
            # ax.plot(x,y,linewidth=0.25,color=overlay_color)

            if i_ter == 0:
                xy_anchor = (CGinfo['ter_ranges'][i_ter][0],CGinfo['ter_ranges'][i_ter][1])
                ax.annotate(text=r''+CGinfo['ter_chromos'][i_ter],xy=xy_anchor,va='bottom',ha='left',fontsize=8,color=overlay_color)
            else:
                xy_anchor = (CGinfo['ter_ranges'][i_ter][0],CGinfo['ter_ranges'][i_ter][1])
                ax.annotate(text=r''+CGinfo['ter_chromos'][i_ter],xy=xy_anchor,va='bottom',ha='right',fontsize=8,color=overlay_color)

    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels,fontsize=8,ha='left')
    ax.set_yticks(ticks)
    ax.set_yticklabels(tick_labels,fontsize=8,rotation=90,va='center')

    ax.set_xlabel(r'Relative DNA content',fontsize=10)
    ax.set_ylabel(r'Relative DNA content',fontsize=10)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    plt.tight_layout()

    fig.savefig(out_dir+out_label+'_mapping.pdf',dpi=300)

    return

def plot_mapping_simple(out_dir,out_label,CGinfo,cb_length,fig_size,CGinfo_flag,overlay_flag):

    fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

    ax = plt.gca()

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    s_expand = 0.015*float(CGinfo['N_CG'])

    ax.set_xlim(xmin=0.0-s_expand,xmax=float(CGinfo['N_CG'])+s_expand)
    ax.set_ylim(ymin=float(CGinfo['N_CG'])+s_expand,ymax=1.0-s_expand)

    map_norm = matplotlib.colors.Normalize(vmin=0.0, vmax=1.0)

    img = np.zeros((CGinfo['N_CG'],CGinfo['N_CG']),dtype=np.single)

    ter_colors = np.zeros((CGinfo['N_ter']),dtype=np.single)

    for i_ter in range(CGinfo['N_ter']):

        CGinfo['ter_ranges'][i_ter] -= 1
        ter_colors[i_ter] = np.mean(CGinfo['ter_ranges'][i_ter])/CGinfo['N_CG']
        print(CGinfo['ter_ranges'][i_ter])

    # fill squares
    for i_ter in range(CGinfo['N_ter']):

        img[CGinfo['ter_ranges'][i_ter][0]:CGinfo['ter_ranges'][i_ter][1]+1,
            CGinfo['ter_ranges'][i_ter][0]:CGinfo['ter_ranges'][i_ter][1]+1] = ter_colors[i_ter]

    # fill checkerboards
    for i_ter in range(CGinfo['N_ter']-1):

        row_ter_range = CGinfo['ter_ranges'][i_ter]
        
        N_rows = CGinfo['ter_ranges'][i_ter][1] - CGinfo['ter_ranges'][i_ter][0] + 1
        N_row_cbs = N_rows//cb_length
        if N_row_cbs*cb_length < N_rows:
            N_row_cbs += 1

        row_cb_ranges = np.zeros((N_row_cbs,2),dtype=np.int32)
        row_cb_ranges[0,0] = row_ter_range[0]

        for i in range(N_row_cbs):

            if i != 0:
                
                row_cb_ranges[i,0] = row_cb_ranges[i-1,1]

            row_cb_ranges[i,1] = min(row_cb_ranges[i,0]+cb_length,row_cb_ranges[0,0]+N_rows+1)
        
        for j_ter in range(i_ter+1,CGinfo['N_ter']):

            col_ter_range = CGinfo['ter_ranges'][j_ter]

            N_cols = CGinfo['ter_ranges'][i_ter][1] - CGinfo['ter_ranges'][i_ter][0] + 1
            N_col_cbs = N_cols//cb_length
            if N_col_cbs*cb_length < N_cols:
                N_col_cbs += 1

            col_cb_ranges = np.zeros((N_col_cbs,2),dtype=np.int32)
            col_cb_ranges[0,0] = col_ter_range[0]

            for i in range(N_col_cbs):

                if i != 0:

                    col_cb_ranges[i,0] = col_cb_ranges[i-1,1]

                col_cb_ranges[i,1] = min(col_cb_ranges[i,0]+cb_length,col_cb_ranges[0,0]+N_cols+1)

            
            # select colors for checkerboard
            if j_ter%2 == 0:
                c1 = ter_colors[i_ter]
                c2 = ter_colors[j_ter]
            else:
                c1 = ter_colors[j_ter]
                c2 = ter_colors[i_ter]

            cb_counter = 0
            for i in range(row_cb_ranges.shape[0]):
                for j in range(col_cb_ranges.shape[0]):

                    if cb_counter%2 == 0:
                        cb_color = c1
                    else:
                        cb_color = c2
                        
                    
                    img[row_cb_ranges[i,0]:row_cb_ranges[i,1],
                        col_cb_ranges[j,0]:col_cb_ranges[j,1]] = cb_color

                    img[col_cb_ranges[j,0]:col_cb_ranges[j,1],
                        row_cb_ranges[i,0]:row_cb_ranges[i,1]] = cb_color

                    cb_counter += 1

            

    im = ax.imshow(img,cmap='viridis',norm=map_norm)
    

    if CGinfo_flag == True:
        ticks = np.arange(0.0,float(CGinfo['N_CG'])/CGinfo['N_base_CG']+0.02,0.5,dtype=np.float32)
        tick_labels = []

        tick_labels.append(r'{:.1f} (ter)'.format(ticks[0]))
        tick_labels.append(r'{:.1f} (ori)'.format(ticks[1]))
        tick_labels.append(r'{:.1f} (ter)'.format(ticks[2]))
        
        if ticks.shape[0] > 3:
            for i in range(3,ticks.shape[0]):
                tick_labels.append(r'{:.1f}'.format(ticks[i]))
        ticks = ticks*CGinfo['N_base_CG']

    if overlay_flag == True:

        overlay_color = 'red'

        for i_ter in range(CGinfo['N_ter']):
        # Create a Rectangle patch
            xy_anchor = (CGinfo['ter_ranges'][i_ter][0],CGinfo['ter_ranges'][i_ter][0])
            w = CGinfo['ter_ranges'][i_ter][1] - CGinfo['ter_ranges'][i_ter][0] + 1
            h = CGinfo['ter_ranges'][i_ter][1] - CGinfo['ter_ranges'][i_ter][0] + 1

            print(xy_anchor)
            print(w)
            print(h)
            rect = patches.Rectangle(xy=xy_anchor,width=w,height=h,linewidth=1.75,edgecolor='white',facecolor='none')
            ax.add_patch(rect)
            rect = patches.Rectangle(xy=xy_anchor,width=w,height=h,linewidth=1.0,edgecolor=overlay_color,facecolor='none')
            ax.add_patch(rect)

            sl = 0.85
            su = 1.0
            x = np.array([float(xy_anchor[0])+sl*float(w)/2.0,float(xy_anchor[0])+su*float(w)/2.0],dtype=np.float32)
            y = np.array([float(xy_anchor[1])+(2.0-sl)*float(h)/2.0,float(xy_anchor[1])+(2.0-su)*float(h)/2.0],dtype=np.float32)

            xy_anchor = (CGinfo['ter_ranges'][i_ter][0],CGinfo['ter_ranges'][i_ter][1])
            ax.annotate(text=r''+CGinfo['ter_chromos'][i_ter],xy=xy_anchor,va='bottom',ha='left',fontsize=6,color=overlay_color)

    print(tick_labels)
    print(ticks)
    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels,fontsize=6,ha='left')
    ax.set_yticks(ticks)
    ax.set_yticklabels(tick_labels,fontsize=6,rotation=90,va='center')

    ax.set_xlabel(r'Total DNA content',fontsize=8)
    ax.set_ylabel(r'Total DNA content',fontsize=8)

    #plt.tight_layout()

    fig.savefig(out_dir+out_label+'_mapping.pdf',dpi=300)

    plt.close()

    

    fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

    ax = plt.gca()

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    s_expand = 0.015*float(CGinfo['N_base_CG'])

    ax.set_xlim(xmin=0.0-s_expand,xmax=float(CGinfo['N_base_CG'])+s_expand)
    ax.set_ylim(ymin=float(CGinfo['N_base_CG'])+s_expand,ymax=1.0-s_expand)

    img_merged = np.zeros((CGinfo['N_base_CG'],CGinfo['N_base_CG']),dtype=np.single)

    print(CGinfo['map'])

    for i_ter in range(CGinfo['N_ter']):

        temp_ter_range = CGinfo['ter_ranges'][i_ter]

        for i in range(temp_ter_range[0],temp_ter_range[1]+1):

            i_m = CGinfo['map'][i,1]

            for j in range(temp_ter_range[0],CGinfo['N_CG']):

                j_m = CGinfo['map'][j,1]

                img_merged[i_m,j_m] = img[i,j]
                # img_merged[i_m,j_m] = img[j,i]

                # if j > temp_ter_range[1]:
                #      img_merged[j_m,i_m] = img[i,j]
                #      img_merged[j_m,i_m] = img[j,i]

    img_merged = np.maximum(img_merged,img_merged.T)

    im_merged = ax.imshow(img_merged,cmap='viridis',norm=map_norm)

    if CGinfo_flag == True:
        ticks = np.arange(0.0,1.0+0.02,0.5,dtype=np.float32)
        tick_labels = []
        tick_labels.append(r'{:.1f} (ter)'.format(ticks[0]))
        tick_labels.append(r'{:.1f} (ori)'.format(ticks[1]))
        tick_labels.append(r'{:.1f} (ter)'.format(ticks[2]))
        ticks = ticks*CGinfo['N_base_CG']

    if overlay_flag == True:

        overlay_color = 'red'

        for i_ter in range(CGinfo['N_ter']):
        # Create a Rectangle patch
            x_mapped = CGinfo['map'][CGinfo['ter_ranges'][i_ter][0]][1]
            y_mapped = CGinfo['map'][CGinfo['ter_ranges'][i_ter][0]][1]
            xy_anchor = (x_mapped,y_mapped)
            w = CGinfo['ter_ranges'][i_ter][1] - CGinfo['ter_ranges'][i_ter][0] + 1
            h = CGinfo['ter_ranges'][i_ter][1] - CGinfo['ter_ranges'][i_ter][0] + 1

            print(xy_anchor)
            print(w)
            print(h)
            rect = patches.Rectangle(xy=xy_anchor,width=w,height=h,linewidth=1.75,edgecolor='white',facecolor='none')
            ax.add_patch(rect)
            rect = patches.Rectangle(xy=xy_anchor,width=w,height=h,linewidth=1.0,edgecolor=overlay_color,facecolor='none')
            ax.add_patch(rect)

            sl = 0.85
            su = 1.0
            x = np.array([float(xy_anchor[0])+sl*float(w)/2.0,float(xy_anchor[0])+su*float(w)/2.0],dtype=np.float32)
            y = np.array([float(xy_anchor[1])+(2.0-sl)*float(h)/2.0,float(xy_anchor[1])+(2.0-su)*float(h)/2.0],dtype=np.float32)

            x_mapped = CGinfo['map'][CGinfo['ter_ranges'][i_ter][0]][1]
            y_mapped = CGinfo['map'][CGinfo['ter_ranges'][i_ter][1]][1]
            xy_anchor = (x_mapped,y_mapped)
            ax.annotate(text=r''+CGinfo['ter_chromos'][i_ter],xy=xy_anchor,va='bottom',ha='left',fontsize=6,color=overlay_color)

    print(tick_labels)
    print(ticks)
    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels,fontsize=6,ha='left')
    ax.set_yticks(ticks)
    ax.set_yticklabels(tick_labels,fontsize=6,rotation=90,va='center')

    ax.set_xlabel(r'Mapped DNA content',fontsize=8)
    ax.set_ylabel(r'Mapped DNA content',fontsize=8)

    #plt.tight_layout()

    fig.savefig(out_dir+out_label+'_merged.pdf',dpi=300)

    plt.close()

    return

def plot_mat_mapped_dev(out_dir,out_label,mat,CGinfo,fig_size,CGinfo_flag,overlay_flag):

    fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

    ax = plt.gca()

    # if norm_flag == 'scaled':
    #     map_norm = matplotlib.colors.Normalize(vmin=0.0, vmax=thresh)
    # elif norm_flag == 'log':
    #     map_norm = matplotlib.colors.LogNorm(vmin=np.min(mat), vmax=thresh)
    # elif norm_flag == 'logmax':
    #     map_norm = matplotlib.colors.LogNorm(vmin=np.min(mat), vmax=np.max(mat))
    # else:
    #     map_norm = matplotlib.colors.Normalize(vmin=0.0, vmax=np.max(mat))

    aligned_mat = contact_mat_alignment(CGinfo['N_base_CG'],mat)
    contact_mean = np.mean(aligned_mat,axis=0)
    print(np.sum(contact_mean))

    contact_stddev = np.std(aligned_mat,axis=0)
    stddev_reciprocal = np.reciprocal(contact_stddev)

    dev_mat = np.zeros((CGinfo['N_base_CG'],CGinfo['N_base_CG']),dtype=np.float32)
    stddev_recip_mat = np.zeros((CGinfo['N_base_CG'],CGinfo['N_base_CG']),dtype=np.float32)
    for i in range(CGinfo['N_base_CG']):
        dev_mat[i,:] = np.roll(contact_mean,i,axis=0)
        stddev_recip_mat[i,:] = np.roll(stddev_reciprocal,i,axis=0)

    alpha_mat = np.copy(dev_mat)/(dev_mat[0,0])
    alpha_mat = 1.0
    dev_mat = (mat - dev_mat)*stddev_recip_mat

    print(dev_mat)

    #map_norm = matplotlib.colors.TwoSlopeNorm(vmin=np.min(dev_mat), vmax=np.max(dev_mat), vcenter=0.0)
    absmax = np.max(np.abs(dev_mat))
    map_norm = matplotlib.colors.TwoSlopeNorm(vmin=-absmax, vmax=absmax, vcenter=0.0)
    # map_norm = matplotlib.colors.TwoSlopeNorm(vmin=-1.0, vmax=1.0, vcenter=0.0)

    im = ax.imshow(dev_mat,cmap='PiYG',norm=map_norm,alpha=alpha_mat)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)

    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label(label=r'Contact Frequency', fontsize=6)
    cbar.ax.tick_params(labelsize=4)


    if CGinfo_flag == True:
        ticks = np.arange(0.0,1.0+0.02,0.5,dtype=np.float32)
        tick_labels = []
        tick_labels.append(r'{:.1f} (ter)'.format(ticks[0]))
        tick_labels.append(r'{:.1f} (ori)'.format(ticks[1]))
        tick_labels.append(r'{:.1f} (ter)'.format(ticks[2]))
        ticks = ticks*CGinfo['N_base_CG']

    if overlay_flag == True:

        overlay_color = 'darkcyan'
        s = 0.9

        for i_ter in range(CGinfo['N_ter']):
            
            if i_ter > 0:

                xmap = float(CGinfo['map'][CGinfo['ter_ranges'][i_ter][0]-1,1])
                ymap = float(CGinfo['map'][CGinfo['ter_ranges'][i_ter][1]-1,1])
                w = float(CGinfo['ter_ranges'][i_ter][1] - CGinfo['ter_ranges'][i_ter][0] + 1)
                h = float(CGinfo['ter_ranges'][i_ter][1] - CGinfo['ter_ranges'][i_ter][0] + 1)

                x = np.array([xmap,xmap,xmap+s*w],dtype=np.float32)
                y = np.array([ymap-s*h,ymap,ymap],dtype=np.float32)

                ax.plot(x,y,linewidth=0.25,color=overlay_color,alpha=0.85)

    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels,fontsize=8,ha='left')
    ax.set_yticks(ticks)
    ax.set_yticklabels(tick_labels,fontsize=8,rotation=90,va='center')

    ax.set_xlabel(r'Combined DNA content',fontsize=10)
    ax.set_ylabel(r'Combined DNA content',fontsize=10)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    plt.tight_layout()

    fig.savefig(out_dir+out_label+'_contacts_dev.pdf',dpi=300)

    return

def plot_mat_diff(out_dir,out_label,mat,CGinfo,fig_size,CGinfo_flag,overlay_flag):

    fig = plt.figure(figsize=(1.15*fig_size[0]*mm,fig_size[1]*mm))

    ax = plt.gca()

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    # if norm_flag == 'scaled':
    #     map_norm = matplotlib.colors.Normalize(vmin=0.0, vmax=thresh)
    # elif norm_flag == 'log':
    #     map_norm = matplotlib.colors.LogNorm(vmin=np.min(mat), vmax=thresh)
    # elif norm_flag == 'logmax':
    #     map_norm = matplotlib.colors.LogNorm(vmin=np.min(mat), vmax=np.max(mat))
    # else:
    #     map_norm = matplotlib.colors.Normalize(vmin=0.0, vmax=np.max(mat))

    # mat_max = np.max(mat)
    # mat_min = np.min(mat)
    mat_max = np.max(np.abs(mat))
    mat_min = -mat_max
    mat_mid = 0.0
    #mat_mid = (mat_min + mat_max)/2.0
    dm = np.mean(np.diagonal(mat))
    print(dm)

    map_norm = matplotlib.colors.TwoSlopeNorm(vmin=mat_min, vcenter=mat_mid, vmax=mat_max)

    im = ax.imshow(mat,cmap='PiYG',norm=map_norm)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.1)

    cbar = fig.colorbar(im, cax=cax, ticks=[-0.1,0,0.1])
    cbar.set_label(label=r'$(\tilde{\mathsf{A}}-\tilde{\mathsf{B}})_{ij}$ - Mapped Contact Frequency Difference', fontsize=7, labelpad=0)
    cbar.ax.tick_params(labelsize=6)

    arrow = dict()
    arrow['width'] = 2
    arrow['headwidth'] = 5
    arrow['headlength'] = 7
    arrow['shrink'] = 0.0
    arrow['ec'] = 'black'
    arrow['fc'] = 'white'

    cax.annotate(r'$\frac{\sum_{i}(\tilde{\mathsf{A}}-\tilde{\mathsf{B}})_{ii}}{\tilde{N}}=$'+' {:.3f}'.format(dm),
                 xy=(-0.11,dm),
                 xytext=(-0.6,dm),
                 xycoords='data',
                 arrowprops=arrow,
                 ha='right',
                 va='center',
                 fontsize=8,
                 annotation_clip=False)

    s_expand = 0.015*float(CGinfo['N_base_CG'])

    ax.set_xlim(xmin=0.0-s_expand,xmax=float(CGinfo['N_base_CG'])+s_expand)
    ax.set_ylim(ymin=float(CGinfo['N_base_CG'])+s_expand,ymax=1.0-s_expand)

    
    if CGinfo_flag == True:
        ticks = np.arange(0.0,1.0+0.02,0.5,dtype=np.float32)
        tick_labels = []
        tick_labels.append(r'{:.1f} (ter)'.format(ticks[0]))
        tick_labels.append(r'{:.1f} (ori)'.format(ticks[1]))
        tick_labels.append(r'{:.1f} (ter)'.format(ticks[2]))
        ticks = ticks*CGinfo['N_base_CG']

    if overlay_flag == True:

        overlay_color = 'red'
        s = 0.95

        for i_ter in range(CGinfo['N_ter']):
            
            if i_ter > 0:

                xmap = float(CGinfo['map'][CGinfo['ter_ranges'][i_ter][0]-1,1])
                ymap = float(CGinfo['map'][CGinfo['ter_ranges'][i_ter][1]-1,1])
                w = float(CGinfo['ter_ranges'][i_ter][1] - CGinfo['ter_ranges'][i_ter][0] + 1)
                h = float(CGinfo['ter_ranges'][i_ter][1] - CGinfo['ter_ranges'][i_ter][0] + 1)

                x = np.array([xmap,xmap,xmap+s*w],dtype=np.float32)
                y = np.array([ymap-s*h,ymap,ymap],dtype=np.float32)

                ax.plot(x,y,linewidth=1.5,color='white',alpha=1.0)
                ax.plot(x,y,linewidth=1.0,color=overlay_color,alpha=1.0)

    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels,fontsize=6,ha='left')
    ax.set_yticks(ticks)
    ax.set_yticklabels(tick_labels,fontsize=6,rotation=90,va='center')

    ax.set_xlabel(r'Mapped DNA content',fontsize=8)
    ax.set_ylabel(r'Mapped DNA content',fontsize=8)

    #plt.tight_layout()

    fig.savefig(out_dir+out_label+'_contacts_diff.png',dpi=300)

    return

def plot_mat_mapped_simple(out_dir,out_label,mat,CGinfo,fig_size,norm_flag,thresh,CGinfo_flag):

    fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

    ax = plt.gca()

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    # if norm_flag == 'scaled':
    #     map_norm = matplotlib.colors.Normalize(vmin=0.0, vmax=thresh)
    # elif norm_flag == 'log':
    #     map_norm = matplotlib.colors.LogNorm(vmin=np.min(mat), vmax=thresh)
    # elif norm_flag == 'logmax':
    #     map_norm = matplotlib.colors.LogNorm(vmin=np.min(mat), vmax=np.max(mat))
    # else:
    #     map_norm = matplotlib.colors.Normalize(vmin=0.0, vmax=np.max(mat))

    if norm_flag == 'scaled':
        map_norm = matplotlib.colors.Normalize(vmin=0.0, vmax=thresh[1])
    elif norm_flag == 'log':
        if np.min(mat) <= 0.0:
            temp_min = mat.flatten()
            temp_min.sort()
            temp_min = temp_min[np.argwhere(temp_min>0.0)[0]]
            mat = np.where(mat==0.0,temp_min,mat)
        else:
            temp_min = np.min(mat)
        temp_max = np.max(mat)
        print('max = '+str(temp_max))
        map_norm = matplotlib.colors.LogNorm(vmin=temp_min, vmax=thresh)
    elif norm_flag == 'logmax':
        if np.min(mat) <= 0.0:
            temp_min = mat.flatten()
            temp_min.sort()
            temp_min = temp_min[np.argwhere(temp_min>0.0)[0]]
            mat = np.where(mat==0.0,temp_min,mat)
        else:
            temp_min = np.min(mat)
        temp_min = np.percentile(mat,10)
        temp_max = np.max(mat)
        print('min = '+str(temp_min))
        print('max = '+str(temp_max))
        map_norm = matplotlib.colors.LogNorm(vmin=temp_min, vmax=temp_max)
    elif norm_flag == 'manual':
        map_norm = matplotlib.colors.LogNorm(vmin=thresh[0], vmax=thresh[1])
    elif norm_flag == 'split':
        mat = np.sqrt(mat)
        mat = -1.0*mat
        map_norm = matplotlib.colors.TwoSlopeNorm(vmin=-thresh[1], vcenter=0.0, vmax=thresh[1])
    else:
        map_norm = matplotlib.colors.Normalize(vmin=0.0, vmax=np.max(mat))

    im = ax.imshow(mat,cmap='PiYG',norm=map_norm)

    s_expand = 0.015*float(CGinfo['N_base_CG'])

    ax.set_xlim(xmin=0.0-s_expand,xmax=float(CGinfo['N_base_CG'])+s_expand)
    ax.set_ylim(ymin=float(CGinfo['N_base_CG'])+s_expand,ymax=1.0-s_expand)

    mat_max = np.max(np.abs(mat))

    print(mat_max)

    
    if CGinfo_flag == True:
        ticks = np.arange(0.0,1.0+0.02,0.5,dtype=np.float32)
        tick_labels = []
        tick_labels.append(r'{:.1f} (ter)'.format(ticks[0]))
        tick_labels.append(r'{:.1f} (ori)'.format(ticks[1]))
        tick_labels.append(r'{:.1f} (ter)'.format(ticks[2]))
        ticks = ticks*CGinfo['N_base_CG']

    ax.set_xticks(ticks)
    ax.set_xticklabels(tick_labels,fontsize=6,ha='left')
    ax.set_yticks(ticks)
    ax.set_yticklabels(tick_labels,fontsize=6,rotation=90,va='center')

    ax.set_xlabel(r'Mapped DNA content',fontsize=8)
    ax.set_ylabel(r'Mapped DNA content',fontsize=8)

    ax.axis('off')
    
    plt.tight_layout()

    #fig.savefig(out_dir+out_label+'_contacts_mapped_simple.pdf',dpi=300)
    fig.savefig(out_dir+out_label+'_contacts_mapped_simple.png',dpi=300)

    plt.close()

    return
