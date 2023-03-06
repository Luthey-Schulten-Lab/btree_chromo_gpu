import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as colormaps
import matplotlib.patches as patches
from matplotlib.lines import Line2D

from mpl_toolkits.axes_grid1 import make_axes_locatable

import os


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

fig_file = '/home/ben/Documents/svn/Minimal_Cell_Chromosome_Organization_2022/tex_figures/replication_topology/raw_plots/topo_replicated_m2_mr1.pdf'
topo_file = '/home/ben/Data/btree_chromo/theta_example/theta_example_scaled_topo_t3.dat'

#fig_file = '/home/ben/Documents/svn/Minimal_Cell_Chromosome_Organization_2022/tex_figures/disentanglement/raw_plots/replication_topology.pdf'
#topo_file = '/home/ben/Data/btree_chromo/disentanglement_testing/wloops_wtopo/wloops_wtopo_topo_rep00001_t20000000.dat'


class topo_single:

    def __init__(self,leaf,size,start_link,start,mid,end,end_link):

        self.leaf = leaf
        self.size = size
        self.start_link = start_link
        self.start = start
        self.mid = mid
        self.end = end
        self.end_link = end_link
        self.circular = False

        if self.end_link == self.start and self.start_link == self.end:

            self.circular = True

    def normalize(self,total_size):

        s = float(total_size)

        self.size = self.size/s
        self.start_link = self.start_link/s
        self.start = self.start/s
        self.mid = self.mid/s
        self.end = self.end/s
        self.end_link = self.end_link/s

        return

class topo_total:

    def __init__(self):

        self.size = 0
        self.N = 0
        self.topos = []

    def parse_topo_file(self,in_topo_file):

        with open(in_topo_file, 'r') as f:

            line = f.readline()
            self.size = int(line.split('=')[1])

            while True:
                line = f.readline()
                if not line:
                    break

                l = line.split(',')
                leaf_info = l[0]

                start_link = int(l[1])
                start = int(l[2])
                mid = int(l[3])
                end = int(l[4])
                end_link = int(l[5])

                leaf_info = leaf_info.split('(')
                leaf = leaf_info[0]
                size = int(leaf_info[1][:-2])

                topo = topo_single(leaf,
                                   size,
                                   start_link,
                                   start,
                                   mid,
                                   end,
                                   end_link)

                self.N += 1
                self.topos.append(topo)

    def get_topo(self,i):

        return self.topos[i]

    def normalize(self):

        for i in range(self.N):

            self.topos[i].normalize(self.size)

        self.size = 1.0

def plot_topology(fig_file,topo_file):

    big_res = 1000
    res = 100

    # cmap = colormaps.get_cmap('plasma')
    # c_space_lower_lim = 0.0
    # c_space_upper_lim = 0.8
    # c_space = np.linspace(c_space_lower_lim,
    #                       c_space_upper_lim,
    #                       N_forks)

    ter_height = 1.0
    fork_height = 0.65
    y_arc_offset = 0.08

    ter_color = 'darkorange'
    ori_color = 'red'
    fork_color = 'magenta'

    ori_marker_size = 8
    ter_marker_size = 5
    fork_marker_size = 5
    
    topo = topo_total()
    topo.parse_topo_file(topo_file)

    topo.normalize()

    fig_size = [87*2.0,87/3.0]

    fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

    ax = plt.gca()

    x = np.linspace(0,topo.size,big_res,dtype=np.double)
    y = np.zeros((big_res),dtype=np.double)

    di = 7
    ax.plot(x[di:-(di+1)],
            y[di:-(di+1)],
            c='white',
            linewidth=7.0,
            zorder=0)


    ax.plot(x[di:-(di+1)],
            y[di:-(di+1)],
            c='black',
            linewidth=6.0,
            zorder=0)

    # needs to be (numlines) x (points per line) x 2 (for x and y)
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    norm = plt.Normalize(0,topo.size)

    lc = matplotlib.collections.LineCollection(segments,
                                               cmap='viridis',
                                               norm=norm)

    lc.set_array(np.linspace(0.0,topo.size,big_res))
    lc.set_linewidth(5.0)
    ax.add_collection(lc)

    x += 1.0/(2.0*big_res)

    # needs to be (numlines) x (points per line) x 2 (for x and y)
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    norm = plt.Normalize(0,topo.size)

    lc = matplotlib.collections.LineCollection(segments,
                                               cmap='viridis',
                                               norm=norm)

    lc.set_array(np.linspace(0.0,topo.size,big_res))
    lc.set_linewidth(5.0)
    ax.add_collection(lc)

    ax.set_xlim(xmin=-0.02,xmax=1.02
                )
    ax.set_ylim(ymin=-0.6,ymax=1.5)
    ax.set_axis_off()

    for i in range(topo.N):

        t = topo.get_topo(i)


        ax.plot([t.mid],
                [0.0],
                marker='o',
                c='white',
                markersize=1.4*ori_marker_size)
        
        ax.plot([t.mid],
                [0.0],
                marker='o',
                c='black',
                markersize=1.2*ori_marker_size)

        ax.plot([t.mid],
                [0.0],
                marker='o',
                c=ori_color,
                markersize=ori_marker_size)

        ax.annotate(text=r''+t.leaf,
                    xy=(t.mid,-0.5),
                    color=ori_color,
                    va='bottom',
                    ha='center')

        if i == 0:

            x_ori_bdry = x=np.array([t.start,t.end])

        else:

            x_ori_bdry = x=np.array([t.end])
            
        ax.vlines(x_ori_bdry,
                  ymin=-0.2,
                  ymax=0.095,
                  lw=1.8,
                  colors='white')

        ax.vlines(x_ori_bdry,
                  ymin=-0.2,
                  ymax=0.095,
                  lw=1.3,
                  colors='black')

        ax.vlines(x=x_ori_bdry,
                  ymin=-0.195,
                  ymax=0.09,
                  lw=0.8,
                  colors=ori_color)

        if t.circular == True:

            s = 0.01
            dx = s*(t.end-t.start)
            x = np.linspace(t.start+dx,t.end-1.8*dx,res,dtype=np.double)

            xm = (x[0]+x[-1])/2.0
            y = ter_height*(x-x[0])*(x-x[-1])/((xm-x[0])*(xm-x[-1])) + y_arc_offset

            ax.plot(x,
                    y,
                    lw=2.5,
                    alpha=1.0,
                    ls='-',
                    c=ter_color,
                    zorder=-5)

            ax.plot([xm],
                    [ter_height + y_arc_offset],
                    marker='s',
                    c='white',
                    markersize=1.4*ter_marker_size)
            
            ax.plot([xm],
                    [ter_height + y_arc_offset],
                    marker='s',
                    c='black',
                    markersize=1.2*ter_marker_size)
            
            ax.plot([xm],
                    [ter_height + y_arc_offset],
                    marker='s',
                    c=ter_color,
                    markersize=ter_marker_size)

            ter_str = t.leaf
            if len(ter_str) > 1:
                ter_str = ter_str[:-1]

            print(ter_str)
                
            ax.annotate(text=r''+ter_str,
                    xy=(xm,ter_height+y_arc_offset+0.4),
                    color=ter_color,
                    va='top',
                    ha='center')

        else:

            s = 0.02
            dx = s*(t.end-t.start)
            
            x = np.linspace(t.start_link,t.start+2.0*dx,res,dtype=np.double)

            xm = (x[0]+x[-1])/2.0
            y = fork_height*(x-x[0])*(x-x[-1])/((xm-x[0])*(xm-x[-1])) + y_arc_offset

            ax.plot(x,
                    y,
                    lw=2.0,
                    alpha=1.0,
                    ls='-',
                    c='white',
                    zorder=-3)
            
            ax.plot(x,
                    y,
                    lw=1.5,
                    alpha=1.0,
                    ls='-',
                    c=fork_color,
                    zorder=-3)

            ax.plot([xm],
                    [fork_height + y_arc_offset],
                    marker='v',
                    c='white',
                    markersize=1.7*fork_marker_size,
                    zorder=-3)
            
            ax.plot([xm],
                    [fork_height + y_arc_offset],
                    marker='v',
                    c='black',
                    markersize=1.5*fork_marker_size,
                    zorder=-3)
            
            ax.plot([xm],
                    [fork_height + y_arc_offset],
                    marker='v',
                    c=fork_color,
                    markersize=fork_marker_size,
                    zorder=-3)


            x = np.linspace(t.end_link,t.end-dx,res,dtype=np.double)

            xm = (x[0]+x[-1])/2.0
            y = fork_height*(x-x[0])*(x-x[-1])/((xm-x[0])*(xm-x[-1])) + y_arc_offset

            ax.plot(x,
                    y,
                    lw=2.5,
                    alpha=1.0,
                    ls='-',
                    c='white',
                    zorder=-3)
            
            ax.plot(x,
                    y,
                    lw=1.75,
                    alpha=1.0,
                    ls='-',
                    c=fork_color,
                    zorder=-3)

            ax.plot([xm],
                    [fork_height + y_arc_offset],
                    marker='v',
                    c='white',
                    markersize=1.7*fork_marker_size,
                    zorder=-3)
            
            ax.plot([xm],
                    [fork_height + y_arc_offset],
                    marker='v',
                    c='black',
                    markersize=1.5*fork_marker_size,
                    zorder=-3)
            
            ax.plot([xm],
                    [fork_height + y_arc_offset],
                    marker='v',
                    c=fork_color,
                    markersize=fork_marker_size,
                    zorder=-3)
            
    plt.tight_layout()

    fig.savefig(fig_file,dpi=300)

    plt.close()

def plot_topology_small(fig_file,topo_file):

    big_res = 1000
    res = 100

    # cmap = colormaps.get_cmap('plasma')
    # c_space_lower_lim = 0.0
    # c_space_upper_lim = 0.8
    # c_space = np.linspace(c_space_lower_lim,
    #                       c_space_upper_lim,
    #                       N_forks)

    ter_height = 1.0
    fork_height = 0.7
    y_arc_offset = 0.08

    ter_color = 'orange'
    ori_color = 'red'
    fork_color = 'magenta'

    ori_marker_size = 8
    ter_marker_size = 6
    fork_marker_size = 5
    
    topo = topo_total()
    topo.parse_topo_file(topo_file)

    topo.normalize()

    fig_size = [87,87/3.0]

    fig = plt.figure(figsize=(fig_size[0]*mm,fig_size[1]*mm))

    ax = plt.gca()

    x = np.linspace(0,topo.size,big_res,dtype=np.double)
    y = np.zeros((big_res),dtype=np.double)


    di = 7
    ax.plot(x[di:-(di+1)],
            y[di:-(di+1)],
            c='white',
            linewidth=7.0,
            zorder=0)


    ax.plot(x[di:-(di+1)],
            y[di:-(di+1)],
            c='black',
            linewidth=6.0,
            zorder=0)

    # needs to be (numlines) x (points per line) x 2 (for x and y)
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    norm = plt.Normalize(0,topo.size)

    lc = matplotlib.collections.LineCollection(segments,
                                               cmap='viridis',
                                               norm=norm)

    lc.set_array(np.linspace(0.0,topo.size,big_res))
    lc.set_linewidth(5.0)
    ax.add_collection(lc)

    x += 1.0/(2.0*big_res)

    #print(x)

    # needs to be (numlines) x (points per line) x 2 (for x and y)
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    norm = plt.Normalize(0,topo.size)

    lc = matplotlib.collections.LineCollection(segments,
                                               cmap='viridis',
                                               norm=norm)

    lc.set_array(np.linspace(0.0,topo.size,big_res))
    lc.set_linewidth(5.0)
    ax.add_collection(lc)

    ax.set_xlim(xmin=-0.02,xmax=1.02
                )
    ax.set_ylim(ymin=-0.5,ymax=1.5) # -0.5,1.3
    ax.set_axis_off()

    for i in range(topo.N):

        t = topo.get_topo(i)


        ax.plot([t.mid],
                [0.0],
                marker='o',
                c='white',
                markersize=1.4*ori_marker_size)
        
        ax.plot([t.mid],
                [0.0],
                marker='o',
                c='black',
                markersize=1.2*ori_marker_size)

        ax.plot([t.mid],
                [0.0],
                marker='o',
                c=ori_color,
                markersize=ori_marker_size)

        ax.annotate(text=r''+t.leaf,
                    xy=(t.mid,-0.45),
                    color=ori_color,
                    va='bottom',
                    ha='center')

        if i == 0:

            x_ori_bdry = x=np.array([t.start,t.end])

        else:

            x_ori_bdry = x=np.array([t.end])

        #print(x_ori_bdry)
            
        ax.vlines(x_ori_bdry,
                  ymin=-0.2,
                  ymax=0.095,
                  lw=2.0,
                  colors='white')

        ax.vlines(x_ori_bdry,
                  ymin=-0.2,
                  ymax=0.095,
                  lw=1.45,
                  colors='black')

        ax.vlines(x=x_ori_bdry,
                  ymin=-0.195,
                  ymax=0.09,
                  lw=0.9,
                  colors=ori_color)

        if t.circular == True:

            s = 0.01
            dx = s*(t.end-t.start)
            x = np.linspace(t.start+dx,t.end-1.8*dx,res,dtype=np.double)

            xm = (x[0]+x[-1])/2.0
            y = ter_height*(x-x[0])*(x-x[-1])/((xm-x[0])*(xm-x[-1])) + y_arc_offset

            ax.plot(x,
                    y,
                    lw=2.5,
                    alpha=1.0,
                    ls='-',
                    c=ter_color,
                    zorder=-5)

            ax.plot([xm],
                    [ter_height + y_arc_offset],
                    marker='s',
                    c='white',
                    markersize=1.4*ter_marker_size)
            
            ax.plot([xm],
                    [ter_height + y_arc_offset],
                    marker='s',
                    c='black',
                    markersize=1.2*ter_marker_size)
            
            ax.plot([xm],
                    [ter_height + y_arc_offset],
                    marker='s',
                    c=ter_color,
                    markersize=ter_marker_size)

            ter_str = t.leaf
            if len(ter_str) > 1:
                ter_str = ter_str[:-1]

            print(ter_str)
                
            ax.annotate(text=r''+ter_str,
                    xy=(xm,ter_height+y_arc_offset+0.4),
                    color=ter_color,
                    va='top',
                    ha='center')

        else:

            s = 0.02
            dx = s*(t.end-t.start)
            
            x = np.linspace(t.start_link,t.start+2.0*dx,res,dtype=np.double)

            xm = (x[0]+x[-1])/2.0
            y = fork_height*(x-x[0])*(x-x[-1])/((xm-x[0])*(xm-x[-1])) + y_arc_offset

            ax.plot(x,
                    y,
                    lw=2.0,
                    alpha=1.0,
                    ls='-',
                    c='white',
                    zorder=-3)
            
            ax.plot(x,
                    y,
                    lw=1.5,
                    alpha=1.0,
                    ls='-',
                    c=fork_color,
                    zorder=-3)

            ax.plot([xm],
                    [fork_height + y_arc_offset],
                    marker='v',
                    c='white',
                    markersize=1.5*fork_marker_size,
                    zorder=-3)
            
            ax.plot([xm],
                    [fork_height + y_arc_offset],
                    marker='v',
                    c='black',
                    markersize=1.3*fork_marker_size,
                    zorder=-3)
            
            ax.plot([xm],
                    [fork_height + y_arc_offset],
                    marker='v',
                    c=fork_color,
                    markersize=fork_marker_size,
                    zorder=-3)


            x = np.linspace(t.end_link,t.end-dx,res,dtype=np.double)

            xm = (x[0]+x[-1])/2.0
            y = fork_height*(x-x[0])*(x-x[-1])/((xm-x[0])*(xm-x[-1])) + y_arc_offset

            ax.plot(x,
                    y,
                    lw=2.5,
                    alpha=1.0,
                    ls='-',
                    c='white',
                    zorder=-3)
            
            ax.plot(x,
                    y,
                    lw=1.75,
                    alpha=1.0,
                    ls='-',
                    c=fork_color,
                    zorder=-3)

            ax.plot([xm],
                    [fork_height + y_arc_offset],
                    marker='v',
                    c='white',
                    markersize=1.5*fork_marker_size,
                    zorder=-3)
            
            ax.plot([xm],
                    [fork_height + y_arc_offset],
                    marker='v',
                    c='black',
                    markersize=1.3*fork_marker_size,
                    zorder=-3)
            
            ax.plot([xm],
                    [fork_height + y_arc_offset],
                    marker='v',
                    c=fork_color,
                    markersize=fork_marker_size,
                    zorder=-3)
            
    plt.tight_layout()

    fig.savefig(fig_file,dpi=300)

    plt.close()

plot_topology_small(fig_file,topo_file)
