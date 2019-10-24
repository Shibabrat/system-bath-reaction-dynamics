#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 17 00:21:35 2018

@author: Shibabrat Naik
"""
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import animation
plt.rcParams['animation.ffmpeg_path'] = '/usr/local/bin/ffmpeg' #point to the ffmpeg install path


import importlib
import ball_rolling
importlib.reload(ball_rolling)
from ball_rolling import g, H_polar, Hr, eSaddle, U, Hx, Hy, \
                        tot_energy, tot_energy_polar, get_traj_energy, \
                        ball_rolling2dof, pr_zero, get_energy_polar_U1, \
                        energy_bndry_polar2cart_U1, get_ptheta, \
                        coords_polar2cart_U1, coords_cart2polar_U1

from pylab import rcParams
    
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'

rcParams['figure.figsize'] = 5, 5 

label_size = 10
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size
mpl.rcParams['axes.labelsize'] = 15

mpl.rcParams['font.weight'] = 'bold'

#"""
#Note: Trajectories on the stable manifold spends a lot of time near the
#periodic orbit, so plotting is made efficient by ignoring some of the initial
#points, skip_pts. Next, to render the surface of the tube, I skip some 
#trajectories, skip_trajs
#"""

skip_pts = int(100)
skip_trajs = int(10)

# line widths for trajectories on the tube
tube_lw = 0.6
po_lw = 0.4

#%% Plot the cylindrical manifold and the intersection with the PSS
from matplotlib import colors as mcolors
from matplotlib.collections import PolyCollection

def cc(arg):
    '''
    Shorthand to convert 'named' colors to rgba format at 60% opacity.
    '''
    return mcolors.to_rgba(arg, alpha=0.6)

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_xlim([0, 15])    
ax.set_ylim([0, 15])  
  
ax.set_xlabel(r'$x \; ({\rm cm})$')
ax.set_ylabel(r'$y \; ({\rm cm})$')    
ax.set_zlabel(r'$v_y \; ({\rm cm/s})$')
file1 = './data/deltaE-100/trajs-500/sosU1p_mfd_params.txt'
file2 = './data/deltaE-100/trajs-500/sosU1p_mfd.txt'
file3m = './data/deltaE-100/trajs-500/sosU1m_xeU1.txt'
file3p = './data/deltaE-100/trajs-500/sosU1p_xeU1.txt'
file4 = './x0po_T_energyPO_DelE100.txt'

# Load the manifold intersection with PSS data
xeU1m = np.loadtxt(file3m)
xeU1p = np.loadtxt(file3p)

# Load periodic orbit data
po = np.loadtxt(file4)
E = po[-1]


# Polar coordinates of the energy boundary
deltaE = 100
numPts = 200 # for the boundary
eBndryPolar = np.empty([numPts, 3]) 
eBndryPolar = get_energy_polar_U1(deltaE, numPts)

eBndryCart = np.empty([numPts, 4])
eBndryCart = energy_bndry_polar2cart_U1(eBndryPolar, deltaE)

# Loading the manifold data
mfd_params = np.loadtxt(file1)
num_trajs = int(mfd_params[-1])
trajs = mfd_params[0:int(num_trajs)].astype(int)
mfd = np.loadtxt(file2)

traj = []
init_pts = np.zeros([int(num_trajs/skip_trajs),4])
inter_U1 = np.zeros([int(num_trajs/skip_trajs),4])

traj = mfd[0: trajs[0]]
init_pts[0] = traj[0,:]
inter_U1[0] = traj[-1,:]

verts_xeU1p = []
verts_xeU1p.append(list(zip(xeU1p[:,1], xeU1p[:,3])))
poly_xeU1p = PolyCollection(verts_xeU1p, facecolors=[cc('blue')])

verts_xeU1m = []
verts_xeU1m.append(list(zip(xeU1m[:,1], xeU1m[:,3])))
poly_xeU1m = PolyCollection(verts_xeU1m, facecolors=[cc('cyan')])

plt.set_cmap('viridis')
#ax.plot(traj[skip_pts:,0], traj[skip_pts:,1], traj[skip_pts:,3], '-', 
#        c = 'cyan', lw = tube_lw, 
#        label = r'$\Delta E = %0.2f$' %deltaE)
ax.plot(traj[skip_pts:,0], traj[skip_pts:,1], traj[skip_pts:,3], '-', 
        c = 'cyan', lw = tube_lw, 
        label = r'$\Delta E = %d$' %deltaE)
ax.plot(traj[skip_pts:,0], traj[skip_pts:,1], '-', zs = -20, zdir = 'z', 
        c = 'cyan', lw  = tube_lw)

### extract individual trajectory
for i in np.arange(0,num_trajs,skip_trajs):
    traj = mfd[np.sum(trajs[0:i], dtype = int) : np.sum(trajs[0:i+1], 
               dtype = int)]    
    init_pts[int(i/skip_trajs)] = traj[0,:]
    inter_U1[int(i/skip_trajs)] = traj[-1,:]
   
    ax.plot(traj[skip_pts:,0], traj[skip_pts:,1],  traj[skip_pts:,3], 
            '-', c = 'cyan', lw  = tube_lw)
    ax.plot(traj[skip_pts:,0], traj[skip_pts:,1], '-', zs = -20, 
            zdir = 'z', c = 'cyan', lw  = tube_lw)
    
ax.plot(init_pts[:,0], init_pts[:,1], init_pts[:,3], '-', c = 'r',
        lw  = 1)
ax.plot(inter_U1[:,0], inter_U1[:,1], inter_U1[:,3], '-', c = 'k',
        lw  = 1)
ax.plot(xeU1m[:,0], xeU1m[:,1], xeU1m[:,3], '-', c = 'k',
        lw  = 1)

ax.plot(init_pts[:,0], init_pts[:,1], '-', zs = -20, zdir = 'z',
        c = 'r', lw  = 1)
ax.plot(inter_U1[:,0], inter_U1[:,1], '-', zs = -20, zdir = 'z',
        c = 'k', lw  = 1)
ax.plot(xeU1m[:,0], xeU1m[:,1], '-', zs = -20, zdir = 'z',
        c = 'k', lw  = 1)

# plot the intersection of the manifolds with the PSS
#ax.plot(xeU1m[:,1], xeU1m[:,3], '-', zs = 12, zdir = 'x', 
#        c = 'black', lw  = 1)
#ax.plot(xeU1p[:,1], xeU1p[:,3], '-', zs = 12, zdir = 'x', 
#        c = 'black', lw  = 1)

#ax.add_collection3d(poly_xeU1m, zs = 12, zdir='x')
#ax.add_collection3d(poly_xeU1p, zs = 12, zdir='x')

#ax.plot(eBndryCart[:,1], eBndryCart[:,3], '-', zs = 8, zdir = 'x', 
#        c = 'magenta', lw  = 2)


ax.set_aspect('equal')

leg = ax.legend(fontsize=15)
for line in leg.get_lines():
    line.set_linewidth(3.0)

#ax.view_init(elev = -25., azim = 15)


#%% 
"""
Make movie of the stable tubes for 3 values of excess energy
"""

energyVals = np.arange(100,700,200)
scaled_cval = [ (float(j)-energyVals[0])/(energyVals[-1] - energyVals[0]) 
                for j in energyVals]


fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_xlim([0, 15])    
ax.set_ylim([0, 15])    
ax.set_xlabel(r'$x \; ({\rm cm})$')
ax.set_ylabel(r'$y \; ({\rm cm})$')    
ax.set_zlabel(r'$v_y \; ({\rm cm/s})$')

def init():
    

    for idx, j in enumerate(energyVals):
        
        # Loading data files of manifolds
        filename1 = './data/deltaE-' + str(j) + '/trajs-1000/mfd_params.txt'
        filename2 = './data/deltaE-' + str(j) + '/trajs-1000/mfd.txt'
#        file3m = './data/deltaE-' + str(j) + '/trajs-500/sosU1m_xeU1.txt'
#        file3p = './data/deltaE-' + str(j) + '/trajs-500/sosU1p_xeU1.txt'
#        file4 = './x0po_T_energyPO_DelE100.txt'
        
#        filename1 = './mfd_params.txt'
#        filename2 = './mfd.txt'
        
        
        # Loading the data files
        #mfd_params = np.loadtxt('./data/deltaE-200/trajs-100/mfd_params.txt')
        mfd_params = np.loadtxt(filename1)
        num_trajs = int(mfd_params[-1])
        trajs = mfd_params[0:int(num_trajs)].astype(int)
        
        #mfd = np.loadtxt('./data/deltaE-200/trajs-100/mfd.txt')
        mfd = np.loadtxt(filename2)
        
        traj = []
        init_pts = np.zeros([int(num_trajs/skip_trajs),4])
        inter_U1 = np.zeros([int(num_trajs/skip_trajs),4])
        
        traj = mfd[0: trajs[0]]
        init_pts[0] = traj[0,:]
        inter_U1[0] = traj[-1,:]
    
        plt.set_cmap('viridis')
        ax.plot(traj[skip_pts:,0], traj[skip_pts:,1], traj[skip_pts:,3], '-', 
                c = cm.viridis(scaled_cval[idx]), 
                lw  = tube_lw)
        
        ### extract individual trajectory
        for i in np.arange(0,num_trajs,skip_trajs):
            traj = mfd[np.sum(trajs[0:i], dtype = int) : np.sum(trajs[0:i+1], dtype = int)]    
            init_pts[int(i/skip_trajs)] = traj[0,:]
            inter_U1[int(i/skip_trajs)] = traj[-1,:]
           
            ax.plot(traj[skip_pts:,0], traj[skip_pts:,1],  traj[skip_pts:,3], 
                    '-', c = cm.viridis(scaled_cval[idx]), 
                    lw  = tube_lw)
            
            
        ax.plot(init_pts[:,0], init_pts[:,1], init_pts[:,3], '-', 
                c = cm.viridis(scaled_cval[idx]), 
                lw  = po_lw, label = r'$\Delta E = %d$'%j)
        ax.plot(inter_U1[:,0], inter_U1[:,1], inter_U1[:,3], '-', 
                c = cm.viridis(scaled_cval[idx]), 
                lw  = po_lw)
    
        leg = ax.legend(fontsize=15)
        for line in leg.get_lines():
            line.set_linewidth(3.0)
        
        ax.set_xlim([0, 12])    
        ax.set_ylim([0, 12])
        ax.set_zlim([-20, 20])        
        ax.set_xticks([0, 4, 8, 12])
        ax.set_yticks([0, 4, 8, 12])
        ax.set_yticks([-20, -10, 0, 10, 20])

        
#        ax.view_init(elev = 30., azim = -115)
    return fig,

def animate(i):
#    for angle in range(0,275,5)
        
    ax.view_init(elev = 30., azim = i)
#    plt.draw()
#    plt.pause(0.1)
#    ax.view_init(elev=10., azim=i)
    return fig,


#init()

# Animate
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=30, interval=15, blit=False)
# Save
#mywriter = animation.FFMpegWriter()
anim.save('tube_animation.mp4', fps=15, dpi=600, bitrate=-1)                     
#plt.show()

#%% 
"""
Make movie of a single trajectory on the stable tubes for 3 values of 
excess energy
"""



energyVals = np.arange(100,700,200)
scaled_cval = [ (float(j)-energyVals[0])/(energyVals[-1] - energyVals[0]) 
                for j in energyVals]


def init():
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    
    for idx, j in enumerate(energyVals):
    
        
        filename1 = './data/deltaE-' + str(j) + '/trajs-1000/mfd_params.txt'
        filename2 = './data/deltaE-' + str(j) + '/trajs-1000/mfd.txt'
        
        # Loading the data files
        #mfd_params = np.loadtxt('./data/deltaE-200/trajs-100/mfd_params.txt')
        mfd_params = np.loadtxt(filename1)
        num_trajs = int(mfd_params[-1])
        trajs = mfd_params[0:int(num_trajs)].astype(int)
        
        #mfd = np.loadtxt('./data/deltaE-200/trajs-100/mfd.txt')
        mfd = np.loadtxt(filename2)
        
        traj = []
        init_pts = np.zeros([int(num_trajs/skip_trajs),4])
        inter_U1 = np.zeros([int(num_trajs/skip_trajs),4])
        
        traj = mfd[0 : trajs[0]]
        init_pts[0] = traj[0,:]
        inter_U1[0] = traj[-1,:]
    
        plt.set_cmap('viridis')
        ax.plot(traj[:,0], traj[:,1], traj[:,3], '-', 
                c = cm.viridis(scaled_cval[idx]), 
                lw  = 2)
        
        ### extract individual trajectory
        for i in np.arange(0,num_trajs,skip_trajs):
            
            traj = mfd[np.sum(trajs[0:i], dtype = int) : np.sum(trajs[0:i+1], dtype = int)]    
            init_pts[int(i/skip_trajs)] = traj[0,:]
            inter_U1[int(i/skip_trajs)] = traj[-1,:]
    #        
    #        ax.plot(traj[skip_pts:,0], traj[skip_pts:,1],  traj[skip_pts:,3], 
    #                '-', c = cm.viridis(scaled_cval[idx]), 
    #                lw  = 0.1)
            
            
        ax.plot(init_pts[:,0], init_pts[:,1], init_pts[:,3], '-', 
                c = cm.viridis(scaled_cval[idx]), 
                lw  = 1, label = r'$\Delta E = %d$'%j)
        ax.plot(inter_U1[:,0], inter_U1[:,1], inter_U1[:,3], '-', 
                c = cm.viridis(scaled_cval[idx]), 
                lw  = 1)
        
    ax.set_xlim([0, 15])    
    ax.set_ylim([0, 15])    
    ax.set_xlabel(r'$x \; ({\rm cm})$')
    ax.set_ylabel(r'$y \; ({\rm cm})$')    
    ax.set_zlabel(r'$v_y \; ({\rm cm/s})$')

    return fig,

def animate(i):

#    for angle in range(0,275,5):
    ax.legend(fontsize=15)
    ax.view_init(elev = 20., azim = i)
#    plt.draw()
#    plt.pause(0.1)
#    ax.view_init(elev=10., azim=i)

    return fig,

# Animate
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=360, interval=20, blit=False)
# Save
#mywriter = animation.FFMpegWriter()
anim.save('traj_animation.mp4', fps=30, dpi=300, bitrate=-1)   #, extra_args=['-vcodec', 'libx264']                 
plt.show()


#%%  Mayavi visualization (work in progress)
from mayavi import mlab

skip_pts = int(150)

traj = []
traj = mfd[0: trajs[0]]
ax.plot(traj[:,0], traj[:,1], traj[:,3], '-c', lw  = 0.1)

### extract individual trajectory
for i in np.arange(0,num_trajs,2):
    traj = mfd[np.sum(trajs[0:i], dtype = int) : np.sum(trajs[0:i+1], dtype = int)]
    s = mlab.plot3d(traj[:,0], traj[:,1], traj[:,3], color=(0,1,1), 
                    line_width=1, extent = [0,1,0,1,0,1])

    
mlab.show()








