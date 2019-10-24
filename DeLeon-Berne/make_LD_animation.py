#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 11:21:02 2019

@author: Shibabrat Naik, shiba@vt.edu
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
plt.rcParams['animation.ffmpeg_path'] = '/usr/local/bin/ffmpeg' #point to the ffmpeg install path

import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'

label_size = 15
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size
mpl.rcParams['axes.labelsize'] = 20

# plt.style.use('seaborn') # use sans-serif fonts

mpl.rcParams['axes.spines.left'] = True   ## display axis spines
mpl.rcParams['axes.spines.bottom'] = True
mpl.rcParams['axes.spines.top'] = True
mpl.rcParams['axes.spines.right'] = True
mpl.rcParams['xtick.top'] = False
mpl.rcParams['ytick.right'] = False
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['xtick.major.size'] = 6
mpl.rcParams['ytick.major.size'] = 6
mpl.rcParams['xtick.major.width'] = 1.0
mpl.rcParams['ytick.major.width'] = 1.0

data_path = "../../data-figures/DeLeon-Berne/LDs/xpx_sosatyw/fig3C2/fixedLD/"
#data_path = "../../data-figures/DeLeon-Berne/LDs/ypy_sosatx0/fig3A1/"
#data_path = "../../data-figures/DeLeon-Berne/LDs/ypy_sosatx0/fig3A2/"
#data_path = "../../data-figures/DeLeon-Berne/LDs/ypy_sosatx0/fig3B2/"
#data_path = "../../data-figures/DeLeon-Berne/LDs/ypy_sosatx0/fig3C2/"
num_contourlines = 50


def find_max_min_nan_mat(input_mat):
    """
    Find maximum and minimum values in a matrix with NaN values
    """
    idx_row, idx_col = np.where(~np.isnan(input_mat))

    input_mat_stripped = input_mat[idx_row, idx_col]

    max_val = np.max(input_mat_stripped)
    min_val = np.min(input_mat_stripped)
                
    return min_val, max_val


def animate_contour():
    
#    seed(1234)
#    x = uniform(-2,2,100)
#    y = uniform(-2,2,100)
#    data = vstack((x*exp(-x**2-y**2),0.5*x*exp(-x**2-y**2),0.2*x*exp(-x**2-y**2)))
#    xi = linspace(min(x), max(x))
#    yi = linspace(min(y), max(y))
#    zi = []
#    numframes = data.shape[0]
#    for ii in range(numframes):
#        zi.append(griddata((x, y), data[ii], (xi[None,:], yi[:,None]), method='cubic'))

    xRes = 500
    yRes = 1000
    tau = 30
    total_energy_vals = np.arange(1.1,9.2,0.1)
    
    # Loading the domain bounds
    params = np.loadtxt(data_path + 'params_deleonberne_M' + 
                        str(yRes) + 'x' + str(xRes) + '_finalT%.6f'%(tau) + \
                        "_E%.6f"%(total_energy_vals[0]) 
                        + '.txt')
    xMin, xMax = params[:2]
    yMin, yMax = params[3:5]
    tau = int(params[7])
    xVec = np.linspace(xMin, xMax, xRes+1)
    yVec = np.linspace(yMin, yMax, yRes+1)
    
    xMesh, yMesh = np.meshgrid(xVec, yVec)
    
#    MMesh = np.genfromtxt(data_path + "deleonberne_forwM" + \
#                          str(yRes) + "x" + str(xRes) + "_finalT%.6f"%(tau) + \
#                          "_E%.6f"%(total_energy_vals[0]) + \
#                          ".txt",delimiter="\t")
    
    MMesh = np.genfromtxt(data_path + "deleonberne_backM" + 
                          str(yRes) + "x" + str(xRes) + \
                          "_finalT%.6f"%(tau) + "_E%.6f"%(total_energy_vals[0]) + \
                          ".txt",delimiter="\t")

#    MMesh = MMesh_back + MMesh_forw
    min_MMesh, max_MMesh = find_max_min_nan_mat(MMesh)
    
#    fig = plt.figure()
#    im = plt.contour(xi, yi, zi[0], 15, linewidths=0.5, colors='k')
#    ax = fig.gca()
    
    fig_cf = plt.figure(figsize=(7,7))
    ax_cf = fig_cf.gca()
    cf = ax_cf.contourf(xMesh, yMesh, MMesh[:,:], num_contourlines, cmap = 'viridis')
    cbar = fig_cf.colorbar(cf, shrink = 0.7, pad = 0.05, 
                           ticks = [min_MMesh, max_MMesh])
    
    ax_cf.set_xlabel(r'$y$')
    ax_cf.set_ylabel(r'$p_y$') 
    cbar.ax.set_yticklabels(['Min', 'Max'])
    cbar.set_label('LD',labelpad = -20)
    
#    cax = fig_cf.colorbar(cf)
#    cb_ax = fig_cf.add_axes([0.83, 0.1, 0.02, 0.8])
#    cax = fig_cf.colorbar(cf, cax=cb_ax)

    ani = animation.FuncAnimation(fig_cf, update_contour_plot, \
                                  frames = range(len(total_energy_vals)), \
                                  fargs=(fig_cf, ax_cf, xRes, yRes, tau, total_energy_vals), \
                                  interval=100, \
                                  blit = False)
    

    ani.save('test.mp4', fps = 5, dpi=300, bitrate=-1)    #, extra_args=['-vcodec', 'libx264']    
#    ani.save('test.mpeg')
#    ani.save('test.gif')

    plt.show()

    return ani


def update_contour_plot(i, fig_cf, ax_cf, xRes, yRes, tau, total_energy_vals):
    
    # Loading the domain bounds
    params = np.loadtxt(data_path + 'params_deleonberne_M' + 
                        str(yRes) + 'x' + str(xRes) + '_finalT%.6f'%(tau) + \
                        "_E%.6f"%(total_energy_vals[i]) 
                        + '.txt')
    xMin, xMax = params[:2]
    yMin, yMax = params[3:5]
    tau = int(params[7])
    xVec = np.linspace(xMin, xMax, xRes+1)
    yVec = np.linspace(yMin, yMax, yRes+1)
    
    xMesh, yMesh = np.meshgrid(xVec, yVec)
    
#    MMesh = np.genfromtxt(data_path + "deleonberne_forwM" + 
#                          str(yRes) + "x" + str(xRes) + 
#                          "_finalT%.6f"%(tau) + "_E%.6f"%(total_energy_vals[i]) + 
#                          ".txt",delimiter="\t")
    
    MMesh = np.genfromtxt(data_path + "deleonberne_backM" + 
                          str(yRes) + "x" + str(xRes) + \
                          "_finalT%.6f"%(tau) + "_E%.6f"%(total_energy_vals[i]) + \
                          ".txt",delimiter="\t")
    
#    MMesh = MMesh_back + MMesh_forw
    min_MMesh, max_MMesh = find_max_min_nan_mat(MMesh)
    
#    fig_cf.clf()
    ax_cf.cla()
    cf = ax_cf.contourf(xMesh, yMesh, MMesh[:,:], num_contourlines, cmap = 'viridis')
    
#    cbar = fig_cf.colorbar(cf, shrink = 0.7, pad = 0.05, 
#                           ticks = [min_MMesh, max_MMesh])
    
#    cax.cla()
#    fig_cf.colorbar(cf)
#    cax = fig_cf.colorbar(cf, cax=cax)
    ax_cf.set_xlabel(r'$y$')
    ax_cf.set_ylabel(r'$p_y$') 
#    cbar.ax_cf.set_yticklabels(['Min', 'Max'])
#    cbar.set_label('LD',labelpad = -20)
    
    plt.title('Total energy:%.3f'%(total_energy_vals[i]))
    
    return cf,



#%%

animate_contour()








