#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 09:23:43 2019

@author: OptimusPrime, shiba@vt.edu
"""


# sol = solve_ivp(lambda t,y: DB2dof.vec_field_DB(t,y,*params), \
#                  [0, final_t], ic, method='RK45', events = cross_stateatyw, \
#                  rtol=rtol_val, atol=atol_val, vectorized = True)

# sol = solve_ivp(lambda t,y: DB2dof.vec_field_DB(t,y,*params), \
#                  [0, final_t], ic, method='RK45', \
#                  dense_output=True, \
#                  rtol=rtol_val, atol=atol_val, vectorized = True)

# if len(sol.t_events) >= 1:
#     num_commit_traj += 1


# vals = np.zeros((len(sol.y[1,:]),1))
# # vals = [1 for y_val in sol.y[1,:] if abs(y_val - (-1/np.sqrt(2))) < 1e-4 ]

# for i, y_val in enumerate(sol.y[1,:]):
#     if abs(y_val - (-1/np.sqrt(2))) < 1e-6:
#         vals[i] = 1
#     else:
#         vals[i] = 0

# # vals
# print(np.sum(vals))
# np.shape(vals)

#%%

import matplotlib as mpl
import matplotlib.pyplot as plt


plt.style.use('default')
mpl.rcParams['mathtext.fontset'] = 'cm'
# mpl.rcParams['mathtext.rm'] = 'serif'


ls_tick = 15
ls_axes = 25
mpl.rcParams['xtick.labelsize'] = ls_tick
mpl.rcParams['ytick.labelsize'] = ls_tick
mpl.rcParams['axes.labelsize'] = ls_axes

# if len(sol.t_events[0]) >= 1:
#     eventSol = sol.sol.__call__(sol.t_events[0][:])

plt.close('all')

#figH = plt.figure(figsize = (7,7))
#ax_xy = figH.subplot(211)

figH, axH = plt.subplots(nrows=2, ncols=1,figsize = (10,10))
ax_xy = axH[0]
ax_ty = axH[1]

ax_xy.plot(sol.y[0,:], sol.y[1,:], '-r', linewidth = 1)
#ax_xy.scatter(eventSol[0,:], eventSol[1,:], 10)
ax_xy.set_xlabel(r'$x$')
ax_xy.set_ylabel(r'$y$') 

ax_ty.plot(sol.t, sol.y[1,:], '-r')
#ax_ty.scatter(sol.t_events[0][0:], eventSol[1,:], 10)
ax_ty.set_xlabel(r'$t$')
ax_ty.set_ylabel(r'$y$') 


  
#print(sol.t_events[0])
#print(eventSol[1,:])



#%%


def characteristic_isomerB(sol):
    """
    The characteristic function for region that defines the isomer B 
    is taken to be the y-coordinate less than -YW. Trajectories that
    satisfy this condition are in the configuration of isomer B.
    
    INPUT
    sol: solve_ivp object containing solution of the vector field
    
    OUTPUT
    commit_ts: time series of the committor probability
    """
    
    commit_ts = np.zeros((len(sol.y[1,:]),1))
    
#     for (i, val) in enumerate(sol.y[1,:]):
#         if (sol.y[1,i] < -1/np.sqrt(2.0)):
#             commit_ts[i,0] = 1

    if len(sol.t_events) >= 1:
        commit_ts[:,-1] = 1
        print(sol.t_events)
        
    return commit_ts







