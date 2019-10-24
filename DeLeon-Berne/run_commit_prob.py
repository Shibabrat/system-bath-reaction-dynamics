#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 11:18:34 2019

@author: Shibabrat Naik, shiba@vt.edu
"""

import numpy as np
from scipy.integrate import solve_ivp

import DeLeonBerne2dof
import importlib
importlib.reload(DeLeonBerne2dof)
import DeLeonBerne2dof as DB2dof



## Integration to compute committor probabilities

# low tolerance
# rtol_val = 1e-8
# atol_val = 1e-10

# high tolerance
rtol_val = 1e-12
atol_val = 1e-14


def cross_stateatyw(t, y):    
    """
    Definition of the section located at y = -yw in the lower well 
    """
    YW  = 1/np.sqrt(2.0)

    return y[1] - (-YW)

cross_stateatyw.terminal = True
cross_stateatyw.direction = 0
    
def get_commit_prob(t_eval, new_ic, params, filename):
    
    commit_ts = np.zeros((len(t_eval)-1,2))
    j = 0
    for (ti,tf) in zip(t_eval[:-1], t_eval[1:]):
        ic = new_ic
        commit_ts[j,0] = tf
        for i in range(len(ic[:,0])):
#         print(ti, tf, ic[i,:])      
            sol = solve_ivp(lambda t,y: DB2dof.vec_field_DB(t,y,*params), \
                            [ti, tf], ic[i,:], method='RK45', \
                            dense_output = True, \
                            events = cross_stateatyw, \
                            rtol=rtol_val, atol=atol_val, vectorized = True)
        
            #update committor probabilities and initial condition 
            new_ic[i,:] = sol.y[:,-1]
            if sol.y[1,-1] - (-1/np.sqrt(2.0)) < 1e-10:
                commit_ts[j,1] += 1
                
        
        j += 1   
        
    np.savetxt(filename,commit_ts)

    return commit_ts



#%%

params = DB2dof.set_parameters('fig3B2/')
initconds_set1 = np.loadtxt('initconds_set1.txt')
initconds_set2 = np.loadtxt('initconds_set2.txt')
initconds_set3 = np.loadtxt('initconds_set3.txt')
initconds_set4 = np.loadtxt('initconds_set4.txt')
initconds_set5 = np.loadtxt('initconds_set5.txt')

num_initconds = np.size(initconds_set1,0)
final_t = -100
delta_t = -1e-2
t_eval = np.arange(0, final_t + delta_t, delta_t)


# new_ic = initconds_set5
commit_ts_set1 = get_commit_prob(t_eval, initconds_set1, params, \
                                'commit_prob_'+ str(num_initconds)+'icset1.txt')
commit_ts_set2 = get_commit_prob(t_eval, initconds_set2, params, \
                                 'commit_prob_'+ str(num_initconds)+'icset2.txt')
commit_ts_set3 = get_commit_prob(t_eval, initconds_set3, params, \
                                 'commit_prob_'+ str(num_initconds)+'icset3.txt')
commit_ts_set4 = get_commit_prob(t_eval, initconds_set4, params, \
                                 'commit_prob_'+ str(num_initconds)+'icset4.txt')
commit_ts_set5 = get_commit_prob(t_eval, initconds_set5, params, \
                                 'commit_prob_'+ str(num_initconds)+'icset5.txt')



#%%









