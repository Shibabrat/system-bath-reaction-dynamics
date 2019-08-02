#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 15 9:02:17 2019

@author: Shibabrat Naik, shiba@vt.edu
"""

import numpy as np
#import matplotlib.pyplot as plt
import matplotlib as mpl

fal = 30 # fontsize axis labels
ftl = 20 # fontsize tick labels
mpl.rcParams['xtick.labelsize'] = ftl
mpl.rcParams['ytick.labelsize'] = ftl
# mpl.rcParams['ztick.labelsize'] = ftl
mpl.rcParams['axes.labelsize'] = fal


#from scipy.optimize import fsolve

def V_DB(x, y, params_pe):
    """
    Potential energy function for the 2 DOF model
    x, y: N x N array of position values as meshgrid data
    params: 1 x 4 array of parameters
    """
    if np.size(x) == 1:
        nX = 1
        nY = 1
    else:
        nX = np.size(x,1)
        nY = np.size(x,0)
        x = np.ravel(x, order = 'C')
        y = np.ravel(y, order = 'C')
    
    epsilon, Dx, alpha, lambd  = params_pe
    
    vy = 4*y**2*(y**2 - 1.0) + epsilon
    vx = Dx*(1.0 - np.exp(-lambd*x))**2
    vyx = 4*y**2*(y**2 - 1)*( np.exp(-alpha*lambd*x) - 1.0 )
    
    pe = np.reshape(vyx + vy + vx, (nY, nX))
#     print(np.size(vy))
#     print(x,y,vyx)
    
    return pe

def get_total_energy(q, p, params_pe, params_ke):
    """
    NEEDS TESTING FOR ARRAY AND ROW INPUTS
    Total energy function 
    q: N x 2 array of position values 
    p: N x 2 array of momenta values
    params_pe: 1 x 4 array of parameters in PE  
    params_ke: 1 x 2 array of parameters in KE
    """
    
    m, mu = params_ke # mass of isomer A and isomer B
    
    total_energy = np.zeros(np.size(q,0))
    for idx in range(np.size(q,0)):
        total_energy[idx] = p[idx,0]**2/(2*m) + p[idx,1]**2/(2*mu) + \
                                V_DB(q[idx,0], q[idx,1], params_pe)
          
    return total_energy


def get_pot_surf_proj(xVec, yVec, params_pe):            

    resX = np.size(xVec)
    resY = np.size(xVec)
    surfProj = np.zeros([resX, resY])
    for i in range(len(xVec)):
        for j in range(len(yVec)):
            surfProj[i,j] = V_DB(xVec[j], yVec[i], params_pe)

    return surfProj


def get_Hills_Region(xVec, yVec, deltaE):

    E = ESADDLE + deltaE
    resX = np.size(xVec)
    eProjCS = np.zeros([resX, resX])

    for i in range(len(yVec)):
        for j in range(len(xVec)):
            eProjCS[i,j] = V_DB(xVec[j], yVec[i], params_pe)
            if ( eProjCS[i,j] > E ):
                eProjCS[i, j] = 0
            else:
                eProjCS[i,j] = np.nan

    return eProjCS



def vec_field_DB(states, *params):
    
    massA, massB, epsilon, Dx, alpha, lambd = params
    
    q1, q2, p1, p2 = states
    q1Dot = p1/massA
    q2Dot = p2/massB
    p1Dot = 2*Dx*lambd*(np.exp(-lambd*q1) - 1)*np.exp(-lambd*q1) + \
            4*alpha*lambd*q2**2*(q2**2 - 1)*np.exp(-alpha*lambd*q1)
    p2Dot = -8*(2*q2**3 - q2)*np.exp(-alpha*lambd*q1)
    
    return np.array([q1Dot, q2Dot, p1Dot, p2Dot])

def get_area_closed_curve(inCurve):
    
    """
    get area enclosed by a closed curve in 2D using Green's theorem
    assumes input curve is not closed
    """
    
    closedCurve = np.append(inCurve, [inCurve[-1,:]], axis = 0)
    
    dC = np.diff(closedCurve, axis = 0)
    
    curveArea = 0
    for i in range(0,np.size(inCurve,0)):
        curveArea = curveArea + 0.5*( inCurve[i,0]*dC[i,1] - \
                                     inCurve[i,1]*dC[i,0])
    
    return np.abs(curveArea)



