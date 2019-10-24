#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 15 9:02:17 2019

@author: Shibabrat Naik, shiba@vt.edu
"""

import numpy as np
import scipy as sp
from scipy.optimize import fsolve
#import matplotlib.pyplot as plt
import matplotlib as mpl

fal = 30 # fontsize axis labels
ftl = 20 # fontsize tick labels
mpl.rcParams['xtick.labelsize'] = ftl
mpl.rcParams['ytick.labelsize'] = ftl
# mpl.rcParams['ztick.labelsize'] = ftl
mpl.rcParams['axes.labelsize'] = fal


#from scipy.optimize import fsolve


# Switch cases for different parameter values

# Fig. 3-A1
# ALPHA = 0.20;
# LAMBDA = 1.00;

# Fig. 3-A2
# ALPHA = 1.00;
# LAMBDA = 1.00;

# Fig. 3-B1
# ALPHA = 1.00;
# LAMBDA = 1.30;

# Fig. 3-B2
# ALPHA = 1.00;
# LAMBDA = 1.5;

# % Fig. 3-C1
# % ALPHA = 1.00;
# % LAMBDA = 2.00;

# % Fig. 3-C2
# ALPHA = 2.30;
# LAMBDA = 1.95;

def set_parameters(label):
    
    # ordering of the parameters: mass A, mass B, epsilon, Dx, alpha (z in paper), lambd
    massA = 8; massB = 8; EPSILON = 1.0; Dx = 10.0;
    
    if label == 'fig3A1/':
        ALPHA = 0.20;
        LAMBDA = 1.00;
    elif label == 'fig3A2/':
        ALPHA = 1.00;
        LAMBDA = 1.00;
    elif label == 'fig3B1/':
        ALPHA = 1.00;
        LAMBDA = 1.30;
    elif label == 'fig3B2/':
        ALPHA = 1.00;
        LAMBDA = 1.50;
    elif label == 'fig3C1/':
        ALPHA = 1.00;
        LAMBDA = 2.00;
    elif label == 'fig3C2/':
        ALPHA = 2.30;
        LAMBDA = 1.95;
    elif label == 'ZETA0.00-LAMBDA1.00/':
        ALPHA = 0.00;
        LAMBDA = 1.00;
    else:
        print('Incompatible arguments, check label string')
        
    params =  [massA, massB, EPSILON, Dx, ALPHA, LAMBDA]
    
    return params

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
    
    if np.size(x) == 1:
        pe = pe[0,0]

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



def vec_field_DB(t, states, *params):
    
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


def px_zero(x, *args):
    """
    function to solve minimum and maximum value of x-coordinate on the 
    surface (x,px) with y = yw, px = py = 0
    """

    Eval, EPSILON, Dx, ALPHA, LAMBDA = args[0:]
    
    H = lambda x: (Dx*(1 - np.exp(-LAMBDA*x))**2 - np.exp(-ALPHA*LAMBDA*x) + EPSILON)
    
    return (Eval - H(x))

def get_energybndry_intersect_sos(params, deltaEnergy, printFlag = False):
    
    massA, massB, EPSILON, Dx, ALPHA, LAMBDA = params 
    saddleEnergy = EPSILON
    
    totalEnergy = saddleEnergy + deltaEnergy
    
    H = lambda x: (Dx*(1 - np.exp(-LAMBDA*x))**2 - np.exp(-ALPHA*LAMBDA*x) + EPSILON)
        
    # minumum and maximum of the energy boundary's x coordinate on the sos
    xMax = fsolve(px_zero,  0.5, args = (totalEnergy,EPSILON, Dx, ALPHA, LAMBDA))
    xMin = fsolve(px_zero,  -1, args = (totalEnergy,EPSILON, Dx, ALPHA, LAMBDA))
#     print(xMin, xMax)
    
    xGrid = np.linspace(xMin,xMax,1000)
    
    pxGridPos = np.zeros(np.size(xGrid))
    pxGridNeg = np.zeros(np.size(xGrid))
    
    for i in range(np.size(xGrid)):
        if (totalEnergy > H(xGrid[i])):
            pxGridPos[i] = sp.sqrt((2*massA)*(totalEnergy - H(xGrid[i])))
            
            pxGridNeg[i] = -sp.sqrt((2*massA)*(totalEnergy - H(xGrid[i])))
            
    
#     print('Excess energy:', deltaEnergy)    
    
    if printFlag:
        fig = plt.figure(figsize=(6, 6))
        ax = fig.gca()

    #     ax = axisH
    #     plt.rcParams['axes.labelweight'] = 'bold'
    #     plt.rcParams['ytick.major.size'] = 5    
    #     plt.rcParams['xtick.major.size'] = 5

    #    ax.xaxis.set_ticks([0, 5, 10, 15])
    #    ax.yaxis.set_ticks(np.linspace(-80,80,9,endpoint=True,dtype=int))

    #     ax.axis([0, rMax + 1, -80, 80])
    #     ax.set_xticklabels(ax.get_xticks(), weight='bold')
    #     ax.set_yticklabels(ax.get_yticks(), weight='bold')

        ax.plot(xGrid, pxGridPos, linewidth = 2, color = 'm')
        ax.plot(xGrid, pxGridNeg, linewidth = 2, color = 'm')

        ax.set_xlabel(r'$x$', fontsize = 25)
        ax.set_ylabel(r'$p_x$', fontsize = 25)
        plt.tick_params(axis = 'both', which = 'major', labelsize = 15)
        
        
    xGrid = np.append(xGrid, np.flip(xGrid,0))
    pxGrid = np.append(pxGridPos, np.flip(pxGridNeg,0))
    
#     print(np.shape(xGrid))
        
    return np.array([xGrid,pxGrid]).T


