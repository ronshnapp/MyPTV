#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat March 19 2022

@author: ron


Implementation of a method for smoothing trajectories.
The method uses polynomial fitting after Luti et al 2005.

"""

from utils import fit_polynomial
from numpy import dot



def smooth_traj_poly(traj, window, polyorder):
    '''
    will smooth the particle position using a mooving polynomial, where
    the velocities and accelerations are calculated by
    differentiating the polynomial coeficcients analytically.
    
    input - 
    traj - a nested list (or numpy array), with shapes (3,N), N being
           the number of position samples.
    windows - integer or tuple with length 3. The window sliding size
              used for smoothing.
    polyorder - integer, the order of the polynomial used for the fitting.
    
    
    returns-
    new_pos - the smoothed particle's position
    new_vel - the particle velocity
    new_acc - the particle acceleration
    '''
    
    if type(window) == int:
        window = (window, window, window)
    
    if type(polyorder) == int:
        polyorder = (polyorder, polyorder, polyorder)
    
    if type(window) != tuple: 
        raise TypeError('window must be either integer or tuple')
    
    if type(polyorder) != tuple: 
        raise TypeError('polyorder must be either a number or tuple')
    
    test = [w%2 != 1 for w in window]
    if sum(test) != 0:
        raise ValueError('windows must be a positive odd integer.')
    
    test = [polyorder[i] > window[i] for i in [0,1,2]]
    if sum(test) != 0:
        raise ValueError('polyorder cannot be larger than window.')
    
    
    N = len(traj[0])
    
    test = [w > N for w in window]
    if sum(test) != 0:
        raise ValueError('window cannot be larger than the trajectory length.')
    
    
    
    new_pos = [[0.0 for i in range(N)] for j in range(3)] 
    new_vel = [[0.0 for i in range(N)] for j in range(3)] 
    new_acc = [[0.0 for i in range(N)] for j in range(3)]
    
    sequence = zip(range(3), window, polyorder)
    for e,win,po in sequence:
        time_ = range(win)
        ev_point = [float(win/2)**i for i in range(po + 1)][::-1]
        
        Deriv_mat = []
        for i in range(po+1):
            a = [0.0 for j in range(po+1)]
            if i!=0:
                a[i-1] = po - (i-1)
            Deriv_mat.append(a)
        # Deriv_mat = np.array(Deriv_mat)
        
        
        for i in range(N):
            if i < win/2:            # get the portion of size widow for
                p_ = traj[e][:win]   # fitting the polynomial
            elif N - (1+i) < win/2:
                p_ = traj[e][-1*win:]
            else:
                p_ = traj[e][i-int(win/2) : i+int(win/2)+1]
            
            C = fit_polynomial(time_, p_, po)      
            C1 = dot(Deriv_mat, C)
            C2 = dot(Deriv_mat, C1)
            
            
            if i < win/2: 
                ev_point = [float(i%(win/2))**k for k in range(po + 1)[::-1]]
            elif N - (1+i) < win/2:
                ev_point = [float(win + i - N)**k for k in range(po + 1)[::-1]]
            else:
                ev_point = [float(win/2-0.5)**k for k in range(po + 1)[::-1]]
            
            new_pos[e][i] = dot(ev_point, C)
            new_vel[e][i] = dot(ev_point, C1)
            new_acc[e][i] = dot(ev_point, C2)

    return new_pos, new_vel, new_acc

