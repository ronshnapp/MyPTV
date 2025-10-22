# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Created on Sun Sep 11 17:52:50 2022

@author: ron


A few useful functions to analyze Lagrangian data.

"""
import numpy as np
import pandas as pd



def load_trajs_as_arrays(fname):
    '''
    given a path, fname, this returns a list of arrays, each corresponding
    to a particlular trajectory.
    '''
    from pandas import read_csv
    from numpy import array
    data = read_csv(fname, sep='\t', header=None)
    trajs = [array(sorted(array(g), key=lambda x: x[-1])) 
             for k,g in data.groupby(0) if k!=-1]
    return trajs




def load_samles_vs_time(fname):
    '''
    given a path, fname, this returns a dictionary where keys are
    frame numbers and values are arrays that contain the samples
    of the flow at these frames.
    '''
    from pandas import read_csv
    from numpy import array
    data = read_csv(fname, sep='\t', header=None)
    samples = [(k,array(g)) for k,g in data.groupby(10) if k!=-1]
    return dict(samples)





def is_inside_a_box(tr, xmin, xmax, ymin, ymax, zmin, zmax):
    '''
    given a trajectory (an array with convention as the result of 
    the smoothing function), this returns a list of boolean variable
    saying whether the coordinate at the corresponding time is inside the
    given box domain.
    
    Note that this can be used with numpy.where on a trajectory
    '''
    res = []
    for xi in tr:
        xin = xi[1] < xmax and xi[1] > xmin
        yin = xi[2] < ymax and xi[2] > ymin
        zin = xi[3] < zmax and xi[3] > zmin
        res.append(xin*yin*zin)
    return res




# ============================================================================
#    Lagrangian Velocity Statistics:


def get_velocity_mean_std(traj_list, kind='x'):
    '''
    For a list of trajectories, this returns the mean and standard deviation 
    of a velocity component. The "kind" parameter defines which component to 
    use: 'x', 'y', 'z' or 'KE' standing for x, y, or z velocity component and 
    KE is the kinetic energy.
    '''
    
    if kind=='x':
        get_component = lambda tr: tr[:,4]
    
    elif kind=='y':
        get_component = lambda tr: tr[:,5]
        
    elif kind=='z':
        get_component = lambda tr: tr[:,6]
        
    elif kind=='KE':
        get_component = lambda tr: 0.5*(np.sum(tr[:,4:7]**2, axis=1))
    
    A = []
    
    for tr in traj_list:
        A += list(get_component(tr))
        
    return np.mean(A), np.std(A)



def get_trajectory_velocity_increments(traj, kind='x'):
    '''
    Given a trajectory, this function returns lists of the temporal increments
    of its velocity. The "kind" parameter defines which component to use:
    'x', 'y', 'z' or 'KE' standing for x, y, or z velocity component and KE is 
    the kinetic energy.
    
    We return a nested list where the second sublist
    contains samples of KE(t+1) - KE(t), and the last sublist contains
    KE(t+T) - KE(t) where KE is the kinetic energy and T is the tracking time
    of the trajectory.
    '''
    if kind=='x':
        A = traj[:,4]
    
    elif kind=='y':
        A = traj[:,5]
        
    elif kind=='z':
        A = traj[:,6]
        
    elif kind=='KE':
        A = 0.5*(np.sum(traj[:,4:7]**2, axis=1))
        
    else:
        raise ValueError('undefined kind "%s"'%kind)
        
    increments = [[0]]
    for i in range(1,len(A)):
        increments.append(list(A[i:] - A[:-i]))
    
    return increments



def get_velocity_increments(traj_list, kind='x'):
    '''
    For a list of trajectories, this returns a nested list of all kinetic energy
    increments samples.
    '''
    increments = []
    
    for tr in traj_list:
        inc = get_trajectory_velocity_increments(tr, kind=kind)
        for i in range(len(inc)):
            if len(increments)==i:
                increments.append(inc[i])
            else:
                increments[i] += inc[i]
    return increments







# ============================================================================
#    Relative data (different trajectories)





def get_relative_samples(data):
    '''
    Returns a list of relative samples, containing 
    relative positions and relative velocities.
    The first index is a tuple of the trajectory
    ids, indexes 1-3 are selative positions, 4-6
    are relative velocitied, and the last is time.
    '''
    # group the data according to frames
    by_time = [np.array(g) for k,g in data.groupby(10)]
    relative_sample = []
    
    for lst in by_time:
        for i in range(len(lst)):
            for j in range(i+1, len(lst)):
                id_ = (lst[i][0], lst[j][0])
                dr = list(lst[i][1:4] - lst[j][1:4])
                dv = list(lst[i][4:7] - lst[j][4:7])
                new_sample = [id_] + dr + dv + [lst[i][-1]]
                relative_sample.append(new_sample)
    
    return relative_sample





def get_binned_relative_velocity_samples(data, r0, rn, n):
    '''
    Returns lists with samples of the relative velocity, binned
    according to their distance.
    
    r0 and rn are the distance limits, 
    n is the number of binnes within the range,
    and data is a DataFrame of the data samples.
    '''
    rel_samps = pd.DataFrame(get_relative_samples(data))
    print(len(rel_samps))
    rel_samps['d'] = np.sum(rel_samps[[1,2,3]]**2, axis=1)**0.5
    rel_samps['dvr'] = np.sum(np.array(rel_samps[[1,2,3]]) * np.array(rel_samps[[4,5,6]]), axis=1) / np.array(rel_samps['d'])
    rel_samps['bin'] = ((np.array(rel_samps['d']) - r0)/(rn-r0)*n).astype(int)
    dv_bins = [list(g['dvr']) for k,g in rel_samps.groupby('bin') if k<n]
    return dv_bins





def get_pairs(traj_list):
    '''
    Returns a list of trajectories of relative position
    and relative velocities at commong time instances.
    
    Input - 
    
    traj_list - a list of arrays that represent single
                particle trajectories. 
    
    Return -
    
    pairs - a list of arrays (n,12). In each array, the 
            indexes 2,3,4 are hte relative position, the
            5,6,7 are relative velocity, the 8,9,10 are
            relative accelerations and 11 is the frame 
            number. Inedx 0 is the id number of particle
            i and 1 is the id of particle j.
    '''
    from numpy import intersect1d, where, hstack, ones
    
    def in_list(arr, lst):
        return [arr[i] in lst for i in range(len(arr))]
        
    pairs = []
    for i in range(len(traj_list)):
        
        for j in range(i+1, len(traj_list)):
            
            common_times = intersect1d(traj_list[i][:,-1], traj_list[j][:,-1])
            
            if len(common_times)>0:
                
                ind_i = where(in_list(traj_list[i][:,-1], common_times))[0]
                ind_j = where(in_list(traj_list[j][:,-1], common_times))[0]
                p_ij = traj_list[i][ind_i, 1:10] - traj_list[j][ind_j, 1:10] 
                
                p_ij = hstack([traj_list[i][ind_i, :1], 
                               traj_list[j][ind_j, :1],
                               p_ij, traj_list[i][ind_i, -1:]])
                pairs.append(p_ij)
            
    return pairs




