# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
Created on Thu Jun 26 17:53:50 2025

@author: ron


A script used to obtain relative trajectories, namely, trajectories
relative to one another.

"""

import numpy as np
import pandas as pd
from tqdm import tqdm



def load_trajs_as_arrays(fname):
    '''
    given a path, fname, this returns a list of arrays, each corresponding
    to a particlular trajectory.
    '''
    from pandas import read_csv
    from numpy import array
    data = read_csv(fname, sep='\t', header=None)
    trajs = [array(sorted(array(g), key=lambda x: x[-1])) 
             for k,g in tqdm(data.groupby(0), desc='reading data') if k!=-1]
    
    time_col = data.shape[1]-1
    Nframes = len(set(data[time_col]))
    
    Nsamples = sum([len(tr) for tr in trajs])
    prnt = (Nsamples/Nframes, Nframes)
    print('\n', 'loaded trajectories, %.1f sample per frame, %d frames\n'%prnt)
    
    return trajs




def get_pairs(traj_list, Np=100):
    '''
    Returns a list of trajectories of relative position
    and relative velocities at commong time instances.
    
    Input - 
    
    traj_list - a list of arrays that represent single
                particle trajectories. 
    
    Np - Integer; the alogrithm divides the trajectories into groups of this
         many trajecotries to make the computation faster, although it means
         not all trajectory pairs will be obtained. 
    
    Return -
    
    pairs - a list of arrays (n,12). In each array, the 
            indexes 2,3,4 are hte relative position, the
            5,6,7 are relative velocity, the 8,9,10 are
            relative accelerations and 11 is the frame 
            number. Inedx 0 is the id number of particle
            i and 1 is the id of particle j.
    '''
    from numpy import intersect1d, where, hstack
    
    def in_list(arr, lst):
        return [arr[i] in lst for i in range(len(arr))]
        
    pairs = []
    
    count = 1
    for k in range(len(traj_list)//Np+1):
        group = traj_list[k*Np:(k+1)*Np]
        
        desc = 'pairing group %d/%d'%(k+1, len(traj_list)//Np+1)
        for i in tqdm(range(len(group)), desc=desc):
            
            for j in range(i+1, len(group)):
                
                common_times = intersect1d(group[i][:,-1], group[j][:,-1])
                
                if len(common_times)>0:
                    
                    ind_i = where(in_list(group[i][:,-1], common_times))[0]
                    ind_j = where(in_list(group[j][:,-1], common_times))[0]
                    p_ij = group[i][ind_i, 1:10] - group[j][ind_j, 1:10] 
                    
                    idcol = np.ones((len(p_ij),1))*count
                    
                    p_ij = hstack([idcol,
                                   group[i][ind_i, :1], 
                                   group[j][ind_j, :1],
                                   p_ij, group[i][ind_i, -1:]])
                    pairs += p_ij.tolist()
                    count+=1
    print('\n', 'finighed pairing: %d pairs\n'%count)
    return pairs








if __name__ == "__main__":
    
    fname = '/home/ron/Desktop/Research/PTV cases/Mani_Benny_Vortex/smoothed_traj_Run26'
    
    # read trajectories sorted by their initial frame number 
    trajs = sorted(load_trajs_as_arrays(fname), key=lambda tr: tr[0,-1])
    
    # pair trajectories in groups
    pairs = get_pairs(trajs[:1000], Np = 50)
        
    # # save the pairs
    saveName = 'pairs'
    fmt = ['%d', '%d', '%d'] + ['%.03f' for i in range(len(pairs[0])-3)]
    np.savetxt(saveName, pairs, delimiter='\t', fmt=fmt)








