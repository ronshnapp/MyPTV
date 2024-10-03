# -*- coding: utf-8 -*-
#!/usr/bin/env python3

"""
Created on Fri May 31 15:01:48 2024

@author: ron
"""

from pandas import read_csv
from numpy import ptp, array, arange, amin, amax
import matplotlib.pyplot as plt



def plot_trajectories(fname, min_length, write_trajID=False, t0=0, te=-1):
    '''
    This function plots trajectories from a given file in 3D.
    
    inputs:
        
    fname - the path to the file that contains the trajectories; the file can 
            be either in trajectories format or in smoothed trajectories format
    
    min_lenth - only trajectories that have more samples than this number will
                be plotted
                
    write_trajID - If True this will desplay the trajectory ID on top of them
    
    t0 and te - used to delineate the time range for which we plot the data. 
                we only plot the samples in the time range starting at frame
                t0 and ending at frame te. Set t0=0 and te=-1 (default) to plot
                trajectories at all times available.
    '''
    
    data = read_csv(fname, header=None, sep='\t')
    trajectories = dict([(g, array(k.values)) 
                         for g,k in data.groupby(0) if g!=-1])
    
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    xm = []
    ym = []
    zm = []
    
    if te==-1:
        te = max(data[data.shape[1]-1])
    
    trajIDs = list(trajectories.keys())
    
    count = 0
    for id_ in trajIDs:
        if len(trajectories[id_][:,1])<min_length: continue
        time = trajectories[id_][:,-1]
        inds = arange(len(trajectories[id_]))
        
        if time[0]>=te or time[-1]<=t0: 
            continue
    
        if time[0]>=t0: 
            i0 = 0
        else:
            i0 = inds[time==t0][0]
        
        if time[-1]<=te: 
            ie = -1
        else:
            ie = inds[time==te][0]    
            
        
        xs = trajectories[id_][i0:ie,1]
        ys = trajectories[id_][i0:ie,2]
        zs = trajectories[id_][i0:ie,3]
        l = ax.plot(xs, zs, ys, 'o-', ms=1, lw=0.5)
        
        xm.append(amin(xs)) ; xm.append(amax(xs))
        ym.append(amin(ys)) ; ym.append(amax(ys))
        zm.append(amin(zs)) ; zm.append(amax(zs))
        
        if write_trajID==True:
            color = l[0].get_color()
            ax.text(xs[0], zs[0], ys[0], str(id_),
                    fontdict={'fontsize': 12, 'color':color})
        
        count += 1
    
    ax.set_box_aspect((ptp(xm), ptp(zm), ptp(ym)))
    
    ax.set_xlabel('x')
    ax.set_zlabel('y')
    ax.set_ylabel('z')
    
    print('plotted %d trajectories'%(count))
    
    plt.show()







def getSamplesFromLongTrajectories(fname, min_len):
    '''
    Reads a trajectory file and returns an array with its samples that
    belong to "long" trajectories, whose length is >= than min_len.
    '''
    data = read_csv(fname, header=None, sep='\t')
    trajectories = dict([(g, array(k.values)) 
                         for g,k in data.groupby(0) if g!=-1])
    
    to_take = []
    for k in trajectories.keys():
        tr = trajectories[k]
        if len(tr)>=min_len:
            for i in range(len(tr)):
                to_take.append(tr[i])
            
    return array(to_take)




    
    
    

def PlotParticlePositionHistogram(fname):
    '''
    This function plots trajectories from a given file in 3D.
    '''
    data = read_csv(fname, header=None, sep='\t')
    
    fig, ax = plt.subplots(1,3)
    
    xm = list(data[0])
    ym = list(data[1])
    zm = list(data[2])
    
    ax[0].hist(xm, bins='auto')
    ax[1].hist(ym, bins='auto')
    ax[2].hist(zm, bins='auto')
    
    return None

