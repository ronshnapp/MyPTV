# -*- coding: utf-8 -*-
#!/usr/bin/env python3

"""
Created on Fri May 31 15:01:48 2024

@author: ron
"""

from pandas import read_csv
from numpy import ptp, array, arange, amin, amax
import matplotlib.pyplot as plt

from moviepy.video.io.bindings import mplfig_to_npimage
import moviepy.editor as mpy






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
    
    xmax = amax(data[1])
    xmin = amin(data[1])
    ymax = amax(data[2])
    ymin = amin(data[2])
    zmax = amax(data[3])
    zmin = amin(data[3])
    
    
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
        c = (1-(xs[0]-xmin)/(xmax-xmin)*0.97, 
             (ys[0]-ymin)/(ymax-ymin)*0.97, 
             (zs[0]-zmin)/(zmax-zmin)*0.97)
        l = ax.plot(xs, zs, ys, 'o-', ms=1, lw=0.5, color=c)
        
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









class animate_trajectories(object):
    
    def __init__(self, fname, min_length, f0=None, fe=None, fps=25,
                 tail_length=4, view_angles = (15,70), rotation_rate=0.1):
        
        
        
        data = read_csv(fname, header=None, sep='\t')
        
        self.trajectories = dict([(g, array(k.values)) 
                             for g,k in data.groupby(0) if g!=-1])
        
        self.longs = [k for k in self.trajectories.keys() 
                      if len(self.trajectories[k])>=min_length]
        
        x_lst, y_lst, z_lst = [], [], []
        for i in self.longs:
            x_lst += list(self.trajectories[i][:,1])
            y_lst += list(self.trajectories[i][:,2])
            z_lst += list(self.trajectories[i][:,3])
        
        self.xmax = amax(x_lst) ; self.xmin = amin(x_lst)
        self.ymax = amax(y_lst) ; self.ymin = amin(y_lst)
        self.zmax = amax(z_lst) ; self.zmin = amin(z_lst)
        
        if f0 is None:
            f0 = int(min(data[data.columns[-1]]))
            
        if fe is None:
            fe = int(max(data[data.columns[-1]]))
        
        self.fps = fps
        self.counter = 0
        self.frames = list(range(f0, fe+1))
        self.duration = (len(self.frames)-1)/self.fps
        self.tl = tail_length
        self.min_length = min_length
        self.angles = view_angles
        self.rotation = rotation_rate
        
        
        
    def update(self, frame):
        frame = self.frames[self.counter]
        cmap = plt.get_cmap('viridis')
        self.ax.clear()
        
        for k in self.longs:
            tr = self.trajectories[k]
            whr = arange(len(tr))[tr[:,-1]==frame]
            if any(whr):
                ind = whr[0]
                x = tr[ind-self.tl:ind+1,1]
                y = tr[ind-self.tl:ind+1,2]
                z = tr[ind-self.tl:ind+1,3]
                if len(x)==0: continue
                dx = ((x[-1]-x[0])**2 + (y[-1]-y[0])**2 + (z[-1]-z[0])**2)**0.5
                c = min([(dx/self.tl) / self.vscale,1])
                self.ax.plot(x, z, y, '-', color = cmap(c*0.9)) #color=(0.1+c*0.9,0,0.8*(1-c)))
        
        self.ax.set_xlim(self.xmin, self.xmax)
        self.ax.set_zlim(self.ymin, self.ymax)
        self.ax.set_ylim(self.zmin, self.zmax)
        
        self.ax.set_xlabel('x')
        self.ax.set_ylabel('z')
        self.ax.set_zlabel('y')

        
        self.ax.w_xaxis.set_pane_color((0.6,0.6,1,0.1))
        self.ax.w_yaxis.set_pane_color((0.6,0.6,1,0.15))
        self.ax.w_zaxis.set_pane_color((0.7,0.6,1,0.2))
        
        self.ax.grid(False)
        
        self.ax.view_init(elev=self.angles[0], 
                          azim=self.angles[1] + self.rotation*self.counter)
        
        
        self.ax.set_box_aspect((self.xmax-self.xmin, 
                                self.zmax-self.zmin, 
                                self.ymax-self.ymin))
        
        plt.tight_layout(0.5)
        
        self.counter += 1
        
        return mplfig_to_npimage(self.fig) # RGB image of the figure

        

    def animate(self):
        '''
        will animate the particle's location, and save the animation
        '''
        #self.prepare_for_animation()
        
        self.vscale = 0
        for i in self.longs:
            tr = self.trajectories[i]
            dt = int(self.min_length/2)
            dx = sum([(tr[dt,j] - tr[0,j])**2 for j in [1,2,3]])**0.5
            if self.vscale<dx/dt:
                self.vscale = dx/dt
        
        self.fig = plt.figure(figsize=(9,9))
        self.ax = self.fig.add_subplot(projection='3d')
        
        animation = mpy.VideoClip(self.update, duration=self.duration)
        animation.write_videofile('animation.mp4',fps = self.fps)
        return animation
    








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

