#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 18:02:07 2018

@author: ron


Contains classes for tracking particles to form trajectories.

"""

from numpy import loadtxt, array


class tracker_four_frames(object):
    '''Implementation of the four-frame 3D particle tracking algorithm'''
    
    
    def __init__(self, fname, mean_flow = 0.0, d_max=1e10):
        '''
        fname - string, path of the particles containing file to which tracking
                should be performed.
                
        mean_flow - a numpy array of the mean flow vector, in units of the 
        calibrations spatial units per frame (e.g. mm per frame). The mean 
        flow is assumed not to change in space and time.
        
        d_max - maximum allowable translation between two frames for the 
                nearest neighbour search, after subtracting the mean flow. 
        '''
        self.fname = fname
        self.U = mean_flow
        self.d_max = d_max
        
        self.particles = {}
        
        data = loadtxt(self.fname)
        self.times = list(set(data[:,-1]))
        
        for tm in self.times:
            self.particles[tm] = []
            p_ = data[data[:,-1]==tm]
            for i in range(p_.shape[0]):
                p = array([-1] + list(p_[i,[0,1,2,-1]]))
                self.particles[tm].append(p)
        
        self.traj_ids = []
        self.traj_lengths = {}
    
    
    def nearest_neighbour(self, frame_num):
        '''For a given frame number, this will attemp to form nearest 
        neighbour trajectories for all the unlinked particles.  
        '''
        p1_lst = self.particles[frame_num]
        p2_lst = self.particles[frame_num+1]
        
        for i in range(len(p1_lst)):
            p1 = p1_lst[i]
            nn_val = self.find_nearest_neighbour(p1, frame_num+1)
            if nn_val[0] < self.d_max:
                p2 = p2_lst[nn_val[1]]
                if i == self.find_nearest_neighbour(p2, frame_num)[1]:
                    
                    if len(self.traj_ids)==0:
                        self.traj_ids.append(0)
                    else:
                        self.traj_ids.append( self.traj_ids[-1] + 1 )
                    
                    id_ = self.traj_ids[-1]
                    self.particles[frame_num][i][0] = id_
                    self.particles[frame_num+1][nn_val[1]][0] = id_
                    self.traj_lengths[id_] = 2
                    
        return None
        
        
    def find_nearest_neighbour(self, particle, frame_num):
        '''For a given particle, this returns the index of its nearest 
        neighbour in the frame number given, adn the distance between them.'''
        dt_particles = (frame_num - particle[4])
        dX = self.U * dt_particles
        dist_particle = lambda p2 : sum((particle[1:4] - (p2[1:4]-dX))**2)**0.5
        values = []
        for i in range(len(self.particles[frame_num])):
            values.append(( dist_particle(self.particles[frame_num][i]), i))
        min_val = min(values, key=lambda x: x[0])
        return min_val








class tracker_nearest_neighbour(object):
    '''A nearest-neighbour 3D particle tracking algorithm'''
    
    
    def __init__(self, fname, mean_flow = 0.0, d_max=1e10):
        '''
        fname - string, path of the particles containing file to which tracking
                should be performed.
                
        mean_flow - a numpy array of the mean flow vector, in units of the 
        calibrations spatial units per frame (e.g. mm per frame). The mean 
        flow is assumed not to change in space and time.
        
        d_max - maximum allowable translation between two frames for the 
                nearest neighbour search, after subtracting the mean flow. 
        '''
        self.fname = fname
        self.U = mean_flow
        self.d_max = d_max
        
        self.particles = {}
        
        data = loadtxt(self.fname)
        self.times = list(set(data[:,-1]))
        
        for tm in self.times:
            self.particles[tm] = []
            p_ = data[data[:,-1]==tm]
            for i in range(p_.shape[0]):
                p = array([-1] + list(p_[i,[0,1,2,-1]]))
                self.particles[tm].append(p)
        
        self.traj_ids = []
        self.traj_lengths = {}
    
    
    
    def track_all_frames(self):
        '''Will perform nearest neighbour tracking over all frames in a loop'''
        for tm in self.times[:-1]:
            self.nearest_neighbour_one_frame(tm)
        print('found %d trajectories'%(len(self.traj_ids)))
            
    
    def return_connected_particles(self):
        '''Will return the list of connected particles. To be used after 
        tracking is complete.'''
        p_list = []
        for tm in self.times:
            for p in self.particles[tm]:
                p_list.append(p)
        return p_list
    
    
    def nearest_neighbour_one_frame(self, frame_num):
        '''For a given frame number, this will attemp to form nearest 
        neighbour trajectories for all the unlinked particles.  
        '''
        p1_lst = self.particles[frame_num]
        p2_lst = self.particles[frame_num+1]
        
        for i in range(len(p1_lst)):
            p1 = p1_lst[i]
            nn_val = self.find_nearest_neighbour(p1, frame_num+1)
            if nn_val[0] < self.d_max:
                p2 = p2_lst[nn_val[1]]
                if i == self.find_nearest_neighbour(p2, frame_num)[1]:
                    
                    if p1[0]==-1:
                        
                        if len(self.traj_ids)==0:
                            self.traj_ids.append(0)
                        else:
                            self.traj_ids.append( self.traj_ids[-1] + 1 )
                        
                        id_ = self.traj_ids[-1]
                        self.particles[frame_num][i][0] = id_
                        self.particles[frame_num+1][nn_val[1]][0] = id_
                        self.traj_lengths[id_] = 2
                    
                    else:
                        id_ = p1[0]
                        self.particles[frame_num+1][nn_val[1]][0] = id_
                        self.traj_lengths[id_] += 1
        return None
        
        
    def find_nearest_neighbour(self, particle, frame_num):
        '''For a given particle, this returns the index of its nearest 
        neighbour in the frame number given, adn the distance between them.'''
        dt_particles = (frame_num - particle[4])
        dX = self.U * dt_particles
        dist_particle = lambda p2 : sum((particle[1:4] - (p2[1:4]-dX))**2)**0.5
        values = []
        for i in range(len(self.particles[frame_num])):
            values.append(( dist_particle(self.particles[frame_num][i]), i))
        min_val = min(values, key=lambda x: x[0])
        return min_val
