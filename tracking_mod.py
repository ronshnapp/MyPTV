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
    
    
    
    def track_all_frames(self):
        '''Will perform tracking over all frames in a loop'''
        for tm in self.times[:-1]:
            self.track_single_frame(tm)
        print('found %d trajectories'%(len(self.traj_ids)))
        N_links = 0
        for k in self.traj_lengths.keys(): N_links += self.traj_lengths[k]
        print('linked %d particles'%(N_links))
            
    
    
    def track_single_frame(self, frame_num):
        '''For a given frame number, this will attempt to link particles and 
        foרm trajectories with particles in the next frame. When possible it 
        will use the 3 frame tracking, and if not it will attempt a nearest 
        neighbour tracking.  
        '''
        p1_lst = self.particles[frame_num]
        p2_lst = self.particles[frame_num+1]
        
        
        # try 3 frame tracking for particles that are connected at least once 
        for i in range(len(p1_lst)):
            p1 = p1_lst[i]
            if p1[0] == -1: continue
            match_val = self.find_velocity_projected_match(p1, frame_num+1)
            # ========================================================
            condition = match_val[0]<1.0               # <--add a test here!
            # ========================================================
            if condition:
                p2 = p2_lst[match_val[1]]
                if p2[0] != -1:
                    id_ = p1[0]
                    self.particles[frame_num+1][match_val[0]][0] = id_
                    self.traj_lengths[id_] += 1
                    
        
        # nearest neighbour tracking on particles that have not yet been
        # connected with any other particle:
        for i in range(len(p1_lst)):
            p1 = p1_lst[i]
            if p1[0] != -1: continue
            nn_val = self.find_nearest_neighbour(p1, frame_num+1)
            if nn_val[0] < self.d_max:
                p2 = p2_lst[nn_val[1]]
                if p2[0] != -1: continue
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
    
    
    def get_particle_by_id(self, id_, frame_num):
        '''Returns the particle with the given id=id_ at the given 
        frame number. If there is no such particle, it returns None.'''
        for p in self.particles[frame_num]:
            if p[0] == id_:
                return p
        return None
        
    
    
    def find_velocity_projected_match(self, particle, frame_num):
        '''
        This looks for the best match according to the three-frame tracking 
        heuristic.
        Namely, given a particle that is linked to at least one more particle, 
        and a given frame number, this will project the expected position 
        based on an estimated constant velocity and search for the nearest
        neighbour to this projection. It returns the index of this nearest 
        neighbour in frame_num, and the distance between it and the 
        projection.'''
        
        id_ = particle[0]
        p_im1 = self.get_particle_by_id(id_, particle[-1]-1)
        
        dt = frame_num - particle[-1]
        v = particle[1:4] - p_im1[1:4]
        x_proj = particle[1:4] + dt * v
        
        dist = lambda p: sum((x_proj - p[1:4])**2)**0.5
        values = []
        for i in range(len(self.particles[frame_num])):
            values.append(( dist(self.particles[frame_num][i]), i))
        min_val = min(values, key=lambda x: x[0])
        return min_val


    def return_connected_particles(self):
        '''Will return the list of connected particles. To be used after 
        tracking is complete.'''
        p_list = []
        for tm in self.times:
            for p in self.particles[tm]:
                p_list.append(p)
        return p_list








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
        N_links = 0
        for k in self.traj_lengths.keys(): N_links += self.traj_lengths[k]
        print('linked %d particles'%(N_links))
            
    
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
                if p2[0] != -1: continue
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
