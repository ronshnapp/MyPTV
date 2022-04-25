# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 18:02:07 2018

@author: ron


Contains classes for tracking particles to form trajectories.

"""

from numpy import loadtxt, array, savetxt
from scipy.spatial import KDTree





class tracker_four_frames(object):
    '''Implementation of a 4-frame 3D particle tracking algorithm using
    the so-called best estimate method in 
    https://doi.org/10.1007/s00348-005-0068-7.
    '''
    
    def __init__(self, fname, mean_flow = 0.0, d_max=1e10, dv_max=1e10):
        '''
        fname - string, path of the particles containing file to which tracking
                should be performed.
                
        mean_flow - a numpy array of the mean flow vector, in units of the 
        calibrations spatial units per frame (e.g. mm per frame). The mean 
        flow is assumed not to change in space and time.
        
        d_max - maximum allowable translation between two frames for the 
                nearest neighbour search, after subtracting the mean flow. 
                
        dv_max - maximum allowable change in velocity for the two-frame 
                 velocity projection search. The radius around the projection
                 is therefore dv_max/dt (where dt = 1 frame^{-1})
        '''
        self.fname = fname
        self.U = mean_flow
        self.d_max = d_max
        self.dv_max = dv_max
        
        self.particles = {}
        
        data = loadtxt(self.fname)
        self.times = list(set(data[:,-1]))
        
        for tm in self.times:
            self.particles[tm] = []
            p_ = data[data[:,-1]==tm]
            for i in range(p_.shape[0]):
                p = array([-1] + list(p_[i,[0,1,2,-1]]))
                self.particles[tm].append(p)
        
        for k in self.particles.keys():
            self.particles[k] = array(self.particles[k])
        
        self.trees = {}
        
        self.traj_ids = []
        self.traj_lengths = {}
        self.N_four_frames = 0
        self.N_nearest_neighbour = 0
    
    
    def track_all_frames(self, frames=None):
        '''
        Will perform tracking over a range of frames in a loop.
        
        input -
        frames - if None (default), will track particles over all the available
                 frames. Else this may be a list of intergers that must be
                 sorted and increasing in increments of one; then, these are
                 the frame numbers used in the tracking.
        '''
        
        if frames == None:
            frames = self.times[:-1]
        else:
            for i in range(len(frames)-1):
                if frames[i+1]-frames[i] != 1 or type(frames[i])!= int:
                    raise ValueError('frame range does not follow the rules.')
        
        for tm in frames:
            print('', end='\r')
            print(' frame: %d'%tm, end='\r')            
            self.track_single_frame(tm)
            
        N_links = 0
        for k in self.traj_lengths.keys(): N_links += self.traj_lengths[k]
        N_p = 0
        for k in self.particles.keys(): N_p += len(self.particles.keys()) 
        NT = len(self.traj_ids)
        print('found %d trajectories (avg. length = %.1f)'%(NT, N_links/NT))
        print('four frame links: %d'%self.N_four_frames)
        print('nearest neighbour links: %d'%self.N_nearest_neighbour)
    
    
    def track_single_frame(self, frame_num):
        '''For a given frame number, this will attempt to link particles and 
        form trajectories with particles in the next frame. When possible it 
        will use the 4 frame tracking. If a particle has no previous links, 
        it will attempt a nearest neighbour tracking.  
        '''
        p1_lst = self.particles[frame_num]
        p2_lst = self.particles[frame_num+1]
        
        
        # try 4 frame tracking for particles that are connected at least once 
        for i in range(len(p1_lst)):
            p1 = p1_lst[i]
            if p1[0] == -1: continue
            best_estimate = self.find_best_estimate_link(p1)
            if best_estimate is None: continue
            # ===================================================
            condition = True             # <-- Make up best estimate threshold
            # ===================================================
            if condition:
                p2 = p2_lst[best_estimate[1]]
                if p2[0] == -1:
                    id_ = p1[0]
                    self.particles[frame_num+1][best_estimate[1]][0] = id_
                    self.traj_lengths[id_] += 1
                    self.N_four_frames += 1
        
        
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
                    self.N_nearest_neighbour += 1
                    
        return None
        
        
    def find_nearest_neighbour(self, particle, frame_num):
        '''For a given particle, this returns the index of its nearest 
        neighbour in the frame number given, adn the distance between them.'''
        dt_particles = (frame_num - particle[4])
        dX = self.U * dt_particles
        p = particle[1:4] + dX
        
        try:
            tree = self.trees[frame_num]
        except:
            tree = KDTree(self.particles[frame_num][:,1:4])
            self.trees[frame_num] = tree
        
        return tree.query(p, k=1)
        

    
    
    def get_particle_by_id(self, id_, frame_num):
        '''Returns the particle with the given id=id_ at the given 
        frame number. If there is no such particle, it returns None.'''
        for p in self.particles[frame_num]:
            if p[0] == id_:
                return p
        return None

        
    
    def find_best_estimate_link(self, particle):
        '''given a particle, this will return a particle from the next frame
        which fulfills the best estimate tracking heuristic condition.
        
        Specifically, the particle is projected assuming constant velocity.
        Then, neighbour of this projection are taken as candidates. To choose
        between the candidates, each is projected again into the n+2 frame, 
        and the candidate whose projection has a nearest neighbour is chosen
        as the "best estimate", and thus returned.
        
        the first projection is determined as:
            
        x_i+1 = x_i + v_i*dt
        with v_i = (x_i - x_i-1)/dt
        
        the second projection is determined as:
        
        x_i+2 = x_i + 2dt*v_i + (2dt)^2 * a_i
        with a_i = (x_i+1 - 2x_i + x_i-1)/(dt)^2
            
        '''
        frame_num = particle[-1] + 1
        
        # 1: find the projection at n+1
        id_ = particle[0]
        p_im1 = self.get_particle_by_id(id_, particle[-1]-1)
        
        dt = frame_num - particle[-1]
        v = (particle[1:4] - p_im1[1:4])/dt
        x_proj = particle[1:4] + dt * v
        
        # 2: find neighbours of the projection at n+1 using KDTree:
        try:
            tree = self.trees[frame_num]
        except:
            tree = KDTree(self.particles[frame_num][:,1:4])
            self.trees[frame_num] = tree
        dist = lambda p: sum((x_proj - p[1:4])**2)**0.5 
        proj_neighbours = [(dist(self.particles[frame_num][i]), i) for i in 
                           tree.query_ball_point(x_proj, self.dv_max, )]

        
        
        # 2.1: if there are no projection neighbours, return None
        if len(proj_neighbours)==0: return None
        # 2.2: if there's only one projection neighbour, return it
        if len(proj_neighbours)==1: return proj_neighbours[0]
        # 2.3: if there is no n+2 frame, return the projection's 
        #      nearest neighbour
        if frame_num+1 > self.times[-1]: 
            return min(proj_neighbours, key=lambda x: x[0])
        
        
        # 3: for each n+1 neighbour, project to frame n+2 and 
        #    locate this projection's nearest neighbour and write it
        
        try:
            tree = self.trees[frame_num+1]
        except:
            tree = KDTree(self.particles[frame_num+1][:,1:4])
        
        for j in range(len(proj_neighbours)):
            d,i = proj_neighbours[j]
            xip1 = self.particles[frame_num][i][1:4]
            xi = particle[1:4]
            xim1 = p_im1[1:4]
            vi = (xip1 - xim1)/(2*dt)
            ai = (xip1 - 2*xi + xim1)/(dt**2)
            xip2 = xi + 2*dt*vi + 2*dt**2*ai     # <-- projection at frame n+2
            
            nnd, dump = tree.query(xip2, k=1)
            proj_neighbours[j] = (nnd, proj_neighbours[j][1])
        
        # 4: return the candidate with the smallest nearest neighbour distance
        return min(proj_neighbours, key=lambda x: x[0])
            

    def return_connected_particles(self):
        '''Will return the list of connected particles. To be used after 
        tracking is complete.'''
        p_list = []
        for tm in self.times:
            for p in self.particles[tm]:
                p_list.append(p)
        return p_list


    def save_results(self, fname):
        '''
        Will save the results after tracking is done.
        '''
        fmt = ['%d', '%.3f', '%.3f', '%.3f', '%.3f']
        savetxt(fname ,self.return_connected_particles(),
                delimiter='\t', fmt=fmt)










class tracker_two_frames(object):
    '''Implementation of a two-frame 3D particle tracking algorithm using
    projection of the particles assuming constant velocity.'''
    
    
    def __init__(self, fname, mean_flow = 0.0, d_max=1e10, dv_max=1e10):
        '''
        fname - string, path of the particles containing file to which tracking
                should be performed.
                
        mean_flow - a numpy array of the mean flow vector, in units of the 
        calibrations spatial units per frame (e.g. mm per frame). The mean 
        flow is assumed not to change in space and time.
        
        d_max - maximum allowable translation between two frames for the 
                nearest neighbour search, after subtracting the mean flow. 
                
        dv_max - maximum allowable change in velocity for the two-frame 
                 velocity projection search. The radius around the projection
                 is therefore dv_max/dt (where dt = 1 frame^{-1})
        '''
        self.fname = fname
        self.U = mean_flow
        self.d_max = d_max
        self.dv_max = dv_max
        
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
        self.N_three_frames = 0
        self.N_nearest_neighbour = 0
    
    
    def track_all_frames(self):
        '''Will perform tracking over all frames in a loop'''
        for tm in self.times[:-1]:
            self.track_single_frame(tm)
        N_links = 0
        for k in self.traj_lengths.keys(): N_links += self.traj_lengths[k]
        N_p = 0
        for k in self.particles.keys(): N_p += len(self.particles.keys()) 
        print('linked %d out of %d particles '%(N_links, N_p))
        NT = len(self.traj_ids)
        print('found %d trajectories (avg. length = %.1f)'%(NT, N_links/NT))
        print('three frame links: %d'%self.N_three_frames)
        print('nearest neighbour links: %d'%self.N_nearest_neighbour)
    
    
    def track_single_frame(self, frame_num):
        '''For a given frame number, this will attempt to link particles and 
        form trajectories with particles in the next frame. When possible it 
        will use the 3 frame tracking, and if not it will attempt a nearest 
        neighbour tracking.  
        '''
        p1_lst = self.particles[frame_num]
        p2_lst = self.particles[frame_num+1]
        
        
        # try 3 frame tracking for particles that are connected at least once 
        D = self.dv_max / 1.0
        for i in range(len(p1_lst)):
            p1 = p1_lst[i]
            if p1[0] == -1: continue
            match_val = self.find_velocity_projected_match(p1, frame_num+1)
            condition = match_val[0] < D  
            if condition:
                p2 = p2_lst[match_val[1]]
                if p2[0] == -1:
                    id_ = p1[0]
                    self.particles[frame_num+1][match_val[1]][0] = id_
                    self.traj_lengths[id_] += 1
                    self.N_three_frames += 1
                    
        
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
                    self.N_nearest_neighbour += 1
                    
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
    
    
    def save_results(self, fname):
        '''
        Will save the results after tracking is done.
        '''
        fmt = ['%d', '%.3f', '%.3f', '%.3f', '%.3f']
        savetxt(fname ,self.return_connected_particles(),
                delimiter='\t', fmt=fmt)







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
    
        for k in self.particles.keys():
            self.particles[k] = array(self.particles[k])
            
        self.trees = {}
        
        self.traj_ids = []
        self.traj_lengths = {}
        
    
    
    def track_all_frames(self, frames=None):
        '''Will perform nearest neighbour tracking over all frames in a loop'''
        
        if frames == None:
            frames = self.times[:-1]
        else:
            for i in range(len(frames)-1):
                if frames[i+1]-frames[i] != 1 or type(frames[i])!= int:
                    raise ValueError('frame range does not follow the rules.')
        
        for tm in frames:
            self.nearest_neighbour_one_frame(tm)
        N_links = 0
        for k in self.traj_lengths.keys(): N_links += self.traj_lengths[k]
        N_p = 0
        for k in self.particles.keys(): N_p += len(self.particles.keys()) 
        print('linked %d out of %d particles '%(N_links, N_p))
        NT = len(self.traj_ids)
        print('found %d trajectories (avg. length = %.1f)'%(NT, N_links/NT))
            
    
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
        p = particle[1:4] + dX
        
        try:
            tree = self.trees[frame_num]
        except:
            tree = KDTree(self.particles[frame_num][:,1:4])
            self.trees[frame_num] = tree
        
        return tree.query(p, k=1)
        
        # dist_particle = lambda p2 : sum((particle[1:4] - (p2[1:4]-dX))**2)**0.5
        # values = []
        # for i in range(len(self.particles[frame_num])):
        #     values.append(( dist_particle(self.particles[frame_num][i]), i))
        # min_val = min(values, key=lambda x: x[0])
        # return min_val
    
    
    def save_results(self, fname):
        '''
        Will save the results after tracking is done.
        '''
        fmt = ['%d', '%.3f', '%.3f', '%.3f', '%.3f']
        savetxt(fname ,self.return_connected_particles(),
                delimiter='\t', fmt=fmt)
