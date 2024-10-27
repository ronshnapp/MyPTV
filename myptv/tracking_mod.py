# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 18:02:07 2018

@author: ron


Contains classes for tracking particles to form trajectories.

"""

from numpy import loadtxt, array, savetxt, hstack, ones, where, mean, zeros
from numpy import poly1d, polyfit
from numpy import sum as npsum
from scipy.spatial import KDTree
from pandas import read_csv
import tqdm




class tracker_four_frames(object):
    '''Implementation of a 4-frame 3D particle tracking algorithm using
    the so-called best estimate method in 
    https://doi.org/10.1007/s00348-005-0068-7.
    '''
    
    def __init__(self, fname, mean_flow = 0.0, d_max=1e10, dv_max=1e10,
                 store_candidates=False):
        '''
        fname - string, path of the particles containing file to which tracking
                should be performed.
                
        mean_flow - a numpy array of the mean flow vector, in units of the 
                    calibrations spatial units per frame (e.g. mm per frame). 
                    The mean flow is assumed not to change in space and time.
        
        d_max - maximum allowable translation between two frames for the 
                nearest neighbour search, after subtracting the mean flow. 
                
        dv_max - maximum allowable change in velocity for the two-frame 
                 velocity projection search. The radius around the projection
                 is therefore dv_max/dt (where dt = 1 frame^{-1})
                 
        store_candidates - boolean indicator, defaul on False. If True, the 
                           tracker will store all the candidate links that were 
                           considered in the tracking process for future 
                           analysis.
        '''
        self.fname = fname
        self.U = mean_flow
        self.d_max = d_max
        self.dv_max = dv_max
        
        # particles are stored in a dictionary. Keys are frame numbers, values
        # are lists of arrays, each array representing a particle
            
        data = read_csv(self.fname, header=None, sep='\t')
        timeIndex = data.shape[1] - 1
        self.particles = dict([(g, hstack([ones((len(k),1))*-1, k.values])) 
                               for g,k in data.groupby(timeIndex)])
        self.times = sorted(list(self.particles.keys()))
        
        self.trees = {}
        
        self.traj_ids = []
        self.traj_lengths = {}
        self.N_four_frames = 0
        self.N_nearest_neighbour = 0
        
        # setting up a dictionary to store the data on all candidate links
        self.store_candidates = store_candidates
        if store_candidates:
            self.candidate_links = dict((tm, []) for tm in self.times)
    
    
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
        try:
            p1_lst = self.particles[frame_num]
            p2_lst = self.particles[frame_num+1]
        except:
            return None
        
        
        # try 4 frame tracking for particles that are connected at least once 
        for i in range(len(p1_lst)):
            p1 = p1_lst[i]
            if p1[0] == -1: continue
            best_estimate = self.find_best_estimate_link(p1)
            if best_estimate is None: continue
            # ===================================================
            cond = True             # <-- Make up best estimate threshold
            # ===================================================
            if cond:
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
                
                # ============================================================
                # optional - check also if p1 is the nearest neighbour of
                # p2 in p1's frame
                
                # cond = i == self.find_nearest_neighbour(p2, frame_num)[1]
                cond = True
                # ============================================================
                
                if cond:
                    
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
        neighbour in the frame number given, and the distance between them.'''
        dt_particles = (frame_num - particle[-1])
        dX = self.U * dt_particles
        p = particle[1:4] + dX
        
        try:
            tree = self.trees[frame_num]
        except:
            tree = KDTree(self.particles[frame_num][:,1:4])
            self.trees[frame_num] = tree
        
        # store all possible candidates if storing option is True
        if self.store_candidates:
            colls = (self.particles[particle[-1]] == particle).all(axis=1)
            particle_index = where(colls)[0][0]
            candidates = tree.query_ball_point(p, self.d_max)
            for cand in candidates:
                cl = [(particle_index, particle[-1]), (cand, frame_num)]
                self.candidate_links[particle[-1]].append(cl)
        
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
        
        
        # if storing option is True, store all possible candidates
        if self.store_candidates:
            colls = (self.particles[particle[-1]] == particle).all(axis=1)
            particle_index = where(colls)[0][0]
            candidates = [pn[1] for pn in proj_neighbours]
            for cand in candidates:
                cl = [(particle_index, particle[-1]), (cand, frame_num)]
                self.candidate_links[particle[-1]].append(cl)
        
        
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
        data_to_save = self.return_connected_particles()
        fmt = ['%d', '%.3f', '%.3f', '%.3f']
        for i in range(len(data_to_save[0])-6):
            fmt.append('%d')
        fmt += ['%.3f', '%.3f']
        savetxt(fname , data_to_save,
                delimiter='\t', fmt=fmt)
        
        
    def plot_candidate_graph(self):
        '''
        After tracking with store_candidates turned on, this function will plot
        a graph showing the network of all possible candidates.
        '''
        if self.store_candidates == False:
            raise ValueError('store_candidate must be True to plot candidate_graph')
        
        from matplotlib.pyplot import figure, plot, show, xlabel, ylabel, title
        fig = figure()
        fig.add_subplot(111)
        
        for tm in self.times:
            x = ones(len(self.particles[tm]))*tm
            y = list(range(len(self.particles[tm])))
            plot(x, y, 'ok', ms=2.5, alpha=0.5)
                
        for k in self.candidate_links.keys():
            for cl in self.candidate_links[k]:
                plot([cl[0][1], cl[1][1]], [cl[0][0], cl[1][0]], '-b',
                     lw=0.7)
                
        for tm in self.times:
            for i in range(len(self.particles[tm])):
                p = self.particles[tm][i]
                if p[0] != -1 and tm < self.times[-1]:
                    whr = where(self.particles[tm+1][:,0] == p[0])[0]
                    if len(whr)>0:
                        plot([tm,tm+1], [i,whr[0]], '-r')
                    
        xlabel('frame number')
        ylabel('particle index')
        title('tracking candidate graph: red = link, blue = regected candidate link')
        show()
        
        










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
        dt_particles = (frame_num - particle[-1])
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
                #p = array([-1] + list(p_[i,[0,1,2,-1]]))
                p = array([-1] + list(p_[i,:]))
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
        dt_particles = (frame_num - particle[-1])
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
        data_to_save = self.return_connected_particles()
        fmt = ['%d', '%.3f', '%.3f', '%.3f']
        for i in range(len(data_to_save[0])-6):
            fmt.append('%d')
        fmt += ['%.3f', '%.3f']
        savetxt(fname , data_to_save,
                delimiter='\t', fmt=fmt)













class tracker_multiframe(object):
    '''
    A muti-frame tracker that connects trajetories over more than one
    frame. This is meant to overcome issues of missed particle images based
    on ideas coming from the 4-frame algorithm best estimate method.
    '''
    
    def __init__(self, fname, max_dt, Ns, mean_flow = 0.0, d_max=1e10, 
                 dv_max=1e10, NSR_th=0.25):
        '''
        fname - string, path of the particles containing file to which tracking
                should be performed.
                
        max_dt - Ihe maximal number of frames over which a skip in tracking
                 is allowed. For example, if max_dt=2 then a trajectory can 
                 be constructed by linking a particle at frame i with another
                 particle at time i+2.
                 
        Ns - The number of frames used to calculate trajectory noise level.
             Should be an odd integer.
                
        mean_flow - a numpy array of the mean flow vector, in units of the 
                    calibrations spatial units per frame (e.g. mm per frame). 
                    The mean flow is assumed not to change in space and time.
        
        d_max - maximum allowable translation between two frames for the 
                nearest neighbour search, after subtracting the mean flow. 
                
        dv_max - maximum allowable change in velocity for the two-frame 
                 velocity projection search. The radius around the projection
                 is therefore dv_max/dt (where dt = 1 frame^{-1})
                 
        NSR_th - Maximum allowed noise level to store a trajectory. We only 
                 store trajectories whose NSR is lower than this number. Note
                 that NSR_th=1 corresponds to a very wiggley trajectory with
                 path length equal to twice its total displacement.
        '''
        self.fname = fname
        self.U = mean_flow
        self.max_dt = max_dt
        self.Ns = Ns
        self.d_max = d_max
        self.dv_max = dv_max
        self.NSR_th = NSR_th
        
        # particles are stored in a dictionary. Keys are frame numbers, values
        # are lists of arrays, each array representing a particle
            
        data = read_csv(self.fname, header=None, sep='\t')
        timeIndex = data.shape[1] - 1
        self.particles = dict([(g, hstack([ones((len(k),1))*-1, k.values])) 
                               for g,k in data.groupby(timeIndex)])
        
        self.times = sorted(list(self.particles.keys()))
        
        # setting up a dictionary to hold KDTrees for candidate searches
        self.trees = {}
        
        # A dictionary to hold the identifiers of particles that were used up
        # already. Keys are frame numbers, items are indexes of used particles
        self.used_particles = dict([(tm, []) for tm in self.times])
        
        # a list to store trajectories
        self.trajs = []



    def build_trajectory(self, particle_index, backwards=False):
        '''
        Given an inital particle this function attempts to construct a
        trajectory out of it by linking it forward in time using the "best 
        estimate" heuritic. The initial particle is given by the frame number
        and its index: particle_index = (time, particle index).
        
        If backwards=True then the tracking is backwards in time.
        '''
        
        t0, ind0 = particle_index
        p0 = self.particles[t0][ind0]
        traj = [p0]
        traj_indexes = [particle_index]
        continue_search = True
        
        while continue_search:

            # 1) get properties at last frame
            tm_i = traj[-1][-1]
            x_i = traj[-1][1:4]
            
            if len(traj)==1: 
                v_i= self.U
                
            else:
                dt_i = traj[-1][-1] - traj[-2][-1]
                dx_i = traj[-1][1:4] - traj[-2][1:4]
                v_i = dx_i / dt_i 
            
            if tm_i == self.times[-1]:
                break
            
            # 2) if this is the first sample we use nearest neighbor search
            if len(traj)==1:
                
                if backwards==False:
                    if tm_i+1 in self.times:
                        cand = self.get_nearest_neighbor(x_i+self.U, tm_i+1)
                    else: 
                        continue_search = False
                        continue
                
                elif backwards==True:
                    if tm_i-1 in self.times:
                        cand = self.get_nearest_neighbor(x_i-self.U, tm_i-1)
                    else:
                        continue_search = False
                        continue
                        
                
                c_tm, c_ind = cand[0]
                test1 = cand[1]<=self.d_max 
                test2 = c_ind not in self.used_particles[c_tm]
                
                if test1 and test2:
                    traj.append(self.particles[c_tm][c_ind])
                    traj_indexes.append(cand[0])
                
                else: continue_search = False
                
                continue
                    
            
            # 3) using "best estimate" heuristic we iterate over all allowable 
            # time separations, j: 1 -> max_dt
            for j in range(1, self.max_dt+1):
                
                # 3.0 ensure that frame ipj exists. Otherwise, stop the search
                if backwards==False:
                    tm_ipj = tm_i + j
                    if tm_ipj > self.times[-1]:
                        continue_search = False
                        break
                
                elif backwards==True:
                    tm_ipj = tm_i - j
                    if tm_ipj < self.times[0]:
                        continue_search = False
                        break
                
                if tm_ipj not in self.times:
                    break
                
                # 3.1 project the particle to frame i+j
                x_ipj = x_i + v_i * (tm_ipj - tm_i)
                
                # 3.2 search for neighbors of the projection within d_max
                candidates = self.search_neighbors(x_ipj, 
                                                   tm_ipj, 
                                                   dist=self.dv_max)
                
                
                # 3.3 remove candidates with dv larger than self.max_dv and 
                # candidates that were used up already
                if len(traj)>1:
                    def get_dv(cand):
                        x_cand = self.particles[cand[0]][cand[1]][1:4]
                        v_cand = (x_cand - x_i) / (cand[0] - tm_i)
                        dv = sum((v_cand - v_i)**2)**0.5
                        return dv
                    
                    candidates[:] = [cand for cand in candidates if 
                                      get_dv(cand)<=self.dv_max]
                
                candidates[:] = [cand for cand in candidates if 
                                 cand[1] not in self.used_particles[cand[0]]]
                
                
                # 3.4 if there are no candidates, continue to next j value
                if len(candidates)==0: continue
            
                # 3.5 if there is no i+j+1 frame, choose the i+j projection's 
                # nearest neighbour given that it is sufficiently close. 
                # Otherwise, stop the search.
                if tm_ipj+1 > self.times[-1]: 
                    NN = self.get_nearest_neighbor(x_ipj, tm_ipj)
                    if NN[1] <= self.d_max:
                        traj_indexes.append(NN[0])
                        traj.append(self.particles[NN[0][0]][NN[0][1]])
                        continue_search = False
                        break
                        
                # 3.6 if there's only one candidate, we choose it and continue
                # to the next frame
                elif len(candidates)==1:
                    traj_indexes.append((tm_ipj, candidates[0][1]))
                    traj.append(self.particles[tm_ipj][candidates[0][1]])
                    break
                
                # 3.7 if there are multiple candidates we choose the one whose
                # projection to frame i+j+1 has the nearest neighbor distance
                cand_nnd_list = []
                for cand in candidates:
                    
                    # 3.7.1 project candidate to i+j+1
                    tm_ipjp1 = tm_i + j + 1
                    x_cand = self.particles[cand[0]][cand[1]][1:4]
                    v_cand = (x_cand - x_i) / (tm_ipj - tm_i)
                    x_ipjp1 = x_cand + v_cand * (tm_ipjp1 - tm_ipj)
                    
                    # 3.7.2 get the i+j+1 projections' nearest neihbor distance
                    nnd = self.get_nearest_neighbor(x_ipjp1, tm_ipjp1)[1]
                    cand_nnd_list.append(nnd)
                # 3.7.3 choose the candidate with minimal nnd
                chosen = candidates[cand_nnd_list.index(min(cand_nnd_list))]
                traj_indexes.append((chosen[0], chosen[1]))
                traj.append(self.particles[chosen[0]][chosen[1]])
                break
            
            
            # 4) if no valid candidates were found then we terminate 
            # the trajectory
            if traj[-1][-1] == tm_i:
                
                continue_search = False
                
        return array(traj), traj_indexes



    
    def get_nearest_neighbor(self, pos, tm):
        '''
        Given a location, pos (= [x,y,z]), this function returns its nearest 
        neighbor at frame tm and the distance to it.
        '''
        # 1) Get the relevant KDTree
        try:
            tree = self.trees[tm]
        except:
            tree = KDTree(self.particles[tm][:,1:4])
            self.trees[tm] = tree
        
        # 2) return particles in the neighborhood
        nn = tree.query(pos) 
        
        return ((tm, nn[1]), nn[0]) 
    
    
    
    
    def search_neighbors(self, pos, tm, dist=None):
        '''
        Given a location, pos (= [x,y,z]), this function returns particles 
        that their distance from pos is up to self.d_max at time tm. The 
        function does not return the particles themselves but their tuple 
        identifires: (time, particle index)
        '''
        # 1) Get the relevant KDTree
        try:
            tree = self.trees[tm]
        except:
            tree = KDTree(self.particles[tm][:,1:4])
            self.trees[tm] = tree
        
        
        # 2) return particles in the neighborhood
        if dist is None: dist=self.d_max
        neighbor_inds = tree.query_ball_point(pos, dist)
        
        return [(tm, ind) for ind in neighbor_inds] 
    
    
    
    
    def build_trajectories_from_frame(self, frame_num, backwards=False, 
                                      p_bar=True):
        '''
        Builds trajectories from particles at a given frame number. We do not
        use particles that have been used up already.
        
        If backwards=True then the tracking is backwards in time.
        '''
        
        # ensure there are particles to track in this frame number
        if frame_num not in list(self.particles.keys()):
            return None
        
        # a list to hold trajectories and trajectory particle identifiers
        trajs = []
        
        # a list to hold the trajectory noise levels
        NSRs = []
        
        # 0) if above 20% of particles in the frame were used, clear them
        if frame_num in list(self.used_particles.keys()):
            Nused = len(self.used_particles[frame_num])
            Nf = len(self.particles[frame_num])
            uf = Nused / Nf
            if uf > 0.2:
                self.clear_used_particles()
        
        # 1) building trajectories from particles in the given frame number
        if backwards==False: 
            msg = 'Forwards tracking frame %d'%frame_num
        else: 
            msg = 'Backwards tracking frame %d'%frame_num
        
        if p_bar==True:
            iter_ = tqdm.tqdm(range(len(self.particles[frame_num])), desc=msg)
        
        else:
            iter_ = range(len(self.particles[frame_num]))
            
        for i in iter_:
            
            # 1.1 check if the particle was used up already
            if i in self.used_particles[frame_num]: continue
            
            # 1.2 build a trajectory, calculate its noise level and store them
            p_id = (frame_num, i)
            trajs.append(self.build_trajectory(p_id, backwards=backwards))
            NSR = traj_NSR(trajs[-1][0], self.Ns)
            if len(trajs[-1][0])<self.Ns: noise_lvl = 1
            else: noise_lvl = mean(NSR[int(self.Ns/2):-int(self.Ns/2)])
            NSRs.append(noise_lvl)
        
        
        # 2) prune the "good" trajectories
        
        # 2.1 sorting the results according to trajectory length
        res = sorted(zip(trajs, NSRs), key=lambda x: len(x[0][0]),reverse=True)
        
        # 2.2 if a trajectory has low NSR we make sure it doesn't have used up
        # particles and if so we add it to self.trajs
        for (tr, tr_ids), NSR in res:
            
            # check noise level
            if NSR > self.NSR_th: continue
                
            # check it wasn't used up
            for p_id in tr_ids:
                if p_id[1] in self.used_particles[p_id[0]]: continue
            
            # add to self.trajs
            if backwards==False:
                self.trajs.append(tr)
            elif backwards==True:
                self.trajs.append(tr[::-1])
            
            # write up the trajectories particles
            for p_id in tr_ids:
                self.used_particles[p_id[0]].append(p_id[1])
    
    
    
    
    def clear_used_particles(self):
        '''
        This function removes the particles that have been used up to make
        good trajectories from self.particles and renews the KDTrees in order
        to make the search for candidates more efficient.
        '''
        
        for tm in self.times:
            used_ind = self.used_particles[tm]
            unused = [i for i in range(len(self.particles[tm])) if i not in used_ind]
            self.particles[tm] = self.particles[tm][unused]
        
        self.trees = {}
        self.used_particles = dict([(tm, []) for tm in self.times])
        
        
        
        
    def track_frames(self, f0=None, fe=None, frame_skips=2):
        '''
        Will build trajectories starting from frames in the range f0 -> fe,
        with given skips. We build trajectories both forward and backwards in 
        time. 
        
        f0 - The frame from which we start. If None then we start from the 
             first available frame.
        
        fe - The frame at which we end. If None then we end at the last 
             available frame.
             
        frame_skips - within the frame range given, we use skips with this 
                      many frames. For example if equal to 2 then we use: 
                      0, 2, 4, 6, ... for tracking forward in time and 
                      -1, -3, -5, -7, ... for tracking backwards in time.
        '''
        
        if f0 is None: f0 = int(self.times[0])
        if fe is None: fe = int(self.times[-1])
        
        msg = 'Tracking forward and backward, round 1'
        
        for frm in tqdm.tqdm(range(f0, fe, frame_skips), desc=msg):
            self.build_trajectories_from_frame(frm, p_bar=False)
            self.build_trajectories_from_frame(fe+f0-frm-1, 
                                              backwards=True, p_bar=False)
            
        msg = 'Tracking forward and backward, round 2'
        
        for frm in tqdm.tqdm(range(f0+int(frame_skips/2), 
                                    fe-int(frame_skips/2), 
                                    frame_skips), desc=msg):
            
            self.build_trajectories_from_frame(frm, p_bar=False)
            self.build_trajectories_from_frame(fe+f0-frm-1, 
                                              backwards=True, p_bar=False)
        
        print('')
        print('')
        print('finished tracking frames')
        print('trajs: ', len(self.trajs))
        print('mean length: ', mean([len(tr) for tr in self.trajs]))
        print('linked particles: ', sum([len(tr) for tr in self.trajs]))
        
        
        
        
        
    def interpolate_trajs(self):
        '''
        Will interpolate trajectories that have skipped frames. We interpolate
        the missing points by using a 3nd order polynomial, namely assuming 
        linear acceleration in the interpolated range.
        '''
        
        N_links_0 = sum([len(tr)+1 for tr in self.trajs])
        
        msg = 'Interpolating skipped frames'
        for i in tqdm.tqdm(range(len(self.trajs)), desc=msg):
            
            tr = self.trajs[i]
            
            # 1) check if the trajectory has skipped frames; if not, we
            #    continue to the next trajectory
            dt = tr[-1,-1] - tr[0,-1]
            l = len(tr)
            if dt == l-1: continue
            
            # 2) interpolate the missing points by interpolation
            interpolated = fill_in_trajectory(tr)
            self.trajs[i] = interpolated
            
        N_links_e = sum([len(tr)+1 for tr in self.trajs])
        N_new = N_links_e-N_links_0
        stats = (N_new, N_new/N_links_e*100)
        print('')
        print('interpolated %d points (%.1f percent)'%stats)
        
        
        
        
    def save_results(self, fname):
        '''
        Will save the tracking results. We save both the trajectories and the
        particles that were not used.
        '''
        to_save = []
        
        # 1) add ids to the trajectories and add them to the list
        for i in range(len(self.trajs)):
            self.trajs[i][:,0] = i+1
            to_save += list(self.trajs[i])
        
        # 2) sort according to the frame number
        to_save = sorted(to_save, key=lambda x: x[-1])
        
        # 3) add the unused particles to the list
        self.clear_used_particles()
        for tm in self.times:
            to_save += list(self.particles[tm])
            
        # 4) save the data in MyPTV format
        fmt = ['%d', '%.3f', '%.3f', '%.3f']
        for i in range(len(to_save[0])-6):
            fmt.append('%d')
        fmt += ['%.3f', '%.3f']
        savetxt(fname , to_save, delimiter='\t', fmt=fmt)
        
    
    
    
    
    # ========================================================================
    # \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    # ========================================================================
    #         Tested and not used functions, kept for legacy.
    # ========================================================================
    # \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    # ========================================================================
        
    
    def get_particle(self, particle_identifier):
        '''
        Given a particle_identifier (=(frame number, particle index)), this
        function returns the particle.
        '''
        tm, ind = particle_identifier
        return self.particles[tm][ind]
    
    
    
    def get_candidate_links(self, traj):
        '''
        Given a trajectory, this function returns a list of particle 
        identifiers that are allowed to be linked to it in the future.
        '''
            
        tm_i = traj[-1][-1] ; x_i = traj[-1][1:4]
        
        if len(traj)==1:
            v_i = self.U
            x_proj = x_i + v_i
            candidates = self.search_neighbors(x_proj, tm_i+1)
            
        else:
            tm_prev = traj[-2][-1] ; x_prev = traj[-2][1:4]
            v_i = (x_i - x_prev) / (tm_i - tm_prev)
            candidates = []
            for dt in range(1, self.max_dt+1):
                x_proj = x_i + v_i * dt
                candidates += self.search_neighbors(x_proj, 
                                                    tm_i+dt, 
                                                    dist=self.dv_max)
                if len(candidates)>0: break
        
        return candidates
    
    
    
    
    def build_candidate_trajectories(self, particle_identifier):
        '''
        Given an initial particle_identifier, this function returns a list of 
        all possible trajectories that could be built from this particle.
        '''
        p0 = self.get_particle(particle_identifier)
        trajs = [[p0]]
        traj_indexes = [[particle_identifier]]
        broken_checks = [True]
        
        test_cands = [True]
        it = 0
        while any(test_cands):
                
            ntraj = len(trajs)
            test_cands = [False for i in range(ntraj)]
            
            for i in range(ntraj):
                
                if broken_checks[i]==False: break
                
                cand_links = self.get_candidate_links(trajs[i])
                
                if len(cand_links)>0:
                    test_cands[i] = True
                
                else:
                    broken_checks[i] = False

                for j in range(len(cand_links)):
                    if j==0:
                        traj_indexes[i].append(cand_links[j])
                        trajs[i].append(self.get_particle(cand_links[j]))

                    else:
                        trajs.append(trajs[i][:-1])
                        trajs[-1].append(self.get_particle(cand_links[j]))
                        traj_indexes.append(traj_indexes[i][:-1])
                        traj_indexes[-1].append(cand_links[j])
                        broken_checks.append(True)

            it+=1
            #if it>3: break
            
        return [array(tr) for tr in trajs], traj_indexes



    def tracking_movie(self, particle_identifier):
        '''
        This animates the tracking of a given particle with all the candidates
        '''
        
        ret = self.build_candidate_trajectories(particle_identifier)[0]
        
        for tr in ret:
            xmin = min(min(tr, key=lambda x: min(x[:,1]))[:,1]) - 2*self.d_max
            xmax = max(max(tr, key=lambda x: max(x[:,1]))[:,1]) + 2*self.d_max
            ymin = min(min(tr, key=lambda x: min(x[:,2]))[:,2]) - 2*self.d_max
            ymax = max(max(tr, key=lambda x: max(x[:,2]))[:,2]) + 2*self.d_max
            t0 = min(min(tr, key=lambda x: min(x[:,-1]))[:,-1]) - 2*self.d_max
            te = max(max(tr, key=lambda x: max(x[:,-1]))[:,-1]) + 2*self.d_max
        
        
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        
        for fe in range(t0, te+1):    
            f0 = max([0, fe-2])
            ax.clear()    
                
            for tr in ret:
                whr = tr[:,-1] <= fe
                ax.plot(tr[whr,1], tr[whr,2], 'o-', alpha=0.3)
            
            for frm in range(f0, fe+1):
                for p in self.particles[frm]:
                    if xmin<p[1]<xmax and ymin<p[2]<ymax:
                        ax.plot(p[1], p[2], 'ko', ms=2)
              
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymin, ymax)
            ax.set_aspect('equal')
            plt.tight_layout()
            fig.savefig('%2d.jpg'%fe)
    
    # ========================================================================
    # \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    # ========================================================================
    # ========================================================================
    # \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    # ========================================================================
        
        
        
        
        
def traj_NSR(traj, Ns):
    '''
    A trajectory is considered "noisy at a scale Ns" if its path is long as
    compared to the displacements it makes at this scale. We thus define the
    noise level (NSR) of a trajectory sample with scale parameter Ns to be
    the path length it makes over a window of size Ns samples around it minus
    the diplacement over this window devided by the displacement:
        
        path_length[i] = sum( dx[i-Ns/2:i+Ns/2] )
        displacement[i] = x[i+Ns/2] - x[i-Ns/2]
        NSR[i] = (path_length[i] - displacement[i]) / displacement[i]
        
    '''
    
    if Ns%2 != 1:
        raise ValueError('Ns must be an odd interger')
        
    w = int(Ns/2)
    
    if len(traj)<=Ns:
        NSR = zeros(len(traj))
        return NSR
    
    NSR = []
    
    for i in range(w):
        NSR.append(0)
        
    for i in range(w, len(traj)-w):
        x_ = traj[i-w:i+w+1, 1:4]
        signal = sum((x_[-1] - x_[0])**2)**0.5
        noise = sum(npsum((x_[1:]-x_[:-1])**2, axis=1)**0.5) - signal
        NSR.append(noise/signal)
    
    for i in range(-w,0):
        NSR.append(0)
    
    return NSR
    
        


def fill_in_trajectory(tr):
    '''
    This function takes in a trajectory and fills in points that were skipped
    in the tracking stage. This is done by interpolating with a 3nd degree
    polynomial.
    '''
    new_points = []
    
    for i in range(len(tr)-1):
        # 1) check if t_i+1 was skipped; if not, we continue
        if tr[i+1,-1] - tr[i,-1] == 1: continue
        
        # 2) get data from 2 frames before and 2 frames after the skipped gap
        # if the missing link is at the second point
        if i==0: 
            tm = tr[:i+4,-1]
            x = tr[:i+4,1]
            y = tr[:i+4,2]
            z = tr[:i+4,3]
        
        # if the missing link is at the second last point
        elif i==len(tr)-2: 
            tm = tr[-4:,-1]
            x = tr[-4:,1]
            y = tr[-4:,2]
            z = tr[-4:,3]
        
        # the normal case of missing point in the middle section
        elif i!=0:
            tm = tr[i-1:i+3,-1]
            x = tr[i-1:i+3,1]
            y = tr[i-1:i+3,2]
            z = tr[i-1:i+3,3]
        
        
        # 3) get a list of the times that have been skipped around i
        tm_interp = [tr[i,-1]+j for j in range(1, int(tr[i+1,-1] - tr[i,-1]))]
        
        # 4) for each axis we calculate the missing data
        t0 = tr[i,-1]
        px = poly1d(polyfit(tm-t0, x, 3)) ; x_interp = px(array(tm_interp)-t0)
        py = poly1d(polyfit(tm-t0, y, 3)) ; y_interp = py(array(tm_interp)-t0)
        pz = poly1d(polyfit(tm-t0, z, 3)) ; z_interp = pz(array(tm_interp)-t0)
        
        # 5) add the new points to the list
        for j in range(len(tm_interp)):
            new_point = [-1 for k in range(len(tr[i]))]
            new_point[0] = tr[i,0]
            new_point[1] = x_interp[j]
            new_point[2] = y_interp[j]
            new_point[3] = z_interp[j]
            new_point[-1] = tm_interp[j]
            new_points.append(array(new_point))
        
    # 6) combine the original and interpolated points to a new trajectory
    new_traj = array(sorted(list(tr) + new_points, key=lambda x: x[-1]))
    
    return new_traj
        
            
        
    
        
        
# =============================================================================
# # Tests:
# if __name__ == "__main__":
#     fname = '/home/ron/Desktop/Research/jetArrayTank/20241020_puffs/Rec18/particles'
#     max_dt = 3 
#     Ns = 11
#     NSR_th = 0.25
#     tmf = tracker_multiframe(fname, max_dt, Ns, d_max=0.5, dv_max=0.5, NSR_th=NSR_th)
#     
#     f0 = None
#     fe = None
#     frame_skips = 3
#     tmf.track_frames(f0=f0, fe=fe, frame_skips=frame_skips)
#     tmf.interpolate_trajs()
#     
#     tmf.save_results('trajectories_multiframe')
# =============================================================================











