# -*- coding: utf-8 -*-
"""
Created on May 2022

@author: ron

This module is used to extend the matching and tracking of particles using
temporal data, using the ideas developed by Schants et al. 2016 in their 
Shake the Box method.

"""




from scipy.spatial import KDTree
from scipy.optimize import minimize
from itertools import combinations, product
from numpy import append as npappend


class extended_linking(object):
    '''
    This class takes in a trajectory that its last point is on frame i, and
    by looking back at the blob files for frame i+1, tries to find whether it
    may have a connection that was missed in the matching stage.
    We do this by 1) estimating the particle's position at frame i+1 using its 
    velocity, 2) projecting this position to camera space coordinates at all
    cameras, 3) taking the blobs that are nearest to the projection in each
    camera, and 4) attemptin to stereo match these nearest neighbouring blobs. 
    If this is successfull, the new particle sample is added to the list of 
    trajectories with the same trajectory number ID as the one given initially.
    '''
    
    def __init__(self, imsys, particles_dic, pixel_neighbourhood_size, dv_max,
                 max_err=1e9,
                 reverse_eta_zeta=False):
        '''
        imsys - an instance of the img_system class
        
        particles_dic - a dictionary with keys being camera names, and values
                        being lists of blobs that belong to the frame i+1 and
                        camera.
                        
        max_err - float, the maximum allowable triangulation error in the 
                  stereo matching
                  
        reverse_eta_zeta - should the order of image space coordinates be
                           reversed.
        '''
        
        self.imsys = imsys
        self.pd = particles_dic
        self.reverse_eta_zeta = reverse_eta_zeta
        self.max_err = max_err
        self.r = pixel_neighbourhood_size
        self.dv_max = dv_max
        
        self.trees = dict([(k, KDTree(self.pd[k][:,:2])) 
                           for k in self.pd.keys()])
    
    
    
    def nearest_blobs_distance(self, x):
        '''
        Given a location x, this returns the sum of the distances between the
        projections of x on all cameras and their nearest neighbours
        '''
        dist_list = []
        for cam in self.imsys.cameras:
            if self.reverse_eta_zeta==False:
                proj = cam.projection(x)[::-1]
            else:
                proj = cam.projection(x)
            dist_list.append( self.trees[cam.name].query([proj])[0][0] )
        
        return dist_list

        
    
    def optimize_particle_position(self, x0, bounds):
        '''
        Given an initial position x0, this function tries to "shake" it
        such that its nearest blobs distance in image space is minimized.
        '''
        cost_func = lambda a: sum(self.nearest_blobs_distance(a))
        
        xmin = minimize(cost_func, x0, bounds=bounds)
        return xmin
    
    
    
    
    def search_for_stereo_matching(self, tr):
        '''
        Given a point in lab space, x, this will search for the blobs that 
        are this points' neighbourhood in image space and will stereo 
        match them.
        '''
        
        # calculate the projected position
        x = tr[-1,1:4] + tr[-1,1:4] - tr[-2,1:4]   
        
        # for each camera, get the neighbouring blobs and put them in lists,
        # according to their camera
        candidates_by_cams = []
        for e,cam in enumerate(self.imsys.cameras):
            
            if self.reverse_eta_zeta==False:
                proj = cam.projection(x)[::-1]
                neighbrs = self.trees[cam.name].query_ball_point(proj, self.r)
                candidates_by_cams.append([(self.pd[cam.name][ind][:2], e, ind)
                                           for ind in neighbrs])
                
            else:
                proj = cam.projection(x)
                neighbrs = self.trees[cam.name].query_ball_point(proj, self.r)
                candidates_by_cams.append([(self.pd[cam.name][ind][:2], e, ind)
                                           for ind in neighbrs])
        
        # make sure there are candidates from at least 2 cameras; else, 
        # return None
        cams_used = sum([len(cand)>0 for cand in candidates_by_cams])
        if cams_used<2: return None
            
        # for all possible group sizes (pairs, triplets, quadruplets, etc.)
        # get all possible blob combinations; triangulate these conbinations
        group_sizes = range(2, len(self.imsys.cameras)+1)[::-1]
        
        # find all possible combinations for ll group sizes:
        blob_combinations = []
        for gs in group_sizes:
            for comb in combinations(candidates_by_cams, gs):
                blob_combinations += product(*comb)
        
        # triangulating the conbinations
        particle_candidates = []
        for comb in blob_combinations:
            dic = {}
            for element in comb:
                dic[element[1]] = element[0]
            triang = self.imsys.stereo_match(dic, 1e9)
            if triang[-1]<= self.max_err:
                particle_candidates.append((triang, comb))
        
        # from all candidates, take the one closest to the original x
        dist = lambda p: sum((p[0][0]-x)**2)**0.5
        sorted_cand = sorted(particle_candidates, key = dist)
        
        # choosing the nearest candidate with admisible velocity change
        for cand in sorted_cand:
            
            #testing velocity change
            v_old = tr[-1,1:4] - tr[-2,1:4] 
            v_new = cand[0][0] - tr[-1,1:4]
            if sum((v_new - v_old)**2)**0.5 / sum(v_old**2)**0.5 > dv_max:
                return cand

        return None
    
    
    
    
    
#%%
if __name__ == '__main__':
    import os
    from myptv.imaging_mod import camera, img_system
    from pandas import read_csv
    from numpy import array
    import numpy as np
    
    folder = '/home/ron/working_PTV_data/1'
    cam_names = ['cam1', 'cam2', 'cam3', 'cam4']
    cam_lst = [camera(cn, (1280, 1024)) for cn in cam_names]
    for cam in cam_lst: 
        cam.load(folder)
        
    imsys = img_system(cam_lst)
    
    blob_fnames = ['blobs_cam1', 'blobs_cam2', 'blobs_cam3', 'blobs_cam4']
    blob_files = [os.path.join(folder, b) for b in blob_fnames]
    
    traj_fname = os.path.join(folder, 'trajectories')
    trajs = [array(g) for k,g in 
             read_csv(traj_fname, sep='\t', header=None).groupby(by=0)]
    
    frames = set([tr[i][-1] for tr in trajs for i in range(len(tr))])
    
    ROI = [[0,70], [0,70], [-20, 10]]
    pixel_neighbourhood_size = 10.0
    dv_max = 0.5
    
    # ============================
    frame = 1.
    count = 0
    
    
    test = True
    iters = 0
    trajs_to_check = []
    
    while test:
        extended_trajs = []
        trajs = sorted(trajs, key = lambda x: x[-1,-1])
        
        for ee,tr in enumerate(trajs[:]):
            
            if iters>0 and tr[0,0] not in trajs_to_check: continue
            
            if tr[0,0]==-1: continue
            
            if tr[-1][-1]+1 != frame:
                frame = tr[-1][-1] + 1
                
                if frame == max(frames): break # <- we are at the last frame
                
                print('iter: ', iters, 'frame: ', frame, end='\r')
                print('')
                
                blobs_data = [read_csv(bf, delimiter='\t', header=None) for bf in blob_files]
                
                pd = {}
                for e, bd in enumerate(blobs_data):
                    cam_name = imsys.cameras[e].name
                    grouped = [(k,array(g)) for k,g in bd.groupby(by=5)]
                    for k,g in grouped:
                        if k == frame:
                            pd[cam_name] = g
                            break
                
                el = extended_linking(imsys, pd, pixel_neighbourhood_size, 
                                      dv_max, max_err=10,reverse_eta_zeta=True)

            res = el.search_for_stereo_matching(tr)
            
            # extending the trajectory with the newly triangulated particle
            if res is not None:
                
                blobs_list = [-1 for i in range(len(imsys.cameras))]
                for blob in res[1]:
                    blobs_list[blob[1]] = blob[2]
                    
                new_p = [tr[0][0]] + list(res[0][0])+blobs_list+[res[0][2]]+[frame]
                tr = npappend(tr, [new_p], axis=0)
                trajs[ee] = tr
                extended_trajs.append(tr[0][0])
                count += 1
                
        test = len(extended_trajs)>0
        trajs_to_check = extended_trajs.copy()
        iters += 1
        
    print(count)

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    