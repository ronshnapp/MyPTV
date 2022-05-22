# -*- coding: utf-8 -*-
"""
Created on May 2022

@author: ron

A module designed for optimizing the stereo-matching and the triangulation 
of particles after the mathcinhg is performed. This is similar to the 
shaking stages used in the iterative-particle-reconstruction and the 
shake-the-box schemes.

"""




from scipy.spatial import KDTree


class particle_optimizer(object):
    '''
    This is an object whose goal is to "correct" the 3D positions of given 
    points in 3D space, such that they will best agree with the given location
    of segmented blobs. More specifically, we want to find particle positions
    that their projection on the cameras minimize the image-space distance to 
    nearest segmented blobs.
    '''
    
    def __init__(self, imsys, particles_dic, search_radius, max_err=1e9, 
                 reverse_eta_zeta=False):
        '''
        imsys - an instance of the img_system class
        
        particles_dic - a dictionary with keys being camera names, and values
                        being lists of blobs that belong to the frame i and
                        camera.
                        
        search_area - the maximum distance in image-space coordinates in which 
                      we look at blobs around a point, x. In the optimization, 
                      the particle will not be moved outside of this region.
                        
        max_err - float, the maximum allowable triangulation error in the 
                  stereo matching
                  
        reverse_eta_zeta - should the order of image space coordinates be
                           reversed.
        '''
        
        self.imsys = imsys
        self.pd = particles_dic
        self.search_radius = search_radius
        self.reverse_eta_zeta = reverse_eta_zeta
        self.max_err = max_err
        
        self.trees = dict([(k, KDTree(self.pd[k][:,:2])) 
                           for k in self.pd.keys()])
    
    
    
    def find_valid_blobs(self, x):
        '''
        given a point in 3D lab-space, x, this function returns the image-
        space coordinates of blobs that are found close to the projection of
        x onto the cameras.
        '''
        
        valid_blobs = []
        for e,k in enumerate(self.trees.keys()):
            
            proj_cam_e = self.imsys.cameras[e].projection(x)
            
            valid_blobs.append([self.pd[k][ind] for ind in
                                self.trees[k].query_ball_point([proj_cam_e], 
                                                          self.search_radius)])
        return valid_blobs
            
    
    
    
    
    
    
    
    
    
    
    
    
    
    