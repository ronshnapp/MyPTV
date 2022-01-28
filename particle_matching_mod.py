#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 12:15:21 2019

@author: ron


Matching module for extracting lab coordinate particles from 
segmented particles' image coordinates.

We are using the Ray Traversal algorithm reported 
by Bourgion and Huisman, 2020 (https://arxiv.org/pdf/2003.12135.pdf).


"""

from math import ceil, floor
from itertools import combinations


class matching(object):
    '''A class for matching particles in images taken simultaniously
    from different cameras.'''
    
    
    def __init__(self, img_system, particles_dic, 
                 RIO, voxel_size):
        '''
        img_system - is an instance of the img_system object with camera 
                     objects. 
                     
        particles_dic - A dictionary with keys that are camera names, and values
                     are lists of particle coordinates segmented in each of 
                     the cameras.
                     
        RIO - A nested list of 3X2 elements. The first holds the minimum and 
              maximum values of x coordinates, the second is same for y, and 
              the third for z coordinates. 
        
        voxel size - the side length of voxel cubes used in the ray traversal
                     algorithm. Given in lab coordinate scales (e.g. mm).
        '''
        
        self.imsys = img_system
        
        # make a list of rays with: 
        # [camera x coord, camera y coord, (cam number, particle number)]
        self.rays = []
        self.ray_camera_indexes = [0]
        for i in range(len(self.imsys.cameras)):
            cam  = self.imsys.cameras[i]
            particles_i = particles_dic[cam.name]
            self.ray_camera_indexes.append(len(particles_i) + 
                                           self.ray_camera_indexes[-1])
            for j in range(len(particles_i)):
                x, y = particles_i[j][0], particles_i[j][1]
                self.rays.append( (x, y, (i,j)) )
        
        self.RIO = RIO
        self.voxel_size = voxel_size
        
        
        # set up lists of voxel centers:
            
        Nx = ceil((RIO[0][1]-RIO[0][0])/voxel_size)
        cx = (RIO[0][1]+RIO[0][0])/2.
        if Nx%2==0: f = (floor(Nx/2)-0.5) * voxel_size
        elif Nx%2!=0: f = floor(Nx/2)*voxel_size
        self.Nx = Nx
        self.x = [i*voxel_size + cx - f for i in range(Nx)]
        
        Ny = ceil((RIO[1][1]-RIO[1][0])/voxel_size)
        cy = (RIO[1][1]+RIO[1][0])/2.
        if Ny%2==0: f = (floor(Ny/2)-0.5) * voxel_size
        elif Ny%2!=0: f = floor(Ny/2)*voxel_size
        self.Ny = Ny
        self.y = [i*voxel_size + cy - f for i in range(Ny)]
        
        Nz = ceil((RIO[2][1]-RIO[2][0])/voxel_size)
        cz = (RIO[2][1]+RIO[2][0])/2.
        if Nz%2==0: f = (floor(Nz/2)-0.5) * voxel_size
        elif Nz%2!=0: f = floor(Nz/2)*voxel_size
        self.Nz = Nz
        self.z = [i*voxel_size + cz - f for i in range(Nz)]
        
        
    def get_traversed_voxels(self):
        '''This will return a list that holds for each ray in self.rays
        the pixels in RIO through which it traverses. '''
        
        traversed_voxels = []
        for ray in self.rays:
            cam  = self.imsys.cameras[ray[2][0]]
            O = cam.O
            r = cam.get_r(ray[0], ray[1])
            r_ = r / sum(r**2)**0.5
            
            for k in range(len(self.z)):
                z_ = self.z[k]
                a = (z_ - O[2])/r_[2]
                x_, y_ = O[0] + r_[0]*a, O[1] + r_[1]*a
                if x_>self.RIO[0][1] or x_<self.RIO[0][0]: continue
                if y_>self.RIO[1][1] or y_<self.RIO[1][0]: continue
                i = int((x_ - self.x[0] - self.voxel_size/2)/self.voxel_size +1)
                j = int((y_ - self.y[0] - self.voxel_size/2)/self.voxel_size +1)
                
                for i_ in range(max([0,i-1]), min([self.Nx-1,i+1])+1):
                    for j_ in range(max([0,j-1]), min([self.Ny-1,j+1])+1):
                        for k_ in range(max([0,k-1]), min([self.Nz-1,k+1])+1):
                            traversed_voxels.append( [(i_,j_,k_), ray[2]] )
        
        return traversed_voxels
        
    
    def get_voxel_dictionary(self):
        '''This generates a dicionary who's keys are voxel indexes and
        who's values are the rays that passed through this voxel.'''
        
        traversed_voxels = self.get_traversed_voxels()
        
        voxel_dic = {}
        for vxl in traversed_voxels:
            if vxl[0] in voxel_dic.keys():
                if vxl[1] not in voxel_dic[vxl[0]]:
                    voxel_dic[vxl[0]].append(vxl[1])
            else:
                voxel_dic[vxl[0]] = [vxl[1]]
        
        self.voxel_dic = voxel_dic
        
    
    def list_candidates(self):
        '''This will make lists of possible candidate rays for
        triangulation, separated for pairs, triplets, quadruplets, etc.
        
        Candidates are based on the voxels of voxel_dic, while calculating the 
        RMS and the maximum distance between the estimated particle location 
        and the epipolar lines.'''
        
        candidate_dic = {}
        for i in range(2, len(self.imsys.cameras)+1):
            candidate_dic[i] = []
        
        for group_size in sorted(list(candidate_dic.keys()), reverse=True):
            for k in self.voxel_dic.keys():
                if len(self.voxel_dic[k]) >= group_size:
                    for comb in combinations(self.voxel_dic[k], group_size):
                        cam_nums = [ray[0] for ray in comb]
                        #ray_nums = [ray[1] for ray in comb]
                        if len(cam_nums) == len(set(cam_nums)):
                            if comb not in candidate_dic[group_size]: 
                                candidate_dic[group_size].append(comb)
                            # dc = {}
                            # for i in range(len(comb)):
                            #     eta, zeta = self.rays[self.ray_camera_indexes[cam_nums[i]]:self.ray_camera_indexes[cam_nums[i]+1]][ray_nums[i]][:2]
                            #     dc[cam_nums[i]] = [eta, zeta]
                            
                            # coord, dump, dist = self.imsys.stereo_match(dc, self.voxel_size*2000)
                            # candidate_dic[group_size].append([k, comb, coord, dist])
        
        self.candidate_dic = candidate_dic
        
        
        
    def triangulate_rays(self, rays):
        '''will return the results of stereo matching of a list of rays'''
        dc = {}
        for ray in rays:
            i = self.ray_camera_indexes[ray[0]] 
            ip1 = self.ray_camera_indexes[ray[0]+1]
            eta, zeta = self.rays[i:ip1][ray[1]][:2]
            dc[ray[0]] = [eta, zeta]
        return self.imsys.stereo_match(dc, 1e19)
        
        
        
    def plot_ray_epipolar_lines(self, ray, zlims, ax):
        '''will plot a ray's epipolar line for a given 3D axis.'''
        import matplotlib.pyplot as plt
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        camera = self.imsys.cameras[ray[0]]
        i = self.ray_camera_indexes[ray[0]] 
        ip1 = self.ray_camera_indexes[ray[0]+1]
        eta, zeta = self.rays[i:ip1][ray[1]][:2]
        camera.plot_3D_epipolar_line(eta, zeta, zlims, ax=ax, color=colors[ray[0]])
        
        
        
        
        
        
        
        
if __name__ == '__main__':
    from imaging_mod import camera, img_system
    import numpy as np
    import os
    
    dirname ='/home/ron/Desktop/Research/plankton_sweeming/experiments/PTV_test3'
    cnames = ['cam1', 'cam3', 'cam4']
    res = 1280, 1024
    cameras = []
    for c in cnames:
        cameras.append(camera(c, res))
    
    for cam in cameras:
        cam.load(dirname)
    
    particles_dic = {}
    for c in cnames:
        blobs = np.loadtxt( os.path.join(dirname, 'blobs_'+c) )
        particles_dic[c] = blobs[:,:2]
    
    RIO = ((0.0, 68.0),
           (0.0, 66.0),
           (-30.0, 20.0))
    
    voxel_size = 2.0
    
    imsys = img_system(cameras)
    
    m = matching(imsys, particles_dic, RIO, voxel_size)
    m.get_voxel_dictionary()
    m.list_candidates()
    

        
        
        
        
        
        
        