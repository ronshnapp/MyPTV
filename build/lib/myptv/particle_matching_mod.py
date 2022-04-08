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
from itertools import combinations, product
from numpy import loadtxt, savetxt





class match_blob_files(object):
    '''A class for obtaining triangulated particles positions from a 
    list of segmented blobs. Use self.get_particles() and after that,
    the particles found are stored in the attribute self.particles.'''
    
    
    def __init__(self, blob_fnames, img_system, RIO, voxel_size, max_err=None,
                 reverse_eta_zeta = False):
        '''
        blob_fname - a list of the file names containing the segmented blob
                     data. The list has to be sorted according the order of
                     cameras in the img_system.
                     
        img_system - an instance of the img_system class with the calibrated
                     cameras.
                     
        RIO - A nested list of 3X2 elements. The first holds the minimum and 
              maximum values of x coordinates, the second is same for y, and 
              the third for z coordinates. 
        
        voxel size - the side length of voxel cubes used in the ray traversal
                     algorithm. Given in lab coordinate scales (e.g. mm).
                     
        max_err - maximum acceptable uncertainty in particle position. If None,
                  (defult), than no bound is used.
                     
        reverse_eta_zeta - Should be false if the eta and zeta coordinates 
                           need to be in reverse order so as to match the
                           calibration. This may be needed if the calibration 
                           data points were given where the x and y coordinates
                           are transposed (as happens, e.g., if using 
                           matplotlib.pyplot.imshow).
        '''
        self.blobs = []
        for fn in blob_fnames:
            self.blobs.append(loadtxt(fn))
        self.imsys = img_system
        self.RIO = RIO
        self.voxel_size = voxel_size
        self.reverse_eta_zeta = reverse_eta_zeta
        self.max_err = max_err
        
        time_lst = []
        for bl in self.blobs:
            for b in bl:
                time_lst.append(b[-1])
        self.time_lst = list(set(time_lst))
        
        
        
    def get_particles(self, frames=None):
        '''
        Use this to match blobs into particlesin 3D.
        
        input - 
        frames - if None, this will match particles at all times. If a list of
                 integers, this will only perform the matching on particles at
                 the given frames.
        '''
        # set up a list of the frames in which particles are matched
        if frames is None:
            frames = self.time_lst
        
        cam_names = [cam.name for cam in self.imsys.cameras]
        
        
        # start matching, one frame at a time
        self.particles = []
        print('')
        for tm in frames:
            print('', end='\r')
            print(' frame: %d'%tm, end='\r')
            
            # set up a blobs dictionary with camera names as key
            pd = {}
            if self.reverse_eta_zeta:
                for i in range(len(self.blobs)):
                    cn = cam_names[i]
                    pd[cn] = self.blobs[i][self.blobs[i][:,-1] == tm][:,1::-1]
            
            else:
                for i in range(len(self.blobs)):
                    cn = cam_names[i]
                    pd[cn] = self.blobs[i][self.blobs[i][:,-1] == tm][:,:2]
            
            # match particles using the matching object
            M = matching(self.imsys, pd, self.RIO, self.voxel_size,
                         max_err = self.max_err)
            #return M  # <-- used for checks
            M.get_voxel_dictionary()
            M.list_candidates()
            M.get_particles()
            
            # extract the matched particles to the list self.particles 
            for p in M.matched_particles:
                self.particles.append(p + [tm])
        
        self.particles = list(filter(lambda p: p[4]<self.max_err,
                                     self.particles))
        
        print('\n','done!')                        
        errors = [p[4] for p in self.particles]
        print('mean error: %.3f'%(sum(errors)/len(errors)))
        Np = len(self.particles)
        times = [p[-1] for p in self.particles]
        Nframes = len(set(times))
        print('avg. particles in frame: %.2f'%(Np/Nframes))

    
    
    def save_results(self, fname):
        '''will save the list of particles obtained'''
        prticles_to_save = []
        Ncams = len(self.imsys.cameras)
        
        for p in self.particles:
            rd = dict(p[3])
            
            p_ = [p[0], p[1], p[2]]
            
            for i in range(Ncams):
                if i in list(rd.keys()):
                    p_.append(rd[i])
                else:
                    p_.append(-1)
            p_.append(p[4])
            p_.append(p[5])
            prticles_to_save.append(p_)
            
        fmt = ['%.3f', '%.3f', '%.3f']
        for i in range(Ncams):
            fmt.append('%d')
        fmt = fmt + ['%.3f', '%.3f']
        savetxt(fname, prticles_to_save, fmt=fmt, delimiter='\t')
        
            
                










class matching(object):
    '''A class for matching particles in images taken simultaniously
    from different cameras.
    
    The relevant functions for use are: 
        1) self.get_voxel_dictionary()
        2) self.list_candidates()
        3) self.get_particles()
    After running these three functions the attribute self.matched_particles
    holds the results of triangulation.
    '''
    
    
    def __init__(self, img_system, particles_dic, 
                 RIO, voxel_size, max_err=None):
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
        max_err - maximum allowable triangulation rms error.
        
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
        self.max_err = max_err
        
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


    
    def ray_traversed_voxels(self, ray):
        '''
        Given a ray, this will add the voxel through which it
        traverses into the list self.traversed_voxels .
        '''
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
            
            
            i0, i1 = max([0,i-1]), min([self.Nx-1,i+1])+1
            j0, j1 = max([0,j-1]), min([self.Ny-1,j+1])+1
            k0, k1 = max([0,k-1]), min([self.Nz-1,k+1])+1
            
            for i_ in range(i0,i1):
                for j_ in range(j0,j1):
                    for k_ in range(k0,k1):
                        self.traversed_voxels.append( [(i_,j_,k_), ray[2]] )
    
    
    def get_voxel_dictionary(self):
        '''This generates a dicionary who's keys are voxel indexes and
        who's values are the rays that passed through this voxel.'''
        
        self.traversed_voxels = []
        for ray in self.rays:
            self.ray_traversed_voxels(ray)
        
        
        voxel_dic = {}
        for vxl in self.traversed_voxels:
            try:
                voxel_dic[vxl[0]].add(vxl[1])
            except:
                voxel_dic[vxl[0]] = set([vxl[1]])
                
        self.voxel_dic = voxel_dic
        
    
    def list_candidates(self):
        '''This will make lists of possible candidate rays for
        triangulation, separated for pairs, triplets, quadruplets, etc.
        
        Candidates are based on the voxels of voxel_dic, while calculating the 
        RMS and the maximum distance between the estimated particle location 
        and the epipolar lines.'''
        
        candidate_dic = {}
        group_sizes = range(2, len(self.imsys.cameras)+1)
        for i in group_sizes:
            candidate_dic[i] = []
        
        for k in self.voxel_dic.keys():
            if len(self.voxel_dic[k]) >= 2:
                # make a nested list of the rays, by their camera number 
                ray_by_cams = [[] for i in range(len(self.imsys.cameras))]
                for ray in self.voxel_dic[k]:
                    ray_by_cams[ray[0]].append(ray)
                
                # find all possible combinations of the rays for various
                # numbers of cameras
                for gs in group_sizes:
                    for comb in combinations(ray_by_cams, gs):
                        candidate_dic[gs] += product(*comb)
        
        for k in candidate_dic.keys():
            candidate_dic[k] = list(set(candidate_dic[k]))
        
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
    
    
    
    def get_particles(self):
        '''Once all candidates are found, this function chooses the "best"
        matches and returns them. The reliability of the matches is considered 
        higher is
        1) they have higher number of cameras participating in the
        triangulation
        2) the RMS of distance between the crossing point and the epipolar 
        lines is smaller
        in this order. Thus, we choose the combinations of rays with highest
        number of camera participating and with the smallest RNS triangulaiton 
        error.
        '''
        matched_particles = []
        used_rays = []
        for k in sorted(self.candidate_dic.keys(), reverse=True):
            
            # triangulate all the candidate rays
            cand_k = self.candidate_dic[k]
            ray_crosses = [self.triangulate_rays(cand) for cand in cand_k]
            
            # zip and sort candidates by RMS error
            key = lambda x: x[1][2]
            dist_sorted_cands = sorted(zip(cand_k, ray_crosses), key = key)
            
            
            # if the rays in the candidate have not been used, adde the
            # particle to the results list 
            for i in range(len(dist_sorted_cands)):
                used_check = False
                for ray in dist_sorted_cands[i][0]:
                    if ray in used_rays: 
                        used_check = True
                        continue
                
                if not used_check:
                    p = dist_sorted_cands[i][1]
                    new_p = [round(p[0][0], ndigits=3), 
                             round(p[0][1], ndigits=3),
                             round(p[0][2], ndigits=3),
                             dist_sorted_cands[i][0],
                             round(p[-1], ndigits=3)]
                    matched_particles.append(new_p)
                    used_rays += dist_sorted_cands[i][0]
        
        
        #d_list = [p[4] for p in matched_particles]
        #print('')
        #print('Found %d particles'%(len(matched_particles)))
        #print('with maximum RMS error of %.2f'%(max(d_list)))
        self.matched_particles = matched_particles
        
        
        
    def plot_ray_epipolar_lines(self, ray, zlims, ax):
        '''will plot a ray's epipolar line for a given 3D axis.'''
        import matplotlib.pyplot as plt
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        camera = self.imsys.cameras[ray[0]]
        i = self.ray_camera_indexes[ray[0]] 
        ip1 = self.ray_camera_indexes[ray[0]+1]
        eta, zeta = self.rays[i:ip1][ray[1]][:2]
        camera.plot_3D_epipolar_line(eta, zeta, zlims, ax=ax, color=colors[ray[0]])
        
        
        



    
# if __name__ == '__main__':
#     from imaging_mod import camera, img_system
#     import numpy as np
#     import os
    
#     dirname ='/home/ron/Desktop/Research/plankton_sweeming/experiments/PTV_test3'
#     cnames = ['cam1', 'cam3', 'cam4']
#     res = 1280, 1024
#     cameras = []
#     for c in cnames:
#         cameras.append(camera(c, res))
    
#     for cam in cameras:
#         cam.load(dirname)
    
#     particles_dic = {}
#     for c in cnames:
#         blobs = np.loadtxt( os.path.join(dirname, 'blobs_'+c) )
#         particles_dic[c] = blobs[:,:2]
    
#     RIO = ((0.0, 68.0),
#            (0.0, 66.0),
#            (-30.0, 20.0))
    
#     voxel_size = 2.0
    
#     imsys = img_system(cameras)
    
#     m = matching(imsys, particles_dic, RIO, voxel_size)
#     m.get_voxel_dictionary()
#     m.list_candidates()
    

        
        
        
        
        
        
        