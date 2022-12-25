# -*- coding: utf-8 -*-
"""
Created on Dec 20 2022

@author: ron


Contains a class for tracking particles to in 2D.

"""

from myptv.tracking_mod import tracker_four_frames

from numpy import loadtxt, array, savetxt
from scipy.spatial import KDTree
from pandas import read_csv



class track_2D(tracker_four_frames):
    '''
    A class used to perform particle tracking in 2D using images from a single
    camera.
    
    The methodology is developed so that we take images from a calibrated 
    camera, and track particles in the x and y coordinates, assuming they all 
    have the same (given as input) z coordinate value. 
    
    To do the tracking, we inherit from the tracker_four_frames class of 3D 
    tracking.
    
    how to use - 
    
    1) initiate the class with the appropriate inputs
    2) run the function self.blobs_to_particles() which transforms the blobs 
       into particles with lab-space coordinates and at the needed format
    3) run the inherited method self.track_all_frames() to track in 2D
    4) save the results on the disk using the inherited method 
       self.save_results(save_name)
    '''
    
    def __init__(self, camera, blob_fnames, z_particles, mean_flow = 0.0, 
                 d_max=1e10, dv_max=1e10, reverse_eta_zeta = False):
        '''
        inputs -
        
        camera - a calibrated camera instance.
        
        blob_fname -  the file name that containes the coordinates of 
                      segmented blobs.
                      
        mean_flow - a numpy array of the mean flow vector, in units of the 
        calibrations spatial units per frame (e.g. mm per frame). The mean 
        flow is assumed not to change in space and time.
        
        d_max - maximum allowable translation between two frames for the 
                nearest neighbour search, after subtracting the mean flow. This
                is in lab-space coordinates (e.g. millimeters)
                
        dv_max - maximum allowable change in velocity for the two-frame 
                 velocity projection search. The radius around the projection
                 is therefore dv_max/dt (where dt = 1 frame^{-1}) in lab-space
                 units.
                 
        reverse_eta_zeta - Should be false if the eta and zeta coordinates 
                           need to be in reverse order so as to match the
                           calibration. This may be needed if the calibration 
                           data points were given where the x, y coordinates
                           are transposed (as happens, e.g., if using 
                           matplotlib.pyplot.imshow).
        '''
        
        self.cam = camera
        self.fname = blob_fnames
        self.z_particles = z_particles
        self.U = mean_flow
        self.d_max = d_max
        self.dv_max = dv_max
        self.reverse_eta_zeta = reverse_eta_zeta
        
        self.blobs = dict([(k, array(g)) for k,g in 
                           read_csv(self.fname,header=None,sep='\t').groupby(5)])
        
        self.times = sorted(list(self.blobs.keys()))
        
        self.trees = {}
        
        self.traj_ids = []
        self.traj_lengths = {}
        self.N_four_frames = 0
        self.N_nearest_neighbour = 0
    
    
    
    def transform_coords(self, eta, zeta):
        '''
        This function transforms an image-space coordinate of a blob and to
        the lab-space coordinate of the particle. For the transofrmation,
        we assume that the particle has the z coordinate given as an input.
        
        input -
        
        eta - pixel "x" coordinate in camera space
        zeta - pixel "y" coordinate in camera space
        '''
        if self.reverse_eta_zeta==False:
            r = self.cam.get_r(eta, zeta)
        
        else:
            r = self.cam.get_r(zeta, eta)
        
        O = self.cam.O
        a = (self.z_particles - O[2])/r[2]
        
        x, y = O[:2]+r[:2]*a
        return x, y
    
    
    
    def blobs_to_particles(self):
        '''
        This function uses the blob data to generate a dicionary of particles 
        in the format used by tracker_four_frames. 
        '''
        self.particles = {}
        
        for k in self.blobs.keys():
            self.particles[k] = []
            
            for i in range(len(self.blobs[k])):
                blob = self.blobs[k][i]
                x, y = self.transform_coords(blob[0], blob[1])
                p = array([-1, x, y, self.z_particles, i, 0.0, blob[-1]])
                self.particles[k].append(p)
                
            self.particles[k] = array(self.particles[k])
        

        
    






    


if __name__=='__main__':
    
    from imaging_mod import camera
    
    fname = '/home/ron/working_PTV_data/blobs_cam1'
    
    cam = camera('cam1', (1024,1280))
    cam.load('/home/ron/working_PTV_data')
    
    z = 10.0
    d_max = 1.0
    dv_max = 0.25
    
    t2d = track_2D(cam, fname, z, d_max=d_max, dv_max = dv_max, 
                   reverse_eta_zeta=True)
    
    t2d.blobs_to_particles()


