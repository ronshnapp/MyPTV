#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 18:02:07 2018

@author: ron

A test function for the calibration module

"""

from myptv.imaging_mod import camera, img_system
from myptv.particle_matching_mod import match_blob_files



def test_matching():
    '''
    A test for the matching class. We attempt to matching a bunch of 
    synthetic particles using synthetic cameras.
    '''
    cam1 = camera('matching_test_cam1', (1280,1024))
    cam1.load('./tests/matching_test_files/')
    
    cam2 = camera('matching_test_cam2', (1280,1024))
    cam2.load('./tests/matching_test_files/')
    
    cam3 = camera('matching_test_cam3', (1280,1024))
    cam3.load('./tests/matching_test_files/')
    
    imsys = img_system([cam1, cam2, cam3])
    
    blob_files = ['./tests/matching_test_files/matching_test_blobs1', 
                  './tests/matching_test_files/matching_test_blobs2', 
                  './tests/matching_test_files/matching_test_blobs3']
    
    ROI = ((-20, 20),
           (-20, 20),
           (-20, 20))
    
    voxel_size = 40.0
    max_blob_dist = 0.0
    
    match = match_blob_files(blob_files, imsys, ROI, voxel_size, max_blob_dist)
    match.get_particles()
    particles = match.particles
    
    test_n_found = len(particles) == 3
    
    x1 = [0.0, 0.0, 0.0]
    p1_x = particles[0][:3]
    err_1 = sum([(x1[i]-p1_x[i])**2 for i in range(3) ])**0.5
    test_p1 = err_1 < 0.1
    assert test_p1 and test_n_found


