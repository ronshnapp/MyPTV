# -*- coding: utf-8 -*-
#!/usr/bin/env python3

"""
Created on April 2026

@author: ron

This script generates synthetic calibraiton images from calibrated
camera files. The images can then be used to calibrate other PTV
softwares. 
""" 

from tqdm import tqdm
from myptv.imaging_mod import camera_wrapper
from PIL import Image, ImageDraw
import numpy as np
from numpy.random import uniform



# ======================================
# Making up lists of calibration points:
# ======================================


x_lst = [0, 10, 30, 40, 50, 60, 70, 90, 100]  
y_lst = list(range(0,101,5))
z_lst = [0,-10,-30,-50,-42,-50,-30,-10, 0]

cal_points = []

for i in range(len(x_lst)):
    for y in y_lst:
        cal_points.append([x_lst[i], y, z_lst[i]])
        
            

# =====================
# Loading cameras:
# =====================

camNames = ['cam1', 'cam2', 'cam3', 'cam4']
cams = [camera_wrapper(cn, '.') for cn in camNames]
for cam in cams:
    cam.load()


# ======================
# Generating the images:
# ======================

resolution = 2048, 2048     # image pixel resolution
n_blobsPerParticle = 1      # number of blobs making up a particle
disparity = 0.00            # random noise of the blobs making up the particles [mm]
s = 4                       # sigma of a gaussian intensity point
I = 250.0                   # intensity of a fiber making point


radius = int(3*s)
x_ = range(-radius, radius+1)
X, Y = np.meshgrid(x_, x_)
blob_image = I * np.exp(-( ((X)**2 + (Y)**2) /2/s**2))


images = [np.zeros((resolution[1], resolution[0])) for cam in cams]
for e, cam in enumerate(cams):
    img = images[e]
    
    for X in cal_points:
        cx, cy = np.array(cam.projection(X)).astype('int')
        i0, i1, j0, j1 = cy-radius, cy+radius+1, cx-radius, cx+radius+1
        
        tst0 = (np.array([i0, i1, j0, j1]) < 0).any()
        tst1 = (np.array([i0, i1]) > resolution[0]).any()
        tst2 = (np.array([j0, j1]) > resolution[1]).any()
        if  tst0 or tst1 or tst2:
            continue
        img[i0:i1, j0:j1] += blob_image
        
       
    for e, cam in enumerate(cams):
        img = images[e]
        img[img>255] = 255
        img = img.astype('int8')

        pil_im = Image.fromarray(img, mode='L')
        pil_im.save('./cam%d_cal_im.tif'%(e+1))
    
    
    np.savetxt('./cal_points', cal_points, delimiter='\t', fmt='%.3f')