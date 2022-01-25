#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 18:02:07 2018

@author: ron


contains a class for segmentation of circular particles
"""

from scipy.signal import convolve2d
from numpy import zeros, arange, exp, gradient, mean
from scipy.ndimage import gaussian_filter




class particle_segmentation(object):
    '''a class for segmenting out particles (blobs) for a given image'''
    
    
    def __init__(self, image, sigma=1.0, threshold=10, mask=1.0,
                 min_xsize=None, max_xsize=None,
                 min_ysize=None, max_ysize=None,
                 min_area=None, max_area=None):
        
        self.im = image
        self.sigma = sigma
        self.th = threshold
        self.mask = mask
        self.bbox_limits = (min_xsize, max_xsize, min_ysize, max_ysize)
        self.area_limits = (min_area, max_area)
        
        
    def get_binary_image(self):
        '''Will mark pixels in the image as background and foreground 
        (particles). We blur the image with a Gaussian
        filter, and look for regions brighter than a global threshold 
        level.'''
        
        if self.sigma is not None:
            blured = gaussian_filter(self.im, self.sigma)
        else:
            blured = self.im
            
        global_filt = blured>self.th
            
        bin_image = 1.0 * global_filt * self.mask
        return bin_image 
    

    
    def blob_labeling(self, image):
        '''Will label connected areas (blobs) in a binary image and return
        these blobs coordinates. The values of image are 0 for background and
        1 for foreground
        
        output - linked: a nested list of connected pixel indexes
        '''
        
        nrow, ncol = image.shape
        
        labeled = zeros((nrow, ncol))
        linked = []
        
        for i in range(1, nrow-1):
            for j in range(1, ncol-1):
                
                if image[i,j]==1 and labeled[i,j]==0:                    
                    linked.append([(i,j)])
                    labeled[i,j] = 1
                    que = [(i,j)]
                    
                    for pixel in que:
                        mn0, mx0 = max([0,pixel[0]-1]), min([nrow, pixel[0]+2])
                        mn1, mx1 = max([0,pixel[1]-1]), min([ncol, pixel[1]+2])
                        for i2 in range(mn0, mx0):
                            for j2 in range(mn1, mx1):
                                if i2==i and j2==j:continue
                                if image[i2,j2]==1 and labeled[i2,j2]==0:
                                    que.append((i2,j2))
                                    linked[-1].append((i2,j2))
                                    labeled[i2,j2]=1
                        
        return linked
    
    
    
    def get_blobs(self):
        '''Returns a list of particle centers, their box size, and area
        
        The center is the weighted mean of the blob coordinates using
        the brightness as weights.
        The box size is the larger lengths than bound the blob in the
        x and y directions.
        The area is the number of pixels belonging to the blob
        
        returns - blobs: a nested list of [ [(center), (box size), area], ... ]
        '''
        
        self.bin_im = self.get_binary_image() 
        self.blob_pixels = self.blob_labeling(self.bin_im)
        
        blobs = []
        
        for i in range(len(self.blob_pixels)):
            
            X = 0.0
            Y = 0.0
            tot = 0.0
            
            for x,y in self.blob_pixels[i]:
                
                if tot == 0.0:
                    xmin, xmax = x, x
                    ymin, ymax = y, y
                else:
                    if x<xmin: xmin=x
                    if x>xmax: xmax=x
                    if y<ymin: ymin=y
                    if y>ymax: ymax=y
                
                bxy = self.im[x,y]
                X += x*bxy
                Y += y*bxy
                tot += bxy
            center = [X/tot, Y/tot]
            box_size = [xmax-xmin+1, ymax-ymin+1]
            area = len(self.blob_pixels[i])
            blobs.append( [center, box_size, area])
            
        self.blobs = blobs
        
        
    def apply_blobs_size_filter(self):
        '''Will filter the list of blobs accoring to their bounding box size 
        and their area.'''
        
        if self.bbox_limits[0] is not None:
            fltr = lambda b: b[1][0] > self.bbox_limits[0]
            self.blobs = list(filter(fltr, self.blobs))
        
        if self.bbox_limits[1] is not None:
            fltr = lambda b: b[1][0] < self.bbox_limits[1]
            self.blobs = list(filter(fltr, self.blobs))
        
        if self.bbox_limits[2] is not None:
            fltr = lambda b: b[1][1] > self.bbox_limits[2]
            self.blobs = list(filter(fltr, self.blobs))
        
        if self.bbox_limits[3] is not None:
            fltr = lambda b: b[1][1] < self.bbox_limits[3]
            self.blobs = list(filter(fltr, self.blobs))
            
        if self.area_limits[0] is not None:
            fltr = lambda b: b[2] > self.area_limits[0]
        
        if self.area_limits[1] is not None:
            fltr = lambda b: b[2] < self.area_limits[1]
            
            
    def plot_blobs(self):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.imshow(self.im, vmin=self.th/3, vmax=3*self.th)
        
        for blb in self.blobs:
            
            ax.errorbar( [blb[0][1]], [blb[0][0]], 
                        xerr=blb[1][1]/2, yerr=blb[1][0]/2,
                        fmt='xr', lw=0.7, capsize=2)
        
    
    
    
        
    
    
    
    



'''
import os
import matplotlib.pyplot as plt
dir_name = '/home/ron/Desktop/Research/plankton_sweeming/experiments/PTV_test3/CoreView_71/Cam1'
images = sorted(os.listdir(dir_name))


av_img = 0
count = 0
for i in range(len(images)):
    av_img += plt.imread(os.path.join(dir_name, images[i]))
    count += 1
av_img = av_img/count
'''



