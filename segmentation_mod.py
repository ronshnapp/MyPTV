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
    
    
    def __init__(self, image, sigma=None, threshold=10, mask=1.0,
                 local_winsize=None, local_threshold=None):
        
        self.im = image
        self.sigma = sigma
        self.th = threshold
        self.lws = local_winsize
        self.lth = local_threshold
        self.bin_im = self.get_binary_image() 
        self.blob_pixels = self.blob_labeling(self.bin_im)
    
        
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

            
        bin_image = 1.0 * global_filt
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
    
    
    
    def local_threshold(self, image):
        '''applies a local threshold to see if each pixel is higher
        than some threshold above the local average.'''
        
        if self.lws is None:
            raise ValueError('local window size not defined')
            
        if self.lth is None:
            raise ValueError('local threshold size not defined')
        
        wi = int(self.lws/2)
        result = zeros(self.im.shape)
        for i in range(wi,self.im.shape[0]-wi):
            for j in range(wi,self.im.shape[1]-wi):
                if image[i,j]>0:
                    local_av = mean(image[i-wi:i+wi+1,j-wi:j+wi+1])
                    if image[i,j] > local_av + self.lth:
                        result[i,j] = 1
        return result


        
    
    
    
    



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



