# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 18:02:07 2018

@author: ron


contains a class for segmentation of circular particles
"""

from numpy import zeros, savetxt
from scipy.ndimage import gaussian_filter
from skimage.io import imread
from numpy import ones
from scipy.signal import convolve2d
from scipy.ndimage.measurements import label, find_objects


class particle_segmentation(object):
    '''a class for segmenting out particles (blobs) for a given image'''
    
    
    def __init__(self, image, sigma=1.0, threshold=10, mask=1.0,
                 local_filter = 15,
                 min_xsize=None, max_xsize=None,
                 min_ysize=None, max_ysize=None,
                 min_area=None, max_area=None):
        
        self.im = image
        self.sigma = sigma
        self.th = threshold
        self.mask = mask
        self.bbox_limits = (min_xsize, max_xsize, min_ysize, max_ysize)
        self.area_limits = (min_area, max_area)
        self.loc_filter = local_filter
        
    
    def local_filter(self, image):
        '''returns a new image where the local mean neighbourhood of
        each pixel is subtracted.'''
        w = self.loc_filter
        window = ones((w, w)) / w**2
        local_mean = convolve2d(image, window, mode='same')
        new_im = image - local_mean
        new_im[new_im<0] = 0
        new_im = new_im.astype('uint8')
        return new_im
        
        
    def get_binary_image(self):
        '''Will mark pixels in the image as background and foreground 
        (particles). We blur the image with a Gaussian
        filter, and look for regions brighter than a global threshold 
        level.'''
        
        # apply a Gussian blur
        if self.sigma is not None:
            blured = gaussian_filter(self.im, self.sigma)
        else:
            blured = self.im
            
        # apply local mean subtraction
        if self.loc_filter is not None:
            filtered = self.local_filter(blured)
        else:
            filtered = blured
        
        # find local maxima and generate a binary image
        global_filt = filtered>self.th
        bin_image = 1.0 * global_filt * self.mask
        
        return bin_image 
            
        
        # find local maxima and generate a binary image
        #from scipy.ndimage import grey_dilation
        #dilated = grey_dilation(filtered, size=self.p_size, mode='constant')
        #bin_image = (filtered==dilated) * (filtered>self.th) * self.mask
        
        #return bin_image
    

    
    def blob_labeling(self, image):
        '''Will label connected areas (blobs) in a binary image and return
        their coordinates. The values of image are 0 for background and
        1 for foreground
        
        output - linked: a nested list of connected pixel indexes
        '''
        # use scipy to label the blobs
        self.labeled = label(image)[0]
        locations = find_objects(self.labeled)
        
        return locations
    
    
    
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
        
        from numpy import meshgrid
        from numpy import sum as npsum
        
        stamp_y, stamp_x = meshgrid(range(self.im.shape[1]), 
                                    range(self.im.shape[0]))
        for loc in self.blob_pixels:
            mask = 1.0*(self.labeled[loc]>0)
            tot = npsum(self.im[loc] * mask)
            X = npsum(stamp_x[loc] * self.im[loc] * mask) / tot
            Y = npsum(stamp_y[loc] * self.im[loc] * mask) / tot
            center = [round(X, ndigits=2), round(Y, ndigits=2)]
            area = npsum(mask)
            box_size = list(mask.shape)
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
            
            
    def plot_blobs(self, vmin=None, vmax=None):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.imshow(self.im, vmin=vmin, vmax=vmax)
        
        for blb in self.blobs:
            
            ax.errorbar( [blb[0][1]], [blb[0][0]], 
                        xerr=blb[1][1]/2, yerr=blb[1][0]/2,
                        fmt='xr', lw=0.7, capsize=2)
        
        
    def save_results(self, fname):
        '''
        This is used to save the blobs found in a text file with 
        the given name fname.
        '''
        blob_list = []
        for blb in self.blobs:
            blob_list.append([blb[0][0], blb[0][1], blb[1][0], blb[1][1],
                              blb[2], 0])
            
        savetxt(fname, blob_list, 
                fmt=['%.02f','%.02f','%d','%d','%d','%d'], delimiter='\t')
        
        
        
        

        
        
        
class loop_segmentation(object):
    
    '''A class for looping over images in a library to segment particles
    and save the results in a file.'''
    
    def __init__(self, dir_name, extension='.tif',
                 N_img = None, sigma=1.0, threshold=10, mask=1.0,
                 local_filter = 15,
                 min_xsize=None, max_xsize=None,
                 min_ysize=None, max_ysize=None,
                 min_area=None, max_area=None):
        '''
        dir_name - string with the name of the directory that holds the 
                   images. Images should have a sequential numbers in their
                   file names. 
        extension - the extension of the images
        
        N_img -     if None, then this will loop over all the images in the 
                    folder. If it is an integer, will loop over the first
                    N images in the folder.
                    
        The rest are parameters for the segmentation class. 
        '''
        self.dir_name = dir_name
        self.extension = extension
        self.N_img = N_img
        self.sigma = sigma
        self.th = threshold
        self.mask = mask
        self.bbox_limits = (min_xsize, max_xsize, min_ysize, max_ysize)
        self.area_limits = (min_area, max_area)
        self.loc_filter = local_filter
    
    
    def get_file_names(self):
        import os
        allfiles = os.listdir(self.dir_name)
        n_ext = len(self.extension)
        fltr = lambda s: s[-n_ext:]==self.extension
        image_files = sorted(list(filter(fltr, allfiles)))
        self.image_files = image_files
    
    
    def segment_folder_images(self):
        '''This loops over the image files in a folder'''
        import os
        
        self.get_file_names()
        
        if self.N_img is None: 
            N = len(self.image_files)
        else:
            N = self.N_img
        
        blob_list = []
        print('Starting loop segmentation.')
        for i in range(N):
            print('', end='\r')
            print(' frame: %d'%i, end='\r')
            im = imread(os.path.join(self.dir_name, self.image_files[i]))
            ps = particle_segmentation(im,
                                       sigma=self.sigma, 
                                       threshold=self.th,
                                       local_filter=self.loc_filter,
                                       mask=self.mask,
                                       max_xsize=self.bbox_limits[1],
                                       min_xsize=self.bbox_limits[0],
                                       max_ysize=self.bbox_limits[3],
                                       min_ysize=self.bbox_limits[2],
                                       min_area=self.area_limits[0],
                                       max_area=self.area_limits[1])
            ps.get_blobs()
            ps.apply_blobs_size_filter()
            for blb in ps.blobs:
                blob_list.append([blb[0][0], blb[0][1], blb[1][0], blb[1][1],
                                  blb[2], i])
        self.blobs = blob_list
        
                                       
    def save_results(self, fname):
        '''
        Will save the extracted blobs. 
        
        The format of the results is
        center_x, center_y, size_x, size_y, area, frame_number
        '''
        savetxt(fname, self.blobs, 
                fmt=['%.02f','%.02f','%d','%d','%d','%d'], delimiter='\t')
        
        
        
    
    



