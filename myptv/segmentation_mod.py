# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 18:02:07 2018

@author: ron


contains a class for segmentation of circular particles
"""

from numpy import ones, savetxt, meshgrid
from numpy import sum as npsum
from numpy import abs as npabs
from numpy import median as npmedian

from skimage.io import imread
from skimage import io

from scipy.signal import convolve2d
from scipy.ndimage import gaussian_filter, median_filter
from scipy.ndimage.measurements import label, find_objects
from scipy.spatial import KDTree



class particle_segmentation(object):
    '''a class for segmenting out particles (blobs) for a given image'''
    
    
    def __init__(self, image, sigma=None, threshold=10, mask=1.0,
                 median = None, local_filter = None, particle_size=3,
                 BG_image = None,
                 min_xsize=None, max_xsize=None,
                 min_ysize=None, max_ysize=None,
                 min_mass=None, max_mass=None,
                 method='labeling'):
        '''
        inputs - 
        
        image - the image (matrix) to be analyzed
        
        particle_size - the expected size of particles; this input is only
                        used if method='dilation'.
        
        sigma - float, standard deviation of the gaussian blur filter. If this 
                is None, then the bluring filter is not applied.
        
        threshold - int, the grey value above which pixels are considered a 
                    particle. Note that the threshold is performed after the
                    blur and local mean subtraction filters are applied.
        
        mask - either array of 0 and 1 with the same shape of imae, or 1. If 
               this is an array of 0 and 1, this is taken to choose regions 
               of interest in the analyzed image (regions with 0 are 
               discarded).
        
        local_filter - int, the window size of the local mean subtraction 
                       filter. In this is None then the filter is not applied.
        
        BG_image - Either a given static image, in which case this image is
                   subtracted from the given image by taking the difference,
                   or this is None, in which case this operation is skipped.
                       
        min/max_( ) - size and area filters for the discovered blobs. If None
                      then filters are not applied.
                      
        method - string. The method used for labeling the blobds. Can be 
                 either 'dilation' or 'labeling'. In dilation, local maxima 
                 are sought in the image and the blobs are considered to be
                 regions of size 'particle_size' around these maxima. In the 
                 'labeling' method, the threshold is performed and blobs are
                 considered to be connected bright pixels. The dilation method
                 is better at distinguishing close particles (saddle points in 
                 the brihtness) if their size is given correctly. The 
                 'labeling' method may be better for non spherical particles.
        '''
        
        self.im = image
        self.p_size = particle_size
        self.sigma = sigma
        self.median = median
        self.th = threshold
        self.mask = mask
        self.bbox_limits = (min_xsize, max_xsize, min_ysize, max_ysize)
        self.mass_limits = (min_mass, max_mass)
        self.loc_filter = local_filter
        self.BG_image = BG_image
        
        if method not in ['dilation', 'labeling']:
            raise ValueError('method "%s" unknown.'%method)
        
        self.method = method
        
        

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
        
    
    
    def process_image(self):
        '''
        Will perform image processing - the gaussian blur and then the loal
        mean subtraction; processing is performed only if the filter parameters
        are not None. Finlly it will apply the mask to cover unwanted regions.
        Results are stored in the self.processed_im.
        '''
        
        # subtract background:
        if self.BG_image is not None:
            imNoBG = npabs(self.im - self.BG_image).astype(self.im.dtype)
        else:
            imNoBG = self.im
        
        # apply a Gussian blur
        if self.sigma is not None:
            blured = gaussian_filter(imNoBG, self.sigma)
        else:
            blured = imNoBG
            
        # apply median filter
        if self.median is not None:
            med = median_filter(blured, size = self.median)
        else:
            med = blured
            
        # apply local mean subtraction
        if self.loc_filter is not None:
            filtered = self.local_filter(med)
        else:
            filtered = med
        
        # apply the mask
        self.processed_im = filtered * self.mask
    
        
    
    
    def get_binary_image(self):
        '''Will mark pixels in the image as background and foreground 
        (particles) using the given method'''
        
        # Do the image processing
        self.process_image()
        
        if self.method=='labeling':
            # find bright regions and generate a binary image
            global_filt = self.processed_im > self.th
            bin_image = 1.0 * global_filt * self.mask
            return bin_image 
                
        elif self.method == 'dilation':
            # dilate, find briht local maxima, and generate a binary image
            from scipy.ndimage import grey_dilation
            dilated = grey_dilation(self.processed_im, 
                                    size=self.p_size, mode='constant')
            bin_im = (self.processed_im==dilated) * (self.processed_im>self.th) 
            
            return bin_im
        
    
    
    def characterize_blob(self, coord, size=None):
        '''
        Returns the characterization of a blob centered around coord in the
        dilation method.
        
        Its' center is defined as the brightness weighted mean of the 
        neighbourhood of size 'size' around the given coordinates 'coord'. Its
        mass is the sum of brighness values inside the blob's pixels. Its bbox
        is the smallest lenghts in x and y that that contain pixels with 
        brightness above the threshold brightness value.
        
        input -
        coord - an tuple or list of size 2 with the x,y coordinates
        size - if None, size is set to be self.particle_size; otherwise it 
               shoud be an integer.
        '''
        
        if size is None:
            size = self.p_size
        else:
            if type(size) != int:
                raise ValueError('Size should be an integer')
        
        # prepare slices to work with the blob neighbourhood
        r = size//2
        
        if coord[0]<r:
            loc_x = slice(0, int(coord[0]+r+1))
        else:
            loc_x = slice(int(coord[0]-r), int(coord[0]+r+1))    
        
        if coord[1]<r:
            loc_y = slice(0, int(coord[1]+r+1))
        else:
            loc_y = slice(int(coord[1]-r), int(coord[1]+r+1))    
            
        # calculate the mass
        mass = npsum(self.processed_im[loc_x,loc_y])
        
        
        if mass == 0: 
            print(coord)
            print(self.processed_im[loc_x,loc_y])
            print(self.im[loc_x,loc_y])
        
        
        # calculate the center of mass
        c = (npsum(self.X[loc_x,loc_y] * self.processed_im[loc_x,loc_y])/mass,
             npsum(self.Y[loc_x,loc_y] * self.processed_im[loc_x,loc_y])/mass)
        
        # calculate the bounding box
        reion_x = self.X[loc_x,loc_y][self.processed_im[loc_x,loc_y] > self.th]
        reion_y = self.Y[loc_x,loc_y][self.processed_im[loc_x,loc_y] > self.th]
        
        try:
            bbox = [max(reion_x) - min(reion_x) + 1, 
                    max(reion_y) - min(reion_y) + 1]
        except:
            # if no pixels above threshold were found, put (0,0). This should
            # raise a red flag.
            bbox = (0, 0)
        
        return c, bbox, mass
        
    
    
    def blob_labeling(self, image):
        '''
        Will label connected areas (blobs) for the labeling method. 
        
        Labels regions with value 1 in a binary image and return their 
        coordinates. The values of image are 0 for background and 1 for 
        foreground.
        
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
        The box size is the largest length that bound the blob in the
        x and y directions.
        The area is the number of pixels belonging to the blob
        
        returns - blobs: a nested list of [ [(center), (box size), area], ...]
        '''    

        if self.method == 'dilation':
                
            # get a list of pixel coordinates that are bright local maxima
            self.bin_im = self.get_binary_image() 
            self.Y, self.X = meshgrid(range(self.im.shape[1]), 
                                      range(self.im.shape[0]))
            coords = list(zip(self.X[self.bin_im>0], self.Y[self.bin_im>0]))
                
            blobs = []
            for coord in coords:
                
                # perform iterations (maximum of 3) to refine particles' position
                for i in range(3) :
                    C, bbox, mass = self.characterize_blob(coord)
                    d = ((C[0]-coord[0])**2 + (C[1]-coord[1])**2)**0.5 
                    coord = C
                    
                    if d < 1.0:
                        break
                
                # round the center of mass
                coord = [round(coord[0], ndigits=2), 
                         round(coord[1], ndigits=2)]
                
                # add final blob to final list
                blobs.append( [coord, bbox, mass] )
                
            
            # search and remove duplicates; duplicates are points that are 
            # closer than self.particle_size/2 away. In this case, we keep the
            # blob with lower mass, 
            if len(blobs) > 0: 
                tree = KDTree([b[0] for b in blobs])
                duplicates = tree.query_pairs(self.p_size/2)
                to_remove = []
                for d in duplicates:
                    if blobs[d[0]][-1] < blobs[d[1]][-1]:
                        to_remove.append(d[1])
                    else:
                        to_remove.append(d[0])
                        
                to_remove = list(set(to_remove))
                        
                for i in sorted(to_remove, reverse=True): 
                    del blobs[i] 
                
            self.blobs = blobs
            
            
            
        elif self.method=='labeling':
            
            self.bin_im = self.get_binary_image() 
            blob_pixels = self.blob_labeling(self.bin_im)
            
            blobs = []
            
            stamp_y, stamp_x = meshgrid(range(self.im.shape[1]), 
                                        range(self.im.shape[0]))
            for loc in blob_pixels:
                mask = 1.0*(self.labeled[loc]>0)
                mass = npsum(self.processed_im[loc] * mask)
                X = npsum(stamp_x[loc] * self.processed_im[loc] * mask) / mass
                Y = npsum(stamp_y[loc] * self.processed_im[loc] * mask) / mass
                center = [round(X, ndigits=2), round(Y, ndigits=2)]
                box_size = list(mask.shape)
                blobs.append( [center, box_size, mass])
                
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
            
        if self.mass_limits[0] is not None:
            fltr = lambda b: b[2] > self.mass_limits[0]
            self.blobs = list(filter(fltr, self.blobs))
        
        if self.mass_limits[1] is not None:
            fltr = lambda b: b[2] < self.mass_limits[1]
            self.blobs = list(filter(fltr, self.blobs))
            
            
            
    def plot_blobs(self, vmin=None, vmax=None):
        import matplotlib.pyplot as plt
        
        if vmax is None:
            vmax = min([self.th*2, max(self.im.ravel())])
        
        fig, ax = plt.subplots()
        ax.imshow(self.processed_im, vmin=vmin, vmax=vmax)
        
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
    
    '''A class for looping over images in a directory to segment particles
    and save the results in a file.'''
    
    def __init__(self, dir_name, extension='.tif',
                 image_start = 0,
                 N_img = None, sigma=1.0, threshold=10, mask=1.0,
                 local_filter = 15, median = None, particle_size=3,
                 remove_ststic_BG = True,
                 min_xsize=None, max_xsize=None,
                 min_ysize=None, max_ysize=None,
                 min_mass=None, max_mass=None,
                 method='labeling'):
        '''
        dir_name - string with the name of the directory that holds the 
                   images. Images should have a sequential numbers in their
                   file names. 
        extension - the extension of the images
        
        image_start - The number of image from which the loop begins. If None, 
                      the loop shall begin from the first image in the folder.
        
        N_img -     if None, then this will loop over all the images in the 
                    folder. If it is an integer, will loop over the first
                    N images in the folder.
                    
        remove_ststic_BG - If true, a "background image" is calculated for the
                           images, defined as the median value of each pixel,
                           and then background_subtraction is done by taking 
                           the difference from the median. If this is False,
                           then this operation is skipped.
                    
        The rest are parameters for the segmentation class. 
        '''
        self.dir_name = dir_name
        self.p_size = particle_size
        self.extension = extension
        self.image_start = image_start
        self.N_img = N_img
        self.sigma = sigma
        self.median = median
        self.th = threshold
        self.mask = mask
        self.bbox_limits = (min_xsize, max_xsize, min_ysize, max_ysize)
        self.mass_limits = (min_mass, max_mass)
        self.loc_filter = local_filter
        self.method = method
        self.BG_remove = remove_ststic_BG
    
    
    def get_file_names(self):
        import os
        allfiles = os.listdir(self.dir_name)
        n_ext = len(self.extension)
        fltr = lambda s: s[-n_ext:]==self.extension
        image_files = sorted(list(filter(fltr, allfiles)))
        
        if self.image_start is not None:
            try:
                image_files = image_files[self.image_start:]
            except:
                raise ValueError('Image start is a positive integer or None')
        
        self.image_files = image_files
    
    
    def calculate_BG(self):
        '''
        Calculates the background image, defined as the median of the given
        images. If there are more than 200 images, then the background is done
        using a subsample of 200 images, to keep the calculation from taking
        too long.
        '''
        import os
        
        print('calculating background...\n')
        # getting the subsample of images for BG calculation
        if len(self.image_files)<=200:
            BG_images = self.image_files
            
        else: 
            BG_images=self.image_files[::int(len(self.image_files)/400+1)][:200]
        
        BG_images = [os.path.join(self.dir_name, im) for im in self.image_files]
        # reading the images for BG subtraction
        ic = io.ImageCollection(BG_images)
        
        self.BG = npmedian(ic, axis=0)
        
        
    def segment_folder_images(self):
        '''This loops over the image files in a folder'''
        import os
        
        self.get_file_names()
        
        if self.N_img is None: 
            N = len(self.image_files)
        else:
            N = self.N_img
            
        if self.BG_remove==True:
            self.calculate_BG()
        else:
            self.BG = None
        
        i0 = (self.image_start is not None) * self.image_start        
        
        blob_list = []
        print('Starting loop segmentation.\n')
        for i in range(N):
            print('', end='\r')
            print(' frame: %d'%(i+i0), end='\r')
            im = imread(os.path.join(self.dir_name, self.image_files[i]))
            ps = particle_segmentation(im,
                                       sigma=self.sigma, 
                                       threshold=self.th,
                                       median=self.median,
                                       local_filter=self.loc_filter,
                                       BG_image=self.BG,
                                       mask=self.mask,
                                       max_xsize=self.bbox_limits[1],
                                       min_xsize=self.bbox_limits[0],
                                       max_ysize=self.bbox_limits[3],
                                       min_ysize=self.bbox_limits[2],
                                       min_mass=self.mass_limits[0],
                                       max_mass=self.mass_limits[1],
                                       method = self.method,
                                       particle_size=self.p_size)
            ps.get_blobs()
            ps.apply_blobs_size_filter()
            for blb in ps.blobs:
                blob_list.append([blb[0][0], blb[0][1], blb[1][0], blb[1][1],
                                  blb[2], i+i0])
        self.blobs = blob_list
        
                                       
    def save_results(self, fname):
        '''
        Will save the extracted blobs. 
        
        The format of the results is
        center_x, center_y, size_x, size_y, area, frame_number
        '''
        savetxt(fname, self.blobs, 
                fmt=['%.02f','%.02f','%d','%d','%d','%d'], delimiter='\t')
        
        
        
    
    



