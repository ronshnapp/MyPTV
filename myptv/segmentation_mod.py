# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 18:02:07 2018

@author: ron


contains a class for segmentation of circular particles
"""

from myptv.tracking_2D_mod import track_2D_multiframe
from myptv.TsaiModel.camera import camera_Tsai
from myptv.tracking_mod import fill_in_trajectory

from pandas import read_csv

from numpy import ones, savetxt, meshgrid, float32
from numpy import sum as npsum
from numpy import abs as npabs
from numpy import median as npmedian
from numpy import array, divide, zeros_like
from numpy import append as npappend

from skimage.io import imread
from skimage import io

from scipy.signal import convolve2d
from scipy.ndimage import gaussian_filter, median_filter
from scipy.ndimage.measurements import label, find_objects
from scipy.spatial import KDTree

import tqdm
import os
from time import time


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
        # w = self.loc_filter
        # window = ones((w, w)) / w**2
        # local_mean = convolve2d(image, window, mode='same')
        # new_im = image - local_mean
        # new_im[new_im<0] = 0
        # new_im = new_im.astype('uint8')
        
        S = self.loc_filter
        flt = image.astype(float32) #/ 255.0
        blur = gaussian_filter(flt, S)
        num = flt - blur
        
        blur = gaussian_filter(num*num, S)
        den = blur**0.5        
        # normed = num / den
        # to ensure no accidental zero division 
        normed = divide(num, den, out = zeros_like(num), where = (den != 0.0))

        return image * (normed>1)
        
    
    
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
            
            # getting the binary image
            self.bin_im = self.get_binary_image() 
            
            # labeling connected foreground pixels to form "blobs"
            blob_pixels = self.blob_labeling(self.bin_im)
            
            blobs = []
            
            stamp_y, stamp_x = meshgrid(range(self.im.shape[1]), 
                                        range(self.im.shape[0]))
            
            for e, loc in enumerate(blob_pixels):
                # extracting blob parameters
                mask = 1.0*(self.labeled[loc]>0)*(self.labeled[loc]==e+1)
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
                 method='labeling',
                 raw_format=False):
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
                           then this operation is skipped. If this is a numpy
                           array, than it is used as the BG image.
                           
        raw_format (default=Flase) - Set to true if the images are in raw 
                                     format (e.g. .dng). Then, we use the
                                     package rawpy to load the images.
        
                    
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
        self.raw_format = raw_format
        
        if self.raw_format == False:
            self.imread_func = lambda x: io.imread(x)
        
        else:
            import rawpy
            self.imread_func = lambda x: rawpy.imread(x).raw_image
        
        # optionally using pre-calculated BG image:
        if type(remove_ststic_BG) == bool:
            self.BG_remove = remove_ststic_BG
            
        else:
            from numpy import ndarray
            if type(remove_ststic_BG) == ndarray:
                self.BG = remove_ststic_BG
            self.BG_remove=None
        
        
        
    
    
    def get_file_names(self):
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
        
        print('calculating background...\n')
        # getting the subsample of images for BG calculation
        if len(self.image_files)<=200:
            BG_images = self.image_files
            
        else: 
            BG_images=self.image_files[::int(len(self.image_files)/400+1)][:200]
        
        BG_images = [os.path.join(self.dir_name, im) for im in BG_images]
        # reading the images for BG subtraction
        
        for i in range(len(BG_images)):
            if i==0:
                #im0 = io.imread(BG_images[i])*1.0
                im0 = self.imread_func(BG_images[i])*1.0
            else:
                #im0 += io.imread(BG_images[i])
                im0 += self.imread_func(BG_images[i])
        
        self.BG = im0 / len(BG_images)
        #ic = io.ImageCollection(BG_images)
        #self.BG = npmedian(ic, axis=0)
        
        
        
        
    def segment_folder_images(self):
        '''This loops over the image files in a folder'''
        
        self.get_file_names()
        
        if self.N_img is None: 
            N = len(self.image_files)
        else:
            N = self.N_img
            
        if type(self.BG_remove)==bool:
            if self.BG_remove==True:
                self.calculate_BG()
            else:
                self.BG = None
        
        i0 = (self.image_start is not None) * self.image_start        
        
        
        params = [self.dir_name, self.image_files, self.sigma, self.th, 
                  self.median, self.loc_filter, self.BG, self.mask, 
                  self.bbox_limits, self.mass_limits, self.method, 
                  self.p_size, i0]
        
        print('Starting loop segmentation.\n')
        try:
            import multiprocessing
            # Running with paralelization:
            print('Running with multiplrocessing...')
            t0 = time()
            args = [(X, 
                     self.imread_func(os.path.join(params[0], 
                                                   params[1][X])),
                     params) for X in range(N)]
            with multiprocessing.Pool() as pool:
                results_list = list(pool.starmap(iter_frame, args))
            print('finished segmentation loop (%.1f sec)'%(time() - t0))
    
        except ImportError as e:
            # Running without paralelization:
            print('Cant import multiprocessing - running on a single core')
            results_list = [iter_frame(i, 
                                       self.imread_func(
                                           os.path.join(params[0], 
                                                        params[1][i])),
                                       params) 
                            for i in tqdm.tqdm(range(N))] 
            
        self.blobs = [b for res_i in results_list for b in res_i]
        
                                       
    def save_results(self, fname):
        '''
        Will save the extracted blobs. 
        
        The format of the results is
        center_x, center_y, size_x, size_y, area, frame_number
        '''
        savetxt(fname, self.blobs, 
                fmt=['%.02f','%.02f','%d','%d','%d','%d'], delimiter='\t')
        
        
        


def iter_frame(i, im, params):
    '''
    Worker function for the multiprocessing
    '''
    
    ps = particle_segmentation(im,
                               sigma=params[2], 
                               threshold=params[3],
                               median=params[4],
                               local_filter=params[5],
                               BG_image=params[6],
                               mask=params[7],
                               max_xsize=params[8][1],
                               min_xsize=params[8][0],
                               max_ysize=params[8][3],
                               min_ysize=params[8][2],
                               min_mass=params[9][0],
                               max_mass=params[9][1],
                               method = params[10],
                               particle_size=params[11])
    ps.get_blobs()
    ps.apply_blobs_size_filter()
    res_i = [[b[0][0], b[0][1], b[1][0], b[1][1], b[2], i+params[12]] for b in ps.blobs]
    print('Frame: %d  ;  Blobs: %d'%(i, len(res_i)))
    return res_i







class tracking_augmented_segmentation(track_2D_multiframe):
    
    def __init__(self, fname, max_dt, Ns, d_max=1e10, dv_max=1e10, NSR_th=0.25):
        '''
        A class that uses 2D tracking of segmented blobs using the multiframe
        algorithm, and then "missing particle" interpolation in attempt to 
        fill in for occlusions.
        
        fname - name of a file with blobs to be tracking augmented
        
        max_dt, Ns, d_max, dv_max, NSR_th - tracking related to the multiframe
        tracking used. See the class over in tracking_mod for details. 
        '''
        
        self.fname = fname
        self.reverse_eta_zeta = False
        self.z_particles = 0.0
        self.U = 0
        
        self.max_dt = max_dt
        self.Ns = Ns
        self.d_max = d_max
        self.dv_max = dv_max
        self.NSR_th = NSR_th
        
        self.blobs = dict([(k, array(g)) for k,g in 
                           read_csv(self.fname,header=None,sep='\t').groupby(5)])
        
        self.times = sorted(list(self.blobs.keys()))
        
        self.trees = {}
        
        self.used_particles = dict([(tm, []) for tm in self.times])
        
        # a list to store trajectories
        self.trajs = []
        
        self.blobs_len = len(self.blobs[self.times[0]][0])
        

    def blobs_to_particles(self):
        '''
        This function uses the blob data to generate a dicionary of particles 
        in the format used by tracker_multiframe. 
        '''
        self.particles = {}
        
        for k in self.blobs.keys():
            self.particles[k] = []
            
            for i in range(len(self.blobs[k])):
                blob = self.blobs[k][i]

                p = array([-1, blob[0], blob[1], 0.0, i, 0.0, blob[-1]])
                self.particles[k].append(p)
                
            self.particles[k] = array(self.particles[k])
            
            
    
    def augment_blobs(self):
        '''
        Will interpolate trajectories that have skipped frames. We interpolate
        the missing points by using a 3nd order polynomial, namely assuming 
        linear acceleration in the interpolated range. Then, we add the 
        interpolated points to the blobs dictionary.
        '''
        
        frame_skips = max([1, int(self.Ns/3)])
        self.track_frames(f0=None, fe=None, frame_skips=frame_skips)
        
        N_links_0 = sum([len(tr)+1 for tr in self.trajs])
        
        msg = 'Interpolating skipped frames'
        for i in tqdm.tqdm(range(len(self.trajs)), desc=msg):
            
            tr = self.trajs[i]
            
            # 1) check if the trajectory has skipped frames; if not, we
            #    continue to the next trajectory
            dt = tr[-1,-1] - tr[0,-1]
            l = len(tr)
            if dt == l-1: continue
            
            # 2) interpolate the missing points by interpolation
            interpolated = fill_in_trajectory(tr)
            self.trajs[i] = interpolated
            
            # 3) get the interpolated points
            tr_times = tr[:,-1]
            interp_points = [interpolated[i] for i in range(len(interpolated))
                             if interpolated[i][-1] not in tr_times]
            
            # 4) add the interpolated points to self.blobs
            for ip in interp_points:
                zeros = [-1 for i in range(self.blobs_len-3)]
                blob = [ip[1], ip[2]] + zeros + [ip[-1]]
                self.blobs[ip[-1]] = npappend(self.blobs[ip[-1]], [blob],
                                              axis=0)
                
            
        N_links_e = sum([len(tr)+1 for tr in self.trajs])
        N_new = N_links_e-N_links_0
        stats = (N_new, N_new/N_links_e*100)
        print('')
        print('interpolated %d points (%.1f percent)'%stats)
        
        
        
    def save_augmented_blobs(self, fname):
        '''
        Will save the blobs with a given file name
        '''
        
        tosave = []
        for k in self.blobs.keys():
            tosave += list(self.blobs[k])
        
        savetxt(fname, tosave, delimiter='\t', 
                fmt=['%.02f','%.02f','%d','%d','%d','%d'])
        


# =============================================================================
# if __name__ == '__main__':
#     fn = '/home/ron/Desktop/Research/jetArrayTank/20241020_puffs/Rec18/blobs_cam4'
#     dt_max = 3 
#     Ns = 9
#     t2d = tracking_augmented_segmentation(fn, dt_max, Ns, d_max=3, dv_max=3)
#     t2d.blobs_to_particles()
#     t2d.augment_blobs()
#     
#     t2d.save_augmented_blobs('blobs_cam4_augmented')
# =============================================================================









