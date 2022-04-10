# -*- coding: utf-8 -*-
"""
Created on Sun 20 March 2022


This script containes a class designed to make using MyPTV a bit easier
for users. 

1) We use a single text file to hold all the parameters used in MyPTV.

2) We havea class that reads the text file and performs the given task:
    segmentaiton, matching, tracking, smoothing and stitching.

"""

from pandas import json_normalize
from pandas import DataFrame as df
from yaml import safe_load






class workflow(object):
    '''
    A class used to run specific MyPTV operations with parameters given 
    in a dedicated text file.
    '''
    
    def __init__(self, param_file, action):
        '''
        input -
        
        param_file - string; the path to a text file with the specified 
                     parameters to be used.
        
        action - string; the name of the PTV action to be performed. Accepted
                 values are: 'segmentaiton', 'matching', 'tracking',
                 'smoothing', and 'stitching'.
        '''
        
        # read the parameter file:
        self.param_file_path = param_file
        self.params = self.read_params_file()
        
        # perform the wanted action:
        if action != None:
            
            if action == 'segmentation':
                self.do_segmentation()
                
            if action == 'matching':
                self.do_matching()
            
            else:
                raise ValueError('given action %s not allowed'%action)
            
    
    def read_params_file(self):
        '''
        Reads the yaml file and formats it as a DataFrame.
        '''
        
        with open(self.param_file_path, 'r') as f:
            params = {}
            sl = safe_load(f)
            for i in range(len(sl)):
                params.update(sl[i])
                    
        as_dict = {'operation':[], 'param': [], 'value': [] }
        for k in params.keys():
            for kk in params[k].keys():
                as_dict['operation'].append(k)
                as_dict['param'].append(kk)
                as_dict['value'].append(params[k][kk])
        
        for i in range(len(as_dict['value'])):
            if as_dict['value'][i] == 'None': 
                as_dict['value'][i] = None
        
        return df(as_dict)
    
    
    
    
    def do_segmentation(self):
        '''
        Will perform segmentation on the images given in the parameters file
        and save the results on the given location.
        '''
        from myptv.segmentation_mod import loop_segmentation
        from myptv.segmentation_mod import particle_segmentation
        from numpy import zeros
        from skimage.io import imread
        import os
        
        par_seg = self.params[self.params['operation']=='segmentation']
        
        dirname = par_seg[par_seg['param']=='images_folder']['value'].iloc[0]
        ext = par_seg[par_seg['param']=='image_extension']['value'].iloc[0]
        N_img = par_seg[par_seg['param']=='Number_of_images']['value'].iloc[0]
        sigma = par_seg[par_seg['param']=='blur_sigma']['value'].iloc[0]
        threshold = par_seg[par_seg['param']=='threshold']['value'].iloc[0]
        local_filter = par_seg[par_seg['param']=='local_filter']['value'].iloc[0]
        max_xsize = par_seg[par_seg['param']=='max_xsize']['value'].iloc[0]
        max_ysize = par_seg[par_seg['param']=='max_ysize']['value'].iloc[0]
        max_area = par_seg[par_seg['param']=='max_area']['value'].iloc[0]
        min_xsize = par_seg[par_seg['param']=='min_xsize']['value'].iloc[0]
        min_ysize = par_seg[par_seg['param']=='min_ysize']['value'].iloc[0]
        min_area = par_seg[par_seg['param']=='min_area']['value'].iloc[0]
        mask = par_seg[par_seg['param']=='mask']['value'].iloc[0]
        plot_res = par_seg[par_seg['param']=='plot_result']['value'].iloc[0]
        save_name = par_seg[par_seg['param']=='save_name']['value'].iloc[0]
        ROI = par_seg[par_seg['param']=='ROI']['value'].iloc[0]
        
        if type(mask)==str:
            mask = imread(mask)
        
        
        # get the shape of the images
        allfiles = os.listdir(dirname)
        n_ext = len(ext)
        image_files = sorted(list(filter(lambda s: s[-n_ext:]==ext, allfiles)))
        image0 = imread(os.path.join(dirname,image_files[0]))
        
        # prepare a mask using the given ROI
        if ROI is not None:
            ROI = [int(val) for val in ROI.split(',')]
            mask_ROI = zeros(image0.shape)
            mask_ROI[ROI[2]:ROI[3]+1, ROI[0]:ROI[1]+1] = 1
            mask = mask * mask_ROI
        
        
        
        if N_img is None or N_img>1:
            loopSegment = loop_segmentation(dirname, 
                                            extension=ext,
                                            N_img=N_img, 
                                            sigma=sigma, 
                                            threshold=threshold, 
                                            local_filter=local_filter, 
                                            max_xsize=max_xsize, 
                                            max_ysize=max_ysize,
                                            max_area=max_area,
                                            min_xsize=min_xsize, 
                                            min_ysize=min_ysize,
                                            min_area=min_area,
                                            mask=mask)
        
            loopSegment.segment_folder_images()
            
            print('\n','blobs found:', len(loopSegment.blobs))
            
            loopSegment.save_results(save_name)
            print('Done, file saved.')
        
        
        
        if N_img == 1:
            print('starting segmentation on a single image.')
            print(os.path.join(dirname,image_files[0]))
            particleSegment = particle_segmentation(image0, 
                                                    sigma=sigma, 
                                                    threshold=threshold, 
                                                    local_filter=local_filter, 
                                                    max_xsize=max_xsize, 
                                                    max_ysize=max_ysize,
                                                    max_area=max_area,
                                                    min_xsize=min_xsize, 
                                                    min_ysize=min_ysize,
                                                    min_area=min_area,
                                                    mask=mask)
            particleSegment.get_blobs()
            particleSegment.apply_blobs_size_filter()
            
            print('blobs found:', len(particleSegment.blobs))
            
            if plot_res:
                from matplotlib.pyplot import show
                particleSegment.plot_blobs()
                show()
                
            if save_name is not None and type(save_name)==str:
                particleSegment.save_results(save_name)
                print('file saved.')
                
            print('done.')
            
            
            
            
    def do_matching(self):
        '''
        Will perform the stereo matching with the file given parameters
        '''
        from myptv.particle_matching_mod import match_blob_files
        from myptv.imaging_mod import camera, img_system
        
        # fetching the parameters
        par_seg = self.params[self.params['operation']=='matching']
        
        blob_fn = par_seg[par_seg['param']=='blob_files']['value'].iloc[0]
        blob_fn = [val.strip() for val in blob_fn.split(',')]
        cam_names = par_seg[par_seg['param']=='camera_names']['value'].iloc[0]
        cam_names = [val.strip() for val in cam_names.split(',')]
        res = par_seg[par_seg['param']=='cam_resolution']['value'].iloc[0]
        res = tuple([float(val) for val in res.split(',')])
        ROI = par_seg[par_seg['param']=='ROI']['value'].iloc[0].split(',')
        ROI = [[float(ROI[2*i]), float(ROI[2*i+1])] for i in [0,1,2]]
        voxel_size = par_seg[par_seg['param']=='voxel_size']['value'].iloc[0]
        max_blob_distance = par_seg[par_seg['param']=='max_blob_distance']['value'].iloc[0]
        max_err = par_seg[par_seg['param']=='max_err']['value'].iloc[0]
        N_frames = par_seg[par_seg['param']=='N_frames']['value'].iloc[0]
        save_name = par_seg[par_seg['param']=='save_name']['value'].iloc[0]
        
        
        # setting up the img_system 
        cams = [camera(cn, res) for cn in cam_names]
        for cam in cams:
            try:
                cam.load('')
            except:
                raise ValueError('camera file %s not found'%cam.name)
        imsys = img_system(cams)
        
        
        mbf = match_blob_files(blob_fn, 
                               imsys, 
                               ROI, 
                               voxel_size, 
                               max_blob_distance,
                               max_err=max_err, 
                               reverse_eta_zeta=True)
        
        # setting the frame range to match
        if N_frames is None:
            frames = None
        else:
            try:
                frames = range(N_frames)
            except:
                tp = type(frames)
                msg = 'N_frames must be an integer of None (given %s).'%tp
                raise TypeError(msg)
                
        
        # mathing
        print('Starting stereo-matching.')
        mbf.get_particles(frames=frames)
        print('particles matched:', len(mbf.particles))
        
        # saving the results
        if save_name is not None:
            print('\n','saving file.')
            mbf.save_results(save_name)
        
        print('\n', 'Finished Matching.')
            





if __name__ == '__main__':
    
    import sys
    fname, action = sys.argv[1], sys.argv[2]
    print('\n','given inputs -')
    print('params file name:', fname)
    print('action:', action, '\n')
    wf = workflow(fname, action)
    
    
    
    
    

