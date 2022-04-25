# -*- coding: utf-8 -*-
"""
Created on Sun 20 March 2022


This script contains a class designed to make using MyPTV a bit easier
for users. 

1) We use a single text file to hold all the parameters used in MyPTV.

2) We have a class that reads the text file and performs the given task:
    segmentation, matching, tracking, smoothing and stitching.

"""


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
                 values are: 'segmentation', 'matching', 'tracking',
                 'smoothing', and 'stitching'.
        '''
        
        # read the parameter file:
        self.param_file_path = param_file
        self.params = self.read_params_file()
        
        # perform the wanted action:
        if action is None:
            print('Started workflow with no particular action.')
            
        elif action != None:
            
            allowed_actions = ['calibration', 'calibration_point_gui', 
                              'match_target_file', 'segmentation', 'matching',
                              'tracking', 'smoothing', 'stitching']
            
            msg1 = 'The given action is unknown.'
            msg2 = 'allowed actions are:'+str(allowed_actions)
            if action not in allowed_actions:
                raise ValueError(msg1+'\n'+msg2)
            
            elif action == 'calibration':
                self.calibration_sequence()
            
            elif action == 'calibration_point_gui':
                self.calibration_point_gui()
            
            elif action == 'match_target_file':
                self.match_target_file()
                
            elif action == 'segmentation':
                self.do_segmentation()
                
            elif action == 'matching':
                self.do_matching()
                
            elif action == 'tracking':
                self.do_tracking()
                
            elif action == 'smoothing':
                self.do_smoothing()
            
            elif action == 'stitching':
                self.do_stitching()
            
            
    
    
    def read_params_file(self):
        '''
        Reads the yaml file and formats it as a DataFrame.
        '''
        
        with open(self.param_file_path, 'r') as f:
            params = {}
            
            try:
                sl = safe_load(f)
            except:
                raise ValueError('Error in loading the parameters file.')
                
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
    
    
    
    
    def get_param(self, act, param):
        '''
        Fetches a parameter value from the self.params DataFrame.
        '''
        if act not in set(self.params['operation']):
            raise ValueError('Cant find action %s in the parameter file.'%act)
        
        par_seg = self.params[self.params['operation']==act]
        
        if param not in set(par_seg['param']):
            msg = 'Cant find the %s -> %s in the parameters file.'%(act,param)
            raise ValueError(msg)
        
        return par_seg[par_seg['param']==param]['value'].iloc[0]
    
    
    
    
    def match_target_file(self):
        '''
        This function is used to match calibration points to a target file.
        '''
        from myptv.imaging_mod import camera
        from myptv.utils import match_calibration_blobs_and_points
        from matplotlib.pyplot import show
        
        # fetch parameters from the file
        cam_name = self.get_param('calibration', 'camera_name')
        points_file = self.get_param('calibration', 'calibration_points_file')
        target_file = self.get_param('calibration', 'target_file')
        segmented_file = self.get_param('calibration', 'segmented_points_file')
        res = self.get_param('calibration', 'resolution').split(',')
        res = (float(res[0]), float(res[1]))
        
        print('Matching target file and segmental calibration points')
        print('for camera: %s'%cam_name)
        # initiate the camera
        cam = camera(cam_name, res)
        cam.load('.')
        
        # match the points
        mtf = match_calibration_blobs_and_points(cam,
                                                 segmented_file,
                                                 target_file)
        mtf.pair_points()
    
        # plot the pairs
        mtf.plot_projections()
        show()
        
        print('Please confirm that the points were pairs correctly')
        print("in the figure, by making sure that the blue points")
        print("and the red x's are close to each other." ,'\n')
        
        print("To confirm and save the calibration point file, enter '1'")
        print("If there are errors, enter '2' and improve the initial calibration.")
        user = input()
        
        if user == '1':
            mtf.save_results(points_file)
            print('\n', 'file saved', '\n')
        
        elif user == '2':
            print('\n', 'file not saved', '\n')
        
        else:
            print('\n','unrecognized command. Leaving the workflow.', '\n')
            
        print('Done.')
        
        
    
    
    
    def calibration_sequence(self):
        '''
        Starts a sequence to guide users through the calibration process.
        '''
        from myptv.imaging_mod import camera
        from myptv.calibrate_mod import calibrate
        from os import listdir
        from os.path import isfile
        from matplotlib.pyplot import subplots, show, imread
        
        # fetch parameters from the file
        cam_name = self.get_param('calibration', 'camera_name')
        blob_file = self.get_param('calibration', 'calibration_points_file')
        cal_image = self.get_param('calibration', 'calibration_image')
        res = self.get_param('calibration', 'resolution').split(',')
        res = (float(res[0]), float(res[1]))
        
        
        # checking that a camera file in the working directory
        ls = listdir('.')                
        
        # if the file is found, start calibration sequence
        if cam_name in ls:
            print('Starting calibration sequence.')
            
            try:
                cam = camera(cam_name, res, cal_points_fname = blob_file)
            except:
                print('\n','Calibration point file (%s) not found!'%blob_file)
                print('\n','Would you like to start the calibration point gui?')
                user = input('1=yes,  else=no : ')
                if user == '1':
                    self.calibration_point_gui()  
                else:
                    print('quitting...')
                return 
            
            cam.load('.')
            print('camera data loaded successfully.')
            cal = calibrate(cam, cam.lab_points, cam.image_points)
            print('initial error: %.3f pixels'%(cal.mean_squared_err()))
            print('')
            
            user = True
            print('Starting calibration sequence:')
            while user != '9':
                print("enter '1' for external parameters calibration")
                print("enter '2' for internal correction ('fine') calibration")
                print("enter '3' to show current camera external parameters")
                print("enter '4' to plot the calibration points' projection")
                print("enter '8' to save the results")
                print("enter '9' to quit")
                user = input('')
                
                if user == '1':
                    print('\n', 'Iterating to minimize external parameters')
                    cal.searchCalibration(maxiter=2000)
                    err = cal.mean_squared_err()
                    print('\n','calibration error: %.3f pixels'%(err),'\n')
                
                if user == '2':
                    print('\n', 'Iterating to minimize correction terms')
                    cal.fineCalibration()
                    err = cal.mean_squared_err()
                    print('\n','calibration error:', err,'\n')
                    
                if user == '3':
                    print('\n', cam, '\n')
                
                if user == '4':
                    img = imread(cal_image)
                    fig, ax = subplots()
                    ax.imshow(img, cmap='gray')
                    cal.plot_proj(ax=ax)
                    show()
                    
                if user == '8':
                    print('\n', 'Saving results')
                    cam.save('.')
                    
            
        # if not, generate an empty file camera file
        else:
            print('')
            print('camera files not detected in the working directory.')
            print('Generating a new empty file and leaving calib. sequence.')
            print('To continue calibration, fill in an initial guess in the')
            print('empty file, and then run again the calibration sequence.')
            cam = camera(cam_name, res)
            cam.save('.')
            print('\n', 'Done.')
    
    
    
    
    def calibration_point_gui(self):
        '''
        This will start the calibration segmentation point gui.
        '''
        from myptv.cal_point_gui import cal_point_gui
        
        # fetch parameters from the file
        blob_file = self.get_param('calibration', 'calibration_points_file')
        cal_image = self.get_param('calibration', 'calibration_image')
        res = self.get_param('calibration', 'resolution').split(',')
        res = (float(res[0]), float(res[1]))
        
        print('\n', 'Starting calibration point segmentation GUI', '\n')
        gui = cal_point_gui(cal_image, blob_file)
    
    
    
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
        
        # fetching parameters
        dirname = self.get_param('segmentation', 'images_folder')
        ext = self.get_param('segmentation', 'image_extension')
        N_img = self.get_param('segmentation', 'Number_of_images')
        sigma = self.get_param('segmentation', 'blur_sigma')
        threshold = self.get_param('segmentation', 'threshold')
        local_filter = self.get_param('segmentation', 'local_filter')
        max_xsize = self.get_param('segmentation', 'max_xsize')
        max_ysize = self.get_param('segmentation', 'max_ysize')
        max_area = self.get_param('segmentation', 'max_area')
        min_xsize = self.get_param('segmentation', 'min_xsize')
        min_ysize = self.get_param('segmentation', 'min_ysize')
        min_area = self.get_param('segmentation', 'min_area')
        mask = self.get_param('segmentation', 'mask')
        plot_res = self.get_param('segmentation', 'plot_result')
        save_name = self.get_param('segmentation', 'save_name')
        ROI = self.get_param('segmentation', 'ROI')
        single_img_name = self.get_param('segmentation', 'single_image_name')
        
        
        # reading preprepared mask
        if type(mask)==str:
            mask = imread(mask)
        
        # get the shape of the images
        allfiles = os.listdir(dirname)
        n_ext = len(ext)
        image_files = sorted(list(filter(lambda s: s[-n_ext:]==ext, allfiles)))
        if single_img_name in image_files:
            image0 = imread(os.path.join(dirname,single_img_name))
        else:
            image0 = imread(os.path.join(dirname,image_files[0]))
        
        # preparing a mask using the given ROI
        if ROI is not None:
            ROI = [int(val) for val in ROI.split(',')]
            mask_ROI = zeros(image0.shape)
            mask_ROI[ROI[2]:ROI[3]+1, ROI[0]:ROI[1]+1] = 1
            mask = mask * mask_ROI
        
        # segmenting the image if there are more than 1 frames
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
            
            if save_name is not None:
                loopSegment.save_results(save_name)
                print('File saved (%s).'%(save_name))
            
            print('Done.')
        
        
        # segmenting the image if there is only 1 frames
        if N_img == 1:
            print('\n','starting segmentation on a single image.')
            if single_img_name not in image_files:
                in_ = os.path.join(dirname,single_img_name)
                msg = 'Image %s not found in the directory.'%in_
                raise ValueError(msg)
            
            print('\n','segmenting image: %s'%single_img_name)
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
        blob_fn = self.get_param('matching', 'blob_files')
        blob_fn = [val.strip() for val in blob_fn.split(',')]
        cam_names = self.get_param('matching', 'camera_names')
        cam_names = [val.strip() for val in cam_names.split(',')]
        res = self.get_param('matching', 'cam_resolution')
        res = tuple([float(val) for val in res.split(',')])
        ROI = self.get_param('matching', 'ROI').split(',')
        ROI = [[float(ROI[2*i]), float(ROI[2*i+1])] for i in [0,1,2]]
        voxel_size = self.get_param('matching', 'voxel_size')
        max_blob_distance = self.get_param('matching', 'max_blob_distance')
        max_err = self.get_param('matching', 'max_err')
        N_frames = self.get_param('matching', 'N_frames')
        save_name = self.get_param('matching', 'save_name')
        
        
        
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
        
        # print matching statistics
        print('particles matched:', len(mbf.particles))
        c4 = sum([1 for p in mbf.particles if len(p[3])==4])
        print('quadruplets:', c4)
        c3 = sum([1 for p in mbf.particles if len(p[3])==3])
        print('triplets:', c3)
        c2 = sum([1 for p in mbf.particles if len(p[3])==2])
        print('pairs:', c2)
        
        
        # save the results
        if save_name is not None:
            print('\n','saving file.')
            mbf.save_results(save_name)
        
        print('\n', 'Finished Matching.')
            
        
        
    def do_tracking(self):
        '''
        Will perform the tracking using the file given parameters.
        '''
        from myptv.tracking_mod import tracker_four_frames
        from numpy import array
        
        # fetching parameters
        particles_fm = self.get_param('tracking', 'particles_file_name')
        N_frames = self.get_param('tracking', 'N_frames')
        d_max = self.get_param('tracking', 'd_max')
        dv_max = self.get_param('tracking', 'dv_max')
        save_name = self.get_param('tracking', 'save_name')
        
        
        
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
        
        
        
        # do the tracking
        t4f = tracker_four_frames(particles_fm, 
                                  d_max=d_max, 
                                  dv_max=dv_max)
        t4f.track_all_frames(frames=frames)
        
        # print some statistics
        tr = array(t4f.return_connected_particles())
        untracked = len(tr[tr[:,0]==-1])
        tot = len(tr)
        print('untracked fraction:', untracked/tot)
        print('tracked per frame:', (tot-untracked)/len(set(tr[:,-1])))

        # save the results
        if save_name is not None:
            print('\n','saving file.')
            t4f.save_results(save_name)
        
        print('\n', 'Finished tracking.')
        
        
        
    def do_smoothing(self):
        '''
        Will smooth the trajectories using the specified file given paramters.
        '''
        from numpy import loadtxt
        from myptv.traj_smoothing_mod import smooth_trajectories
        
        # fetching the smoothing parameters
        trajectory_file = self.get_param('smoothing', 'trajectory_file')
        window = self.get_param('smoothing', 'window_size')
        polyorder = self.get_param('smoothing', 'polynom_order')
        save_name = self.get_param('smoothing', 'save_name')
        
        traj_list = loadtxt(trajectory_file)
        
        
        # smoothing the trajectories     
        print('Starting to smooth trajectories.')
        sm = smooth_trajectories(traj_list, window, polyorder)
        sm.smooth()
        
        # saving the data
        if save_name is not None:
            print('\n', 'Saving the smoothed data (%s).'%save_name)
            sm.save_results(save_name)
        
        print('\n', 'Done.')
        
        
    
    
    def do_stitching(self):
        '''
        Will perfrom trajectory stitching using the file given parameters.
        '''
        from numpy import loadtxt
        from myptv.traj_stitching_mod import traj_stitching
        
        # fetchhing the stitching parameters
        trajectory_file = self.get_param('stitching', 'trajectory_file')
        Ts = self.get_param('stitching', 'max_time_separation')
        dm = self.get_param('stitching', 'max_distance')
        save_name = self.get_param('stitching', 'save_name')
        
        traj_list = loadtxt(trajectory_file)
        
        # stitch the trajectories
        ts = traj_stitching(traj_list, Ts, dm)
        ts.stitch_trajectories()
        
        # saveing the data
        if save_name is not None:
            print('\n', 'Saveing the data.')    
            ts.save_results(save_name)
        
        print('\n', 'Done.')
        
        
#%%
        
        
        


if __name__ == '__main__':
    
    import sys
    fname, action = sys.argv[1], sys.argv[2]
    print('\n','given inputs -')
    print('params file name:', fname)
    print('action:', action, '\n')
    wf = workflow(fname, action)
    
    
    
    
    

