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
        
        
        self.allowed_actions = ['help', 'initial_calibration', 
                                'final_calibration',
                                'calibration_with_particles', 'matching', 
                                'segmentation',
                                'smoothing', 'stitching', 'tracking', 
                                'calibration', 'calibration_point_gui', 
                                'match_target_file', '2D_tracking', 
                                'manual_matching',
                                'run_extention']
        
        
        # perform the wanted action:
        if action is None:
            print('Started workflow with no particular action.')
            
        elif action != None:
            
            msg1 = 'The given action is unknown.'
            msg2 = 'allowed actions are:'+str(self.allowed_actions)
            if action not in self.allowed_actions:
                raise ValueError(msg1+'\n'+msg2)
            
            elif action == 'initial_calibration':
                self.initial_calibration()
                
            elif action == 'final_calibration':
                self.final_calibration()
                
            elif action == 'calibration_with_particles':
                self.calibration_with_particles()
                
            elif action == 'segmentation':
                self.do_segmentation()
                
            elif action == 'matching':
                self.do_matching()
                
            elif action == 'tracking':
                self.do_tracking()
                
            elif action == '2D_tracking':
                self.do_2d_tracking()
                
            elif action == 'smoothing':
                self.do_smoothing()
            
            elif action == 'stitching':
                self.do_stitching()
            
            elif action == 'manual_matching':
                self.do_manual_matching()
            
            elif action == 'run_extention':
                self.do_run_extention()    
            
            elif action == 'help':
                self.help_me()
                
                
            # legacy functions:
            elif action == 'calibration':
                print('Note: you are running an outdated action!')
                print('consider using the initial_calibration and')
                print('final_calibration actions instead.')
                self.calibration_sequence()
            
            elif action == 'calibration_point_gui':
                print('Note: you are running an outdated action!')
                print('consider using the initial_calibration and')
                print('final_calibration actions instead.')
                self.calibration_point_gui()
            
            elif action == 'match_target_file':
                print('Note: you are running an outdated action!')
                print('consider using the initial_calibration and')
                print('final_calibration actions instead.')
                self.match_target_file()
            
            
    
    def help_me(self):
        '''
        Prints a message that might help users with the allowable commands.
        '''
        print('\nThe workflow script is intended to help users utilize MyPTVs')
        print('capabilities in their 3D particle tracking experiments. \n')
        
        print('To use the workflow script, run it with Python, usin one of')
        print('the following actions that you would like to perform:\n')
        
        for e, act in enumerate(self.allowed_actions):
            print('%d) "%s"'%(e, act))
            
        print('\nThe script will now close, so the wanted action could be run.')
        
        print('\nGood luck!')
        
        print('\nP.S. - try using the user manual that is found on the main')
        print('Github repository.')
            
    
            
    
    
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
    
    
    
    
    def initial_calibration(self):
        '''
        Starts the initial calibration GUI
        '''
        from myptv.gui_intial_cal import initial_cal_gui
        from matplotlib.pyplot import imread
        
        # fetch parameters from the file
        cam_name = self.get_param('calibration', 'camera_name')
        cal_image = self.get_param('calibration', 'calibration_image')
        target_file = self.get_param('calibration', 'target_file')
        res = self.get_param('calibration', 'resolution').split(',')
        res = (float(res[0]), float(res[1]))
        
        image = imread(cal_image)
        if image.shape[1] != res[0] or image.shape[0] != res[1]:
            msg = 'The given resolution doesnt match the image size'
            raise ValueError(msg)
        
        gui = initial_cal_gui(cam_name, cal_image, target_file)
        
        
        
    def final_calibration(self):
        '''
        Starts the initial calibration GUI
        '''
        from myptv.gui_final_cal import cal_gui
        from myptv.imaging_mod import camera
        from myptv.calibrate_mod import calibrate
        from matplotlib.pyplot import imread
        import os
        
        # fetch parameters from the file
        cam_name = self.get_param('calibration', 'camera_name')
        cal_image = self.get_param('calibration', 'calibration_image')
        res = self.get_param('calibration', 'resolution').split(',')
        res = (float(res[0]), float(res[1]))
        
        
        # checking that a camera file in the working directory
        ls = os.listdir('.')                
        
        # make sure camera file exists
        if cam_name not in ls:
            msg = 'No camera file was found. Start with initial calibration.'
            raise ValueError(msg)
            
        # detect the calibration folder
        cal_folder = '.'
        for fname in ls:
            if fname in ['calibration', 'Calibration', 'cal', 'Cal']:
                if os.path.isdir(os.path.join('.', fname)):
                    cal_folder = os.path.join('.', fname)
                    break
        
        # get the blob file and setup the camera instance
        blob_file = os.path.join(cal_folder, cam_name+'_cal_points')
        
        try:
            cam = camera(cam_name, res, cal_points_fname = blob_file)
        except:
            msg = 'Calibration point file (%s) not found!'%blob_file
            msg2 = 'Make sure the initial calibration was completed fully.'
            raise ValueError(msg+msg2)
            
        
        # load the camera
        cam.load('.')
        print('camera data loaded successfully.')
        cal = calibrate(cam, cam.lab_points, cam.image_points)
        print('initial error: %.3f pixels\n'%(cal.mean_squared_err()))
        
        # run the final calibration gui
        print('starting calibration GIU\n')
        gui = cal_gui(cal, cal_image)                                   
            
        
    
    
    
    # ========================================================================
    # /\/\//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    #  Legacy functions that are no longer needed due to the cal_gui's
    # /\/\//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    # ========================================================================
    
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
        #from matplotlib.pyplot import subplots, show, imread
        
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
            print('initial error: %.3f pixels\n'%(cal.mean_squared_err()))
            
            
            print('starting calibratino GIU\n')
            from myptv.gui_final_cal import cal_gui
            gui = cal_gui(cal, cal_image)
                    
            
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
        
        
    # ========================================================================
    # /\/\//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    #                         End of legacy functions
    # /\/\//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    # ========================================================================
    
    
    
    def calibration_with_particles(self):
        '''
        This starts the calibrate with particles sequence
        '''
        from myptv.gui_final_cal import cal_gui
        from myptv.imaging_mod import camera
        from myptv.calibrate_mod import calibrate_with_particles
        from  matplotlib.pyplot import subplots, show
        
        # fetch parameters from the file
        camera_name =  self.get_param('calibration_with_particles',
                                      'camera_name')
        resolution = self.get_param('calibration_with_particles',
                             'resolution').split(',')
        resolution = (float(resolution[0]), float(resolution[1]))
        traj_filename = self.get_param('calibration_with_particles',
                                      'traj_filename')
        cam_number = self.get_param('calibration_with_particles',
                                      'cam_number') 
        blobs_fname = self.get_param('calibration_with_particles',
                                      'blobs_fname')
        min_traj_len = self.get_param('calibration_with_particles',
                                      'min_traj_len')
        max_point_number = self.get_param('calibration_with_particles',
                                      'max_point_number')
        cal_image = self.get_param('calibration_with_particles', 
                                   'calibration_image')
        print('\n', 'starting calibration with particles')
        
        # setting up a camera instance            
        cam = camera(camera_name, resolution)
        cam.load('./')
        
        
        # set up the calibration object
        cal_with_p = calibrate_with_particles(traj_filename, cam, cam_number, 
                                      blobs_fname, min_traj_len=min_traj_len,
                                      max_point_number=max_point_number)
        
        
        cal = cal_with_p.get_calibrate_instance()
        
        # run the final calibration gui
        print('starting calibration GIU using calibration with particles\n')
        gui = cal_gui(cal, cal_image) 
        
        
        # print('\n', 'ready to calibrate')
        # print('initial error: %.3f pixels'%(cal.mean_squared_err()))
        # print('')
        
        # user = True
        # print('Starting calibration sequence:')
        # while user != '9':
        #     print("enter '1' for external parameters calibration")
        #     print("enter '2' for internal correction ('fine') calibration")
        #     print("enter '3' to show current camera external parameters")
        #     print("enter '4' to plot the calibration points' projection")
        #     print("enter '8' to save the results")
        #     print("enter '9' to quit")
        #     user = input('')
            
        #     if user == '1':
        #         print('\n', 'Iterating to minimize external parameters')
        #         cal.searchCalibration(maxiter=2000)
        #         err = cal.mean_squared_err()
        #         print('\n','calibration error: %.3f pixels'%(err),'\n')
            
        #     if user == '2':
        #         print('\n', 'Iterating to minimize correction terms')
        #         cal.fineCalibration()
        #         err = cal.mean_squared_err()
        #         print('\n','calibration error:', err,'\n')
                
        #     if user == '3':
        #         print('\n', cam, '\n')
            
        #     if user == '4':
        #         fig, ax = subplots()
        #         cal.plot_proj(ax=ax)
        #         show()
        #     if user == '8':
        #         print('\n', 'Saving results')
        #         cam.save('.')
        
    
    
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
        image_start = self.get_param('segmentation', 'image_start')
        sigma = self.get_param('segmentation', 'blur_sigma')
        threshold = self.get_param('segmentation', 'threshold')
        median = self.get_param('segmentation', 'median')
        local_filter = self.get_param('segmentation', 'local_filter')
        max_xsize = self.get_param('segmentation', 'max_xsize')
        max_ysize = self.get_param('segmentation', 'max_ysize')
        max_mass = self.get_param('segmentation', 'max_mass')
        min_xsize = self.get_param('segmentation', 'min_xsize')
        min_ysize = self.get_param('segmentation', 'min_ysize')
        min_mass = self.get_param('segmentation', 'min_mass')
        mask = self.get_param('segmentation', 'mask')
        plot_res = self.get_param('segmentation', 'plot_result')
        save_name = self.get_param('segmentation', 'save_name')
        ROI = self.get_param('segmentation', 'ROI')
        single_img_name = self.get_param('segmentation', 'single_image_name')
        method = self.get_param('segmentation', 'method')
        p_size = self.get_param('segmentation', 'particle_size')
        
        # reading preprepared mask
        if type(mask)==str:
            mask = imread(mask)
        
        if method not in ['dilation', 'labeling']:
            raise ValueError('Method can be only "dilation" or "labeling"')
        
        if method=='dilation' and type(p_size) != int:
            raise ValueError('In dilation, particle_size can only be integer')
        
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
                                            particle_size=p_size,
                                            extension=ext,
                                            image_start=image_start,
                                            N_img=N_img, 
                                            sigma=sigma, 
                                            median=median,
                                            threshold=threshold, 
                                            local_filter=local_filter, 
                                            max_xsize=max_xsize, 
                                            max_ysize=max_ysize,
                                            max_mass=max_mass,
                                            min_xsize=min_xsize, 
                                            min_ysize=min_ysize,
                                            min_mass=min_mass,
                                            mask=mask,
                                            method=method)
        
            loopSegment.segment_folder_images()
            
            print('\n','blobs found:', len(loopSegment.blobs))
            
            # saving the semented blobs:
            if save_name is not None and type(save_name)==str:
                cwd_ls = os.listdir(os.getcwd())
                if save_name in cwd_ls or os.path.exists(save_name):
                    print('\n The file name "%s" already exists in'%save_name)
                    print(' the working directory. Should I save anyways?')
                    usr = input('(1=yes, else=no)')
                    if usr == '1':
                        loopSegment.save_results(save_name)
                        print('\nfile saved.')
                    else:
                        print('\nskipped saving')
                    
                else:
                    loopSegment.save_results(save_name)
                    print('\nfile saved.')    
            print('\nDone.\n')
            
        
        
        # segmenting the image if there is only 1 frames
        if N_img == 1:
            print('\n','starting segmentation on a single image.')
            if single_img_name not in image_files:
                in_ = os.path.join(dirname,single_img_name)
                msg = 'Image %s not found in the directory.'%in_
                raise ValueError(msg)
            
            print('\n','segmenting image: %s'%single_img_name)
            particleSegment = particle_segmentation(image0, 
                                                    particle_size=p_size,
                                                    sigma=sigma, 
                                                    median=median,
                                                    threshold=threshold, 
                                                    local_filter=local_filter, 
                                                    max_xsize=max_xsize, 
                                                    max_ysize=max_ysize,
                                                    max_mass=max_mass,
                                                    min_xsize=min_xsize, 
                                                    min_ysize=min_ysize,
                                                    min_mass=min_mass,
                                                    mask=mask,
                                                    method=method)
            particleSegment.get_blobs()
            particleSegment.apply_blobs_size_filter()
            
            print('blobs found:', len(particleSegment.blobs))
            
            if plot_res:
                from matplotlib.pyplot import show
                particleSegment.plot_blobs()
                show()
                
                
            # Saving the segmented blobs:
            if save_name is not None and type(save_name)==str:
                cwd_ls = os.listdir(os.getcwd())
                if save_name in cwd_ls or os.path.exists(save_name):
                    print('\n The file name "%s" already exists in'%save_name)
                    print(' the working directory. Should I save anyways?')
                    usr = input('(1=yes, else=no)')
                    if usr == '1':
                        particleSegment.save_results(save_name)
                        print('\nfile saved.')
                    else:
                        print('\nskipped saving')
                else:
                    particleSegment.save_results(save_name)
                    print('\nfile saved.')
                
            print('\nDone.\n')
            
            
            
            
    def do_matching(self):
        '''
        Will perform the stereo matching with the file given parameters
        '''
        from myptv.particle_matching_mod import match_blob_files
        from myptv.imaging_mod import camera, img_system
        from os import getcwd, listdir
        from os.path import exists as pathExists
        
        
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
        frame_start = self.get_param('matching', 'frame_start')
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
        ts = int(mbf.time_lst[0])
        te = int(mbf.time_lst[-1])
        print('segmented particles time range: %d -> %d'%(ts,te),'\n')
        
        if frame_start is not None:
            if frame_start>=ts and frame_start <=te:
                ts = frame_start
            else: 
                raise ValueError('frame_start outside the available frame range')
        
        if N_frames is None:
            frames = range(ts, te+1)
        else:
            try:
                frames = range(ts, ts+N_frames)
            except:
                tp = type(frames)
                msg = 'N_frames must be an integer or None (given %s).'%tp
                raise TypeError(msg)
                
                
        # mathing
        print('Starting stereo-matching.')
        print('Matching in the frame range: %d -> %d'%(frames[0], frames[-1]))
        mbf.get_particles(frames=frames)
        
        # print matching statistics
        print('particles matched: %d \n'%(len(mbf.particles)))
        
        Nframes = len(frames)
        c4 = sum([1 for p in mbf.particles if len(p[3])==4]) / Nframes
        print('quadruplets: %.1f per frame'%c4)
        c3 = sum([1 for p in mbf.particles if len(p[3])==3]) / Nframes
        print('triplets: %.1f per frame'%c3)
        c2 = sum([1 for p in mbf.particles if len(p[3])==2]) / Nframes
        print('pairs: %.1f per frame \n'%c2)
        
        
        # save the results
        if save_name is not None:
            cwd_ls = listdir(getcwd())
            if save_name in cwd_ls or pathExists(save_name):
                print('\n The file name "%s" already exists in'%save_name)
                print(' the working directory. Should I save anyways?')
                usr = input('(1=yes, else=no)')
                if usr == '1':
                    print('\n','saving file.')
                    mbf.save_results(save_name)
                else:
                    print('\n','skiped saving.')
                
            else:
                print('\n','saving file.')
                mbf.save_results(save_name)
        
        print('\n', 'Finished Matching.\n')
            
        
        
    def do_tracking(self):
        '''
        Will perform the tracking using the file given parameters.
        '''
        from myptv.tracking_mod import tracker_four_frames
        from numpy import array
        from os import getcwd, listdir
        from os.path import exists as pathExists
        
        # fetching parameters
        particles_fm = self.get_param('tracking', 'particles_file_name')
        frame_start = self.get_param('tracking', 'frame_start')
        N_frames = self.get_param('tracking', 'N_frames')
        d_max = self.get_param('tracking', 'd_max')
        dv_max = self.get_param('tracking', 'dv_max')
        mean_flow = self.get_param('tracking', 'mean_flow')
        candidate_graph = self.get_param('tracking', 'plot_candidate_graph')
        save_name = self.get_param('tracking', 'save_name')
        
        
        # initiate the tracker
        t4f = tracker_four_frames(particles_fm, 
                                  d_max=d_max, 
                                  dv_max=dv_max,
                                  mean_flow=array(mean_flow),
                                  store_candidates = candidate_graph)
        
        #setting up the frame range
        ts = int(t4f.times[0])
        te = int(t4f.times[-1])
        
        print('available particles time range: %d -> %d'%(ts,te),'\n')
        
        if candidate_graph and (te-ts)>100:
            print('Warning: you are about to plot a candidate graph with')
            print('more than 100 frames.')
            ans = input('Do you wish to proceed (1 = Yes , else = No)?  ')
            
            if ans=='1':
                pass
            
            else:
                print('quitting ')
                return None
            
        
        if frame_start is not None:
            if frame_start>=ts and frame_start <=te:
                ts = frame_start
            else: 
                print('Warning: frame_start outside the available frame range')
                #raise ValueError('frame_start outside the available frame range')
        
        if N_frames is None:
            frames = range(ts, te)
        else:
            try:
                frames = range(ts, ts+N_frames)
            except:
                tp = type(frames)
                msg = 'N_frames must be an integer or None (given %s).'%tp
                raise TypeError(msg)
        
        
        # do the tracking
        t4f.track_all_frames(frames=frames)
        
        # print some statistics
        tr = array(t4f.return_connected_particles())
        untracked = len(tr[tr[:,0]==-1])
        tot = len(tr)
        print('untracked fraction:', untracked/tot)
        print('tracked per frame:', (tot-untracked)/len(set(tr[:,-1])))
        
        if candidate_graph:
            t4f.plot_candidate_graph()
        
        # save the results
        if save_name is not None:
            cwd_ls = listdir(getcwd())
            if save_name in cwd_ls or pathExists(save_name):
                print('\n The file name "%s" already exists in'%save_name)
                print(' the working directory. Should I save anyways?')
                usr = input('(1=yes, else=no)')
                if usr == '1':
                    print('\n','saving file.')
                    t4f.save_results(save_name)
                else:
                    print('\n', 'skipped saving.')
            
            else:
                print('\n','saving file.')
                t4f.save_results(save_name)
        
        print('\n', 'Finished tracking.')
        
        
        
    def do_smoothing(self):
        '''
        Will smooth the trajectories using the specified file given paramters.
        '''
        from numpy import loadtxt
        from myptv.traj_smoothing_mod import smooth_trajectories
        from os import getcwd, listdir
        from os.path import exists as pathExists
        
        # fetching the smoothing parameters
        trajectory_file = self.get_param('smoothing', 'trajectory_file')
        window = self.get_param('smoothing', 'window_size')
        polyorder = self.get_param('smoothing', 'polynom_order')
        min_traj_length = self.get_param('smoothing', 'min_traj_length')
        repetitions = self.get_param('smoothing', 'repetitions')
        save_name = self.get_param('smoothing', 'save_name')
        
        if min_traj_length <= polyorder:
            raise ValueError('min_traj_length must be larger than polyorder')

        traj_list = loadtxt(trajectory_file)
        
        
        # smoothing the trajectories     
        print('Starting to smooth trajectories.')
        sm = smooth_trajectories(traj_list, 
                                 window, 
                                 polyorder,
                                 repetitions=repetitions,
                                 min_traj_length=min_traj_length)
        sm.smooth()
        
        # saving the data
        if save_name is not None:
            cwd_ls = listdir(getcwd())
            if save_name in cwd_ls or pathExists(save_name):
                print('\n The file name "%s" already exists in'%save_name)
                print(' the working directory. Should I save anyways?')
                usr = input('(1=yes, else=no)')
                if usr == '1':
                    print('\n', 'Saving the smoothed data (%s).'%save_name)
                    sm.save_results(save_name)
                else:
                    print('\n', 'Skipped saving file.')
            
            else:
                print('\n', 'Saving the smoothed data (%s).'%save_name)
                sm.save_results(save_name)
        
        print('\n', 'Done.')
        
        
    
    
    def do_stitching(self):
        '''
        Will perfrom trajectory stitching using the file given parameters.
        '''
        from numpy import loadtxt
        from myptv.traj_stitching_mod import traj_stitching
        from os import getcwd, listdir
        from os.path import exists as pathExists
        
        # fetchhing the stitching parameters
        trajectory_file = self.get_param('stitching', 'trajectory_file')
        Ts = self.get_param('stitching', 'max_time_separation')
        dm = self.get_param('stitching', 'max_distance')
        save_name = self.get_param('stitching', 'save_name')
        
        traj_list = loadtxt(trajectory_file)
        
        # stitch the trajectories
        ts = traj_stitching(traj_list, Ts, dm)
        ts.stitch_trajectories()
        
        # saving the data
        if save_name is not None:
            cwd_ls = listdir(getcwd())
            if save_name in cwd_ls or pathExists(save_name):
                print('\n The file name "%s" already exists in'%save_name)
                print(' the working directory. Should I save anyways?')
                usr = input('(1=yes, else=no)')
                if usr == '1':
                    print('\n', 'Saveing the data.')    
                    ts.save_results(save_name)
                else:
                    print('\n', 'Skipped saving file.')
            
            else:
                print('\n', 'Saveing the data.')    
                ts.save_results(save_name)
        
        print('\n', 'Done.')
        
        
        
    def do_2d_tracking(self):
        '''
        Will perform 2D tracking of segmented blobs using give data.
        '''
            
        from myptv.imaging_mod import camera
        from myptv.tracking_2D_mod import track_2D
        
        # fetchhing the stitching parameters
        fname = self.get_param('2D_tracking', 'blob_file')
        frame_start = self.get_param('2D_tracking', 'frame_start')
        N_frames = self.get_param('2D_tracking', 'N_frames')
        cam_name = self.get_param('2D_tracking', 'camera_name')
        res = self.get_param('2D_tracking', 'camera_resolution')
        res = tuple([float(val) for val in res.split(',')])
        z_particles = self.get_param('2D_tracking', 'z_particles')
        d_max = self.get_param('2D_tracking', 'd_max')
        dv_max = self.get_param('2D_tracking', 'dv_max')
        save_name = self.get_param('2D_tracking', 'save_name')

        print('\ninitiating 2D tracking...')

        cam = camera(cam_name, res)
        cam.load('')
        
        print('\nloading blobs and transforming to lab-space coordinates')
        t2d = track_2D(cam, fname, z_particles, d_max=d_max, dv_max = dv_max, 
                       reverse_eta_zeta=True)
        
        t2d.blobs_to_particles()
        
        
        #setting up the frame range
        ts = int(t2d.times[0])
        te = int(t2d.times[-1])
        
        print('\navailable particles time range: %d -> %d'%(ts,te),'\n')
        
        if frame_start is not None:
            if frame_start>=ts and frame_start <=te:
                ts = frame_start
            else: 
                print('Warning: frame_start outside the available frame range')
                #raise ValueError('frame_start outside the available frame range')
        
        if N_frames is None:
            frames = range(ts, te)
        else:
            try:
                frames = range(ts, ts+N_frames)
            except:
                tp = type(frames)
                msg = 'N_frames must be an integer or None (given %s).'%tp
                raise TypeError(msg)
        
        
        print('\ntrackin particles...')
        
        t2d.track_all_frames(frames=frames)
        
        print('\nsaving results...')
        
        t2d.save_results(save_name)
        
        print('\nDone!')
        
        
    
    def do_manual_matching(self):
        '''
        Runs a GUI that helps performing manual stereo-matching of
        points from images. You simply click on the images from different
        cameras and the GUI gives back the 3D coordinates of this point. 
        '''
        from myptv.gui_manual_matching import man_match_gui
        
        # fetchhing the stitching parameters
        camera_names = self.get_param('manual_matching_GUI', 'cameras')
        im_fname = self.get_param('manual_matching_GUI', 'images')
        
        print(camera_names)
        print(im_fname)
        
        gui = man_match_gui(camera_names, im_fname, cameras_folder='.')
    
        
        
    def do_run_extention(self):
        '''
        This is an option to load extrenal extentions to MyPTV. Get it done
        by setting the propper parameters in the params_file.
        '''
        
        # fetchhing the stitching parameters
        path_to_extention = self.get_param('run_extention', 'path_to_extention')
        action_name = self.get_param('run_extention', 'action_name')
        extention_params_file = self.get_param('run_extention', 'extention_params_file')
        
        # 1) import the script  "path_to_extention"
        
        # 2) load the extensions' parameter from extention_params_file
        
        # 3) run the class given as action_name, with the parameter given
        
        
        
        return None
        
#%%
        
        
        


if __name__ == '__main__':
    
    import sys
    fname, action = sys.argv[1], sys.argv[2]
    print('\n','given inputs -')
    print('params file name:', fname)
    print('action:', action, '\n')
    wf = workflow(fname, action)
    
    
    
    
    

