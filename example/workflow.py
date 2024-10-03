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
                                'analyze_calibration_error',
                                'calibration_with_particles', 'matching', 
                                'segmentation',
                                'smoothing', 'stitching', 'tracking', 
                                'calibration', 'calibration_point_gui', 
                                'match_target_file', '2D_tracking', 
                                'manual_matching',
                                'fiber_orientations',
                                'plot_trajectories',
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
                
            elif action == 'analyze_calibration_error':
                self.calibration_error_estimation()
                
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
                
            elif action == 'fiber_orientations':
                self.do_orientations()
                
            elif action == 'plot_trajectories':
                self.do_plot_trajectories()
            
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
        from matplotlib.pyplot import imread
        
        # fetch parameters from the file
        model_name = self.get_param('calibration', '3D_model')
        cam_name = self.get_param('calibration', 'camera_name')
        cal_image = self.get_param('calibration', 'calibration_image')
        target_file = self.get_param('calibration', 'target_file')
        res = self.get_param('calibration', 'resolution').split(',')
        res = (float(res[0]), float(res[1]))
        
        
        if model_name == 'Tsai':
            from myptv.TsaiModel.gui_intial_cal import initial_cal_gui
            image = imread(cal_image)
            if image.shape[1] != res[0] or image.shape[0] != res[1]:
                msg = 'The given resolution doesnt match the image size'
                raise ValueError(msg)
        
            gui = initial_cal_gui(cam_name, cal_image, target_file)
        
        
        elif model_name == 'extendedZolof':
            from myptv.extendedZolof.gui_intial_cal import initial_cal_gui
            image = imread(cal_image)
            gui = initial_cal_gui(cam_name, cal_image, target_file)
            
        
        else:
            models = str(['Tsai', 'extendedZolof'])[1:-1]
            msg = 'Unknown 3D model; permisible model names are: '
            raise ValueError(msg + models)
        
        
        
    def final_calibration(self):
        '''
        Starts the initial calibration GUI
        '''
        import os
        
        # fetch parameters from the file
        model_name = self.get_param('calibration', '3D_model')
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
        
        if model_name == 'Tsai':
            from myptv.TsaiModel.gui_final_cal import cal_gui
            from myptv.TsaiModel.camera import camera_Tsai
            from myptv.TsaiModel.calibrate import calibrate_Tsai
            
            try:
                cam = camera_Tsai(cam_name, cal_points_fname = blob_file)
            except:
                msg = 'Calibration point file (%s) is not right!'%blob_file
                msg2 = 'check that the file exists and that it has no errors.'
                raise ValueError(msg+msg2)
                
            
            # load the camera
            cam.load('.')
            print('camera data loaded successfully.')
            cal = calibrate_Tsai(cam, cam.lab_points, cam.image_points)
            print('initial error: %.3f pixels\n'%(cal.mean_squared_err()))
            
            # run the final calibration gui
            print('starting calibration GIU\n')
            gui = cal_gui(cal, cal_image)   

        
        elif model_name == 'extendedZolof':
            from myptv.extendedZolof.camera import camera_extendedZolof
            from myptv.extendedZolof.calibrate import calibrate_extendedZolof
            from myptv.extendedZolof.gui_final_cal import cal_gui
            
            try:
                cam = camera_extendedZolof(cam_name, cal_points_fname = blob_file)
            except:
                msg = 'Calibration point file (%s) is not right!'%blob_file
                msg2 = 'check that the file exists and that it has no errors.'
                raise ValueError(msg+msg2)
            
            cam.load('.')
            print('camera data loaded successfully.')
            cal = calibrate_extendedZolof(cam, 
                                          cam.image_points, 
                                          cam.lab_points)
            print('Starting calibration gui')
            gui = cal_gui(calibrate_obj=cal)
            #cal.calibrate()
            err = cal.mean_squared_err()
            print('Calibration finished. The calibration error is: %.3e'%err)
            #cam.save('.')
            
        
        else:
            models = str(['Tsai', 'extendedZolof'])[1:-1]
            msg = 'Unknown 3D model; permisible model names are: '
            raise ValueError(msg + models)                                
            
    
    
    
    
    def calibration_error_estimation(self):
        '''
        Performs stereo matching of the calibration points and compares
        then with the ground truth. 
        '''
        from numpy import loadtxt, array, mean, median, savez
        from myptv.imaging_mod import camera_wrapper, img_system
        from pandas import DataFrame
        
        cam_names = self.get_param('analyze_calibration_error', 'camera_names')
        cam_names = [val.strip() for val in cam_names.split(',')]
        plot = self.get_param('analyze_calibration_error', 'plot_histogram')
        
        
        # setting up the img_system 
        cams = [camera_wrapper(cn,'./') for cn in cam_names]
        for cam in cams:
            try:
                cam.load()
            except:
                raise ValueError('camera file %s not found'%cam.name)
        imsys = img_system(cams)
        
        
        # read calibration point files and organize in a dictionary
        point_dic = {}
        for e, cn in enumerate(cam_names):
            filename = './Calibration/%s_cal_points'%cn
            data = loadtxt(filename)
            for i in range(len(data)):
                try:
                    point_dic[tuple(data[i][2:])][e] = data[i][:2]
                except:
                    point_dic[tuple(data[i][2:])] = {e: data[i][:2]}
        
        
        # for each point in the dictionary, get the calibration error
        errors = []
        errsX, errsY, errsZ = [], [] ,[]
        x, y, z = [], [], []
        for k in point_dic.keys():
            if len(point_dic[k])!=len(cam_names): continue
            ground_truth = array(k)
            triangulation = imsys.stereo_match(point_dic[k], 1e20)[0]
            diff = triangulation - ground_truth
            err = sum((diff)**2)**0.5
            
            errors.append(err)
            errsX.append(diff[0])
            errsY.append(diff[1])
            errsZ.append(diff[2])
            x.append(ground_truth[0])
            y.append(ground_truth[1])
            z.append(ground_truth[2])
            
        
        print('Calibration error in in lab-space units:')
        print('RMS of full error: %.3e'%(mean(errors)))
        print('median of full error: %.3e'%(median(errors)))
        print('x error: %.3e'%(mean(abs(array(errsX)))))
        print('y error: %.3e'%(mean(abs(array(errsY)))))
        print('z error: %.3e'%(mean(abs(array(errsZ)))))
        print('')
        
        if plot == True:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(1,4)
            
            fig.suptitle('Errors in lab-space unites')
            titles = ['total error', 'x error', 'y error', 'z error']
            for e, lst in enumerate([errors, errsX, errsY, errsZ]):
                h = ax[e].hist(lst, bins='auto')
                ax[e].set_title(titles[e])
            plt.show()
        
        
        errs_df = DataFrame({'x':x, 'y':y, 'z':z,
                             'x_err':errsX, 'y_err':errsY, 'z_err':errsZ})
        
        errs_df.to_csv('calibration_errors_data', 
                       sep='\t', 
                       index=False,
                       float_format='%.4e')
        
        print('Calibration errors data saved as "calibration_errors_data" \n')
            
            
        
    
    
    def calibration_with_particles(self):
        '''
        This starts the calibrate with particles sequence
        '''
        from  matplotlib.pyplot import subplots, show
        
        # fetch parameters from the file
        camera_name =  self.get_param('calibration_with_particles',
                                      'camera_name')
        #resolution = self.get_param('calibration_with_particles',
        #                     'resolution').split(',')
        #resolution = (float(resolution[0]), float(resolution[1]))
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
        
        
        cam_file = open('./'+camera_name, 'r')
        model = cam_file.readline().split()[0]
        
        print('\n Camera model: %s'%model)
        
        
        if model == 'Tsai':
            from myptv.TsaiModel.gui_final_cal import cal_gui
            from myptv.TsaiModel.camera import camera_Tsai
            from myptv.TsaiModel.calibrate import calibrate_with_particles_Tsai
            
            # setting up a camera instance            
            cam = camera_Tsai(camera_name)
            cam.load('./')
            
            
            # set up the calibration object
            cal_with_p = calibrate_with_particles_Tsai(traj_filename, cam, 
                                                       cam_number, 
                                                       blobs_fname, 
                                                       min_traj_len=min_traj_len,
                                                       max_point_number=max_point_number)
            
            cal = cal_with_p.get_calibrate_instance()
            
            # run the final calibration gui
            print('starting calibration GIU using calibration with particles\n')
            gui = cal_gui(cal, cal_image) 
            
            
        
        if model == 'extendedZolof':
            from myptv.extendedZolof.gui_final_cal import cal_gui
            from myptv.extendedZolof.camera import camera_extendedZolof
            from myptv.extendedZolof.calibrate import calibrate_with_particles_EZ
            
            # setting up a camera instance            
            cam = camera_extendedZolof(camera_name)
            cam.load('./')
            
            
            # set up the calibration object
            cal_with_p = calibrate_with_particles_EZ(traj_filename, cam, 
                                                       cam_number, 
                                                       blobs_fname, 
                                                       min_traj_len=min_traj_len,
                                                       max_point_number=max_point_number)
            
            cal = cal_with_p.get_calibrate_instance()
            
            # run the final calibration gui
            print('starting calibration GIU using calibration with particles\n')
            gui = cal_gui(cal, cal_image) 
        
    
    
    def do_segmentation(self):
        '''
        Will perform segmentation on the images given in the parameters file
        and save the results on the given location.
        '''
        
        from numpy import zeros, amax
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
        shape = self.get_param('segmentation', 'shape')
        remove_BG = self.get_param('segmentation', 'remove_background')
        raw_format = self.get_param('segmentation', 'raw_format')
        
        
        # reading preprepared mask
        if type(mask)==str:
            mask = imread(mask)
        
        if method not in ['dilation', 'labeling']:
            raise ValueError('Method can be only "dilation" or "labeling"')
        
        if method=='dilation' and type(p_size) != int:
            raise ValueError('In dilation, particle_size can only be integer')
        
        if shape not in ['particles', 'fibers']:
            raise ValueError('Shape can be only "particles" or "fibers"')
        
        if raw_format==False:
            imread_func = lambda x: imread(x)
        else:
            import rawpy
            imread_func = lambda x: rawpy.imread(x).raw_image
        
        
        # get the shape of the images
        allfiles = os.listdir(dirname)
        n_ext = len(ext)
        image_files = sorted(list(filter(lambda s: s[-n_ext:]==ext, allfiles)))
        if single_img_name in image_files:
            image0 = imread_func(os.path.join(dirname,single_img_name))
        else:
            image0 = imread_func(os.path.join(dirname,image_files[0]))
        
        # preparing a mask using the given ROI
        if ROI is not None:
            ROI = [int(val) for val in ROI.split(',')]
            mask_ROI = zeros(image0.shape)
            mask_ROI[ROI[2]:ROI[3]+1, ROI[0]:ROI[1]+1] = 1
            mask = mask * mask_ROI
            mask = (mask / amax(mask)).astype('uint')
            
            
        def calculate_BG_image(dirname, extension):
            '''
            Calculates the background of images, defined as the median over a
            subsample of 200 images from the image folder.
            '''
            import os
            from skimage import io
            from numpy import median
            
            print('\ncalculating background...')
            
            allfiles = os.listdir(dirname)
            n_ext = len(extension)
            fltr = lambda s: s[-n_ext:]==extension
            image_files = sorted(list(filter(fltr, allfiles)))
            image_files = [os.path.join(dirname, fn) for fn in image_files]
            
            if len(image_files)<=200:
                ic = io.ImageCollection(image_files)
                
            else:
                ic = io.ImageCollection(
                            image_files[::int(len(image_files)/400+1)][:200])
            
            BG = median(ic, axis=0)
            
            return BG
            
        
        if shape=='particles':
            
            from myptv.segmentation_mod import loop_segmentation
            from myptv.segmentation_mod import particle_segmentation
        
            # segmenting the image if there are more than 1 frames
            if N_img is None or N_img>1:
                loopSegment = loop_segmentation(dirname, 
                                                particle_size=p_size,
                                                extension=ext,
                                                image_start=image_start,
                                                N_img=N_img,
                                                remove_ststic_BG=remove_BG,
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
                                                method=method,
                                                raw_format=raw_format)
            
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
                
                if type(remove_BG)==str:
                    print('\n','using given background image')
                    BG = imread(remove_BG)*1.0
                elif remove_BG==True:
                    print('\n','calculating background image')
                    BG = calculate_BG_image(dirname, ext)
                else:
                    BG=None
                    
                print('\n','segmenting image: %s'%single_img_name)
                particleSegment = particle_segmentation(image0, 
                                                        particle_size=p_size,
                                                        sigma=sigma, 
                                                        median=median,
                                                        BG_image=BG,
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
                
            
        elif shape=='fibers':
            
            from myptv.fibers.fiber_segmentation_mod import fiber_segmentation
            from myptv.fibers.fiber_segmentation_mod import loop_fiber_segmentation
            
            # segmenting the image if there are more than 1 frames
            if N_img is None or N_img>1:
                loopSegment = loop_fiber_segmentation(dirname, 
                                                particle_size=p_size,
                                                extension=ext,
                                                image_start=image_start,
                                                N_img=N_img, 
                                                sigma=sigma,
                                                remove_ststic_BG=remove_BG,
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
                                                method=method,
                                                raw_format=raw_format,
                                                pca_limit=1.0)
                
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
                            loopSegment.save_results_direction(save_name+'_directions')
                            print('\nfile saved.')
                        else:
                            print('\nskipped saving')
                        
                    else:
                        loopSegment.save_results(save_name)
                        loopSegment.save_results_direction(save_name+'_directions')
                        print('\nfile saved.')    
                print('\nDone.\n')
                
            
            
            # segmenting the image if there is only 1 frames
            if N_img == 1:
                print('\n','starting segmentation on a single image.')
                if single_img_name not in image_files:
                    in_ = os.path.join(dirname,single_img_name)
                    msg = 'Image %s not found in the directory.'%in_
                    raise ValueError(msg)
                
                if type(remove_BG)==str:
                    print('\n','using given background image')
                    BG = imread(remove_BG)*1.0
                elif remove_BG==True:
                    print('\n','calculating background image')
                    BG = calculate_BG_image(dirname, ext)
                else:
                    BG=None
                
                print('\n','segmenting image: %s'%single_img_name)
                particleSegment = fiber_segmentation(image0, 
                                                        particle_size=p_size,
                                                        sigma=sigma, 
                                                        median=median,
                                                        BG_image=BG,
                                                        threshold=threshold, 
                                                        local_filter=local_filter, 
                                                        max_xsize=max_xsize, 
                                                        max_ysize=max_ysize,
                                                        max_mass=max_mass,
                                                        min_xsize=min_xsize, 
                                                        min_ysize=min_ysize,
                                                        min_mass=min_mass,
                                                        mask=mask,
                                                        method=method,
                                                        pca_limit=1.0)
                
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
                            particleSegment.save_results_direction(save_name+'_directions')
                            print('\nfile saved.')
                        else:
                            print('\nskipped saving')
                    else:
                        particleSegment.save_results(save_name)
                        particleSegment.save_results_direction(save_name+'_directions')
                        print('\nfile saved.')
                    
                print('\nDone.\n')
            
            
            
            
    def do_matching(self):
        '''
        Will perform the stereo matching with the file given parameters
        '''
        from myptv.particle_matching_mod import matching_with_marching_particles_algorithm
        from myptv.imaging_mod import camera_wrapper, img_system
        from os import getcwd, listdir
        from os.path import exists as pathExists
        from time import localtime, strftime
        
        
        # fetching the parameters
        blob_fn = self.get_param('matching', 'blob_files')
        blob_fn = [val.strip() for val in blob_fn.split(',')]
        cam_names = self.get_param('matching', 'camera_names')
        cam_names = [val.strip() for val in cam_names.split(',')]
        # res = self.get_param('matching', 'cam_resolution')
        # res = tuple([float(val) for val in res.split(',')])
        ROI = self.get_param('matching', 'ROI').split(',')
        ROI = [float(ROI[i]) for i in range(6)]
        voxel_size = self.get_param('matching', 'voxel_size')
        N0 = self.get_param('matching', 'N0')
        max_err = self.get_param('matching', 'max_err')
        min_cam_match = self.get_param('matching', 'min_cam_match')
        frame_start = self.get_param('matching', 'frame_start')
        N_frames = self.get_param('matching', 'N_frames')
        march_forwards = self.get_param('matching', 'march_forwards')
        march_backwards = self.get_param('matching', 'march_backwards')
        save_name = self.get_param('matching', 'save_name')
        
        
        if N0==0 and voxel_size==None:
            raise ValueError('No initial guess method given (N0=0, voxel_size=None)')
        
        if min_cam_match<2:
            raise ValueError('min_cam_match needs to be at least 2.')
        
        # setting up the img_system 
        cams = [camera_wrapper(cn,'./') for cn in cam_names]
        for cam in cams:
            try:
                cam.load()
            except:
                raise ValueError('camera file %s not found'%cam.name)
        imsys = img_system(cams)
        
        
        mps = matching_with_marching_particles_algorithm(imsys, 
                                               blob_fn, 
                                               max_err, 
                                               ROI,
                                               N0,
                                               voxel_size,
                                               min_cam_match=min_cam_match,
                                               reverse_eta_zeta=True)

        
        
        # setting the frame range to match
        ts = int(mps.frames[0])
        te = int(mps.frames[-1])
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
        print('Starting stereo-matching at: ', strftime("%H:%M:%S", localtime()))
        
        if march_forwards==True:
            print('Matching forwards. Frames: %d -> %d'%(frames[0], frames[-1]))
            for f in frames:
                mps.match_frame(f)
                
        if march_backwards==True:
            print('\n','Matching backwards. Frames: %d -> %d'%(frames[-1], frames[0]))
            for f in frames[::-1]:
                mps.match_frame(f)
        
        
        
        # print matching statistics
        print('')
        print('Finished! \n')
        print('particles matched: %d \n'%(len(mps.matches)))
        
        Nframes = len(frames)
        c4 = sum([1 for p in mps.matches if len(p[1])==4]) / Nframes
        print('quadruplets: %.1f per frame'%c4)
        c3 = sum([1 for p in mps.matches if len(p[1])==3]) / Nframes
        print('triplets: %.1f per frame'%c3)
        c2 = sum([1 for p in mps.matches if len(p[1])==2]) / Nframes
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
                    mps.save_particles(save_name)
                else:
                    print('\n','skiped saving.')
                
            else:
                print('\n','saving file.')
                mps.save_particles(save_name)
        
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
    
        
    
    
    
    
    def do_orientations(self):
        '''
        A part of Eric Aschari's Fiber tracking extension (MyFTV):
            
        Will perform a fiber orientation analysis
        '''
        from numpy import loadtxt, empty, array, zeros, pi, sign, savetxt
        from numpy import abs as absnp
        from numpy.linalg import norm
        from pandas import read_csv
        from myptv.fibers.fiber_orientation_mod import FiberOrientation
        from myptv.imaging_mod import camera, img_system
        # from myptv.particle_matching_mod import match_blob_files
    
        
        # fetching the parameters
        blob_fn = self.get_param('fiber_orientations', 'blob_files')
        blob_fn = [val.strip() for val in blob_fn.split(',')]
        cam_names = self.get_param('fiber_orientations', 'camera_names')
        cam_names = [val.strip() for val in cam_names.split(',')]
        res = self.get_param('fiber_orientations', 'cam_resolution')
        res = tuple([float(val) for val in res.split(',')])
        trajectory_file = self.get_param('fiber_orientations','trajectory_file')
        
        save_name = self.get_param('fiber_orientations', 'save_name')
        print(save_name)
        # setting up the img_system 
        cams = [camera(cn, res) for cn in cam_names]
        for cam in cams:
            try:
                cam.load('')
            except:
                raise ValueError('camera file %s not found'%cam.name)
        imsys = img_system(cams)
        # fibers = loadtxt(fibers_file)
        trajectories = loadtxt(trajectory_file)
        
        oris = empty((len(trajectories[:,0]),8))
        shape = (3,1)
        blobs = []
        for fn in blob_fn:
            blobs.append(array(read_csv(fn, sep='\t', header=None)))

        count = 0
        for frame in range(int(trajectories[-1,-1])+1):
            
            frame_trajectories = [line for line in trajectories if line[7] == frame]

            for line in range(len(frame_trajectories)):    
                X = array([zeros(shape),zeros(shape)])
                B = array([zeros(shape),zeros(shape)])
                
                n_x_prev = 0
                n_y_prev = 0
                n_z_prev = 0
                
                for cam in range(2):
                    
                    frame_blobs = [line for line in blobs[cam] if line[5] == frame]
                    reso = cams[cam].resolution[0]
                    index = int(frame_trajectories[line][4+cam])
                    
                    x_corr = cams[cam].xh
                    y_corr = cams[cam].yh
                    
                    temp1_x = (frame_blobs[index][0]) - (reso/2. + x_corr)
                    temp1_b = temp1_x + frame_blobs[index][6]*100
                    
                    temp2_x = (frame_blobs[index][1]) - (reso/2. + y_corr)
                    temp2_b = temp2_x + frame_blobs[index][7]*100
                    
                    # coordinate transformation
                    X[cam][0,0] = temp2_x
                    B[cam][0,0] = temp2_b
                    
                    X[cam][1,0] = temp1_x
                    B[cam][1,0] = temp1_b
                    
                # get orientations
                o = FiberOrientation(X, B)
                c,u,ori = o.image2fiber(imsys.cameras)
                ori = ori/pi*180
                u /= norm(u)
                
                if line == 0:
                    n_x_prev = u[0]
                    n_y_prev = u[1]
                    n_z_prev = u[2]
                
                # correct sign
                value = 0.2
                if absnp(absnp(u[0]) - absnp(n_x_prev)) < value and sign(u[0]) != sign(n_x_prev):
                    u[0] *= -1
                if absnp(absnp(u[1]) - absnp(n_y_prev)) < value and sign(u[1]) != sign(n_y_prev):
                    u[1] *= -1
                if absnp(absnp(u[2]) - absnp(n_z_prev)) < value and sign(u[2]) != sign(n_z_prev):
                    u[2] *= -1
                
                n_x_prev = u[0]
                n_y_prev = u[1]
                n_z_prev = u[2]
            
                temp = frame_trajectories[line]

                temp[1] = u[0]
                temp[2] = u[1]
                temp[3] = u[2]

                oris[count] = temp
                count += 1

                
                
        # saving the data
        if save_name is not None:
            print('\n', 'Saving the data.')    
            # o.save_results(save_name, oris)
            savetxt(save_name, oris, fmt = 
                       ['%d', '%.3f', '%.3f', '%.3f', '%d', '%d', '%.2f', '%.2f'], delimiter='\t')
        print('\n', 'Done.')
    
    
    
    
    
    def do_plot_trajectories(self):
        '''
        This function is used to generate a 3D plot of the trajectories in 
        a given file.
        '''
        from myptv.makePlots.plot_trajectories import plot_trajectories
        
        # fetching the parameters
        file_name = self.get_param('plot_trajectories', 'file_name')
        min_length = self.get_param('plot_trajectories', 'min_length')
        write_trajID = self.get_param('plot_trajectories', 'write_trajID')
        t0 = self.get_param('plot_trajectories', 't0')
        te = self.get_param('plot_trajectories', 'te')
        
        plot_trajectories(file_name, 
                          min_length, 
                          write_trajID=write_trajID, 
                          t0=t0, 
                          te=te)
    
    
    
    
        
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
        
    
    
    
    # ========================================================================
    # /\/\//\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
    #  Legacy functions that are no longer needed due to the cal_gui
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
    
    
#%%
        
        
        


if __name__ == '__main__':
    
    import sys
    fname, action = sys.argv[1], sys.argv[2]
    print('\n','given inputs -')
    print('params file name:', fname)
    print('action:', action, '\n')
    wf = workflow(fname, action)
    
    
    
    
    

