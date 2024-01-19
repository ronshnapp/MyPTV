# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 12:15:21 2019

@author: ron


calibration module for the cameras

"""

from numpy import mean, sum, hstack, array, loadtxt, zeros, ones, arange
from pandas import read_csv



class calibrate_Tsai(object):
    '''
    This object is used to calibrate cameras against a given list
    of lab and camera point coordinates. 
    The main problem is to minimize the distance between the projected
    lab coordinates and the given image coordinates. The minimization 
    procedure uses scipy minimize method (through searchCalibration()).
    '''
    
    def __init__(self, camera, lab_coords, img_coords, random_sampling=20):
        '''
        input -
        camera - the camera instance to be calibrated
        lab_coords - a list of 1d arrays of length 3, giving the calibration 
                     points coordinates in lab space
        lab_coords - a list of 1d arrays of length 2, giving the calibration 
                     points coordinates in the camera's pixels
        random_sampling - an integer, using self.stochastic_search_calibration
                          in each iteration of the minimization the cost
                          function is calculated with a given subset of 
                          calibration points with this number of samples. 
        
        '''
        self.camera = camera
        self.img_coords = img_coords
        self.lab_coords = lab_coords
        self.D_lst = [self.mean_squared_err()]
        self.sep = sum((array(self.img_coords[1])-array(self.img_coords[0]))**2)**0.5
        self.random_sampling = random_sampling
        
        
        
    def mean_squared_err(self, correction=True, points=None):
        '''
        This calculates the mean squared distance between the 
        projection and the given coordinates (in units of pixel).

        (in the calibration we want to minimize this D)
        
        if points=None, then the calculation is using all the calibratino 
        points. Otherwise, it can be a tuple of two lists, the first being a 
        list of lab_points and the second is a list of img_points. 
        '''
        
        if points is None:
            lp, imp = self.lab_coords, self.img_coords
        else:
            lp, imp = points
            
        z_lst = []        
        for x in lp:
            z_lst.append(self.camera.projection(x, correction=correction))
        e = array(z_lst) - array(imp)
        D = mean( sum(e**2, axis=1)**0.5 )
        
        return D
        
    
    
    def searchCalibration(self, maxiter=5000, fix_f=True, points=None):
        '''
        using scipy minimize function to obtain calibration
        parameters for the camera.
        '''
        from scipy.optimize import minimize
        
        
        # 1) find center and rotation angles
        
        def func(X):
            self.camera.O = X[:3]
            self.camera.theta = X[3:6]
            self.camera.xh = X[6]
            self.camera.yh = X[7]
            
            if not fix_f:
                self.camera.f = X[-1]
                
            self.camera.calc_R()
            
            meanSquaredErr = self.mean_squared_err(points=points)
            self.D_lst.append( meanSquaredErr )
            return meanSquaredErr
        
        c = self.camera
        
        if fix_f:
            X0 = hstack([c.O, c.theta, c.xh, c.yh])
        
        else:
            X0 = hstack([c.O, c.theta, c.xh, c.yh, c.f])
            
        if (self.camera.E == 0).all():
            res = minimize(func, X0, method='BFGS',
                            options={'maxiter': maxiter},
                            jac = '2-point')
            
        else:
            res = minimize(func, X0, method='nelder-mead',
                            options={'maxiter': maxiter})
        
        return res.message, res.x
    
    
    
    def stochastic_searchCalibration(self, iterSteps=2000):
        '''
        This function divides the calibration points into subsets with 
        self.random_sampling points each, and then performs iterSteps 
        minimization steps using searchCalibration for each point subset.
        '''
        from random import shuffle
        
        print('\nRunning stochastic self-calibration\n')
        
        nPointsTot = len(self.lab_coords)
        
        
        # in case there are not enough calibraiton points
        if nPointsTot < self.random_sampling + 1:
            print('Number of calibration points too small for stochastic ')
            print('calibration. Falling back to regular minimization.')
            self.searchCalibration()
            return
            
        # if there are sufficient amount of points -     
        # We do the minimization over many different subsets
        # of increasing sizes
        
        subset_sizes = [self.random_sampling]
        while subset_sizes[-1]*2<nPointsTot:
            subset_sizes.append(int(subset_sizes[-1]*2))
        
        
        err_0 = self.mean_squared_err()
        ee = 0
        for subset_size in subset_sizes:
            print('\niterating at subset size %d'%subset_size)
            
            itermax = int((iterSteps * self.random_sampling)/subset_size)
            
            # get random subsets of the points -
            nSubsets = int(nPointsTot / subset_size)
            subsets = [([],[]) for i in range(nSubsets)]
            
            shuffledIndexes = list(range(nPointsTot))
            shuffle(shuffledIndexes)
            
            for e, ind in enumerate(shuffledIndexes):
                subsets[e%nSubsets][0].append(self.lab_coords[ind])
                subsets[e%nSubsets][1].append(self.img_coords[ind])
            
            # Run minimization steps using the subsets:
            e=0
            err_i = self.mean_squared_err()
            Xi = hstack([self.camera.O, self.camera.theta,
                          self.camera.xh, self.camera.yh])
            convergenceTracker = [err_i]
            for subset in subsets:
                e+=1
                print('starting subset %d/%d; err=%.3f'%(e, len(subsets), 
                                                self.mean_squared_err()))
                self.searchCalibration(maxiter=itermax, points=subset)
                err_ip1 = self.mean_squared_err()
                Xip1 = hstack([self.camera.O, self.camera.theta,
                                self.camera.xh, self.camera.yh])
                
                # if the error increased, go back to previous state
                if err_ip1 > err_i:
                    self.camera.O = Xi[:3]
                    self.camera.theta = Xi[3:6]
                    self.camera.xh = Xi[6]
                    self.camera.yh = Xi[7]
                    self.camera.calc_R()
                    convergenceTracker.append(0.0)
                    
                else:
                    convergenceTracker.append(abs((err_i - err_ip1)/err_i))
                    err_i = err_ip1
                    Xi = Xip1
                
                # check for converence of the error
                if ee>0:
                    if all([conv < 0.005 for conv in convergenceTracker[-3:]]):
                        print('\nreached dead end. Increasing subset size.')
                        break                
                ee+=1
                
        err = self.mean_squared_err()
        print('\n','Error reduced by: %.2f%%'%( (err_0-err)/err_0*100 ))
        print('stochastic iteration done; err=%.3f'%(err))
        
    
    
    def fineCalibration(self, maxiter=500, points=None):
        '''
        Calibration for the nonlinear error term. 
        This function attempts to find the 27 parameters that minimize the
        calibration error using scipy.minimize.
        '''
        from scipy.optimize import minimize
        
        def func(X):
            shape = (2, self.camera.E.shape[1])
            X_ = X.reshape(shape)
            self.camera.E[0,:] = X_[0,:]
            self.camera.E[1,:] = X_[1,:]
            meanSquaredErr = self.mean_squared_err(correction=True, 
                                                   points=points )
            self.D_lst.append( meanSquaredErr)
            return meanSquaredErr
        
        c = self.camera
        X0 = c.E[:2,:].flatten()
        res = minimize(func, X0, method='nelder-mead', 
                       options={'maxiter': maxiter})
        
        return res.message, res.x
    
    
    def stochastic_fineCalibration(self, iterSteps=2000):
        
        from random import shuffle
        
        nPointsTot = len(self.lab_coords)
        subset_sizes = [self.random_sampling]
        while subset_sizes[-1]*2<nPointsTot:
            subset_sizes.append(int(subset_sizes[-1]*2))
        
        err_0 = self.mean_squared_err()
        ee = 0
        for subset_size in subset_sizes:
            
            itermax = int((iterSteps * self.random_sampling)/subset_size)
            
            # get random subsets of the points -
            nSubsets = int(nPointsTot / subset_size)
            subsets = [([],[]) for i in range(nSubsets)]
            
            shuffledIndexes = list(range(nPointsTot))
            shuffle(shuffledIndexes)
            
            for e, ind in enumerate(shuffledIndexes):
                subsets[e%nSubsets][0].append(self.lab_coords[ind])
                subsets[e%nSubsets][1].append(self.img_coords[ind])
            
            # Run minimization steps using the subsets:
            e=0
            err_i = self.mean_squared_err()
            Xi = self.camera.E[:2,:].copy()
            convergenceTracker = [err_i]
            for subset in subsets:
                e+=1
                print('starting subset %d/%d; err=%.3f'%(e, len(subsets), 
                                                self.mean_squared_err()))

                self.fineCalibration(maxiter=itermax, points=subset)
                err_ip1 = self.mean_squared_err()
                Xip1 = self.camera.E[:2,:].copy()
                
                # if the error increased, go back to previous state
                if err_ip1 > err_i:
                    self.camera.E[0,:] = Xi[0,:]
                    self.camera.E[1,:] = Xi[1,:]
                    convergenceTracker.append(0.0)
                    
                else:
                    convergenceTracker.append(abs((err_i - err_ip1)/err_i))
                    err_i = err_ip1
                    Xi = Xip1
                    
                # check for converence of the error
                if ee>0:
                    if all([conv < 0.005 for conv in convergenceTracker[-4:]]):
                        print('\nreached dead end. Increasing subset size.')
                        break
            ee+=1
                
        err = self.mean_squared_err()
        print('\n','Error reduced by: %.2f%%'%( (err_0-err)/err_0*100 ))
        print('\n','stochastic iteration done; err=%.3f'%(self.mean_squared_err()))
                
        
        
# =============================================================================
#    An attempt to calculate the error gradient numerically. This turns out
#    to be somwhat of a challange, and the speed does not increase in the
#    minimization. So, we will drop this for now.     
#
#     def errorGradient(self, X):
#         '''
#         returns the gradient of the mean squared error. 
#         '''
#         from numpy import sin, cos, dot
#         
#         # the inverse of the rotation matrix is its transpose
#         R_inv = self.camera.R.T
#         
#         # calculating the derivatives of the rotation matrix
#         tx,ty,tz = self.camera.theta
#         Rx = array([[1,0,0],
#                     [0,cos(tx),-sin(tx)],
#                     [0,sin(tx),cos(tx)]])
#         Ry = array([[cos(ty),0,sin(ty)],
#                      [0,1,0],
#                      [-sin(ty),0,cos(ty)]])
#         Rz = array([[cos(tz),-sin(tz),0],
#                     [sin(tz),cos(tz),0],
#                     [0,0,1]])
#         
#         dxRx = array([[1,0,0],
#                       [0,-sin(tx),-cos(tx)],
#                       [0,cos(tx),-sin(tx)]])
#         dyRy = array([[-sin(ty),0,cos(ty)],
#                       [0,1,0],
#                       [-cos(ty),0,-sin(ty)]])
#         dzRz = array([[-sin(tz),-cos(tz),0],
#                       [cos(tz),-sin(tz),0],
#                       [0,0,1]])
#         
#         dxR = dot(dot(dxRx,Ry), Rz).T
#         dyR = dot(dot(Rx,dyRy), Rz).T
#         dzR = dot(dot(Rx,Ry), dzRz).T
#         
#         a = dot(-self.camera.O, R_inv)[2] / self.camera.f
#         
#         
#         # calculating the gradients by summing over all calibration points
#         grad_E = zeros(8)
#         for i in range(len(self.lab_coords)):
#             eta, zeta = self.img_coords[i]
#             eta_, zeta_ = self.camera.projection(self.lab_coords[i])
#             ei = ((eta-eta_)**2 + (zeta-zeta_)**2)**0.5
#             
#             # derivative of O
#             grad_E[0] += -((eta_-eta)*R_inv[0,0] + (zeta_-zeta)*R_inv[0,1])/ei/a
#             grad_E[1] += -((eta_-eta)*R_inv[1,0] + (zeta_-zeta)*R_inv[1,1])/ei/a
#             #grad_E[2] += -((eta_-eta)*R_inv[2,0] + (zeta_-zeta)*R_inv[2,1])/ei/a
#             
#             dx = self.lab_coords[i] - self.camera.O
#             deta_dOz = (R_inv[0,2]*dot(dx, R_inv[:,2])-R_inv[2,2]*dot(dx, R_inv[:,0]))/dot(dx, R_inv[:,2])**2
#             dzeta_dOz = (R_inv[1,2]*dot(dx, R_inv[:,2])-R_inv[2,2]*dot(dx, R_inv[:,1]))/dot(dx, R_inv[:,2])**2
#             grad_E[2] += -((eta_-eta)*deta_dOz + (zeta_-zeta)*dzeta_dOz)/ei*self.camera.f
#             
#             # derivative of theta
#             Ax = (self.lab_coords[i] - self.camera.O).dot(dxR)
#             grad_E[3] += -((eta-eta_)*Ax[0] + (zeta-zeta_)*Ax[1])/ei/a
#             Ay = (self.lab_coords[i] - self.camera.O).dot(dyR)
#             grad_E[4] += -((eta-eta_)*Ay[0] + (zeta-zeta_)*Ay[1])/ei/a
#             Az = (self.lab_coords[i] - self.camera.O).dot(dzR)
#             grad_E[5] += -((eta-eta_)*Az[0] + (zeta-zeta_)*Az[1])/ei/a
#         
#             # derivative of xh and yh
#             grad_E[6] += ((eta-eta_) + (zeta-zeta_))/ei
#             grad_E[7] += ((eta-eta_) + (zeta-zeta_))/ei
#             
#         #print(grad_E / len(self.lab_coords))
#         return grad_E / len(self.lab_coords)
# =============================================================================
    
    
    
    def plot_proj(self, ax = None):
        import matplotlib.pyplot as plt
        
        if ax == None:
            fig, ax = plt.subplots()
        
        imc = array(self.img_coords)
        ax.plot(imc[:,0], imc[:,1], 'ob')
        for i in range(imc.shape[0]):
            ax.text(imc[i,0], imc[i,1], '%d'%i, color = 'b')
        
        z_lst = array([self.camera.projection(x) for x in self.lab_coords])
        ax.plot( z_lst[:,0], z_lst[:,1], 'xr' )
        for i in range(z_lst.shape[0]):
            ax.text(z_lst[i,0], z_lst[i,1], '%d'%i, color = 'r')
            
        ax.set_aspect('equal')
        
        
    def plot_err_distribution(self, ax = None):
        import matplotlib.pyplot as plt
        
        if ax == None:
            fig, ax = plt.subplots()
        
        imc = array(self.img_coords)
        z_lst = array([self.camera.projection(x) for x in self.lab_coords])
        err = sum((imc-z_lst)**2, axis=1)**0.5
        
        h = ax.hist( err, bins='auto')
        ax.set_xlabel('Camera projection err [px]')
        ax.set_ylabel('Counts')
        
        
    def manual_calibration(self):
        '''
        A manual calibration app - manually change calibration
        parameters and see the error changing live.
        
        ===================================================
        IN DEVELOPMENT: Currently it only supports actively 
        changing camera center and angle. 
        To add: live plotting.
        ===================================================
        
        keys:
            a - change Ox
            s - change Oy
            d - change Oz
            
            z - change Theta_x            
            x - change Theta_y
            c - change Theta_z
        '''       
        print('manual calibration; enter q to quit.')
        cmd = 'n'
        known_command = ['a','z','s','x','d','c']
        
        while cmd != 'q':
            cmd = input('enter calibration command:')
            if cmd in known_command:
                increment = float(input('choose increments: '))
                if cmd == 'a':
                    self.camera.O[0] += increment
                if cmd == 'z':
                    self.camera.theta[0] += increment
                if cmd == 's':
                    self.camera.O[1] += increment
                if cmd == 'x':
                    self.camera.theta[1] += increment
                if cmd == 'd':
                    self.camera.O[2] += increment
                if cmd == 'c':
                    self.camera.theta[2] += increment
                
                D = self.mean_squared_err()
                self.D_lst.append( D )
                print('D = %f'%D)
                
            elif cmd == 'q':
                print('quitting')
                break
                
            else:
                print('unknown command \n')
                
        self.plot_proj()
        
        
        







class calibrate_with_particles_Tsai(object):
    '''
    A class used to refine the calibration using particles data. In short,
    after the primary clibration is done, matching and tracking can be used
    to obtain trajectories from the experimental data. Here, we can leverage
    the trajectories obtained to minimize further the calibration error.
    The assumption is that longer trajectories are considered more reliable 
    as compared to shorter trajectories, so we use only "long" trajectories 
    in this process.
    '''
    
    def __init__(self, traj_filename, camera, cam_number, blobs_fname, 
                 min_traj_len = 15, max_point_number = 400):
        '''
        input -
        
        traj_filename - the name of the file containing the trajectories from
                        which the calibration points are taken.
        
        camera - an instance of the camera we wish to try and re-calibrate
        
        cam_number - int, >= 1; the number (index) of the camera to be
                     calibrated. For example, if this is 1 then it takes
                     the blobs from the 4th column of the trajectory file.
        
        blobs_fname - The name of the file that contains the segmented 
                      particles' data.
        
        min_traj_len - only trajectories longer then this number will be used
                       in the calibration
                       
        max_point_number - the maximum number of points that shall be taken
                           to re-calibrate the camera. Note that too many 
                           points might lead to long calculation times.
        '''
        
        self.traj_fname = traj_filename
        self.camera = camera
        self.cam_number = cam_number
        self.blobs_fname = blobs_fname
        self.min_traj_len = min_traj_len
        self.max_point_number = max_point_number
        
        # gathering points from trajectory file and matching them to blobs
        self.fetch_points()
        
        
    def fetch_points(self):
        '''
        This will fetch the calibration points from the trajectory file
        '''
        print('fetching data from trajectories and blobs file...')
        
        # load the trajectories and sort them in a dictionary according to id
        fltr = lambda k,g: k!=-1 and len(g)>=self.min_traj_len
        tr_data = read_csv(self.traj_fname, delimiter='\t', header=None)
        self.trajs = dict([(k,array(g)) for k, g in tr_data.groupby(by=0) 
                                                                if fltr(k,g)])

        # load blobs an arrange in lists according to their frame
        blob_data = read_csv(self.blobs_fname, delimiter='\t', header=None)
        self.blobs = dict([(k,array(g)) for k, g in blob_data.groupby(by=5)])

        # the index of the blobs column in the trajectories file
        ind = self.cam_number + 3
        
        skip = int(self.min_traj_len/4)
        
        # extract the data from blobs and trajs dictionaries
        all_valid_points = []
        for k in self.trajs.keys():
            if len(self.trajs[k])>=self.min_traj_len:
                for p in self.trajs[k][::skip]:     # <--- taking only once some 
                    blob_num = int(p[ind])     #      data points in each traj
                    if blob_num==-1:continue
                    frame = p[-1]
                    err = p[-2]
                    blob_coords = self.blobs[frame][blob_num][:2]
                    
                    point = p[1:4]
                    
                    #self.cal_points.append((point, blob_coords[::-1]))
                    all_valid_points.append( (err,(point, blob_coords[::-1])) )
        
        
        if len(all_valid_points)==0:
            msg1 = 'No trajectories were found. min_traj_len may be too high.'
            raise ValueError(msg1)
        
        
        # sort the points according to their triangulation error
        all_valid_points = sorted(all_valid_points, key=lambda x: x[0])
        
        # get the best N points into a final list
        self.cal_points = [p[1] 
                           for p in all_valid_points[:self.max_point_number]]

        

    def get_calibrate_instance(self):
        
        # initiating a calibrate object using this data
        self.cal = calibrate(self.camera, 
                             [p[0] for p in self.cal_points], 
                             [p[1] for p in self.cal_points])
        return self.cal
        
        
        