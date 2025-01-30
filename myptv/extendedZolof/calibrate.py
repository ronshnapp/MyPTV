# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 12:15:21 2019

@author: ron


calibration module for Extended Zolof camera instances. This is used to 
obtain the [A], [B], and O Extended Zolof model parameters.

"""

from numpy import array, dot
from numpy import sum as npsum
from numpy.linalg import lstsq, norm
from scipy.optimize import minimize

from myptv.extendedZolof.camera import camera_extendedZolof
from myptv.utils import line, get_nearest_line_crossing

from pandas import read_csv



class calibrate_extendedZolof(camera_extendedZolof):
    '''
    This object is used to calibrate cameras against a given list
    of lab and camera point coordinates. 
    '''
    
    
    def __init__(self, camera, x_list, X_list, quadratic=False):
        '''
        Given a list of 2D points, x=(x,y), and a list of 3D point X=(X,Y,Z), 
        we assume that given a point X, we can compute x by a polynomial 
        of degree 3, as - 
        
        x = A0 + A1*X + A2*Y + A3*Z +
            A4*X^2 + A5*Y^2 + A6*Z^2 + A7*XY + A8*YZ + A9*ZX + A108XYZ
            A11*XY^2 + A12*XZ^2 + A13*YX^2 + A14*YZ^2 + A15*ZX^2 + A16*ZY^2  
            
        if quadratic==True, then only the quadratic terms are used. This is 
        used in the initial calibration.
        '''
        self.cam = camera
        self.A = [[0.0 for i in range(17)] for j in [0,1]]
        self.B = [[0.0 for i in range(10)] for j in [0, 1, 2]]
        self.x_list = x_list
        self.X_list = X_list
        self.quadratic = quadratic
        
        
        
    
    def calibrate(self):
        '''
        Given a list of points, x and X, this function attempts to determine
        the A coefficients. 
        '''
        # 1) finding the A coefficients - 
        if self.quadratic==False:
            XColumns = [self.cam.get_XCol(Xi) for Xi in self.X_list]
        
        elif self.quadratic==True:
            XColumns = []
            for Xi in self.X_list:
                Xcol_i = self.cam.get_XCol(Xi)
                # (Here, -7 is for quadratic, and -13 is linear)
                for i in range(-7,0): Xcol_i[i] = 0  
                XColumns.append(Xcol_i)
        
        res = lstsq(XColumns, self.x_list, rcond=None)
        self.A = res[0]
        
        # 2) finding the best camera center -
        line_list = []
        for i in range(0, len(self.X_list)):
            O, e = self.get_ray_from_x(self.x_list[i], X0=self.X_list[i])
            line_list.append(line(O, e)) 
        self.O = get_nearest_line_crossing(line_list)
        
        # 3) finding the unit vector for each X -
        r_list = []
        for Xi in self.X_list:
            r = (Xi - self.O)/norm(Xi - self.O)
            r_list.append(r)
        
        # 4) finding the B coefficients -
        xColumns = [self.cam.get_xCol(xi) for xi in self.x_list]
        res = lstsq(xColumns, r_list, rcond=None)
        self.B = res[0]
        
        self.cam.O = self.O
        self.cam.A = self.A
        self.cam.B = self.B
        
        
            
    def get_ray_from_x(self, x, X0=None):
        '''
        Given a point in 2D image space, this function returns a line in 3D
        that passes through this point. The line is represented with six 
        parameters: one point in 3D, O, and one unit vector in 3D, e.
        '''
        
        func = lambda X: sum((array(self.projection(X)) - array(x))**2)
        
        if X0 is None:
            X0 = array([0,0,0])
        
        X02 = array(X0) + array([1,1,1])
            
        O = minimize(func, X0).x
        dX = minimize(func, X02).x
        e = (O-dX)/sum((O-dX)**2)**0.5
        
        return O, e
        
    
    
    def mean_squared_err(self):
        '''
        Calculates and returns the mean square of the deviations in camera
        space.
        '''
        errorsSquard = []
        
        for i in range(len(self.X_list)):
            xProj = dot(self.get_XCol(self.X_list[i]), self.cam.A)
            errorsSquard.append( norm(array(xProj)-array(self.x_list[i]))**2 )
        
        return (sum(errorsSquard)/len(errorsSquard))**0.5
        
        
        
        
    def plot_err_distribution(self, ax = None):
        import matplotlib.pyplot as plt
        from numpy import sum as npsum
        
        if ax == None:
            fig, ax = plt.subplots()
        
        imc = array(self.x_list)
        z_lst = array([self.cam.projection(x) for x in self.X_list])
        err = npsum((imc-z_lst)**2, axis=1)**0.5
        
        h = ax.hist( err, bins='auto')
        ax.set_xlabel('Camera projection err [px]')
        ax.set_ylabel('Counts')
        
        
        
        
    def plot_proj(self, ax = None):
        import matplotlib.pyplot as plt
        
        if ax == None:
            fig, ax = plt.subplots()
        
        imc = array(self.x_list)
        ax.plot(imc[:,0], imc[:,1], 'ob')
        for i in range(imc.shape[0]):
            ax.text(imc[i,0], imc[i,1], '%d'%i, color = 'b')
        
        z_lst = array([self.cam.projection(x) for x in self.X_list])
        ax.plot( z_lst[:,0], z_lst[:,1], 'xr' )
        for i in range(z_lst.shape[0]):
            ax.text(z_lst[i,0], z_lst[i,1], '%d'%i, color = 'r')
            
        ax.set_aspect('equal')
        
        
        
        
        
        
        
        
        
        
        
        
        

class calibrate_with_particles_EZ(object):
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
                 min_traj_len = 15, max_point_number = 1000):
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
                    
                    all_valid_points.append( (err,(point, blob_coords[::-1])) )
        
        
        if len(all_valid_points)==0:
            msg1 = 'No trajectories were found. min_traj_len may be too high.'
            raise ValueError(msg1)
        
        
        # sort the points according to their triangulation error
        all_valid_points = sorted(all_valid_points, key=lambda x: x[0])
        
        # get the best N points into a final list
        self.cal_points = [p[1] 
                           for p in all_valid_points[:self.max_point_number]]


    def get_particle_disparity(self):
        '''
        Returns a list with the discrepancies between the segmented blob 
        positions and the projections of their 3D positions onto the camera.
        '''
        cam = self.camera
        disparities = [cam.projection(p[0]) - p[1] for p in self.cal_points]
        return disparities
        

    def get_calibrate_instance(self):
        
        # initiating a calibrate object using this data
        self.cal = calibrate_extendedZolof(self.camera.camera, 
                                           [p[1] for p in self.cal_points], 
                                           [p[0] for p in self.cal_points])
        return self.cal
        
        