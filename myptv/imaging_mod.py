# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 18:02:07 2018

@author: ron


Imaging Module:

This model contains the img_system and camera_wrapper classes used to 
handle the tasks of 2D to 3D transformations.

"""

import os

from numpy import array, eye
from numpy.linalg import solve, LinAlgError
from myptv.utils import line_dist
from myptv.TsaiModel.calibrate import calibrate_Tsai
from myptv.TsaiModel.camera import camera_Tsai
from myptv.extendedZolof.camera import camera_extendedZolof




class img_system(object):
    '''
    an object that holds a number of camera wrappers and can be used to 
    perform stereo matching operations.
    '''
    
    def __init__(self, camera_list):
        '''
        input - camera_list is a list of loaded camera wrappers.
        '''
        self.cameras = camera_list
    
    
    def stereo_match(self, coords, d_max, strict_match=False):
        '''
        given n particle images [(eta, zeta) coords in camera space], this will
        determine whereather there is a good candidate point for the intersection
        of epipolar lines, and if so returns it. Here it is assumed that
        all images correspond to the same "real world" particle.
        This point is estimated as the average of the crossing point of the
        epipolar lines, that cross at distances smaller than a maximum value.
        
        The math: for lines with origin O_i and unit direction r_i, X
        minimizing sum_i ||(I - r_i r_i^T)(X - O_i)||^2 satisfies
            A X = b
            A = sum_i (I - r_i r_i^T) = N*I - R^T R
            b = sum_i (I - r_i r_i^T) O_i
        where R is the (N,3) matrix of stacked direction vectors.
        
        input - 
        coords (dic) - keys are camera number, values are the image space 
                       coordinates of each point. Must have at least 2 entries
        d_max (float) - maximum allowable distance separating two lines
        strict_match (bolean) - If this is True, returns non-None answer only
                                if the distance from the point to all epipolar
                                lines is at most d_max.
        
        output - either -
        X (numpy array, 3) - lab space coordinates of the sought point 
        cams (list) - list of camera indexes for which the point was found
        dist - average distance to crossing points
        
        or - (if all epipolar lines )
        None
        '''

        keys = list(coords.keys())
        N = len(keys)
        if N < 2:
            return None
 
        epilines = {k: self.cameras[k].get_epipolarline(coords[k][0], coords[k][1])
                    for k in keys}
 
        # ---- N == 2: lean direct two-line path (no O(N^2) to remove here) ----
        if N == 2:
            O1, r1 = epilines[keys[0]]
            O2, r2 = epilines[keys[1]]
 
            D, X = line_dist(O1, r1, O2, r2)
            if D > d_max:
                return None
 
            return X, set(keys), D
 
        # ---- N >= 3: closed-form multi-line least-squares solve ----
        O = array([epilines[k][0] for k in keys])
        R = array([epilines[k][1] for k in keys])
        R = R / ((R**2).sum(axis=1, keepdims=True))**0.5  # defensive re-normalize, O(N)
 
        A = N*eye(3) - R.T.dot(R)
        c = (R*O).sum(axis=1)
        b = O.sum(axis=0) - R.T.dot(c)
 
        try:
            X = solve(A, b)
        except LinAlgError:
            return None
 
        diff = X[None, :] - O
        perp = diff - (diff*R).sum(axis=1, keepdims=True)*R
        dists = (perp**2).sum(axis=1)**0.5
 
        if strict_match:
            if (dists > d_max).any():
                return None
            return X, set(keys), dists.mean()
 
        else:
            passed = dists <= d_max
            if passed.sum() < 2:
                return None
            cams = set(array(keys)[passed])
            return X, cams, dists[passed].mean()
        
        





class camera_wrapper(object):
    '''
    A camera is an object that can transform epipolar lines to pixels
    and pixel coordinates into epipolar lines. There are several methods that
    could be used to obtain these transformations depending on the 3D model used.
    The only requirement is that the camera is "calibrated".
    Operationally, different 3D models are used via individual camera classes.
    The "calibration" of a camera is represented by a file saved on the disk. 
    The camera wrapper is a class that wraps around various types of cameras
    and can handle the tasks of epipolar lines <-> pixel transformations. 
    '''
    
    def __init__(self, fileName, dirPath):
        '''
        Input:
            
        fileName (string) - this is the name (path) to a file that contains
                            the calibration parameters related to its 3D model.
        
        dirPath (string) - path of the directory in which fileName is found.
        '''
        self.fileName = fileName
        self.dir = dirPath
        self.ListOfModels = ['Tsai', 'extendedZolof']
        self.camera = None
        
    
    def __repr__(self):
        msg1 = 'Camera Wrapper instance'
        
        if self.camera is None:
            return msg1 + '; no loaded camera.'
    
        else:
            return msg1 + '; camera loaded:\n\n' + self.camera.__repr__()
    
    
    def load(self):
        '''
        Here we read the first line of fileName, and from that we infer
        the 3D model and the class that should be used. Then, the camera
        is loaded into the wrapper, which later allows to perform the 3D
        epipolar line <-> pixel transformations. 
        
        dirPath - the path of the directory in which the file is stored.
        '''
        fullPath = os.path.join(self.dir, self.fileName)
        f = open(fullPath, 'r')
        firstLine = f.readline()
        f.close()
        
        self.modelName = firstLine.split()[0]
        if self.modelName not in self.ListOfModels:
            raise ValueError('model "%s" not identified'%(self.modelName))
            
        # =====================================================================
        
        if self.modelName == 'Tsai':
            self.camera = camera_Tsai(self.fileName)
            self.camera.load(self.dir)
        
        if self.modelName =='extendedZolof':
            self.camera = camera_extendedZolof(self.fileName)
            self.camera.load(self.dir)
            
    
    
    def save(self):
        
        self.camera.save(dir_path = self.dirPath)
    
        
    
    def getCalibrator(self, lab_coords, img_coords):
        '''
        Based on the camera model, this method initiates an instance of the
        matching calibration class for it and returns is.  
        '''
        # 1) identify which class should be used for the calibration
         
        if self.modelName == 'Tsai':
            calibratorClass = calibrate_Tsai
             
        elif self.modelName == 'extendedZolof':
            calibratorClass = ...
            
        calibrator = calibratorClass(self.camera, lab_coords, img_coords)
        
        return calibrator
    
    
    
    def projection(self, x):
        '''
        Return the image pixel coordinates (eta, zeta) of a lab-space point, x.
        
        input - x (array,3) - 3D lab-space coordinates
        
        output - (array,2) - camera coordinates of the projection of x 
                             (eta, zeta) 
        '''
        
        return self.camera.projection(x)
        
        # if self.modelName == 'Tsai':
        #     return self.camera.projection(x)
        
        # elif self.modelName == 'extendedZolof':
        #     return self.camera.projection(x)
            
    
    
    def get_epipolarline(self, eta, zeta):
        '''
        This method takes in camera pixel coordinates and retuns the 
        direction and origin of the associated epipolar line.
        
        input (two floats) - pixel coordinates (eta, zeta) seen by the camera
        
        output:
        O (array, 3) - the origin of the epipolar line
        r (array, 3) - the direction vector of the epipolar line
        '''
        
        return (self.camera.O, self.camera.get_r(eta, zeta))
        
        # if self.modelName == 'Tsai':
        #     return (self.camera.O, self.camera.get_r(eta, zeta))
        
        # elif self.modelName == 'extendedZolof':
        #     return (self.camera.O, self.camera.get_r(eta, zeta))
    
    
    
    def get_r(self, eta, zeta):
        '''
        input - pixel coordinates (eta, zeta) seen by the camera
        output - direction vector in real space
        '''
        
        return self.camera.get_r(eta, zeta)
    
        # if self.modelName == 'Tsai':
        #     return self.camera.get_r(eta, zeta)
        
        # elif self.modelName == 'extendedZolof':
        #     return self.camera.get_r(eta, zeta)



    def get_r_ori(self, u):
        '''
        A function used for the orientation of fibers, written 
        by Eric Aschari.
        
        input - pixel coordinates (eta, zeta) seen by the camera
        output - direction vector in real space
        '''
        if self.modelName == 'Tsai':
            return self.camera.get_r_ori(u)
        
        elif self.modelName == 'extendedZolof':
            raise ValueError('extendedZolof model does not yet include the '+
                             'orientation feature')
            
    
    @property
    def O(self):
        '''
        Returns the camera origin, O. This method will work only for 3D models
        in which a cantral origin point is used, and otherwise it will raise 
        an error. It is basically made to help in backward compatibility.
        '''
        
        return self.camera.O
        
        # if self.modelName == 'Tsai':
        #     return self.camera.O
        
        # elif self.modelName == 'extendedZolof':
        #     return self.camera.O
        
        
    @property
    def name(self):
        '''
        Returns the camera name
        '''
        if self.camera is not None:
            return self.camera.name
        else:
            print('no loaded camera; returning fileName instead.')
            return self.fileName
        


