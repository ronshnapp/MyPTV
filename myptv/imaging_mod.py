# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 18:02:07 2018

@author: ron


Imaging Module:

This model contains the img_system and camera_wrapper classes used to 
handle the tasks of 2D to 3D transformations.

"""

import os
from math import sin, cos
from numpy import zeros, array, dot, transpose
from numpy.linalg import inv
from myptv.utils import line_dist, point_line_dist

from myptv.TsaiModel.calibrate import calibrate_Tsai
from myptv.TsaiModel.camera import camera_Tsai

from myptv.extendedZolof.camera import camera_extendedZolof
from myptv.extendedZolof.calibrate import calibrate_extendedZolof




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
        N = len(coords)
        x = []
        cams = []
        d = []
        keys = list(coords.keys())
        for i in range(N):
            for j in range(i+1,N):
                ki = keys[i]
                #O1 = self.cameras[ki].O
                #r1 = self.cameras[ki].get_r(coords[ki][0], coords[ki][1])
                eta, zeta = (coords[ki][0], coords[ki][1])
                O1, r1 = self.cameras[ki].get_epipolarline(eta, zeta)
                
                kj = keys[j]
                #O2 = self.cameras[kj].O
                #r2 = self.cameras[kj].get_r(coords[kj][0], coords[kj][1])
                eta, zeta = (coords[kj][0], coords[kj][1])
                O2, r2 = self.cameras[kj].get_epipolarline(eta, zeta)
                
                D, x_ij = line_dist(O1, r1, O2, r2)
                
                if D <= d_max:
                    x.append(x_ij)
                    cams.append(ki)
                    cams.append(kj)
                    d.append(D)
                   
                    
        if len(x)==0:
            return None


        if strict_match==True:
            
            # calculate the mean of all crossing points
            X = sum(x)/len(x)
            
            # Check the mean crossing point is close to all blobs' lines of 
            # sight. If not, return None.  
            for ki in keys:
                ri = self.cameras[ki].get_r(coords[ki][0],coords[ki][1])
                di = point_line_dist(self.cameras[ki].O, ri, X)
                if di>d_max:
                    return None
                
            return X, set(cams), sum(d)/len(x)
            
# =============================================================================
#             if len(x)==len(coords):
#                 X = sum(x)/1.0/len(x)
#                 
#                 for ki in keys:
#                     ri = self.cameras[ki].get_r(coords[ki][0],coords[ki][1])
#                     di = point_line_dist(self.cameras[ki].O, ri, X)
#                     if di>d_max:
#                         return None
#                     
#                 return sum(x)/1.0/len(x), set(cams), sum(d)/1.0/len(x)
# =============================================================================
                
            
        else:
            return sum(x)/len(x), set(cams), sum(d)/len(x)






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
        
        if self.modelName == 'Tsai':
            return self.camera.projection(x)
        
        elif self.modelName == 'extendedZolof':
            return self.camera.projection(x)
            
    
    
    def get_epipolarline(self, eta, zeta):
        '''
        This method takes in camera pixel coordinates and retuns the 
        direction and origin of the associated epipolar line.
        
        input (two floats) - pixel coordinates (eta, zeta) seen by the camera
        
        output:
        O (array, 3) - the origin of the epipolar line
        r (array, 3) - the direction vector of the epipolar line
        '''
        
        if self.modelName == 'Tsai':
            return (self.camera.O, self.camera.get_r(eta, zeta))
        
        elif self.modelName == 'extendedZolof':
            return (self.camera.O, self.camera.get_r(eta, zeta))
    
    
    
    def get_r(self, eta, zeta):
        '''
        input - pixel coordinates (eta, zeta) seen by the camera
        output - direction vector in real space
        '''
        
        if self.modelName == 'Tsai':
            return self.camera.get_r(eta, zeta)
        
        elif self.modelName == 'extendedZolof':
            return self.camera.get_r(eta, zeta)



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
        if self.modelName == 'Tsai':
            return self.camera.O
        
        elif self.modelName == 'extendedZolof':
            return self.camera.O
        


