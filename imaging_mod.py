#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 18:02:07 2018

@author: ron


Imaging Module:
    
containts the camera and imaging system classes that handle 
the transformation for camera space coordinates to lab coordinates.

"""

import sys, os
from numpy import zeros, array, cos, sin, dot
from numpy.linalg import inv
from utils import *




class img_system(object):
    '''
    an object that holds a number of cameras.
    '''
    
    def __init__(self):
        self.cameras = []
    
    def add_cameras(self,n):
        for i in range(n):
            self.cameras.append( camera() )
    
    def stereo_match(self, coords, d_max):
        '''
        given n particle images [(eta, zeta) coords in camera space], this will
        determine whereather there is a good candidate point for the intersection
        of these epipolar lines, and if so returns it. Here it is assumed that
        all images correspond to the same "real world" particle.
        This point is estimated as the average of the crossing point of the
        epipolar lines, that cross at distances smaller than a maximum value.
        
        input - 
        coords (dic) - keys are camera number, values are the image space 
                       coordinates of each point. Must have at least 2 entries
        d_max (float) - maximum allowable distance separating two lines
        
        
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
        for i in range(N):
            for j in range(i+1,N):
                ki = coords.keys()[i]
                O1 = self.cameras[ki].O
                r1 = self.cameras[ki].get_r(coords[ki][0], coords[ki][1])
                
                kj = coords.keys()[j]
                O2 = self.cameras[kj].O
                r2 = self.cameras[kj].get_r(coords[kj][0], coords[kj][1])
                
                D, x_ij = line_dist(O1, r1, O2, r2)
                
                if D <= d_max:
                    x.append(x_ij)
                    cams.append(ki)
                    cams.append(kj)
                    d.append(D)

        if len(x)>=1:
            return sum(x)/1.0/len(x), set(cams), sum(d)/1.0/len(x)
        else:
            return None





class camera(object):
    '''
    an object that holds the calibration information for
    each camera. It can be used to:
    1) obtain image coordinates from given lab coordinates. 
    2) vice versa if used together with other cameras at 
       other locations (img_system.stereo_match).
      
    input:
    name - string name for the camera
    resolution - tuple (2) two integers for the camera pixels
    '''
    
    def __init__(self, name, resolution):    
        self.O = zeros(3)     # camera location
        self.theta = zeros(3) # rotation angles
        self.f = 1.0          # focal depth
        self.calc_R()
        self.resolution = resolution
        self.give_name(name)
        
        
    def calc_R(self):
        '''
        calculates the rotation matrix for the camera's angles
        '''
        tx,ty,tz = self.theta
        Rx = array([[1,0,0],
                    [0,cos(tx),-sin(tx)],
                    [0,sin(tx),cos(tx)]])
        Ry = array([[cos(ty),0,sin(ty)],
                     [0,1,0],
                     [-sin(ty),0,cos(ty)]])
        Rz = array([[cos(tz),-sin(tz),0],
                    [sin(tz),cos(tz),0],
                    [0,0,1]])
        self.R = dot(dot(Rx,Ry), Rz)
    
    
    def get_r(self, eta, zeta):
        '''
        input - pixel coordinates (eta, zeta) seen by the camera
        output - direction vector in real space
        '''
        self.calc_R()
        r = dot(array([-eta, -zeta, -self.f]), self.R)
        return r
    
    
    def projection(self,x):
        '''
        will return the image coordinate (eta, zeta) of a real point x.
        
        input - x (array,3) - real world coordinates
        output - (eta, zeta) (array,2) - camera coordinates of the projection 
                                         of x
        '''
        v = dot(x - self.O, inv(self.R))
        a = -1.0 * v[2] / self.f
        eta = (-1.0 * v[0]) / a  + self.resolution[0]/2
        zeta = (-1.0 * v[1]) / a  + self.resolution[1]/2
        return array([eta, zeta])
        
    
    def give_name(self, name):
        '''
        adds a name for the camera
        '''
        if type(name) == str:
            self.name = name
        else:
            raise TypeError('name must be string')
    
    
    def save(self, dir_path = ''):
        '''
        will save the camera on the hard drive
        '''
        full_path = os.path.join(dir_path, self.name)
        
        f = file(full_path, 'w')
        f.write(self.name+'\n')
        
        S = ''
        for s in self.O:
            S+= str(s)+' '
        f.write(S+'\n')
        
        S = ''
        for s in self.theta:
            S+= str(s)+' '
        f.write(S+'\n')
        
        f.write(str(self.f))
        f.close()
        
    def load(self, dir_path):
        '''
        will load camera data from the jard disk
        '''
        full_path = os.path.join(dir_path, self.name)
        
        f = file(full_path)
        self.name = f.readline()[:-2]
        
        S = f.readline()[:-2]
        self.O = array([float(s) for s in S.split()])
        
        S = f.readline()[:-2]
        self.theta = array([float(s) for s in S.split()])
        
        self.f = float(f.readline()[:-2])
        f.close()
        
        self.calc_R()
        
    
    
    
    
    
    
if __name__ == '__main__':
    from numpy import pi 
    c1 = camera('1', (10,10))
    c2 = camera('2', (10,10))
    c3 = camera('3', (10,10))
    
    c1.O = array([10.0 ,0,0])
    c2.O = array([0,10.0 ,0])
    c3.O = array([10.0,10.0 ,0])
    
    c2.theta[0] = -pi / 2.0
    c1.theta[1] = pi / 2.0
    
    c1.calc_R()
    c2.calc_R()
    c3.calc_R()
    
    x = array([0.0,0.1,0.1])        
    
    imgsys = img_system()
    imgsys.cameras.append(c1)
    imgsys.cameras.append(c2)
    imgsys.cameras.append(c3)
    
    coords = {0: c1.projection(x),
              1: c2.projection(x)*1.01,
              2: c3.projection(x)*0.99}
        
    print imgsys.stereo_match(coords, 0.5)
        
        
        
        
        
        
        
        
        
