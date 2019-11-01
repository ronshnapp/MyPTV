#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 12:15:21 2019

@author: ron


calibration module for the cameras

"""

from imaging_mod import * 
from utils import *
from numpy import mean




class calibrate(object):
    '''
    This object is used to calibrate cameras against a given list
    of lab and camera point coordinates
    '''
    
    def __init__(self, camera, lab_coords, img_coords):
        
        self.camera = camera
        self.img_coords = img_coords
        self.lab_coords = lab_coords
        self.D_lst = [self.mean_squared_err()]
        self.sep = sum((array(self.img_coords[1])-array(self.img_coords[0]))**2)**0.5
        
    def mean_squared_err(self):
        '''
        for the current parameters of the camera, this 
        calculaes the mean squared distance between the 
        projection and the given coordinates.
        
        (in the calibration we want to minimize this D)
        '''
        z_lst = [self.camera.projection(x) for x in self.lab_coords]        
        D = mean((array(z_lst) - array(self.img_coords))**2)**0.5
        return D
        
    
    def gradient(self):
        '''
        calculate the gradient vector of the three origins, angles
        and the focal depth. used in the minimization.
        '''
        dO = .001*self.sep
        dth = 0.0001
        df = 0.0005
        
        dDPls = []
        D0 = self.D_lst[-1]
        
        for i in range(len(self.camera.O)):
            self.camera.O[i] = self.camera.O[i] + dO
            dDPls.append(self.mean_squared_err()-D0)
            self.camera.O[i] = self.camera.O[i] - dO
        
        for i in range(len(self.camera.theta)):
            self.camera.theta[i] = self.camera.theta[i] + dth
            self.camera.calc_R()
            dDPls.append(self.mean_squared_err()-D0)
            self.camera.theta[i] = self.camera.theta[i] - dth
            self.camera.calc_R()
            
        self.camera.f += df
        dDPls.append(self.mean_squared_err()-D0)
        self.camera.f -= df
        
        
        dDMns = []
        D0 = self.D_lst[-1]
        
        for i in range(len(self.camera.O)):
            self.camera.O[i] = self.camera.O[i] - dO
            dDMns.append(self.mean_squared_err()-D0)
            self.camera.O[i] = self.camera.O[i] + dO
        
        for i in range(len(self.camera.theta)):
            self.camera.theta[i] = self.camera.theta[i] - dth
            self.camera.calc_R()
            dDMns.append(self.mean_squared_err()-D0)
            self.camera.theta[i] = self.camera.theta[i] + dth
            self.camera.calc_R()
            
        self.camera.f -= df
        dDMns.append(self.mean_squared_err()-D0)
        self.camera.f += df
        
        
        dD = array(dDPls) - array(dDMns)
        dDdO = dD / array([dO,dO,dO,dth,dth,dth,df]) / 2.0
        
        self.dDdO = dDdO
        return dDdO
    
    
    def iterate(self):
        '''
        makes a single gradient decent iteration step
        '''
        dO = 0.2*self.sep
        dth = 0.001
        df = 0.002
        
        dDdO = self.gradient()
        self.camera.O[0] -= dDdO[0] * dO
        self.camera.O[1] -= dDdO[1] * dO
        self.camera.O[2] -= dDdO[2] * dO
        self.camera.theta[0] -= dDdO[3] * dth
        self.camera.theta[1] -= dDdO[4] * dth
        self.camera.theta[2] -= dDdO[5] * dth
        self.camera.f -= dDdO[6] * df
        
        self.camera.calc_R()
        
        self.D_lst.append( self.mean_squared_err() )
    
    
    def plot_proj(self, ax = None):
        
        if ax == None:
            fig, ax = plt.subplots()
        
        ax.plot(array(self.img_coords)[:,0], array(self.img_coords)[:,1], 'ob')
        
        z_lst = [self.camera.projection(x) for x in self.lab_coords]        
        ax.plot( array(z_lst)[:,0], array(z_lst)[:,1], 'xr' )
        ax.set_aspect('equal')
        
        
        
        

if __name__ == '__main__':
    from numpy import pi , array
    import matplotlib.pyplot as plt
    
    # set up a camera:
    c1 = camera('1')
    c1.O = array([20.0 ,0, 0])
    c1.theta = array( [ 0.0, pi / 2.0, 0.0] )
    c1.f = 10.0
    c1.calc_R()
    
    # make up calibration points and projections:
    x_lst = [array([2.0, 0.0, 5.0]),
             array([0.0, 0.0, 0.0]),
             array([0.0, 0.0, -5.0]),
             array([0.0, 5.0, 5.0]),
             array([0.0, 5.0, 0.0]),
             array([2.0, 5.0, -5.0]),
             array([0.0, -5.0, 5.0]),
             array([0.0, -5.0, 0.0]),
             array([0.0, -5.0, -5.0]),
             array([2.0, -5.0, -5.0]),
             array([4.0, -1.0, -1.0]),
             array([-4.0, -1.0, 1.0]),
             array([4.0, 1.0, -1.0]),
             array([-4.0, 1.0, 1.0]),]
    z_lst = [c1.projection(x) for x in x_lst]

    # move the camera to another position:
    c1.O = array([20.0 , 5.0 , 0.0])
    c1.theta = np.array( [ 0.5, pi / 2.0 - .1, 0.2] )
    c1.calc_R()
    c1.f = 10.0
    
    
    cal = calibrate(c1, array(x_lst), array(z_lst))
    #cal.plot_proj()
    
    
    fig, ax = plt.subplots()
    
    for i in range(2000):
        cal.iterate()
        if i%2 == 0 :
            cal.plot_proj(ax)
        if (cal.D_lst[-1] / cal.sep) < 0.005:
            break
        
    print c1.O
    print c1.theta
    print c1.f
    
    fig, ax = plt.subplots()
    ax.semilogy( array(cal.D_lst)/cal.sep ) 

        
        
        
        
        
        
        
        