#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 12:15:21 2019

@author: ron


calibration module for the cameras

"""

from imaging_mod import * 
from utils import *
from numpy import mean, sum, hstack, array
import matplotlib.pyplot as plt



class calibrate(object):
    '''
    This object is used to calibrate cameras against a given list
    of lab and camera point coordinates. 
    The main problem is to minimize the distance between the projected
    lab coordinates and the given image coordinates. The minimalization 
    procedure uses scipy's minimize method (through searchCalibration()).
    '''
    
    def __init__(self, camera, lab_coords, img_coords):
        
        self.camera = camera
        self.img_coords = img_coords
        self.lab_coords = lab_coords
        self.D_lst = [self.mean_squared_err()]
        self.sep = sum((array(self.img_coords[1])-array(self.img_coords[0]))**2)**0.5
        
    def mean_squared_err(self):
        '''
        This calculaes the mean squared distance between the 
        projection and the given coordinates  (in units of pixel).
        
        (in the calibration we want to minimize this D)
        '''
        z_lst = [self.camera.projection(x) for x in self.lab_coords]        
        e = array(z_lst) - array(self.img_coords)
        D = mean( sum(e**2, axis=1)**0.5 )
        return D
        
    
    
    def searchCalibration(self, maxiter=5000, fix_f=True):
        '''
        using scipy's minimize function to obtain calibration
        parameters for the camera.
        '''
        from scipy.optimize import minimize
        
        def func(X):
            self.camera.O = X[:3]
            self.camera.theta = X[3:6]
            if not fix_f:
                self.camera.f = X[-1]
                
            self.camera.calc_R()
            
            meanSquaredErr = self.mean_squared_err()
            self.D_lst.append( meanSquaredErr )
            return meanSquaredErr
        
        c = self.camera
        
        if fix_f:
            X0 = hstack([c.O, c.theta])
        
        else:
            X0 = hstack([c.O, c.theta, c.f])
            
            
        res = minimize(func, X0, method='nelder-mead', 
                       options={'disp': True, 'maxiter': maxiter})
        return res
    
    
    
    def plot_proj(self, ax = None):
        
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
        
        
    
        

        
        
        
        
        