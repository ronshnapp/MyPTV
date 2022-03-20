#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 12:15:21 2019

@author: ron


A class for reading calibration image coordinate files.

The image coordinate files have tab separated columns with coordiantes as:

x_image, y_image, x_lab, y_lab, z_lab

"""



from numpy import array




class Cal_image_coord(object):
    
    
    def __init__(self, fname):
        
        self.image_coords = []
        self.lab_coords = []
        self.fname = fname
        self.read_file()
        
        
    def read_file(self):
        
        with open(self.fname) as f:
            lines = f.readlines()
            self.N_points = len(lines)
            
            for ln in lines:
                ln_ = ln.split()
                self.image_coords.append([float(ln_[0]), float(ln_[1])])
                self.lab_coords.append( [float(ln_[2]), 
                                         float(ln_[3]),
                                         float(ln_[4])] )
                f.close()
        
        
    
        

        
        
        
        
        