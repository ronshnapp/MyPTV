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
        
        
        
        
        
        
        