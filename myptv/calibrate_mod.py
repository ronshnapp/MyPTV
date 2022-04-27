# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 12:15:21 2019

@author: ron


calibration module for the cameras

"""

from numpy import mean, sum, hstack, array




class calibrate(object):
    '''
    This object is used to calibrate cameras against a given list
    of lab and camera point coordinates. 
    The main problem is to minimize the distance between the projected
    lab coordinates and the given image coordinates. The minimization 
    procedure uses scipy minimize method (through searchCalibration()).
    '''
    
    def __init__(self, camera, lab_coords, img_coords):
        
        self.camera = camera
        self.img_coords = img_coords
        self.lab_coords = lab_coords
        self.D_lst = [self.mean_squared_err()]
        self.sep = sum((array(self.img_coords[1])-array(self.img_coords[0]))**2)**0.5
        
        
    def mean_squared_err(self, correction=True):
        '''
        This calculates the mean squared distance between the 
        projection and the given coordinates  (in units of pixel).
        
        (in the calibration we want to minimize this D)
        '''
        z_lst = []        
        for x in self.lab_coords:
            z_lst.append(self.camera.projection(x, correction=correction))
        e = array(z_lst) - array(self.img_coords)
        D = mean( sum(e**2, axis=1)**0.5 )
        return D
        
    
    def searchCalibration(self, maxiter=5000, fix_f=True):
        '''
        using scipy minimize function to obtain calibration
        parameters for the camera.
        '''
        from scipy.optimize import minimize
        
        def func(X):
            self.camera.O = X[:3]
            self.camera.theta = X[3:6]
            self.camera.xh = X[6]
            self.camera.yh = X[7]
            
            if not fix_f:
                self.camera.f = X[-1]
                
            self.camera.calc_R()
            
            meanSquaredErr = self.mean_squared_err()
            self.D_lst.append( meanSquaredErr )
            return meanSquaredErr
        
        c = self.camera
        
        if fix_f:
            X0 = hstack([c.O, c.theta, c.xh, c.yh])
        
        else:
            X0 = hstack([c.O, c.theta, c.xh, c.yh, c.f])
                        
        res = minimize(func, X0, method='nelder-mead', 
                       options={'disp': True, 'maxiter': maxiter})
        return res
    
    
    
    def fineCalibration(self, maxiter=500):
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
            meanSquaredErr = self.mean_squared_err(correction=True)
            self.D_lst.append( meanSquaredErr )
            return meanSquaredErr
        
        c = self.camera
        X0 = c.E[:2,:]
        res = minimize(func, X0, method='nelder-mead', 
                       options={'disp': True, 'maxiter': maxiter})
        return res
    
    
    
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
        
    
        

        
        
        
        
        