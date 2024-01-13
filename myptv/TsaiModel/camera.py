# -*- coding: utf-8 -*-
#!/usr/bin/env python3

"""
Created on Thu Dec 28 07:41:56 2023

@author: ron


containts the camera class that represents the Tsai 3D camera model.
The class handles the transformation from camera space coordinates to 
lab coordinates using the Tsai, pin-hole model, with a non-linear 
polynomial correction term.


For the 3D-model we use a pin-hole camera model. Each camera has an
imaging center (O) and an orientation vector (theta), that describe 
its location and rotation in 3D space. Each point in image space, whose 
coordinates are (eta, zeta), is related to a ray going from the imaging 
center, and has a direction vector b. The vector b can be calculated as:
    
    (1)    b = [eta, zeta, f] * [R]

where [R] is the rotation matrix that corresponds to theta, and f is the
focal length of the camera optical system.


Equation (1) is a linear model. It might fit in ideal cases with no
image distortion or multimedia problems, however in realistic cases
it might not be sufficiently accurate. To overcome this difficulty,
we add a correction term, redefining (1) as:
    
    (2)    b = ([eta, zeta, f] + e) * [R]
    
The correction term is given as a cubic polynomial in eta and zeta as follows:

    (3)    e = [E] * Z3
    
           Z3 = [eta, zeta, eta^2, zeta^2, eta*zeta, eta^3, eta^2*zeta, zeta^2*eta, zeta^3]  

so Z3 are the polymer terms and [E] is a (3X9) matrix with a total 
of 27 coefficients.  
"""


from myptv.utils import Cal_image_coord

import os
from math import sin, cos
from numpy import zeros, array, dot, transpose





class camera_Tsai(object):
    '''
    an object that holds the calibration information for each camera using the
    Tsai model with polynomial non-linear term.
    It can be used to:
    1) obtain image coordinates from given lab coordinates. 
    2) vice versa if used together with other cameras at 
       other locations (img_system.stereo_match).
      
    input:
    name - string name for the camera
    resolution - tuple (2) two integers for the camera pixels
    cal_points_fname - path to a file with calibration coordinates for the cam
    '''
    
    def __init__(self, name, cal_points_fname = None):    
        self.O = zeros(3) + 1.     # camera location
        self.theta = zeros(3) + 1. # rotation angles
        self.f = 1.0               # focal depth / magnification
        self.xh = 0.0              # image center correction, x
        self.yh = 0.0              # image center correction, y
        self.calc_R()
        self.resolution = (1, 1)
        self.give_name(name)
        
        self.E = zeros((3,5))     # correction coefficients matrix
        #self.E = zeros((3,9))     # correction coefficients matrix
        
        if cal_points_fname is not None:
            cic = Cal_image_coord(cal_points_fname)
            self.image_points = cic.image_coords
            self.lab_points = cic.lab_coords
    

    
        
    def __repr__(self):
        
        ret = ('Tsai model camera instace. ' + 
               self.name +
               '\n O: ' + str(self.O) +
               '\n theta:' + str(self.theta) +
               '\n f:' + str(self.f))
        return ret
    
    
    
    def give_name(self, name):
        '''
        adds a name for the camera
        '''
        if type(name) == str:
            self.name = name
        else:
            raise TypeError('name must be string')    
    
    
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
        
        #bR_inv = dot(-self.O, self.R.T)
        #self.a = bR_inv[2] / self.f
    
    
    def get_r(self, eta, zeta):
        '''
        r = ([eta, zeta, f] + e) * [R]
        
        input - pixel coordinates (eta, zeta) seen by the camera
        output - direction vector in real space
        '''
        # self.calc_R()
        eta_ = eta - self.resolution[0]/2.0 - self.xh
        zeta_ = zeta - self.resolution[1]/2.0  - self.yh
        
        
        Z3 = [eta, zeta, eta**2, zeta**2, eta * zeta]
        #Z3 = [eta, zeta, eta**2, zeta**2, eta * zeta, 
        #      eta**3, eta**2*zeta, eta*zeta**2, zeta**3]
        
        e = dot(self.E, Z3)
        
        r = dot(array([-eta_, -zeta_, -self.f]) - e, self.R)
        r = r / (r[0]**2 + r[1]**2 + r[2]**2)**0.5
        return r
    
    
    def get_r_ori(self, u): # from Eric
        '''
        A function used for the orientation of fibers, written 
        by Eric Aschari.
        
        input - pixel coordinates (eta, zeta) seen by the camera
        output - direction vector in real space
        '''
         
        eta = u[0,0]
        zeta = u[1,0]
         
        Z3 = [eta, zeta, eta**2, zeta**2, eta * zeta]
         
        e = dot(self.E, Z3)
        
        r = dot(array([-eta, -zeta, -self.f]) - e, self.R) # previously with - before eta_ and zeta_
        r = r / (r[0]**2 + r[1]**2 + r[2]**2)**0.5
        
        return transpose(array([self.O + r])) 
    
    
    def projection(self, x, correction=True):
        '''
        will return the image coordinate (eta, zeta) of a real point x.
        
        input - x (array,3) - real world coordinates
                correction - if True, will return the coordinates after
                the non-linear error correction. If False, we not do the
                correction.
        output - (eta, zeta) (array,2) - camera coordinates of the projection 
                                         of x
        '''
        b = x - self.O
        v = dot(b, self.R.T) #inv(self.R))
        a =  v[2] / self.f
        eta_ = v[0] / a  + self.resolution[0]/2 + self.xh
        zeta_ = v[1] / a + self.resolution[1]/2 + self.yh
        
        # if we add the error correction term.
        if correction:
            eta, zeta = self.eta_zeta_from_bRinv(eta_, zeta_)
            return array([eta, zeta])
        
        # if we don't use the error correction term.
        else:
            return array([eta_, zeta_])
    
    
    def eta_zeta_from_bRinv(self, eta_, zeta_):
        '''
        the projection equation is 
        [eta, zeta, f] = b * [R]^-1 + e(eta, zeta)
        This function returns (eta, zeta) for an input of b*[R]^-1. 
        To make inverting the non-linear correction solvable we 
        linearize the error term with a Taylor series expantion.
        '''
        
        Z3 = [eta_, zeta_, eta_**2, zeta_**2, eta_ * zeta_]
        #Z3 = [eta_, zeta_, eta_**2, zeta_**2, eta_ * zeta_,
        #      eta_**3, eta_**2*zeta_, eta_*zeta_**2, zeta_**3]
        e_ = dot(self.E, Z3)
        
        # calculating the derivatives of the error term:
        e_0 = e_[0]
        a, b, c, d, ee = self.E[0,:]
        e_eta_0 = a + 2*c*eta_ + ee*zeta_
        e_zeta_0 = b + 2*d*zeta_ + ee*eta_
        #a, b, c, d, ee, f, g, h, i = self.E[0,:]
        #e_eta_0 = a + 2*c*eta_ + ee*zeta_ + 3*f*eta_**2 + 2*g*eta_*zeta_ + h*zeta_**2
        #e_zeta_0 = b + 2*d*zeta_ + ee*eta_ + g*eta_**2 + 2*h*eta_*zeta_ + 3*i*zeta_**2
        
        e_1 = e_[1]
        a, b, c, d, ee = self.E[1,:]
        e_eta_1 = a + 2*c*eta_ + ee*zeta_
        e_zeta_1 = b + 2*d*zeta_ + ee*eta_
        #a, b, c, d, ee, f, g, h, i = self.E[1,:]
        #e_eta_1 = a + 2*c*eta_ + ee*zeta_ + 3*f*eta_**2 + 2*g*eta_*zeta_ + h*zeta_**2
        #e_zeta_1 = b + 2*d*zeta_ + ee*eta_ + g*eta_**2 + 2*h*eta_*zeta_ + 3*i*zeta_**2
        
        A11 = 1.0 + e_eta_0
        A12 = e_zeta_0
        A21 = e_eta_1
        A22 = 1.0 + e_zeta_1
        
        rhs1 = eta_*(1.0 + e_eta_0) + zeta_*e_zeta_0 - e_0
        rhs2 = zeta_*(1.0 + e_zeta_1) + eta_*e_eta_1 - e_1
        
        Ainv = array([[A22, -A12],[-A21, A11]]) / (A11*A22 - A12*A21)
        eta, zeta = dot(Ainv, [rhs1, rhs2])
        
        return eta, zeta
    
    
    
    def save(self, dir_path = ''):
        '''
        will save the camera on the hard drive
        '''
        full_path = os.path.join(dir_path, self.name)
        
        f = open(full_path, 'w')
        f.write('Tsai model camera\n')
        f.write(self.name+'\n')
        
        S = ''
        for s in self.resolution:
            S+= str(s)+' '
        f.write(S+'\n')
        
        S = ''
        for s in self.O:
            S+= str(s)+' '
        f.write(S+'\n')
        
        S = ''
        for s in self.theta:
            S+= str(s)+' '
        f.write(S+'\n')
        
        f.write(str(self.f)+'\n')
        
        S = ''
        for s in [self.xh, self.yh]:
            S+= str(s)+' '
        f.write(S+'\n')
        
        for i in range(3):
            S = ''
            for s in self.E[i,:]:
                S+= str(s)+' '
            f.write(S+'\n')
        
        f.close()
        
        
        
    def load(self, dir_path):
        '''
        will load camera data from the hard disk
        '''
        full_path = os.path.join(dir_path, self.name)
        
        f = open(full_path, 'r')
        modelName = f.readline()
        name = f.readline()
        
        S = f.readline().strip()
        self.resolution = array([int(s) for s in S.split()])
        
        S = f.readline()[:-2]
        self.O = array([float(s) for s in S.split()])
        
        S = f.readline()[:-2]
        self.theta = array([float(s) for s in S.split()])
        
        self.f = float(f.readline()[:-2])
        
        S = f.readline()[:-2]
        self.xh, self.yh = array([float(s) for s in S.split()])
        
        for i in range(3):
            S = f.readline()[:-2]
            self.E[i,:] = array([float(s) for s in S.split()])
        
        f.close()
        
        self.calc_R()
        
    
    def plot_3D_epipolar_line(self, eta, zeta, zlims, ax=None, color=None):
        '''Will plot a 3D epipolar line for the image point (eta,zeta) in 
        between two given z values.
        
        This requires matplotlib.'''
        
        z0, z1 = zlims
        r = self.get_r(eta, zeta)
        a0 = (z0 - self.O[2])/r[2]
        a1 = (z1 - self.O[2])/r[2]
        x0, x1 = self.O[0]+a0*r[0], self.O[0]+a1*r[0]
        y0, y1 = self.O[1]+a0*r[1], self.O[1]+a1*r[1]
        
        if ax is None:
            from mpl_toolkits import mplot3d
            import matplotlib.pyplot as plt
            fig = plt.figure()
            ax = plt.axes(projection='3d')
            if color is None:
                ax.plot3D([x0,x1], [z0,z1], [y0,y1])
                ax.plot3D([self.O[0]], [self.O[2]], [self.O[1]], 'ko')
            else:
                ax.plot3D([x0,x1], [z0,z1], [y0,y1], c=color)
                ax.plot3D([self.O[0]], [self.O[2]], [self.O[1]], 'o', c=color)
        
            return fig, ax
       
        else:
            if color is None:
                ax.plot3D([x0,x1], [z0,z1], [y0,y1])
                ax.plot3D([self.O[0]], [self.O[2]], [self.O[1]], 'ko')
            else:
                ax.plot3D([x0,x1], [z0,z1], [y0,y1], c=color)
                ax.plot3D([self.O[0]], [self.O[2]], [self.O[1]], 'o', c=color)
            
    
    
    







