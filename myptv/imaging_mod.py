# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 18:02:07 2018

@author: ron


Imaging Module:
    
containts the camera and imaging system classes that handle 
the transformation from camera space coordinates to lab coordinates.


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

import os
from math import sin, cos
from numpy import zeros, array, dot
from numpy.linalg import inv
from myptv.utils import line_dist





class img_system(object):
    '''
    an object that holds a number of cameras.
    '''
    
    def __init__(self, camera_list):
        self.cameras = camera_list
    
    
    def stereo_match(self, coords, d_max):
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
                O1 = self.cameras[ki].O
                r1 = self.cameras[ki].get_r(coords[ki][0], coords[ki][1])
                
                kj = keys[j]
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
    cal_points_fname - path to a file with calibration coordinates for the cam
    '''
    
    def __init__(self, name, resolution, cal_points_fname = None):    
        self.O = zeros(3) + 1.     # camera location
        self.theta = zeros(3) + 1. # rotation angles
        self.f = 1.0               # focal depth / magnification
        self.xh = 0.0              # image center correction, x
        self.yh = 0.0              # image center correction, y
        self.calc_R()
        self.resolution = resolution
        self.give_name(name)
        
        self.E = zeros((3,5))     # correction coefficients matrix
        #self.E = zeros((3,9))     # correction coefficients matrix
        
        if cal_points_fname is not None:
            cic = Cal_image_coord(cal_points_fname)
            self.image_points = cic.image_coords
            self.lab_points = cic.lab_coords
    

    
        
    def __repr__(self):
        
        ret = (self.name +
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
        B = x - self.O
        b = B / (B[0]**2 + B[1]**2 + B[2]**2)**0.5
        v = dot(b, inv(self.R))
        a =  v[2] / self.f
        eta_ = v[0] / a  + self.resolution[0]/2 + self.xh
        zeta_ = v[1] / a + self.resolution[1]/2 + self.yh
        
        # add the error correction term.
        if correction:
            eta, zeta = self.eta_zeta_from_bRinv(eta_, zeta_)
            return array([eta, zeta])
        
        # do not use the error correction term.
        else:
            return array([eta_, zeta_])
    
    
    def eta_zeta_from_bRinv(self, eta_, zeta_):
        '''
        the projection equation is 
        [eta, zeta, f] = b * [R]^-1 + e(eta, zeta)
        This function returns (eta, zeta) for an input of b*[R]^-1
        by solving a least squares equation
        '''
        
        Z3 = [eta_, zeta_, eta_**2, zeta_**2, eta_ * zeta_]
        #Z3 = [eta_, zeta_, eta_**2, zeta_**2, eta_ * zeta_,
        #      eta_**3, eta_**2*zeta_, eta_*zeta_**2, zeta_**3]
        e_ = dot(self.E, Z3)
        
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
        f.write(self.name+'\n')
        
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
        name = f.readline()
        
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
            
    






class Cal_image_coord(object):
    '''
    A class used for reading the calibration image files. This is called
    by the camera class if given a filename with calibration points. 
    '''
    
    def __init__(self, fname):
        '''
        input - 
        fname - String, the path to your calibration point file. The file is
                holds tab separated values with the meaning of: 
                    [x_image, y_image, x_lab, y_lab, z_lab]
        '''
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






    
# =============================================================================
# A live testing of the imaging module:    
# 
# if __name__ == '__main__':
#     from numpy import pi 
#     c1 = camera('1', (1000.,1000.))
#     c2 = camera('2', (1000.,1000.))
#     c3 = camera('3', (1000.,1000.))
#     
#     c1.O = array([400.0 , 0, 1])
#     c2.O = array([0, 400.0, -1])
#     c3.O = array([200.0, 400.0 ,400])
#     
#     c1.f = 4000
#     c2.f = 4000
#     c3.f = 4000
#     
#     c1.theta = array([0.0, -1*pi / 2.0, 0.0])
#     c2.theta = array([pi / 2.0, 0., 0.])
#     c3.theta = array([0.8, -0.4, 0.0])
#     
#     c1.calc_R()
#     c2.calc_R()
#     c3.calc_R()
#     
#     c1.xh = 1.0
#     c1.yh = -1.0
#     
#     x = array([0.1,0.1,0.1])        
#     
#     imgsys = img_system([c1,c2,c3])
#     
#     E = 0.0
#     c1.E[0,:] = E
#     c1.E[1,:] = E
#     c2.E[0,:] = E
#     c2.E[1,:] = E
#     c3.E[0,:] = E
#     c3.E[1,:] = E
#     
#     proj1 = c1.projection(x)
#     proj2 = c2.projection(x)
#     proj3 = c3.projection(x)
#     
#     err1 = array([1.0, 1.0]) * 0
#     err2 = array([1.0, -2.0]) * 0
#     err3 = array([-1.0, 1.0]) * 0
#     
#     coords = {0: proj1 + err1,
#               1: proj2 + err2,
#               2: proj3 + err3}
#     
#     print('\n', proj1)
#     print('',proj2)
#     print('',proj3)
#     
#     res = imgsys.stereo_match(coords, 1e9)
#     print('\n', res)
# 
#    
# 
#     e,z = coords[0]
#     fig, ax = c1.plot_3D_epipolar_line(e, z, (-2,4), color='b')
#    
#     e,z = coords[1]
#     c2.plot_3D_epipolar_line(e, z, (-1,1), ax=ax, color='r')
#         
#     e,z = coords[2]
#     c3.plot_3D_epipolar_line(e, z, (-50,200), ax=ax, color='g')     
# 
#     ax.plot3D(x[0], x[1], x[2], 'ko')
#     x_ = res[0]
#     ax.plot3D(x_[0], x_[1], x_[2], 'yo', fillstyle='none')
#         
#     ax.set_xlim(-25,25)
#     ax.set_ylim(-25,25)
#     ax.set_zlim(-25,25)
# 
# =============================================================================
