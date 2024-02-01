# -*- coding: utf-8 -*-
#!/usr/bin/env python3

"""
Created on Jan 13 2024

@author: ron


containts the camera class that represents the Extended Zolof 3D camera model.
The class handles the transformation from camera space coordinates to 
lab coordinates using the extended Zolof model, with a third-order back-and-
forth polynomial transformation.



The original forward transformation model:
    
The Extended Zolof 3D model is an extension of the existing Zolof polynomial
method of 3D transformation. In the original Zolof model, each 3D point in 
lab-space, X, is related to a 2D camera-space point, x, through a third degree
poynomial. Specifically, let X = (X1, X2, X3) and x = (eta, zeta). The 
transformation from X to x is: 
    
    x = [A] * P(X);

In MyPTV, P(X) is the following third order polynomial with 17 terms:
    
    P(X) = 1 + X1 + X2 + X3 + X1^2 + X2^2 + X3^2 + X1*X2 + X2*X3 + X3*X1 + 
           X1*X2^2 + X1*X3^2 + X2*X1^2 + X2*X3^2 + X3*X1^2 + X3*X2^2 + X1*X2*X3
    
and [A] is a vector of 17 coefficients:
    
    [A] = A1, A2, ..., A17 .
    
Note that unlike other models, the Zolof model does not assume anything about
the physics of the camera, and simply applied a least-squares fitting to find
the coefficients [A] that reduce the transformaiton error. The polynomial used
in MyPTV was tested to give good results for standard PTV experiments.



The extension to backward transformation:
    
The original Zolof model can only perform a lab-to-camera transformation. 
However, to complete a PTV analysis, we need to be able to transform camera
space coordinates back into lab space coordinates, and this does not exist in 
the original Zolof model. Therefore, in MyPTV we extended this model to 
include also the inverse, camera-to-lab transformation. 

To do so we use an approach that combines the pin-hole model with the Zolof 
polynomial. Specifically, we assume that each point in camera-space (e.g. 
each camera pixel) corresponds to a stright physical line in lab-space; these 
lines are called epipolar lines. Furthermore, we assume that all the epipolar 
lines for a certain camera cross each other in a single point. These two 
assumptions were tested with two separate experiments, which gave good results.
With the two assumptions given above, it is possible to define the backwards, 
namely, camera-to-lab transformation as a function that given coordinates
in camera space, returns an epipolar line.   

With the second assumption above, epipolar lines belonging to a camera can 
be defined using a direction vector. Thus, the backwards, camera-to-lab 
transformation is made as follows. Let O = (O1, O2, O3) denote the point in 
which all epipolar lines cross, and let r = (r1, r2, r3) denote a direction
vector in lab-space. Points laying along an epipolar line can thus be written 
as:
    
    l(x) = O + r(x) * a

where a is a free parameter running from -infinity to +infinity. The 
camera-to-lab transformation is achieved by:
    
    r(x) = [B] G(x)
    
where G(x) is a third degree polynomial with the following 10 terms

    G(x) = 1 + x1 + x2 + x1^2 + x2^2 + x1*x2 + x1^2*x2 + x2^2*x1 + x1^3 + x2^3
    
and where [B] is a vector of 10 polynomial coeficients:
    
    [B] = B1, B2, ... , B10 .
    
    
    
    
The advantage of the Extended Zolof model:
    
The advantage of the Extended Zolof model is that 1) it is fully defined and
and perfectly fits to be used in PTV experiments because it has the forward 
and backwards transformations; 2) it is straightforward to determine the 
model parameters in the calibration process. Specifically, the problem can be
formulated as a normal least-squares minimization, which gives out fast and 
unique results without resorting to numerial minimization. This is true for
all the model parameters, [A], [B], and O.  

"""

from myptv.utils import Cal_image_coord

from numpy import dot, array, zeros
import os






class camera_extendedZolof(object):
    '''
    an object that holds the calibration information for each camera using the
    Extended Zolof model.
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
        '''
        inputs:
        name - (string) name of the camera and name of the file in which 
                        it is stored on the disk.
        
        cal_points_fname (string) - path of a file that stores the calibration
                                    points data. Optional. 
        '''
        
        self.name = name
        self.O = zeros(3) + 1.     # camera location
        self.A = array([[0.0 for i in range(17)] for j in [0,1]]).T
        self.B = array([[0.0 for i in range(10)] for j in [0, 1, 2]]).T
        
        if cal_points_fname is not None:
            cic = Cal_image_coord(cal_points_fname)
            self.image_points = cic.image_coords
            self.lab_points = cic.lab_coords
        
    

    def __repr__(self):
        
        ret = ('Extended Zolof model camera instace. ' + 
               self.name + '\n O: ' + str(self.O))
        return ret
    
    
    
    def get_XCol(self, X):
        '''
        Given a point in 3D lab-space, this method returns its 17 P(X) 
        polynomial terms.
        '''
        X1,X2,X3 = X[0],X[1],X[2]
        XColumn = [1.0, X1, X2, X3,
                   X1**2, X2**2, X3**2, X1*X2, X2*X3, X3*X1, X1*X2*X3,
                   X1*X2**2, X1*X3**2, X2*X1**2, X2*X3**2, X3*X1**2, X3*X2**2]
        return XColumn
    


    def get_xCol(self, x):
        '''
        Given a point in 2D camera-space, this method returns its 10 G(x) 
        polynomial terms
        '''
        x1, x2 = x[0], x[1]
        xColumn = [1.0, x1, x2, x1**2, x2**2, x1*x2,
                   x1**2*x2, x2**2*x1, x1**3, x2**3]
        return xColumn
    
    
    
    def projection(self, X):
        '''
        Given a point in 3D, X, this method returns its 2D camera-space 
        projection.
        '''
        XColumn = self.get_XCol(X)
        res = dot(XColumn, self.A)
        return [res[0], res[1]]
    
    
    
    def get_r(self, eta, zeta):
        '''
        Given a point in 2D, x, this method returns its 3D direction vector
        '''
        x = [eta, zeta]
        xColumn = self.get_xCol(x)
        res = dot(xColumn, self.B)
        return array([res[0], res[1], res[2]])
    
    
    
    def get_r_ori():
        msg1 = 'The extended Zolof model is not yet capable of returning'
        msg2 = 'particle orientation results. Use the Tsai model instead.'
        raise TypeError(msg1 + msg2)
    
    
    
    def save(self, dir_path = ''):
        full_path = os.path.join(dir_path, self.name)
        
        f = open(full_path, 'w')
        f.write('extendedZolof model camera\n')
        f.write(self.name+'\n')
        
        S = ''
        for i in range(len(self.O)):
            S += 'O %s \t %s \n'%(i, self.O[i])
        
        for i in range(len(self.A)):
            for j in range(len(self.A[i])):
                S += 'A %s %s \t %s \n'%(i, j, self.A[i][j])
        
        for i in range(len(self.B)):
            for j in range(len(self.B[i])):
                S += 'B %s %s \t %s \n'%(i, j, self.B[i][j])
        
        f.write(S)
        f.close()
        
    

    def load(self, dir_path):
        '''
        will load camera data from the hard disk
        '''
        full_path = os.path.join(dir_path, self.name)
        
        f = open(full_path, 'r')
        model = f.readline()
        modelName = model.split()[0]
        if modelName != 'extendedZolof':
            msg = 'Camera file is not an extendedZolof camera (%s)'%modelName
            raise ValueError(msg)
        name = f.readline()
        
        lines = f.readlines()
        f.close()
        
        for i in range(len(lines)):
            lines[i] = lines[i].strip().split()
        
        O_vals = list(filter(lambda l: l[0]=='O', lines))
        for v in O_vals:
            ind = int(v[1])
            self.O[ind] = float(v[-1])
        
        A_vals = list(filter(lambda l: l[0]=='A', lines))
        Ashape = (int(max(A_vals, key=lambda x: float(x[1]))[1])+1,
                  int(max(A_vals, key=lambda x: float(x[2]))[2])+1)
        self.A = array([[0.0 for i in range(Ashape[0])] for j in range(Ashape[1])]).T
        for v in A_vals:
            i, j = int(v[1]), int(v[2])
            self.A[i][j] = float(v[-1])
        
        B_vals = list(filter(lambda l: l[0]=='B', lines))
        Bshape = (int(max(B_vals, key=lambda x: float(x[1]))[1])+1,
                  int(max(B_vals, key=lambda x: float(x[2]))[2])+1)
        self.B = array([[0.0 for i in range(Bshape[0])] for j in range(Bshape[1])]).T
        for v in B_vals:
            i, j = int(v[1]), int(v[2])
            self.B[i][j] = float(v[-1])
        
    
            
    






