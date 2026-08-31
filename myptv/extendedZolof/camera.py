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
           X1*X2^2 + X1*X3^2 + X2*X1^2 + X2*X3^2 + X3*X1^2 + X3*X2^2 + 
           X1*X2*X3 + X1^3 + X2^3
    
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
lines are called epipolar lines. 

In the ORIGINAL model, we further assumed that all the epipolar lines for a
certain camera cross each other in a single point, O = (O1, O2, O3), i.e. 
that the camera behaves like an ideal pinhole. Points laying along an
epipolar line are then written as:
    
    l(x) = O + r(x) * a

where a is a free parameter running from -infinity to +infinity, and

    r(x) = [B] G(x)

is the ray direction, with G(x) the following third degree, 10-term
polynomial in camera-space:

    G(x) = 1 + x1 + x2 + x1^2 + x2^2 + x1*x2 + x1^2*x2 + x2^2*x1 + x1^3 + x2^3

and [B] a vector of 10 polynomial coefficients.


The extension to a variable (per-pixel) origin:

Real cameras/lenses are never perfect pinholes, so forcing every epipolar
line through a single point O is itself an approximation, and it can be a
significant source of error for wide-angle lenses or cameras close to the
calibration volume. This class optionally supports relaxing that assumption:
instead of a single fixed O, the origin of the epipolar line is itself
allowed to be a function of the pixel, using the same 10-term basis G(x):

    O(x) = [C] G(x)
    l(x) = O(x) + r(x) * a ,      r(x) = [B] G(x)

[C] is a vector of 10 polynomial coefficients (one set per lab-space
dimension), fit alongside [B]. This is controlled by the `variable_origin`
flag: when True, `get_ray()` returns the per-pixel origin O(x); when False
(the default, and the behavior of all previously-calibrated camera files),
it returns the single fixed self.O, exactly as before.



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
from numpy.linalg import norm
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
        self.O = zeros(3) + 1.     # camera location (fixed-origin model)
        self.A = array([[0.0 for i in range(19)] for j in [0,1]]).T
        self.B = array([[0.0 for i in range(10)] for j in [0, 1, 2]]).T
        self.C = array([[0.0 for i in range(10)] for j in [0, 1, 2]]).T

        # Whether get_ray() should use the per-pixel origin O(x) = [C]G(x)
        # (True) 
        self.variable_origin = False

        # Pixel-coordinate normalization used ONLY when variable_origin is
        # True: get_xCol_scaled() centers/scales (eta,zeta) before building
        # the polynomial basis, so that fitting [C] is numerically 
        # well-conditioned.
        # Set by fit_variable_origin_model(); identity (no-op) by default.
        self.px_center = array([0.0, 0.0])
        self.px_scale = array([1.0, 1.0])

        # How the per-pixel origin is parameterized when variable_origin
        # is True:
        #   'free'  - O(x) = [C] G(x), with [C] of shape (n_c, 3), i.e.
        #             the origin floats freely in 3D. Note this has a gauge
        #             degeneracy: sliding an origin ALONG its own ray does
        #             not change the epipolar line at all, so one of the
        #             three components is unidentifiable from the fit and
        #             is driven purely by noise.
        #   'plane' - O(x) is constrained to a fixed plane that passes
        #             through self.O with normal self.plane_n:
        #                 O(x) = self.O + plane_u*a(x) + plane_v*b(x)
        #             where (a, b) = [C] G(x) with [C] of shape (n_c, 2),
        #             and plane_u/plane_v are orthonormal in-plane axes.
        #             This removes the degenerate direction above and uses
        #             2 (rather than 3) coefficient sets per basis term.
        # 'free' is the default so that previously-saved variable-origin
        # cameras keep their original behavior.
        self.origin_mode = 'free'
        self.plane_u = array([1.0, 0.0, 0.0])
        self.plane_v = array([0.0, 1.0, 0.0])
        self.plane_n = array([0.0, 0.0, 1.0])
        
        if cal_points_fname is not None:
            cic = Cal_image_coord(cal_points_fname)
            self.image_points = cic.image_coords
            self.lab_points = cic.lab_coords
        
    

    def __repr__(self):
        
        ret = ('Extended Zolof model camera instace. ' + 
               self.name + '\n O: ' + str(self.O) +
               '\n variable_origin: ' + str(self.variable_origin))
        return ret
    
    
    
    def get_XCol(self, X):
        '''
        Given a point in 3D lab-space, this method returns its 17 P(X) 
        polynomial terms.
        '''
        X1,X2,X3 = X[0],X[1],X[2]
        mx = max(array(self.A).shape)
        if mx==17: # compatibility with v1.3.5 and lower
            XColumn = [1.0, X1, X2, X3,
                       X1**2, X2**2, X3**2, X1*X2, X2*X3, X3*X1, X1*X2*X3,
                       X1*X2**2, X1*X3**2, X2*X1**2, X2*X3**2, X3*X1**2, X3*X2**2]
        elif mx==19: # for v1.3.6 and above
            XColumn = [1.0, X1, X2, X3,
                       X1**2, X2**2, X3**2, X1*X2, X2*X3, X3*X1, X1*X2*X3,
                       X1*X2**2, X1*X3**2, X2*X1**2, X2*X3**2, X3*X1**2, X3*X2**2,
                       X1**3, X2**3]
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



    def get_xCol_scaled(self, x):
        '''
        Same 10-term polynomial basis as get_xCol(), but built from
        pixel coordinates normalized by self.px_center / self.px_scale
        first. Used for the variable-origin extension (both fitting [C]
        and [B], and evaluating them in get_ray()) to keep the polynomial
        well-conditioned; with the default px_center=[0,0], px_scale=[1,1]
        this is identical to get_xCol().
        '''
        x1 = (x[0] - self.px_center[0]) / self.px_scale[0]
        x2 = (x[1] - self.px_center[1]) / self.px_scale[1]
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
        Given a point in 2D, x, this method returns its 3D direction vector.

        Note: kept for backwards compatibility. This does NOT normalize the
        vector and does NOT account for a variable origin. New code should
        generally prefer get_ray(), which handles both models correctly.
        '''
        x = [eta, zeta]
        xColumn = self.get_xCol(x)
        res = dot(xColumn, self.B)
        return array([res[0], res[1], res[2]])



    def get_ray(self, eta, zeta):
        '''
        Given a point in 2D camera-space, x = (eta, zeta), this method
        returns the epipolar line associated with it in lab-space, as an
        origin point, O, and a unit direction vector, e:

            l(x) = O + e * a

        If self.variable_origin is True (the camera was calibrated with the
        variable-origin extension), O is computed per-pixel from [C] using
        NORMALIZED pixel coordinates (get_xCol_scaled, via self.px_center /
        self.px_scale), in one of two parameterizations set by
        self.origin_mode:

            'free'  : O(x) = [C] G_c(x),  [C] shape (n_c, 3)
            'plane' : O(x) = self.O + plane_u*a(x) + plane_v*b(x),
                      where (a,b) = [C] G_c(x),  [C] shape (n_c, 2)

        [C] may also use a REDUCED polynomial order relative to [B] - see
        calibrate_variable_origin.fit_variable_origin_model(c_order=...).
        The basis terms are ordered by increasing degree
        (1, x, y, x^2, y^2, xy, x^2y, xy^2, x^3, y^3), so the first
        self.C.shape[0] terms of the full basis are exactly the reduced-
        order basis [C] was fit against - no separate stored "order" is
        needed, it's inferred from C's shape.

        Otherwise (self.variable_origin is False), the single fixed camera
        center self.O is used with the raw-pixel basis (get_xCol), exactly
        as in the original Extended Zolof model.
        '''
        x = [eta, zeta]

        if self.variable_origin:
            xColumn = array(self.get_xCol_scaled(x))
        else:
            xColumn = array(self.get_xCol(x))

        r = dot(xColumn, self.B)
        e = array([r[0], r[1], r[2]])
        nrm = norm(e)
        if nrm > 0:
            e = e / nrm

        if self.variable_origin:
            n_c = array(self.C).shape[0]   # C may use fewer basis terms than B
            coef = dot(xColumn[:n_c], self.C)

            if self.origin_mode == 'plane':
                # O(x) constrained to the calibrated plane through self.O
                O = (array(self.O)
                     + array(self.plane_u) * coef[0]
                     + array(self.plane_v) * coef[1])
            else:
                O = array([coef[0], coef[1], coef[2]])
        else:
            O = array(self.O)

        return O, e
    
    
    
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
        S += 'variable_origin %s \n' % (self.variable_origin,)

        for i in range(len(self.O)):
            S += 'O %s \t %s \n'%(i, self.O[i])
        
        for i in range(len(self.A)):
            for j in range(len(self.A[i])):
                S += 'A %s %s \t %s \n'%(i, j, self.A[i][j])
        
        for i in range(len(self.B)):
            for j in range(len(self.B[i])):
                S += 'B %s %s \t %s \n'%(i, j, self.B[i][j])

        if self.variable_origin:
            for i in range(len(self.C)):
                for j in range(len(self.C[i])):
                    S += 'C %s %s \t %s \n'%(i, j, self.C[i][j])

            for i in range(len(self.px_center)):
                S += 'PXC %s \t %s \n'%(i, self.px_center[i])
            for i in range(len(self.px_scale)):
                S += 'PXS %s \t %s \n'%(i, self.px_scale[i])

            S += 'ORIGIN_MODE %s \n'%(self.origin_mode,)
            if self.origin_mode == 'plane':
                for i in range(3):
                    S += 'PU %s \t %s \n'%(i, self.plane_u[i])
                for i in range(3):
                    S += 'PV %s \t %s \n'%(i, self.plane_v[i])
                for i in range(3):
                    S += 'PN %s \t %s \n'%(i, self.plane_n[i])
        
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

        # variable_origin flag - absent in files saved before this feature
        # existed, so default to False (fixed-origin) for backwards
        # compatibility.
        vo_vals = list(filter(lambda l: len(l)>0 and l[0]=='variable_origin',
                               lines))
        if len(vo_vals) > 0:
            self.variable_origin = (vo_vals[0][-1] == 'True')
        else:
            self.variable_origin = False
        
        O_vals = list(filter(lambda l: len(l)>0 and l[0]=='O', lines))
        for v in O_vals:
            ind = int(v[1])
            self.O[ind] = float(v[-1])
        
        A_vals = list(filter(lambda l: len(l)>0 and l[0]=='A', lines))
        Ashape = (int(max(A_vals, key=lambda x: float(x[1]))[1])+1,
                  int(max(A_vals, key=lambda x: float(x[2]))[2])+1)
        self.A = array([[0.0 for i in range(Ashape[0])] for j in range(Ashape[1])]).T
        for v in A_vals:
            i, j = int(v[1]), int(v[2])
            self.A[i][j] = float(v[-1])
        
        B_vals = list(filter(lambda l: len(l)>0 and l[0]=='B', lines))
        Bshape = (int(max(B_vals, key=lambda x: float(x[1]))[1])+1,
                  int(max(B_vals, key=lambda x: float(x[2]))[2])+1)
        self.B = array([[0.0 for i in range(Bshape[0])] for j in range(Bshape[1])]).T
        for v in B_vals:
            i, j = int(v[1]), int(v[2])
            self.B[i][j] = float(v[-1])

        # C coefficients - only for variable-origin
        C_vals = list(filter(lambda l: len(l)>0 and l[0]=='C', lines))
        if len(C_vals) > 0:
            Cshape = (int(max(C_vals, key=lambda x: float(x[1]))[1])+1,
                      int(max(C_vals, key=lambda x: float(x[2]))[2])+1)
            self.C = array([[0.0 for i in range(Cshape[0])] for j in range(Cshape[1])]).T
            for v in C_vals:
                i, j = int(v[1]), int(v[2])
                self.C[i][j] = float(v[-1])

        PXC_vals = list(filter(lambda l: len(l)>0 and l[0]=='PXC', lines))
        if len(PXC_vals) > 0:
            for v in PXC_vals:
                self.px_center[int(v[1])] = float(v[-1])

        PXS_vals = list(filter(lambda l: len(l)>0 and l[0]=='PXS', lines))
        if len(PXS_vals) > 0:
            for v in PXS_vals:
                self.px_scale[int(v[1])] = float(v[-1])

        OM_vals = list(filter(lambda l: len(l)>0 and l[0]=='ORIGIN_MODE',
                               lines))
        if len(OM_vals) > 0:
            self.origin_mode = OM_vals[0][-1]
        else:
            self.origin_mode = 'free'

        for tag, attr in [('PU', 'plane_u'), ('PV', 'plane_v'),
                          ('PN', 'plane_n')]:
            vals = list(filter(lambda l: len(l)>0 and l[0]==tag, lines))
            if len(vals) > 0:
                vec = array([0.0, 0.0, 0.0])
                for v in vals:
                    vec[int(v[1])] = float(v[-1])
                setattr(self, attr, vec)
        
    
            
    
