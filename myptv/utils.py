# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 13:33:57 2018

@author: ron


Utility code to use for the MyPTV package.

"""



from numpy import dot, array, loadtxt, savetxt
from numpy import append as NPappend
from numpy.linalg import inv, norm






def line_dist(O1, r1, O2, r2):
    '''
    2 lines are defined as (O1 + a r1) and  (O2 + b r2), where O are origins,
    r are direction vectors, and a,b are variables (in n dimensions). 
    This utility calculates the minimal distance between these 2 lines.
    
    input - 
    O1,O2,r1,r2 (arrays, n) - line parameters
    
    output - 
    dist (float) -the minimum distance between the lines
    x (array, n)- the point that is nearest to the two points crossing
    '''
    
    # find the a,b that minimize the distance:
    # A = array([[dot(r1,r1), -dot(r1,r2)],
    #            [dot(r1,r2), -dot(r2,r2)]])
    
    r1r2 = r1[0]*r2[0] + r1[1]*r2[1] + r1[2]*r2[2]
    r12 = r1[0]**2 + r1[1]**2 + r1[2]**2
    r22 = r2[0]**2 + r2[1]**2 + r2[2]**2

    #Ainv = [[-r22, r1r2],
    #        [-r1r2, r12]]
    
    dO = O2-O1
    B = [r1[0]*dO[0] + r1[1]*dO[1] + r1[2]*dO[2],
         r2[0]*dO[0] + r2[1]*dO[1] + r2[2]*dO[2]]
    
    try:
        #a,b = dot(Ainv, B)
        a = (-r22*B[0] + r1r2*B[1])/(r1r2**2 - r12 * r22)
        b = (-r1r2*B[0] + r12*B[1])/(r1r2**2 - r12 * r22)
    except:
        a,b = 0.0, 0.0
    
    # use the a,b to calc the minimum distance:
    l1,l2 = O1 + a*r1 , O2 + b*r2
    dist = sum((l1 - l2)**2)**0.5
    x = (l1+l2)*0.5
    
    return dist, x



def point_line_dist(O,r,P):
    '''
    for a line (O + a r) and a point P, this returns the distance between
    the line and the point.
    '''
     
    anum = sum([ r[i]*(P[i]-O[i]) for i in range(3)])
    adenum = sum([ r[i]*r[i] for i in range(3)])
    a = anum / adenum
    l = [O[i] + a*r[i] for i in range(3)]
    d = ((l[0]-P[0])**2 + (l[1]-P[1])**2 + (l[2]-P[2])**2)**0.5
    return d






class line(object):
    
    def __init__(self, O, r):
        '''
        A helper class for calculating line and point distances.
        O - a point through which the line passes
        r - a vector of the line direction
        '''
        self.O = O
        self.r = r
        
    def distance_to_point(self, P):
        '''
        Returns the least distance between the line and a given point P.
        '''
        s1 = sum([self.r[i]*(P[i]-self.O[i]) for i in [0,1,2]])
        s2 = sum([ri**2 for ri in self.r])
        amin = s1 / s2
        
        dmin = sum([(self.O[i]-P[i]+amin*self.r[i])**2 for i in [0,1,2]])**0.5
        
        return dmin, amin, self.O + amin*self.r
        
    


def get_nearest_line_crossing(line_list):
    '''
    Given a list of line objects, this function find a point that minimizes
    the sum of the distances to it.
    '''
    from scipy.optimize import minimize
    func = lambda P: sum([l.distance_to_point(P)[0] for l in line_list])
    P = minimize(func, array([0,0,0])).x
    return P









def fit_polynomial(x, y, n):
    '''
    A polynomial of degree n is written as:
        
        p_n(x) = an*x^n + ... + a0
        
    This function solves the least squares equation minimize  
    the coefficients an ... a0 that minimize the residuals to 
    the data points (x_i, y_i):  
        R = sum_i | p_n(x_i) - y_i |^2
        
    input - 
    x - a list of observations for the dependent variable 
    y - a list of observations for the independent variable 
    n - degree of polynomial to fit
    
    output - 
    an - a list of n+1 polynomial coefficients.
    '''
    N = len(y)
    X = []
    
    for i in range(N):
        X.append([])
        for j in range(n, -1, -1):
            X[-1].append(x[i]**j)
    
    X = array(X)
    A = dot( inv(dot(X.T, X)) , X.T )
    an = dot(A, y)
    return an
    







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






class match_calibration_blobs_and_points(object):
    '''
    This is a tool used to matxh between a calibration "target file" and
    the blobs extracted from the image. It takes in a pre-calibrated 
    camera object, the target file data and the blob file data.
    
    To do the pairing, use 
    match_calibration_blobs_and_points.pair_points() 
    
    To save the results, use
    match_calibration_blobs_and_points.save_results(fname)
    
    '''
    
    def __init__(self, camera, blob_fname, target_file_fname):
        '''
        input - 
        camera - an instance of the camera class, already calibrated with at
                 least a rough calibration
        blob_fname - string, path of the file containing the segmented 
                     calibration points
        target_file_fname - string, path to the calibration target file.
        '''
        self.cam = camera
        self.blobs = loadtxt(blob_fname)
        self.targets = loadtxt(target_file_fname)
        self.get_neighbor_separation()
        
    
    def pair_points(self):
        '''
        Does the actual pairing of blobs and target points.
        The paired points are stored in the attribute self.point_pairs.
        '''
        
        target_points_eta_zeta = [self.cam.projection(p) for p in self.targets]
        point_pairs = []
        
        for i in range(len(self.blobs)):
            j_min = 10000
            d_min = 1000
            for j in range(len(target_points_eta_zeta)):
                d = norm(target_points_eta_zeta[j] - self.blobs[i][1::-1])
                if d<d_min:
                    d_min = d
                    j_min = j
            
            if d_min < self.dmin/2:
                blob_point_pair = NPappend(self.blobs[i][1::-1],
                                            self.targets[j_min])
                point_pairs.append(blob_point_pair)
        
        self.point_pairs = point_pairs
        
        
    def get_neighbor_separation(self):
        '''
        returns an estimate for the minimum separation distance between blobs.
        '''
        bxy = self.blobs[:,1::-1]
        dmin = 1e20
        N = len(bxy)
        i=0
        
        for i in range(min([N,20])):
            for j in range(N):
                if j==i: continue
                d = sum((bxy[i] - bxy[j])**2)**0.5
                if d<dmin:
                    dmin = d
        print('dmin: ', dmin)
        self.dmin = dmin
        
        
        
    def save_results(self, fname):
        '''
        Saves the pairs of blobs and target points in a given file name.
        '''
        fmt = ['%.2f', '%.2f', '%.2f', '%.2f', '%.2f']
        savetxt(fname, self.point_pairs, fmt=fmt, delimiter='\t')
        
        
    def plot_projections(self):
        '''
        Plots the projection of the target points on the camera, and the 
        blobs that were segmented.
        
        Requires matplotlib.
        '''
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        
        for p in self.targets:
            eta, zeta = self.cam.projection(p)
            ax.plot(eta, zeta, 'ob')
            
        ax.plot(self.blobs[:,1], self.blobs[:,0], 'rx') 
        
        ax.plot([],[],'ob', label='target points')
        ax.plot([],[],'rx', label='segmented blobs')
        ax.legend()
        return fig, ax
    
    
    


import pandas as pd
import numpy as np

class get_residual_blobs(object):
    '''
    A class to retriev blobs that have not been used up in the matching step.
    '''
    
    def __init__(self, blob_filess, particles_file):
        
        self.particles = dict([(k, g.values.tolist()) for k,g in 
                               pd.read_csv(particles_file, 
                                           sep='\t', 
                                           header=None).groupby(7)])
        self.blobfiles = blob_filess
        
    
    def get_residual_blobs(self, frame_lst=None):
        
        self.residual_blobs = []
        
        for i in range(len(self.blobfiles)):
            blob_column = i+3
            blobs_i = dict([(k, g.values.tolist()) for k,g in 
                               pd.read_csv(self.blobfiles[i], 
                                           sep='\t', header=None).groupby(5)])
            
            unused_blobs_i = []
            
            if frame_lst is None:
                lst = blobs_i.keys()
            
            else:
                lst = frame_lst
            
            for k in lst:
                
                if k not in list(self.particles.keys()):
                    unused_blobs_i += blobs_i[k]
                
                else:
                    used_blob_inds = [blb[blob_column] for blb in self.particles[k]]
                    all_ind = list(range(len(blobs_i[k])))
                    unused_inds = list(set(all_ind).difference(set(used_blob_inds)))
                    unused_blobs_i += [blobs_i[k][j] for j in unused_inds]
                
            self.residual_blobs.append(unused_blobs_i)
    
    
    
    def save_residual_blobs(self):
        
        fmt = ['%.3f', '%.3f', '%d', '%d', '%d', '%d']
        
        for i in range(len(self.residual_blobs)):
            fname = 'residual_blobs_%d'%i
            np.savetxt(fname, self.residual_blobs[i],delimiter='\t', fmt=fmt)
         



