#!/usr/bin/env python2
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
    A = array([[dot(r1,r1), -dot(r1,r2)],
               [dot(r1,r2), -dot(r2,r2)]])
    
    B = array([[dot(r1,O2-O1)],[dot(r2,O2-O1)]])
    
    try:
        a,b = dot(inv(A), B)
    except:
        a,b = 0.0, 0.0
    
    # use the a,b to calc the minimum distance:
    l1,l2 = O1 + a*r1 , O2 + b*r2
    dist = norm(l1 - l2)
    x = (l1+l2)*0.5
    
    return dist, x



def point_line_dist(O,r,P):
    '''
    for a line (O + a r) and a point P, this returns the distance between
    the line and the point.
    '''
    a = dot(r, P - O) / dot(r, r)
    l = O + a*r
    d = norm(l-P)
    return d



def fit_polynomial(x, y, n):
    '''
    A polynomial of degree n is written as:
        
        p_n(x) = an*x^n + ... + a0
        
    This function solves the least squares equation minimize  
    the coefficients an ... a0 that minimize the resudials to 
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
            
            if d_min < 30:
                blob_point_pair = NPappend(self.blobs[i][1::-1],
                                            self.targets[j_min])
                point_pairs.append(blob_point_pair)
        
        self.point_pairs = point_pairs
        
        
    def save_results(self, fname):
        '''
        Saves the pairs of blobs and target points in a given file name.
        '''
        fmt = ['%.3f', '%.3f', '%.3f', '%.3f', '%.3f']
        savetxt(self.point_pairs, fmt=fmt, delimeter='\t')
        
        
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
        ax.plot([],[],'rx', label='segmented plobs')
        ax.legend()
        return fig, ax