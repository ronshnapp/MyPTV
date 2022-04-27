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
    
    
    





#=============================================================================
#Line intersects. Not in use at the moment - only, testing.


def find_point_nearest_to_lines(line_list):
    '''
    Will use least squares in order to solve the problem of finding the
    point that is the closest to several lines. (Note that for two lines there
    is an analytical solution that is used in the function line_dist()).
    
    This function uses scipy.optimize.least_squares() !
    
    input -
    line_list - a list of tuples (O, r) where O and r are the origin and 
                direction vectors (numpy arrays, 3) that define a line in 
                3d space.
                
    returns -
    P - the point in 3d space that minimizes the distance to the lines in the 
        list.
    '''
    from scipy.optimize import least_squares
    
    def dist_sum(P):
        dist = 0
        for O,r in line_list:
            dist += point_line_dist(O, r, P)
        return dist
    
    x0 = line_dist(line_list[0][0], line_list[0][1], 
                   line_list[1][0], line_list[1][1])[1]
    P = least_squares(dist_sum, x0, method='dogbox', ftol=1e-8)
    return P
    




def nearest_intersect(lines):
    import numpy as np
    I = np.eye(3)
    
    S1 = 0
    S2 = 0
    for i in range(len(lines)):
        Oi = np.array([lines[i][0]])
        ri = lines[i][1]
        S1 += I - np.dot(Oi.T, Oi)
        S2 += np.dot(I - np.dot(Oi.T, Oi), ri)
    return np.linalg.lstsq(S1, S2, rcond=None)
        
#=============================================================================



