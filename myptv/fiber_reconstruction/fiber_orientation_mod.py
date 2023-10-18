# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 15:59:35 2022

@author: Eric Aschari

"""

import numpy as np
import math
from scipy import linalg

class FiberOrientation(object):
    '''A class to obtain the 3D fiber orientation from two fiber points on every segmented image'''
    
    def __init__(self, X: np.ndarray,
                 B: np.ndarray):
        self.X = X # Center of fiber
        self.B = B # Endpoints of fiber
 

    def image2fiber(self, cams):
        '''
        input:
            cams: array of cameras
        output:
            cAvg: point on 3D line passing through fiber
            uAvg: 3D orientation of fiber
        '''
        
        cs = self.X
        bs = self.B
        # alphas = self.get_alphas()
        s = np.shape(cams)[0]
        c = np.zeros((3, math.comb(s,2)));
        u = np.zeros((3, math.comb(s,2)));
        r = 0 # non-nan result counter
        
        # for all camera pairs calculate the intersecting line
        for i in range(0,s-1):
            for j in range((i + 1),s):
                # get_r_ori???
                p1_n, p1_m = self.getPlane(np.transpose(np.array([cams[i].O])),
                              cams[i].get_r_ori(cs[i]),
                              cams[i].get_r_ori(bs[i]))
                p2_n, p2_m = self.getPlane(np.transpose(np.array([cams[j].O])),
                              cams[j].get_r_ori(cs[j]),
                              cams[j].get_r_ori(bs[j]))
                ctemp,utemp = self.intersectPlanes(p1_n,p1_m,p2_n,p2_m)

                if type(utemp[0]) != 'NaN' and type(utemp[1]) != 'NaN' and type(utemp[2]) != 'NaN':
                    c[:,[r]] = ctemp
                    u[:,[r]] = utemp
                    r += 1
        
        cAvg,uAvg = self.averageLine(c,u)
        
        ori = self.get_ori(uAvg)
        
        return cAvg,uAvg,ori
    
    
    ### helper functions
    
    def solve_svd(self, A, b):
        '''
        input:
            A: matrix
            b: column vector
        output:
            x: vector so that A*x = b
        '''
        # compute svd of A
        U,s,Vh = linalg.svd(A)
        # U diag(s) Vh x = b <=> diag(s) Vh x = U.T b = c
        c = np.dot(U.T,b)
        # diag(s) Vh x = c <=> Vh x = diag(1/s) c = w (trivial inversion of a diagonal matrix)
        w = np.dot(np.diag(1/s),c)
        # Vh x = w <=> x = Vh.H w (where .H stands for hermitian = conjugate transpose)
        x = np.dot(Vh.conj().T,w)
        return x
    

    def getPlane(self, P1, P2, P3): 
        '''
        input:
            P123: three 3D points
        output:
            p_n: normal vector to plane spun by P123
            p_m: ax + by + cz = *m*
        '''
        A = np.array([[P1[0,0], P1[1,0], -1.0],
                      [P2[0,0], P2[1,0], -1.0],
                      [P3[0,0], P3[1,0], -1.0]])
        z = np.array([[-P1[2,0]], [-P2[2,0]], [-P3[2,0]]])
        # b = np.linalg.solve(A,z)
        # print('A___',A)
        # print('z___',z)
        b = self.solve_svd(A,z)
        
        p_n = np.array([[b[0,0]], [b[1,0]], [1]])
        p_m = b[2] / np.linalg.norm(p_n)
        p_n = p_n / np.linalg.norm(p_n)
        return p_n, p_m[0]


    def intersectPlanes(self, p1_n, p1_m, p2_n, p2_m): 
        '''
        input:
            p12_n: normal vectors of two planes
            p12_m: plane constant
        output:
            c,u: f(k) = c + k * u 
        '''
        n = np.concatenate([np.transpose(p1_n),np.transpose(p2_n)])
        if np.dot(np.transpose(p1_n),p2_n) > 0.999:  # acos(0.999) = 2.5° between planes
            c = float("nan")
            u = float("nan")
            return c,u

        c = np.linalg.lstsq(n,np.array([[p1_m],[p2_m]]),rcond=-1)[0]
        u = linalg.null_space(n)
        return c,u
    

    def averageLine(self, c, u): 
        '''
        input:
            c,u: f(k) = c + k * u 
        output:
            cAvg: average of c
            uAvg: average of u
        '''
        A = np.zeros(np.shape(u))
        for i in range(0,np.shape(A)[1]):
            A[:,i] = u[:,0]
        
        s = np.sign(np.dot(np.transpose(A),u))[0]
        u = np.multiply(s,u)
        uAvg = np.mean(u,1)
        cAvg = np.mean(c,1)
        return cAvg,uAvg
    

    def get_alphas(self):
        '''
        input:
            self: centers and endpoints of fibers
        output:
            alphas: array of angles on camera planes
        '''
        alpha1 = np.arctan2((self.B[0,1] - self.X[0,1]), (self.B[0,0] - self.X[0,0]))[0]
        alpha2 = np.arctan2((self.B[1,1] - self.X[1,1]), (self.B[1,0] - self.X[1,0]))[0]
        # alpha3 = np.arctan2((self.B[2,1] - self.X[2,1]), (self.B[2,0] - self.X[2,0]))[0]
        # alpha4 = np.arctan2((self.B[3,1] - self.X[3,1]), (self.B[3,0] - self.X[3,0]))[0]
        
        # alphas = np.array([alpha1, alpha2, alpha3, alpha4])
        alphas = np.array([alpha1, alpha2])

        return alphas
    
    
    def get_ori(self,u):
        '''
        input:
            u: direction vector
        output:
            ori: angle array []
        '''
        xy_ori = np.arctan2(u[1],u[0])
        xz_ori = np.arctan2(u[2],u[0])
        ori = np.array([xy_ori,xz_ori])
        
        return ori
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    