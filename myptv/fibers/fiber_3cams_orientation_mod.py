# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 15:59:35 2022

@author: Alessandro Gambino

"""

from numpy import linalg
from tqdm import tqdm
from PIL import Image, ImageTk
from tkinter import Label, Canvas, LabelFrame, Entry, Tk, Scrollbar, Button
import os
import math
import numpy as np
from skimage.io import imsave, imread #save and read
from numpy import array, dot
import math
import matplotlib.pyplot as plt
from matplotlib.pyplot import imread
import matplotlib.pyplot as plt #draw fibers
from matplotlib.patches import FancyArrowPatch #draw fibers

# Modules from the main_cam_reproj
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Modules to save
from os import getcwd, listdir
from os.path import exists as pathExists

class Camera:
    def __init__(self, name, origin):
        self.name = name
        self.O = None
        self.theta = None
        self.f = None
        self.xh = None
        self.yh = None
        self.E = np.zeros((5, 3))
        self.resolution = None
        self.R = None
        self.pos2D = None
        self.ori2D = None
        self.frames2D = None
        self.pos3D = None
        self.ori3D = None
        self.blob3D = None  #the blob index in the given camera
        self.frames3D = None
        self.error3D = None
        self.target3D = None
        self.indtraj = None

        self.indline3D = None
        self.badline3D = None #flag for badlines to consider in the orientation output


def read_camera_data(folder_path, cam_file):
    cam_path = os.path.join(folder_path, cam_file)
    res = None
    O = None
    theta = None
    f = None
    xh = None
    yh = None
    E = []

    with open(cam_path, 'r') as file:
        lines = file.readlines()
        res = list(map(int, lines[2].strip().split()))
        O = list(map(float, lines[3].strip().split()))
        theta = list(map(float, lines[4].strip().split()))
        f = list(map(float, lines[5].strip().split()))
        xh, yh = map(float, lines[6].strip().split())
        for line in lines[7:]:
            row = list(map(float, line.strip().split()))
            E.append(row)

    #print('ecco E in componenti:', E[0][:])
    return res, O, theta, f, xh, yh, E


def initialize_camera(camera, name, res, O, theta, f, xh, yh, E):
    from math import sin, cos
    from numpy import dot

    camera.name = name
    camera.O = np.array(O)
    camera.theta = np.array(theta)
    camera.f = f[0]
    camera.xh = xh
    camera.yh = yh
    camera.E = np.array([E[0][:], E[1][:], [0, 0, 0, 0, 0]]) 
    camera.resolution = res

    tx,ty,tz = camera.theta
    Rx = array([[1,0,0],
                [0,cos(tx),-sin(tx)],
                [0,sin(tx),cos(tx)]])
    Ry = array([[cos(ty),0,sin(ty)],
                [0,1,0],
                [-sin(ty),0,cos(ty)]])
    Rz = array([[cos(tz),-sin(tz),0],
                [sin(tz),cos(tz),0],
                [0,0,1]])
    camera.R = dot(dot(Rx,Ry), Rz)


def assign_2d_positions_and_orientations(camera, file_path, fiber=False):
    positions = []
    orientations = []
    frames = []

    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            data = list(map(float, line.strip().split()))
            x, y = data[:2]
            positions.append((x, y))
            fr = data[5]
            frames.append(fr)

            if fiber:
                ori_x, ori_y = data[6:8]
                orientations.append((ori_x, ori_y))

    camera.pos2D = np.array(positions)
    camera.frames2D = np.array(frames)
    #print('cam pos 2D: ', camera.pos2D)
    if fiber:
        camera.ori2D = np.array(orientations)
    #print('cam ori 2D: ', camera.ori2D)



def assign_3d_targetfile(camera, file_path):
    positions = []

    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            data = list(map(float, line.strip().split()))
            x, y, z = data[:3] 
            positions.append((x, y, z))

    camera.target3D = np.array(positions)
    #print('cam target 3D: ', camera.target3D)




def assign_3d_positions(camera, camind, camnum, file_pos_path):
    indtraj = []
    positions = []
    frames = []
    error_t = []
    blob_ind = []
    indline = []
    badline = [] #lines which contain faulty particle positions

    with open(file_pos_path, 'r') as file:
        lines = file.readlines()
        for iline, line in enumerate(lines):
            data = list(map(float, line.strip().split()))
            ind = data[0] #"trajectories"
            x, y, z = data[1:4] #"trajectories"
            bl = data[3+camind] #"trajectories"
            err = data[2*camnum+1] #"trajectories"
            fr = data[2*camnum+2] #"trajectories"

            if err == 100: # "missing frame flag"
                badline.append(iline)

            indtraj.append(ind)
            positions.append((x, y, z))
            blob_ind.append(bl)
            error_t.append(err)
            frames.append(fr)

            indline.append(iline)

    #print(np.argwhere(np.array(error_t)>100))
    camera.indtraj = np.array(indtraj)
    camera.pos3D = np.array(positions)
    camera.frames3D = np.array(frames)
    camera.blob3D = np.array(blob_ind)
    camera.indline3D = np.array(indline)
    camera.badline3D = np.array(badline)
    camera.error3D = np.array(error_t)

# =============================================================================
# Functions for computing orientations.  

def get_r(self, eta, zeta):
    '''
    r = ([eta, zeta, f] + e) * [R]
        
    input - pixel coordinates (eta, zeta) seen by the camera
    output - direction vector in real space
    '''
    eta_ = eta - self.resolution[0]/2.0 - self.xh
    zeta_ = zeta - self.resolution[1]/2.0  - self.yh
    Z3 = [eta, zeta, eta**2, zeta**2, eta * zeta]
        
    e = dot(self.E, Z3)
        
    r = dot(array([-eta_, -zeta_, -self.f]) - e, self.R)
    r = r / (r[0]**2 + r[1]**2 + r[2]**2)**0.5 #no normalization
    rout = np.array(r)

    return rout


def getPlane(P1, P2, P3): 
    '''
    input:
        P123: three 3D points
    output:
        p_n: normal vector to plane spun by P123
        p_m: ax + by + cz = *m*
    '''
    A = np.array([[P1[0], P1[1], -1.0],
                [P2[0], P2[1], -1.0],
                [P3[0], P3[1], -1.0]])
    z = np.array([[-P1[2]], [-P2[2]], [-P3[2]]])
    b = solve_svd(A,z)
        
    p_n = np.array([b[0,0], b[1,0], 1])
    nrm = np.linalg.norm(p_n)
    p_n = p_n / nrm
    #p_m = b[2] / np.linalg.norm(p_n)
    #return p_m[0]
    return p_n

def solve_svd(A, b):
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

def intersectPlanes(p1, p2): 
    '''
    input:
        p12_n: normal vectors of two planes
        p12_m: plane constant
    output:
        c,u: f(k) = c + k * u 
    '''
    angle = np.abs(np.dot(p1,p2))
    if angle > 0.995:  # acos(0.995) = 5.73196797° between planes 
        a = np.ones((3))
        a = np.nan
        #print('size a: ', np.size(a))
        return a
    u = np.cross(p1, p2)
    #print('size u: ', np.size(u))
    return u

def is_invalid_vector(vector):
    return np.all(np.isnan(vector)) or np.all(vector == 0)

def do_orientations(self1, self2, self3, min_err): #main function to compute orientations
    
    ori12 = np.zeros((self1.pos3D.shape[0], 3)) #8 
    ori13 = np.zeros((self1.pos3D.shape[0], 3)) #8
    ori23 = np.zeros((self1.pos3D.shape[0], 3)) #8
    error12 = np.zeros((self1.pos3D.shape[0], 1)) #8 
    error13 = np.zeros((self1.pos3D.shape[0], 1)) #8
    error23 = np.zeros((self1.pos3D.shape[0], 1)) #8

    for f in tqdm(self1.frames3D): # frame number
        #print('ecco f: ', f, ' out of ', cam1.frames3D.shape[0])
        # indices3D = indices for the trajectories at frame i
        indices3D = [index for index, value in enumerate(self1.frames3D) if value == f]
        #print('ecco indices3D: ', indices3D)
        #print('ecco f: ', f)
        u12, u13, u23, err12, err13, err23 = do_ori_cams(self1, self2, self3, indices3D,f)
   
        ori12[indices3D ,:] = u12
        ori13[indices3D ,:] = u13
        ori23[indices3D ,:] = u23
        error12[indices3D ,:] = err12
        error13[indices3D ,:] = err13
        error23[indices3D ,:] = err23

    # Create the final orientations based on reprojection errors
    errtot = np.hstack((error12, error13, error23)) 
    min_indices = np.argmin(np.nan_to_num(errtot, nan=1), axis=1)  
    errmin = errtot[np.arange(errtot.shape[0]), min_indices].reshape(-1, 1)
    orimin = np.zeros_like(ori12)
    for i in range(len(min_indices)):
        if min_indices[i] == 0:
            orimin[i, :] = ori12[i, :]
        elif min_indices[i] == 1:
            orimin[i, :] = ori13[i, :]
        elif min_indices[i] == 2:
            orimin[i, :] = ori23[i, :]

    
    # Prepare the output file
    res12 = np.zeros((self1.pos3D.shape[0], 9)) 
    res13 = np.zeros((self1.pos3D.shape[0], 9)) 
    res23 = np.zeros((self1.pos3D.shape[0], 9)) 
    resmin = np.zeros((self1.pos3D.shape[0], 9)) 

    indd12 = self1.indtraj.reshape(-1, 1)
    indd13 = self1.indtraj.reshape(-1, 1)
    indd23 = self1.indtraj.reshape(-1, 1)
    inddmin = self1.indtraj.reshape(-1, 1)

    blob1 = self1.blob3D.reshape(-1, 1)
    blob2 = self2.blob3D.reshape(-1, 1)
    blob3 = self3.blob3D.reshape(-1, 1)
    frames = self1.frames3D.reshape(-1, 1)

    error12 = error12.reshape(-1, 1)
    error13 = error13.reshape(-1, 1)
    error23 = error23.reshape(-1, 1)
    # =====================================
    # Fiber pre-processing. 
    # =====================================

    #1) clean high error
    mask12 = error12 > min_err
    mask13 = error13 > min_err
    mask23 = error23 > min_err
    maskmin = errmin > min_err
    ori12[mask12.flatten()] = np.nan
    ori13[mask13.flatten()] = np.nan
    ori23[mask23.flatten()] = np.nan
    orimin[maskmin.flatten()] = np.nan
    print('Removed data with angle error higher than ', min_err, ': ', sum(np.logical_and(errmin > min_err, errmin < 1))[0], 'out of ', orimin.shape[0])
    
    #2) sign switch
    ori12 = restore_sign_switches_multiple_trajectories(indd12, ori12)
    ori13 = restore_sign_switches_multiple_trajectories(indd13, ori13)
    ori23 = restore_sign_switches_multiple_trajectories(indd23, ori23)
    orimin = restore_sign_switches_multiple_trajectories(inddmin, orimin)

    #3) normalization of orientation vectors
    ori12 = ori12 / np.linalg.norm(ori12, axis=1, keepdims=True)
    ori13 = ori13 / np.linalg.norm(ori13, axis=1, keepdims=True)
    ori23 = ori23 / np.linalg.norm(ori23, axis=1, keepdims=True)
    orimin = orimin / np.linalg.norm(orimin, axis=1, keepdims=True)

    #4) interpolate NaN values
    #orimin = polyFitNaNs(orimin, inddmin, maskmin, 3)

    res12 = np.concatenate((indd12,ori12, blob1, blob2, blob3, error12, frames), axis=1)
    res13 = np.concatenate((indd13,ori13, blob1, blob2, blob3, error13, frames), axis=1)
    res23 = np.concatenate((indd23,ori23, blob1, blob2, blob3, error23, frames), axis=1)
    resmin = np.concatenate((inddmin,orimin, blob1, blob2, blob3, errmin, frames), axis=1)

    return res12, res13, res23, resmin

def do_ori_cams(self1, self2, self3, indices3D, f):

    N = len(indices3D)

    # Extract blob indices for each camera at the given 3D indices
    blob3DFrame_1 = self1.blob3D[indices3D]
    blob3DFrame_2 = self2.blob3D[indices3D]
    blob3DFrame_3 = self3.blob3D[indices3D]

    pos3DFrame = self1.pos3D[indices3D]

    # Find corresponding frames in each segmented blob file
    indices2D_1 = [index for index, value in enumerate(self1.frames2D) if value == f]
    indices2D_2 = [index for index, value in enumerate(self2.frames2D) if value == f]
    indices2D_3 = [index for index, value in enumerate(self3.frames2D) if value == f]

    pos2DFrame_1 = self1.pos2D[indices2D_1]
    pos2DFrame_2 = self2.pos2D[indices2D_2]
    pos2DFrame_3 = self3.pos2D[indices2D_3]
    ori2DFrame_1 = self1.ori2D[indices2D_1]
    ori2DFrame_2 = self2.ori2D[indices2D_2]
    ori2DFrame_3 = self3.ori2D[indices2D_3]

    # Reverse the order of columns (y, x) to (x, y)
    pos2DFrame_1 = pos2DFrame_1[:, ::-1]
    pos2DFrame_2 = pos2DFrame_2[:, ::-1]
    pos2DFrame_3 = pos2DFrame_3[:, ::-1]
    ori2DFrame_1 = ori2DFrame_1[:, ::-1]
    ori2DFrame_2 = ori2DFrame_2[:, ::-1]
    ori2DFrame_3 = ori2DFrame_3[:, ::-1]

    # Initialize result arrays
    u12 = np.zeros((N, 3))
    u13 = np.zeros((N, 3))
    u23 = np.zeros((N, 3))
    error12 = np.zeros((N, 1))
    error13 = np.zeros((N, 1))
    error23 = np.zeros((N, 1))

    # Convert cam1.badline3D to a set for efficient look-up
    bad_lines_set = set(self1.badline3D)

    for i, index in enumerate(indices3D):
        if index in bad_lines_set:
            # Initialize u12, u13, and u23 to nan and ones for bad lines
            u12[i, :] = np.full((1, 3), np.nan)
            u13[i, :] = np.full((1, 3), np.nan)
            u23[i, :] = np.full((1, 3), np.nan)
            error12[i] = 1
            error13[i] = 1
            error23[i] = 1
            continue
        
        # Proceed with normal computation for non-bad lines
        mask_1 = blob3DFrame_1[i] >= 0
        mask_2 = blob3DFrame_2[i] >= 0
        mask_3 = blob3DFrame_3[i] >= 0

        if mask_1:
            idx_1 = int(blob3DFrame_1[i])
            XFrame_1 = pos2DFrame_1[idx_1]
            BFrame_1 = XFrame_1 + 200 * ori2DFrame_1[idx_1]
            rX_1 = get_r(self1, XFrame_1[0], XFrame_1[1])
            rB_1 = get_r(self1, BFrame_1[0], BFrame_1[1])
            p_1 = getPlane(self1.O, self1.O + rX_1, self1.O + rB_1)
        if mask_2:
            idx_2 = int(blob3DFrame_2[i])
            XFrame_2 = pos2DFrame_2[idx_2]
            BFrame_2 = XFrame_2 + 200 * ori2DFrame_2[idx_2]
            rX_2 = get_r(self2, XFrame_2[0], XFrame_2[1])
            rB_2 = get_r(self2, BFrame_2[0], BFrame_2[1])
            p_2 = getPlane(self2.O, self2.O + rX_2, self2.O + rB_2)
        if mask_3:
            idx_3 = int(blob3DFrame_3[i])
            XFrame_3 = pos2DFrame_3[idx_3]
            BFrame_3 = XFrame_3 + 200 * ori2DFrame_3[idx_3]
            rX_3 = get_r(self3, XFrame_3[0], XFrame_3[1])
            rB_3 = get_r(self3, BFrame_3[0], BFrame_3[1])
            p_3 = getPlane(self3.O, self3.O + rX_3, self3.O + rB_3)

        if mask_1 and mask_2 and not mask_3:
            u12[i, :] = intersectPlanes(p_1, p_2)
            error12[i] = 0
            error13[i] = 1
            error23[i] = 1
        if mask_1 and mask_3 and not mask_2:
            u13[i, :] = intersectPlanes(p_1, p_3)
            error12[i] = 1
            error13[i] = 0
            error23[i] = 1
        if mask_2 and mask_3 and not mask_1:
            u23[i, :] = intersectPlanes(p_2, p_3)
            error12[i] = 1
            error13[i] = 1
            error23[i] = 0
        if mask_1 and mask_2 and mask_3:
            u12[i, :] = intersectPlanes(p_1, p_2)
            u13[i, :] = intersectPlanes(p_1, p_3)
            u23[i, :] = intersectPlanes(p_2, p_3)
            error12[i] = oriError(self3,pos3DFrame[i,:],u12[i, :],ori2DFrame_3[idx_3])
            error13[i] = oriError(self2,pos3DFrame[i,:],u13[i, :],ori2DFrame_2[idx_2])
            error23[i] = oriError(self1,pos3DFrame[i,:],u23[i, :],ori2DFrame_1[idx_1])

    '''
    error12 = np.nan_to_num(error12, nan=1)
    error13 = np.nan_to_num(error13, nan=1)
    error23 = np.nan_to_num(error23, nan=1)
    '''

    return u12, u13, u23, error12, error13, error23


def oriError(self,x,p,pframe):
    proj = projectionFiber(self, x, p, 100, correction=True) 
    dx = proj[4] - proj[2]
    dy = proj[5] - proj[3]
    d_unit = np.array([dx, dy])/np.sqrt(dx**2 + dy**2)
    ori = np.array(pframe)
    ori = ori/np.sqrt(ori[0]**2 + ori[1]**2)
    cos = np.dot(d_unit, ori)
    #print('ecco proj: ', proj)
    return 1-np.abs(cos)


def restore_sign_switches_multiple_trajectories(indd, ori):

    corrected_orientations = np.zeros_like(ori)
    unique_indices = np.unique(indd)
    
    for index in unique_indices:
        trajectory_indices = np.where(indd == index)[0]
        trajectory_orientations = ori[trajectory_indices,:]
        corrected_trajectory = np.zeros_like(trajectory_orientations)
        corrected_trajectory[0, :] = trajectory_orientations[0, :]
        for i in range(1, trajectory_orientations.shape[0]):
            dot_product = np.dot(corrected_trajectory[i-1,:], trajectory_orientations[i,:])
            
            if dot_product < 0:
                corrected_trajectory[i,:] = -trajectory_orientations[i,:]
            else:
                corrected_trajectory[i,:] = trajectory_orientations[i,:]
        
        corrected_orientations[trajectory_indices,:] = corrected_trajectory
    
    return corrected_orientations



def polyFitNaNs(ori, indd, mask,deg):
    return interp_ori

# =============================================================================
# Functions for computing back reprojections for errors.  

def projectionFiber(camera, x, p, L, correction=True):
        '''
        will return the image coordinate (eta, zeta) of a real point x.
        
        input - x (array,3) - real world coordinates
                correction - if True, will return the coordinates after
                the non-linear error correction. If False, we not do the
                correction.
        output - (eta, zeta) (array,2) - camera coordinates of the projection 
                                         of x
        '''

        scaleF = L
        #print('ecco x:', x)
        #print('ecco p:', p)
        xp1 = x + scaleF*p
        xp2 = x - scaleF*p
        #print('ecco xp1:', xp1)
        #print('ecco xp2:', xp2)
        b = x - camera.O
        bp1 = xp1 - camera.O
        bp2 = xp2 - camera.O

        v = dot(b, camera.R.T) #inv(camera.R))
        vp1 = dot(bp1, camera.R.T) 
        vp2 = dot(bp2, camera.R.T) 

        a =  v[2] / camera.f
        ap1 =  vp1[2] / camera.f
        ap2 =  vp2[2] / camera.f

        eta_ = v[0] / a  + camera.resolution[0]/2 + camera.xh
        eta_p1 = vp1[0] / ap1  + camera.resolution[0]/2 + camera.xh
        eta_p2 = vp2[0] / ap2  + camera.resolution[0]/2 + camera.xh

        zeta_ = v[1] / a + camera.resolution[1]/2 + camera.yh
        zeta_p1 = vp1[1] / ap1 + camera.resolution[1]/2 + camera.yh
        zeta_p2 = vp2[1] / ap2 + camera.resolution[1]/2 + camera.yh

        # if we add the error correction term.
        if correction:
            eta, zeta = eta_zeta_from_bRinv(camera, eta_, zeta_)
            etap1, zetap1 = eta_zeta_from_bRinv(camera, eta_p1, zeta_p1)
            etap2, zetap2 = eta_zeta_from_bRinv(camera, eta_p2, zeta_p2)
            return array([eta, zeta, etap1, zetap1, etap2, zetap2])
        
        # if we don't use the error correction term.
        else:
            return array([eta_, zeta_, eta_p1, zeta_p1, eta_p2, zeta_p2])

        
def eta_zeta_from_bRinv(camera, eta_, zeta_):
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
        e_ = dot(camera.E, Z3)
        
        # calculating the derivatives of the error term:
        e_0 = e_[0]
        a, b, c, d, ee = camera.E[0,:]
        e_eta_0 = a + 2*c*eta_ + ee*zeta_
        e_zeta_0 = b + 2*d*zeta_ + ee*eta_
        #a, b, c, d, ee, f, g, h, i = camera.E[0,:]
        #e_eta_0 = a + 2*c*eta_ + ee*zeta_ + 3*f*eta_**2 + 2*g*eta_*zeta_ + h*zeta_**2
        #e_zeta_0 = b + 2*d*zeta_ + ee*eta_ + g*eta_**2 + 2*h*eta_*zeta_ + 3*i*zeta_**2
        
        e_1 = e_[1]
        a, b, c, d, ee = camera.E[1,:]
        e_eta_1 = a + 2*c*eta_ + ee*zeta_
        e_zeta_1 = b + 2*d*zeta_ + ee*eta_
        #a, b, c, d, ee, f, g, h, i = camera.E[1,:]
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

def run_3cams_orientation(ori_lim,cam_names,blob_fn,trajectory_file,save_name):
        local_path = r'./'
        camn = np.shape(cam_names)[0]

        cam1 = Camera('0', (0., 0.))
        cam2 = Camera('1', (1., 1.))
        cam3 = Camera('2', (2., 2.))

        initialize_camera(cam1, '1', *read_camera_data(local_path, cam_names[0]))
        initialize_camera(cam2, '2', *read_camera_data(local_path, cam_names[1]))
        initialize_camera(cam3, '3', *read_camera_data(local_path, cam_names[2]))

        segm1_path = os.path.join(blob_fn[0])
        segm2_path = os.path.join(blob_fn[1])
        segm3_path = os.path.join(blob_fn[2])
        pos3D_path = os.path.join(trajectory_file) 
        save_path = os.path.join(local_path, save_name)

        assign_2d_positions_and_orientations(cam1, segm1_path, fiber=True)
        assign_2d_positions_and_orientations(cam2, segm2_path, fiber=True)
        assign_2d_positions_and_orientations(cam3, segm3_path, fiber=True)
        assign_3d_positions(cam1, 1, camn, pos3D_path)
        assign_3d_positions(cam2, 2, camn, pos3D_path)
        assign_3d_positions(cam3, 3, camn, pos3D_path)

        res12, res13, res23, res = do_orientations(cam1, cam2, cam3, ori_lim)
    
    
        if save_name is not None:
            cwd_ls = listdir(getcwd())
            if save_name in cwd_ls or pathExists(save_name):
                print('\n The file name "%s" already exists in'%save_name)
                print(' the working directory. Should I save anyways?')
                usr = input('(1=yes, else=no)')
                if usr == '1':
                    print('\n', 'Saving the fiber orientation data (%s).'%save_name)
                    np.savetxt(save_path, res, fmt = ['%d','%.3f', '%.3f', '%.3f','%d','%d','%d','%.8f','%d'], delimiter='\t')
                else:
                    print('\n', 'Skipped saving file.')
        
            else:
                print('\n', 'Saving the fiber orientation data (%s).'%save_name)
                np.savetxt(save_path, res, fmt = ['%d','%.3f', '%.3f', '%.3f','%d','%d','%d','%.8f','%d'], delimiter='\t')
        
        print('\n', 'Done.')


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    