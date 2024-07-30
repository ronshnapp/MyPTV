# -*- coding: utf-8 -*-
"""
Created on May 2024

@author: ron

Code for 3D matching of particle silouettes 
"""

from numpy import loadtxt, array, savetxt, hstack, ones, where, argmax, mean 
from scipy.optimize import minimize
from numpy import arctan2
from pandas import read_csv, DataFrame
from random import sample
from itertools import product




class get_silhouettes(object):
    
    def __init__(self, imsys, blob_files, edge_files, traj_file, maxMatchErr, 
                 Np=3, maxCandidatesCheck=5):
        '''
        
        input:
        -----
        
        imsys - A calibrated myptv img_system instance.
        
        blob_files - a list of paths, pointing to the files in which the 
                     particles' blobs are stored.
        
        edge_files - a list of paths, pointing to the files in which the 
                     particles' edges are stored.
        
        traj_file - string, the path to which the particles' trajectories 
                    are stored.
                    
        maxMatchErr - the heighest acceptable matching error, e.g. in mm.
                    
        Np - The maximum number of points to stereo match along the particles'
             silhouettes. 
             
        maxCandidatesCheck - When stereo matching edge pixels along the edges,
                             for each camera we will only consider the 
                             maxCandidatesCheck number of best pixel candidates
        '''
        
        self.imsys = imsys
        self.blob_files = blob_files
        self.edge_files = edge_files
        self.traj_file = traj_file
        self.nCams = len(self.imsys.cameras)
        self.maxMatchErr = maxMatchErr
        self.Np = Np
        self.maxCandidatesCheck = maxCandidatesCheck
        
        # read the trajectory file
        data = read_csv(self.traj_file, header=None, sep='\t')
        timeIndex = data.shape[1] - 1
        self.trajs = dict([(g, k.values) for g,k in data.groupby(timeIndex)])
        self.times = sorted(list(self.trajs.keys()))
        
        
        # read the edge files; After this block we have a dataset with the
        # following structure:
        #   self.edge_pixels[camera number][frame number] -> array with 5 columns:
        #    [traj number, blobNumber, xpixel, ypixel, framenumber]
        self.edge_pixels = {}
        for e, fn in enumerate(self.edge_files):
            edged_data = read_csv(fn, header=None, sep='\t')
            timeIndex = edged_data.shape[1] - 1
            d = dict([(g, hstack([ones((len(k),1))*-1, k.values])) 
                      for g,k in edged_data.groupby(timeIndex)])
            self.edge_pixels[e] = d
            
        for frameNum in self.trajs.keys():
            for tr in self.trajs[frameNum]:
                trajNum = tr[0]
                blobNumColumns = tr[-2-self.nCams:-2]
                for camNum, blobNum in enumerate(blobNumColumns):
                    whr = self.edge_pixels[camNum][frameNum][:,1]==blobNum
                    self.edge_pixels[camNum][frameNum][whr,0] = trajNum
                    
        # a list that will contain the result of finding edges. The format
        # is [traj_id, x, y, z, dist, framenum] where x, y, z are coordinates
        # of a stereomatched point along the silhouette of particle number 
        # traj_id, at time framenum, and dist is the stereo matching error in 
        # e.g. mm.
        self.particleEdges = []
                    
                    
                    
                    
                    
    def get_particle_edge_pixels(self, frameNum, trajNum):
        '''
        Given trajectory and frame numbers, this function returns the 
        pixels that form its silhouette. 
        '''
        particleEdges = {}
        for camNum in range(self.nCams):
            whr = self.edge_pixels[camNum][frameNum][:,0]==trajNum
            particleEdges[camNum] = self.edge_pixels[camNum][frameNum][whr,2:4]
        return particleEdges
    
    
    def stereo_match_particle_edge(self, frameNum, trajNum):
        '''
        Given trajectory and frame numbers, this function attempts to stereo 
        match points along its silhouette. 
        '''
        
        # get the pixels of the edges for all cameras
        pixels = self.get_particle_edge_pixels(frameNum, trajNum)
        
        # form list of cameras
        camNums = list(pixels.keys())
        
        # sort the list accoring to the number of pixels seen by each camera
        pxCount = [len(pixels[cn]) for cn in camNums]
        camNums = sorted(camNums, reverse=True, key=pxCount.__getitem__)
        
        # get the index of the camera with most pixels
        indMax = camNums[0]
        
        # The reference camera is the camera that has the 
        # most pixels. For this camera we choose Np pixels.
        
        # Choosing particles evenly, sorted according to their polar angle in 
        # image space
        center = mean(pixels[indMax], axis=0)
        getAng = lambda point: arctan2(point[1]-center[1], point[0]-center[0])
        sortedPixels = sorted(list(pixels[indMax]), key=getAng)
        pixelsToMatch = sortedPixels[::int(len(sortedPixels)/self.Np)+1]
        
        # --------------------------------------------------------------------
        # Choosing randomly - abandoned this option
        #pixelsToMatch = sample(list(pixels[indMax]), 
        #                       min([len(pixels[indMax]), self.Np]))
        # --------------------------------------------------------------------
        
        # for each pixel chosen from the reference camera, we search for 
        # matches in the other cameras
        edgePoints = []
        for px in pixelsToMatch:
            # 1. We match the chosen pixel with the particle edges along each
            # of the other cameras to locate candiates for the stereo matching
            candidates = dict([(cn, []) for cn in camNums[1:]])
            for cn in camNums[1:]:
                camPixels = pixels[cn]
                for i in range(len(camPixels)):
                    cdic = {indMax: px[::-1], cn: camPixels[i][::-1]}
                    X, C, dist = self.imsys.stereo_match(cdic, 1e10)
                    if dist <= self.maxMatchErr:
                        candidates[cn].append((i, dist))
                    
                    # sorting the candidates of each camera according to the 
                    # epipolar line distance
                    candidates[cn] = sorted(candidates[cn], key=lambda x:x[1])
            
            # Using the best candidates from all cameras we form all the 
            # possible combination of candidates to find the best combination
            a = [candidates[cn][:self.maxCandidatesCheck] for cn in camNums[1:]]
            candidateCombinations = product(*a)
            candidatesAllCameras = []
            for comb in candidateCombinations:
                cdic = dict([(camNums[1:][i], pixels[camNums[1:][i]][comb[i][0]][::-1])
                            for i in range(len(comb))] + [(indMax, px[::-1])])
                
                X, C, dist = self.imsys.stereo_match(cdic, 1e10)
                
                if dist<=self.maxMatchErr:
                    candidatesAllCameras.append((X, dist))
            
            # add the best stereo matched point to a list
            if len(candidatesAllCameras)>0:
                edgePoints.append(min(candidatesAllCameras, key=lambda x: x[1]))
        
        # here we store the results in self.particleEdges
        for ep in edgePoints:
            self.particleEdges.append([trajNum, ep[0][0], ep[0][1], ep[0][2], 
                                       ep[1], frameNum])
        
        return edgePoints
    
    
    
    
    def characterizePoints(self, points):
        '''
        Given a list of points in 3D, this function returns:
            
            1) c -    the mean position of the points
            2) dmaj - the maximum of the distance between c and all points
            3) dmin - the minimum of the distance between c and all points
            4) n -    a unit vector normal to the plane that minimizes the squared
                      distance to all the points.
        '''
        cx, cy, cz = mean(points, axis=0) 
        c = array([cx, cy, cz])
        
        
        # bounds for a and b
        distPC = []
        for p in points:
            distPC.append( sum((p - c)**2)**0.5 )
        dmaj = max(distPC)
        dmin = min(distPC)
        
        n = fitPlaneToPoints(points)
        
        return c, dmaj, dmin, n
    
    
    
    
    
    def save_results(self, fname):
        '''
        Here we save two files:
            1) the list of edgePoints
            2) we find a unit vector normal to the plane of all the edge 
               points, their center and their largest and minimal distances
               of points from the center.
        '''
        # save the points
        savetxt(fname+'_particleEdges', self.particleEdges, fmt='%.3f', 
                delimiter='\t')
        
        # characterize the points
        dic = {}
        for frameNum in self.trajs.keys():
            dic[frameNum] = {}
            
        for pe in self.particleEdges:
            id_, x, y, z, err, fnum = pe
            try:
                dic[fnum][id_].append([x,y,z])
            except:
                dic[fnum][id_] = [[x,y,z]]
                
        characterizedPoints = []
        for fnum in dic.keys():
            for id_ in dic[fnum].keys():
                c, dmaj, dmin, n = self.characterizePoints(dic[fnum][id_])
                characterizedPoints.append([id_, c[0], c[1], c[2], dmaj, dmin,
                                            n[0], n[1], n[2], fnum])
        
        fmt = ['%d', '%.3f', '%.3f', '%.3f', '%.3f', '%.3f',
               '%.5f', '%.5f', '%.5', '%d']
        
        print(characterizedPoints)
        
        savetxt(fname+'_sheetTrajectories', characterizedPoints, fmt=fmt, 
                delimiter='\t')
            
            
            
            
        
        
        
        
        
        
        
        
# ==========================================================
# two helper functions for fitting planes to the edge points

def distancePointPlane(A, B, z0, point):
    '''
    A 2D plane and a point in 3D are defined as
    
        z(x,y) = Ax + By + z0   and   (px, py, pz).
    
    This function returns the distance, d, between the point and the plane.
    '''
    
    # calculate the normal vector
    r = (A**2 + B**2 + 1)**0.5
    n = array([-A/r, -B/r, 1/r])
    
    # calculate the distance
    d = (point[0]*A + point[1]*B + z0 - point[2]) / (n[2] - A*n[0] - B*n[1])
    
    return abs(d)



def fitPlaneToPoints(points):
    '''
    This function searches for the 2D plane whose sum of the distances from a 
    given list of points is the minimum. The plane is defined as 
    
        z(x,y) = Ax + By + z0 .
        
    This function returns the vector normal to this plane 
    
        n = [-A/r, -B/r, 1/r], where r = (A**2 + B**2 + 1)**0.5
    '''
    
    def err(params):
        A, B, z0 = params
        errs = 0
        for p in points:
            errs += distancePointPlane(A, B, z0, p)**2
            
        return errs
        
    params0 = [0.1, 0.1, 0.1]
    sol = minimize(err, params0)
    
    A, B, z0 = sol.x
    r = (A**2 + B**2 + 1)**0.5
    n = array([-A/r, -B/r, 1/r])
    
    return n

# ==========================================================










if __name__ == '__main__':
    
    from myptv.imaging_mod import camera_wrapper, img_system
    import os
    import time
    
    cam_names = ['cam1', 'cam2', 'cam3', 'cam4']
    cam_dir = '/home/ron/Desktop/Research/myptv_tests/sheet_segmentation'
    cam_lst = [camera_wrapper(cn, cam_dir) for cn in cam_names]
    for cam in cam_lst:
        cam.load()
        
    imsys = img_system(cam_lst)
    
    traj_file = '/home/ron/Desktop/Research/myptv_tests/sheet_segmentation/trajectories'
    
    blob_files = [os.path.join(cam_dir, 'blobs_cam1'),
                  os.path.join(cam_dir, 'blobs_cam2'),
                  os.path.join(cam_dir, 'blobs_cam3'),
                  os.path.join(cam_dir, 'blobs_cam4')]
    
    edge_files = [os.path.join(cam_dir, 'blobs_cam1_edges'),
                  os.path.join(cam_dir, 'blobs_cam2_edges'),
                  os.path.join(cam_dir, 'blobs_cam3_edges'),
                  os.path.join(cam_dir, 'blobs_cam4_edges')]
    
    maxMatchErr = 0.2
    Np = 50
    maxCandidatesCheck=5
    
    gs = get_silhouettes(imsys, blob_files, edge_files, traj_file, maxMatchErr,
                         Np=Np, maxCandidatesCheck=maxCandidatesCheck)
    
    t0 = time.time()
    edgePoints = gs.stereo_match_particle_edge(0,1)
    print('time: ', time.time() - t0)
    print('Found points: ', len(edgePoints))
    
    
    import matplotlib.pyplot as plt
    import numpy as np
    
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    
    edgePoints = np.array([ep[0] for ep in edgePoints])
    ax.plot(edgePoints[:,0], edgePoints[:,1], edgePoints[:,2],
            'bo-', ms=5, lw=0.5)
    
    
    theta = np.linspace(0,2*np.pi,num=1000)
    x_ = 5*np.sin(theta) +20
    y_ = 5*np.cos(theta)
    z_ = 0.5*y_ - x_
    ax.plot(x_, y_, z_, 'ok', ms=2)
    
    
    
    fig, ax = plt.subplots(1,3)
    
    # stereo matching errors
    smerrors = [pe[4] for pe in gs.particleEdges]
    ax[0].hist(smerrors, bins='auto')
    ax[0].set_xlabel('stere matching err')
    
    
    # real errors
    groundTruth = np.array([x_, y_, z_])
    errs = [min(np.linalg.norm(np.array(groundTruth).T - np.array(pe[1:4]), axis=1))
            for pe in gs.particleEdges]
    ax[1].hist(errs, bins='auto')
    ax[1].set_xlabel('real err')
    
    
    ax[2].plot(smerrors, errs, 'o')
    ax[2].set_xlabel('stere matching err')
    ax[2].set_ylabel('real err')
    ax[2].set_xlim(0, 2*max(smerrors))
    ax[2].set_ylim(0, 1.05*max(errs))
    
    
    
    
    
    
    
    
    
    