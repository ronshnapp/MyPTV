# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 12:15:21 2019

@author: ron


Matching module for extracting lab coordinate particles from 
segmented particles' image coordinates.

"""


from myptv.utils import line_dist
from math import ceil, floor
from itertools import combinations, product
from numpy import savetxt, array, inf
from numpy.random import uniform
from scipy.spatial import KDTree
from math import isinf

from pandas import read_csv





class matching_with_marching_particles_algorithm(object):
    '''
    This is a class used to perform stereo matching with the 
    "matching particles" algorithm. The algorithm uses the follwing principle:
    
    1) Given an initial guess for a particle position in lab space coordinates, 
       we can project it onto camera space coordinates of our camera system. 
       The results of the projection can be denoted [(x1, y1), ..., (xn, yn)],
       where n is the number of cameras in our imaging system.
       
    2) The segmentation results from the experiment will be called blobs. 
       For each point (xi, yi) from (1), we choose the blob
       that has the smallest distance to it. If the nearest neighbor blob in
       a certain camera was already used up previously, then we take the 2nd
       nearest blob (or the 3rd,...). We denote the nearest neighbor blobs we
       found as [(xNN1, yNN1), ..., (xNNn, yNNn)]. 
       
    3) Takingg the NN blobs, we try to see if the epipolar lines associated
       with them cross each other. If more than a threshold number of lines
       cross, then the point at which they crossed is considdered the position
       of a real particle. Otherwise, the initial guess is deemed bad and
       we dump it.
       
    The thing that is left it to have a method for choosing initial guesses for
    the lab space particles. We are using for that a combination of three 
    methods - 
    
    1) If particles were found in the previous frame (say frame number i-1), 
       then these positions are used as initial guesses for frame i. This
       allows to have contious trajectories already from the matching process
       and really speeds up the process. It is also possible to do the same
       process going backwards, i.e. using particles from frame i+1 as the 
       initial guesses. From this concept comes the name of this algorithm.
      
    2) Random guesses across the measurement volume. This is not so efficient
       but it does work and it allows starting up the marching process.
       
    3) Using the ray traversal algorithm over pairs of cameras to get initial
       guesses. This works quite good, however when there are many blobs in 
       each frame it can take up a lot of time. This also allows starting up 
       the marching process.
    '''
    
    
    def __init__(self, imsys, blob_files, max_d_err, ROI, N0, voxel_size,
                 min_cam_match=3, reverse_eta_zeta=False):
        '''
        inputs 
        
        imsys - An instance of the myptv img_system with loaded cameras inside
        
        blob_files - a list of paths pointing to blob files. The list shoud
                     be orders according to the order of the camera in imsys.
        
        max_d_err - the maximal uncertainty in the matching, namely, the max 
                    distance between any two epipolar line pairs. This is in
                    lab space units, e.g. millimeters.
                    
        ROI (tuple) - Reion of interest in which we are searching for particles 
              in lab space. This is a tuple of length 6 with the following
              format: (xmin, xmax, ymin, ymax, zmin, zmax)
        
        N0 (integer) - The number of random initial points to try and match in 
                       the random initial guess method.
             
        voxel_size (float of None) - Voxel size when getting initial guesses 
                                     with the ray traversal algorithm. If this
                                     is None then the ray traversal search is 
                                     skipped.
        
        min_cam_match (integer) - The minimum number of cameras allowed when 
                        matching a point, e.g. 2 for pairs, 3 for quadruplets,
                        etc...
        '''
        
        print('\ninitializin matcher:\n')
        
        self.imsys = imsys
        self.blob_files = blob_files
        self.max_d_err = max_d_err
        self.min_cam_match = min_cam_match
        self.ROI = ROI
        self.N0 = N0
        self.matches = []
        self.Ncams = len(self.imsys.cameras)
        
        self.voxel_size= voxel_size
        
        # k is the k nearest neighbour blobs out of which we search
        self.max_k = 2
        
        # a dictionary that holds identifires for the blobs that have been 
        # matched, with keys that are frame numbers; identifiers are
        # each a tuple with (cam number, frame number, blob index)
        self.matchedBlobs = {0.0: set([])}
        
        print('loading blob data...')
        # extract the blob data - each blobfile is a dictionay in a list, where 
        # keys are frame numbers and values are the blob data as arrays. 
        self.frames = set([])
        self.blobs = []
        for fn in blob_files:
            bd = read_csv(fn, sep='\t', header=None)
            
            if reverse_eta_zeta:
                ncols = bd.shape[1]
                ind = list(range(ncols))
                ind[0]=1 ; ind[1]=0
                bd = bd[ind]
            
            self.blobs.append(dict([(k,array(g)) for k,g in bd.groupby(5)]))
            self.frames.update(self.blobs[-1].keys())
            
        self.frames = sorted(list(self.frames))
        
        # a dicionary that is used to hold kd trees of blob coordinates
        self.B_ik_trees = {'frame': None}
 
    
 
        
    def match_nearest_blobs(self, x, frame):
        '''
        Given a point in lab space, x, and a frame number, this function 
        finds in each camera the blobs nearest to this point's projection, 
        it stereo matches them, and returns the results. 
        
        Results are returned only if the stereo match is considered
        successfull, which means the max_d_err and min_cam_match tests 
        were successfull.
        '''
        
        # (1) if the KDTrees are setup with the wromg frame, we fix this
        if self.B_ik_trees['frame'] != frame:
            self.generate_blob_trees(frame)
            #for camNum in range(self.Ncams):
            #    self.B_ik_trees[camNum] = KDTree(self.blobs[camNum][frame][:,:2])
        
        # (2) prepare a "coords" dictionary with the items being the NN blobs 
        coords = {}
        matchBlobs = {}
        for camNum in range(self.Ncams):
            if frame not in self.blobs[camNum].keys(): continue
            cam = self.imsys.cameras[camNum]
            projection = cam.projection(x)
            kNN = self.B_ik_trees[camNum].query([projection], k=self.max_k)  
            ind = kNN[1][0]
            dist = kNN[0][0]
            for i in range(len(ind)):
                

                if dist[i]==inf:
                    continue

                identifier = (camNum, frame, ind[i])
                if identifier in self.matchedBlobs[frame]: # this blob has been used
                    continue
                
                else:
                    blob = self.blobs[camNum][frame][ind[i]]
                    coords[camNum] = blob[:2]
                    matchBlobs[camNum] = (blob[:2], ind[i])
                    break
                
        
        
        # (3) perform the stereo matching; If it fails, return None.
        # res = self.imsys.stereo_match(coords, self.max_d_err*self.Ncams)
        res = self.imsys.stereo_match(coords, self.max_d_err, 
                                      strict_match=True)
        if res is None: return None
        else: xNew, pairedCams, err = res 
        
        # if the min_cam_match passes, return the results; else, return None
        if len(pairedCams)<self.min_cam_match: return None
        
        # if the error is too big, return None
        elif err>self.max_d_err: return None
        
        else:
            for camNum in list(matchBlobs.keys()):
                if camNum not in pairedCams:
                    del(matchBlobs[camNum])
                else:
                    self.matchedBlobs[frame].add((camNum, frame, matchBlobs[camNum][1]))
        
        return xNew, matchBlobs, err, frame
        
        
        
    
    def generate_blob_trees(self, frame):
        '''
        Will create new KDTrees of the blob coordinates for a given frame
        number. Blobs that have been used up already do not appear in the 
        trees' dataset.
        '''
            
        for camNum in range(self.Ncams):
            # whr = [i for i in range(self.blobs[camNum][frame].shape[0]) 
            #                              if i not in used_blob_indexes[camNum]]
            # self.B_ik_trees[camNum] = KDTree(self.blobs[camNum][frame][whr,:2])

            if frame not in self.blobs[camNum].keys():
                self.B_ik_trees[camNum] = None
                
            else:
                self.B_ik_trees[camNum] = KDTree(self.blobs[camNum][frame][:,:2])
            
        self.B_ik_trees['frame'] = frame
    
    
    
    
    def stereo_match_frame_with_random_initial_points(self, frame, message=False):
        '''
        This function iterates over the points in the initial point list
        and attempt to match each of then with the nearest neighbout bolb
        projections.
        '''
        
        if frame not in self.matchedBlobs.keys():
            self.matchedBlobs[frame] = set([])
            
        self.generate_blob_trees(frame)
        
        count = 0
        for i in range(self.N0):
            x0 = [uniform(self.ROI[0], self.ROI[1]),
                  uniform(self.ROI[2], self.ROI[3]),
                  uniform(self.ROI[4], self.ROI[5])]
            res = self.match_nearest_blobs(x0, frame)
            if res is not None:
                self.matches.append(res)
                count += 1
                
        if message:
            print('', 'matches using random guesses: %d'%count)
            
        
        
        
    def stereo_match_frame_with_given_points(self, points, frame, message=False):
        '''
        Given a list of points in lab space, this function iterates over them
        and attepms to match each of them with their nearest neighbout blob
        projections.
        '''
        
        if frame not in self.matchedBlobs.keys():
            self.matchedBlobs[frame] = set([])
            
        self.generate_blob_trees(frame)
        
        count = 0
        for x0 in points:
            res = self.match_nearest_blobs(x0, frame)
            if res is not None:
                self.matches.append(res)
                count += 1
                
        if message:
            print('', 'matches using given points: %d'%count)
        
        
        
        
    def stereo_match_frame_with_previous_matches(self, frame, message=False,
                                                 backwards=False):
        '''
        This function iterates over the points that were found in frame i-1
        (if any exist), and attemps to stereo match their nearest 
        neighbout blob projections in the given frame (i).
        
        if backwards==True, then we try to use particles from i+1 (is exist)
        to find matches at frame i.
        '''
        
        if frame not in self.matchedBlobs.keys():
            self.matchedBlobs[frame] = set([])
            
        self.generate_blob_trees(frame)
        
        if backwards==False:
            pointToMatchOn = [m[0] for m in self.matches if m[3]==frame-1]
            
        elif backwards==True:
            pointToMatchOn = [m[0] for m in self.matches if m[3]==frame+1]
        
        count = 0
        for x0 in pointToMatchOn:
            res = self.match_nearest_blobs(x0, frame)
            if res is not None:
                self.matches.append(res)
                count += 1
                
        if message:
            print('matches using prvious results: %d'%count)
        
        
        
        
    def find_candidates_with_two_cameras(self, camNum1, camNum2, frame):
        '''
        This function uses epipolar voxel traversal search to stereo match 
        blob pairs in two given cameras. The points that are found are then
        returned.
        '''
        # ensure the cameras have blobs in this frame
        if frame not in list(self.blobs[camNum1].keys()): return []
        if frame not in list(self.blobs[camNum2].keys()): return []
        
        # fetching the cameras and the blobs
        cam1 = self.imsys.cameras[camNum1] ; cam2 = self.imsys.cameras[camNum2]
        # O1 = cam1.O ; O2 = cam2.O
        blobs1 = self.blobs[camNum1][frame]
        blobs2 = self.blobs[camNum2][frame]
        
        # getting the center point of the ROI (for epipolar line seach)
        O_ROI = [(self.ROI[1]+self.ROI[0])/2, 
                 (self.ROI[3]+self.ROI[2])/2, 
                 (self.ROI[5]+self.ROI[4])/2]
        
        # getting the size of the ROI diagonal
        a_range = sum([(self.ROI[2*i+1]-self.ROI[2*i])**2 for i in range(3)])**0.5
        da = self.voxel_size/4.0
        
        # the number of voxel in each direction
        nx = int((self.ROI[1]-self.ROI[0])/self.voxel_size)+1
        ny = int((self.ROI[3]-self.ROI[2])/self.voxel_size)+1
        nz = int((self.ROI[5]-self.ROI[4])/self.voxel_size)+1
        
        
        # a dicionary that holds the blob numbers traversed in each voxel
        voxel_dic = {}
        
        # listing the traversed volxels for blobs in camera 2
        for e,b in enumerate(blobs2):
            
            identifier = (camNum2, frame, e)
            if identifier in self.matchedBlobs[frame]: # this blob has been used
                continue
            
            #r = cam2.get_r(b[0], b[1])
            O2, r = cam2.get_epipolarline(b[0], b[1])
            a_center = sum([r[i]*(O_ROI[i]-O2[i]) for i in range(3)]) 
            a1, a2 = a_center - a_range/2 , a_center + a_range/2 
            
            # traversing the blob from O2+a1*r to O2+a2*r to list the voxels
            blob_voxels = set([])
            a = a1
            while a<=a2:
                x, y, z = O2[0] + r[0]*a, O2[1] + r[1]*a, O2[2] + r[2]*a
                
                if self.ROI[0] < x < self.ROI[1]:
                    if self.ROI[2] < y < self.ROI[3]:
                        if self.ROI[4] < z < self.ROI[5]:
                            i = int((x-self.ROI[0])/(self.ROI[1]-self.ROI[0])*nx)
                            j = int((y-self.ROI[2])/(self.ROI[3]-self.ROI[2])*ny)
                            k = int((z-self.ROI[4])/(self.ROI[5]-self.ROI[4])*nz)
                            blob_voxels.add((i,j,k))
                a += da
            
            for voxel in blob_voxels:
                try:
                    voxel_dic[voxel].append(e)
                except:
                    voxel_dic[voxel] = [e]
            
        candidate_pairs = set([])
        # traversing the blobs in camera1 to obtain candidates pairs
        for e,b in enumerate(blobs1):
            
            identifier = (camNum1, frame, e)
            if identifier in self.matchedBlobs[frame]: # this blob has been used
                continue
            
            # r = cam1.get_r(b[0], b[1])
            O1, r = cam1.get_epipolarline(b[0], b[1])
            a_center = sum([r[i]*(O_ROI[i]-O1[i]) for i in range(3)]) 
            a1, a2 = a_center - a_range/2 , a_center + a_range/2 
            
            # traversing the blob from O2+a1*r to O2+a2*r to list the candidates
            a = a1
            while a<=a2:
                x, y, z = O1[0] + r[0]*a, O1[1] + r[1]*a, O1[2] + r[2]*a
                
                if self.ROI[0] < x < self.ROI[1]:
                    if self.ROI[2] < y < self.ROI[3]:
                        if self.ROI[4] < z < self.ROI[5]:
                            i = int((x-self.ROI[0])/(self.ROI[1]-self.ROI[0])*nx)
                            j = int((y-self.ROI[2])/(self.ROI[3]-self.ROI[2])*ny)
                            k = int((z-self.ROI[4])/(self.ROI[5]-self.ROI[4])*nz)
                            try:
                                candidates = voxel_dic[(i,j,k)]
                                for cnd in candidates:
                                    candidate_pairs.add((e, cnd))
                            except:
                                pass
                a += da
        
        
        candidate_points = []
        for e1, e2 in candidate_pairs:
            coords = {camNum1:blobs1[e1], camNum2:blobs2[e2]}
            # stereoMatch = self.imsys.stereo_match(coords, self.max_d_err)
            stereoMatch = self.imsys.stereo_match(coords, self.max_d_err, 
                                                  strict_match=True)
            if stereoMatch is not None:
                candidate_points.append(stereoMatch[0])
        
        return candidate_points
        
        
        
        
        
        
    def match_frame(self, frame, backwards=False):
        '''
        Will stereo match particles in the given frame number.
        '''
        print('\n')
        print('frame: %d'%frame)
        
        if frame not in self.matchedBlobs.keys():
            self.matchedBlobs[frame] = set([])
        
        m0 = len(self.matches)
        self.stereo_match_frame_with_previous_matches(frame, backwards=backwards)
        newPrevFrame = len(self.matches) - m0
        
        # Matching with random points
        self.stereo_match_frame_with_random_initial_points(frame)
        
        # matching with pair candidates
        if self.voxel_size is not None:
            for i in range(self.Ncams-1):
                camNum1=i ; camNum2=i+1
                cands = self.find_candidates_with_two_cameras(camNum1, 
                                                              camNum2, 
                                                              frame)
                self.stereo_match_frame_with_given_points(cands, frame)
        
        newEpipolarCands = len(self.matches) - m0 - newPrevFrame
        
        newTot = newEpipolarCands + newPrevFrame
        prnt = (newTot, newPrevFrame, newEpipolarCands)
        print('Found %d matches: %d from prev. frame + %d new'%prnt)
        
        
        
       
        
    def plot_disparity_map(self, camNum):
        '''
        A disparity map is a graph that hels in assessing the uncertainty
        of the calibration. It is done by taking the results of the stereo
        matching, projecting them back on a cameras' sensor, and plotting the
        difference between the projected point and the blob with which this 
        point was matched. The disparity map is drawn separately for each
        cameras, and hence the camNum parameter.
        '''
        import matplotlib.pyplot as plt
        disparities_x = []
        disparities_y = []
        cam = self.imsys.cameras[camNum]
        for m in self.matches:
            
            try: blob = m[1][camNum][0]
            except: continue  # < cameras did not participate in this match
                
            x = m[0]
            proj = cam.projection(x)
            disparities_x.append(blob[0] - proj[0])
            disparities_y.append(blob[1] - proj[1])
        
        fig, ax = plt.subplots()
        ax.hist2d(disparities_x, disparities_y, bins=20)
        ax.scatter(disparities_x, disparities_y, s=2)
        ax.set_aspect('equal')
        return fig, ax
        
        
        
        
    def save_particles(self, saveName):
        '''
        Will save the matched particles on the disk with the given
        file name.
        '''
        
        toSave = []
        for m in self.matches:
            p = [m[0][0], m[0][1], m[0][2]]
            for cn in range(self.Ncams):
                try: 
                    p.append(m[1][cn][1])
                except:
                    p.append(-1)
            p.append(m[2])
            p.append(m[3])
            toSave.append(p)
        
        fmt = ['%.3f', '%.3f', '%.3f']
        for i in range(self.Ncams):
            fmt.append('%d')
        fmt = fmt + ['%.3f', '%.3f']
        savetxt(saveName, toSave, fmt=fmt, delimiter='\t')
        
        
        

































# =============================================================================
# =============================================================================
#    Legacy: classes for stereo matching using the Ray travesal algorithm
# =============================================================================
# =============================================================================



'''
Below there are classes that use two matching algorithms in conjunction:
    
1) the Ray Traversal algorithm reported 
by Bourgion and Huisman, 2020 (https://arxiv.org/pdf/2003.12135.pdf).

2) a novel algorithm which prioritizes the matching of blobs who are trackable 
   in the 2D images.
   
The two algorithms are stitched to work together in the match_blob_files()
class, first running the trackable blobs search and then the Ray Traversal
algorithm.
'''






class match_blob_files_Ray_Traversal(object):
    '''A class for obtaining triangulated particles positions from a 
    list of segmented blobs. Use self.get_particles() and after that,
    the particles found are stored in the attribute self.particles.'''
    
    
    def __init__(self, blob_fnames, img_system, RIO, voxel_size, max_blob_dist,
                 max_err=1e9, reverse_eta_zeta = False, 
                 travered_voxel_rep = True):
        '''
        blob_fname - a list of the file names containing the segmented blob
                     data. The list has to be sorted according the order of
                     cameras in the img_system.
                     
        img_system - an instance of the img_system class with the calibrated
                     cameras.
                     
        RIO - A nested list of 3X2 elements. The first holds the minimum and 
              maximum values of x coordinates, the second is same for y, and 
              the third for z coordinates. 
        
        voxel size - the side length of voxel cubes used in the ray traversal
                     algorithm. Given in lab coordinate scales (e.g. mm).
                     
        max_blob_dist - the largest distance for which blobs are concidered 
                        neighbours in the image space coordinates (namely, the
                        highest permissible particle displacement in pixels).
                     
        max_err - maximum acceptable uncertainty in particle position. If None,
                  (defult), than no bound is used.
                     
        reverse_eta_zeta - Should be false if the eta and zeta coordinates 
                           need to be in reverse order so as to match the
                           calibration. This may be needed if the calibration 
                           data points were given where the x, y coordinates
                           are transposed (as happens, e.g., if using 
                           matplotlib.pyplot.imshow).
                           
        travered_voxel_rep - If true, this will repeat the segmentation with 
                             translated voxels to make sure no blobs were 
                             missed due to an alliasing.
        '''
        self.blobs = []
        for fn in blob_fnames:
            #self.blobs.append(loadtxt(fn))
            self.blobs.append(array(read_csv(fn, sep='\t', header=None)))
                     
        self.imsys = img_system
        self.RIO = RIO
        self.voxel_size = voxel_size
        self.reverse_eta_zeta = reverse_eta_zeta
        self.max_blob_dist = max_blob_dist
        self.max_err = max_err
        
        time_lst = []
        for bl in self.blobs:
            for b in bl:
                time_lst.append(b[-1])
        self.time_lst = sorted(list(set(time_lst)))
        self.cam_names = [cam.name for cam in self.imsys.cameras]
        self.travered_voxel_rep = travered_voxel_rep
        
    def get_particles_dic(self, frame):
        '''
        returns a particles dictionary (with camera names as keys and blob
        values) for particles in the given frame.
        '''
        pd = {}
        if self.reverse_eta_zeta:
            for i in range(len(self.blobs)):
                cn = self.cam_names[i]
                arr = self.blobs[i][self.blobs[i][:,-1] == frame][:,1::-1]
                pd[cn] = arr.tolist()
        
        else:
            for i in range(len(self.blobs)):
                cn = self.cam_names[i]
                arr = self.blobs[i][self.blobs[i][:,-1] == frame][:,:2]
                pd[cn] = arr.tolist()
                
        return pd
        
        
        
    def get_particles(self, frames=None):
        '''
        Use this to match blobs into particlesin 3D.
        
        input - 
        frames - if None, this will match particles at all times. If a list of
                 integers, this will only perform the matching on particles at
                 the given frames.
        '''
        # set up a list of the frames in which particles are matched
        if frames is None:
            frames = self.time_lst
        
        
        # start matching, one frame at a time
        self.particles = []
        previous_particles = []
        print('')
        
        for e,tm in enumerate(frames):
            
            countParticlesInThisFrame = 0
            
            # set up a blobs dictionary with camera names as key
            pd = self.get_particles_dic(tm)
            nb = sum([len(pd[k]) for k in pd.keys()]) /len(pd.keys())
            #print('', end='\r')
            
            
            # (1) for iterations after the first run, use the time
            # augmented matching
            useTimeMatching = True
            if e>0 and useTimeMatching:
                mut = matching_using_time_Ray_Traversal(self.imsys, pd, previous_particles,
                                          max_err = self.max_err)
                #return mut  # <-- used for checks
                mut.triangulate_candidates()
                for p in mut.matched_particles:
                    self.particles.append(p + [tm])
                    countParticlesInThisFrame += 1
                
                # update the particles dictionary
                pd = mut.return_updated_particle_dict()
                
                
# =============================================================================
#            The initiation of matching with time has a bug. For now
#            We comment it out, and will fix this in the future.
#
#             # for the first iteration, initiate search using the neighbouring
#             # blobs paradigm
#             if e==0 and len(frames)>1:
#                 pd1 = self.get_particles_dic(frames[e+1])
#                 itm = initiate_time_matching_Ray_Traversal(self.imsys, pd, pd1, 
#                                              self.max_blob_dist, self.RIO, 
#                                              self.voxel_size, 
#                                              max_err = self.max_err)
#                 itm.choose_blobs_with_neighbours()
#                 itm.match_blobs_with_neighbours()
#                 for p in itm.matched_particles:
#                     self.particles.append(p + [tm])
#                 pd = itm.return_updated_particle_dict()
# =============================================================================
                                
            # (2) match particles using the voxel method
            M = matching_Ray_Traversal(self.imsys, pd, self.RIO, self.voxel_size,
                         max_err = self.max_err)
            #return M  # <-- used for checks
            M.get_voxel_dictionary()
            M.list_candidates()
            M.get_particles()
            
            # extract the matched particles to the list self.particles 
            for p in M.matched_particles:
                self.particles.append(p + [tm])
                countParticlesInThisFrame += 1
            
            
            
            # (3) match remaining particles using traversed voxels
            if self.travered_voxel_rep:
                dv = self.voxel_size/2.0
                new_ROI = tuple([(self.RIO[i][0]+dv,self.RIO[i][1]-dv) 
                                  for i in range(3)])
                # updating the particle dictionary 
                new_pd = pd.copy()
                for p in M.matched_particles:
                    blob_info = p[3] 
                    for ci, (rayNum, xy) in blob_info:
                        cn = self.cam_names[ci]
                        new_pd[cn][rayNum] = [-1,-1]
                for k in new_pd.keys():
                    i=0
                    while i<len(new_pd[k]):
                        if new_pd[k][i]==[-1,-1]:
                            del new_pd[k][i]
                        else:
                            i+=1
    
                M2 = matching_Ray_Traversal(self.imsys, new_pd, new_ROI, self.voxel_size,
                              max_err = self.max_err)
                # return M2  # <-- used for checks
                M2.get_voxel_dictionary()
                M2.list_candidates()
                M2.get_particles()
                
                # extract the matched particles to the list self.particles 
                for p in M2.matched_particles:
                    self.particles.append(p + [tm])
                    countParticlesInThisFrame += 1
            
            
                
            # (4) list the particles found in this frame from the next  
            # iteration of the time matching, and print statistics
            previous_particles = self.particles[-countParticlesInThisFrame:]
            
            c4 = sum([1 for p in previous_particles if len(p[3])==4])
            c3 = sum([1 for p in previous_particles if len(p[3])==3])
            c2 = sum([1 for p in previous_particles if len(p[3])==2])
            pc = ' quads. trips. pairs. = (%d, %d, %d)'%(c4, c3, c2)
            print(' frame: %d ; %.1f blobs/cam ;'%(tm, nb) + pc)

            
        
        # filter particles with large triangulation error
        self.particles = list(filter(lambda p: p[4]<self.max_err,
                                     self.particles))
        print('\n')
        print('done!\n')                        
        errors = [p[4] for p in self.particles]
        print('mean error: %.3f'%(sum(errors)/len(errors)))
        Np = len(self.particles)
        times = [p[-1] for p in self.particles]
        Nframes = len(set(times))
        print('avg. particles in frame: %.2f'%(Np/Nframes))

    
    
    def save_results(self, fname):
        '''will save the list of particles obtained'''
        particles_to_save = []
        Ncams = len(self.imsys.cameras)
        
        for p in self.particles:
            rd = dict(p[3])
            
            p_ = [p[0], p[1], p[2]]
            
            for i in range(Ncams):
                if i in list(rd.keys()):
                    p_.append(rd[i][0])
                else:
                    p_.append(-1)
            p_.append(p[4])
            p_.append(p[5])
            particles_to_save.append(p_)
            
        fmt = ['%.3f', '%.3f', '%.3f']
        for i in range(Ncams):
            fmt.append('%d')
        fmt = fmt + ['%.3f', '%.3f']
        savetxt(fname, particles_to_save, fmt=fmt, delimiter='\t')
        
            
                










class matching_Ray_Traversal(object):
    '''A class for matching particles in images taken simultaneously
    from different cameras.
    
    The relevant functions for use are: 
        1) self.get_voxel_dictionary()
        2) self.list_candidates()
        3) self.get_particles()
    After running these three functions the attribute self.matched_particles
    holds the results of triangulation.
    '''
    
    
    def __init__(self, img_system, particles_dic, 
                 RIO, voxel_size, max_err=None):
        '''
        img_system - is an instance of the img_system object with camera 
                     objects. 
                     
        particles_dic - A dictionary, keys are camera names, and values
                     are lists of particle coordinates segmented in each of 
                     the cameras.
                     
        RIO - A nested list of 3X2 elements. The first holds the minimum and 
              maximum values of x coordinates, the second is same for y, and 
              the third for z coordinates. 
        
        voxel size - the side length of voxel cubes used in the ray traversal
                     algorithm. Given in lab coordinate scales (e.g. mm).
        max_err - maximum allowable triangulation rms error.
        
        '''
        
        self.imsys = img_system
        
        # make a list of rays with: 
        # [camera x coord, camera y coord, (cam number, particle number)]
        self.rays = []
        self.ray_camera_indexes = [0]
        for i in range(len(self.imsys.cameras)):
            cam  = self.imsys.cameras[i]
            particles_i = particles_dic[cam.name]
            self.ray_camera_indexes.append(len(particles_i) + 
                                           self.ray_camera_indexes[-1])
            for j in range(len(particles_i)):
                x, y = particles_i[j][0], particles_i[j][1]
                r_ij = cam.get_r(x, y)
                self.rays.append( (x, y, (i,j), r_ij) )
        
        self.RIO = RIO
        self.voxel_size = voxel_size
        self.max_err = max_err
        
        # set up lists of voxel centers:
            
        Nx = ceil((RIO[0][1]-RIO[0][0])/voxel_size)
        cx = (RIO[0][1]+RIO[0][0])/2.
        if Nx%2==0: f = (floor(Nx/2)-0.5) * voxel_size
        elif Nx%2!=0: f = floor(Nx/2)*voxel_size
        self.Nx = Nx
        self.x = [i*voxel_size + cx - f for i in range(Nx)]
        
        Ny = ceil((RIO[1][1]-RIO[1][0])/voxel_size)
        cy = (RIO[1][1]+RIO[1][0])/2.
        if Ny%2==0: f = (floor(Ny/2)-0.5) * voxel_size
        elif Ny%2!=0: f = floor(Ny/2)*voxel_size
        self.Ny = Ny
        self.y = [i*voxel_size + cy - f for i in range(Ny)]
        
        Nz = ceil((RIO[2][1]-RIO[2][0])/voxel_size)
        cz = (RIO[2][1]+RIO[2][0])/2.
        if Nz%2==0: f = (floor(Nz/2)-0.5) * voxel_size
        elif Nz%2!=0: f = floor(Nz/2)*voxel_size
        self.Nz = Nz
        self.z = [i*voxel_size + cz - f for i in range(Nz)]



    def ray_traversed_voxels(self, ray):
        '''
        Given a ray, this will add the voxel through which it
        traverses into the list self.traversed_voxels .
        '''
        if [ray[0], ray[1]] == [-1,-1]:
            return
        
        cam  = self.imsys.cameras[ray[2][0]]
        # O = cam.O
        # r = cam.get_r(ray[0], ray[1])
        O, r = cam.get_epipolarline(ray[0], ray[1])
        r_ = r / sum(r**2)**0.5

        a1, a2 = (self.RIO[2][0] - O[2])/r_[2], (self.RIO[2][1] - O[2])/r_[2]
        if a2<a1:
            a2, a1 = a1, a2
        
        da = self.voxel_size/4.0
        
        ray_voxels = set([])
        a = a1
        while a<=a2:
            x_, y_, z_ = O[0] + r_[0]*a, O[1] + r_[1]*a, O[2] + r_[2]*a
            if x_>self.RIO[0][1] or x_<self.RIO[0][0]: 
                a += self.voxel_size/4.0
                continue
            
            if y_>self.RIO[1][1] or y_<self.RIO[1][0]: 
                a += self.voxel_size/4.0
                continue
            
            i = int((x_ - self.x[0] - self.voxel_size/2)/self.voxel_size +1)
            j = int((y_ - self.y[0] - self.voxel_size/2)/self.voxel_size +1)
            k = int((z_ - self.z[0] - self.voxel_size/2)/self.voxel_size +1)
            ray_voxels.add(((i, j, k), ray[2]))
            
            #--------------------------------------------------------------
            # An attempt to fix voxel aliasing that did not improve results 
            # di = x_ - self.x[i]
            # dj = y_ - self.y[j]
            # i_ = i + int(di / abs(di))
            # j_ = j + int(dj / abs(dj))
            # ray_voxels.add(((i_, j, k), ray[2]))
            # ray_voxels.add(((i, j_, k), ray[2]))
            # ray_voxels.add(((i_, j_, k), ray[2]))
            #--------------------------------------------------------------
            
            a += da
        
        self.traversed_voxels += list(ray_voxels)


    
    
    def get_voxel_dictionary(self):
        '''This generates a dictionary who's keys are voxel indexes and
        who's values are the rays that passed through this voxel.'''

        self.traversed_voxels = []
        for ray in self.rays:
            self.ray_traversed_voxels(ray)

        
        voxel_dic = {}
        for vxl in self.traversed_voxels:
            try:
                voxel_dic[vxl[0]].add(vxl[1])
            except:
                voxel_dic[vxl[0]] = set([vxl[1]])
                
        self.voxel_dic = voxel_dic
        
    
    def list_candidates(self):
        '''This will make lists of possible candidate rays for
        triangulation, separated for pairs, triplets, quadruplets, etc.
        
        Candidates are based on the voxels of voxel_dic, while calculating the 
        RMS and the maximum distance between the estimated particle location 
        and the epipolar lines.'''
        
        self.candidate_dic = {}
        group_sizes = range(2, len(self.imsys.cameras)+1)
        for i in group_sizes:
            self.candidate_dic[i] = []
        
        for k in self.voxel_dic.keys():
            if len(self.voxel_dic[k]) >= 2:
                # make a nested list of the rays, by their camera number 
                ray_by_cams = [[] for i in range(len(self.imsys.cameras))]
                for ray in self.voxel_dic[k]:
                    ray_by_cams[ray[0]].append(ray)
                
                # find all possible combinations of the rays for various
                # numbers of cameras
                for gs in group_sizes:
                    for comb in combinations(ray_by_cams, gs):
                        self.candidate_dic[gs] += product(*comb)
        
        for k in self.candidate_dic.keys():
            self.candidate_dic[k] = list(set(self.candidate_dic[k]))
        


    def triangulate_rays(self, rays):
        '''will return the results of stereo matching of a list of rays'''
        
        dc = {}
        cams = []
        for ray in rays:
            i = self.ray_camera_indexes[ray[0]] 
            ip1 = self.ray_camera_indexes[ray[0]+1]
            ri = self.rays[i:ip1][ray[1]][3]
            Oi = self.imsys.cameras[self.rays[i:ip1][ray[1]][2][0]].O
            dc[ray[0]] = Oi, ri
            cams.append(ray[0])
        
        n = len(rays)
        d, x = [], []
        
        for i in range(n):
            Oi, ri = dc[cams[i]]
            for j in range(i+1, n):
                Oj, rj = dc[cams[j]]
                D, x_ij = line_dist(Oi, ri, Oj, rj)
                #d.append(D)
                #x.append(x_ij)
                
                if D<4*self.max_err:
                    d.append(D)
                    x.append(x_ij)
                else:
                    return x_ij, cams, 1e9
        
        return sum(x)/1.0/len(x), cams, sum(d)/1.0/len(x)

    
    
    def is_used(self, cand):
        '''
        Returns True if al least one of the rays in the candidate had been
        used up, and False if all of the rays have not been used yet.
        '''
        for ray in cand:
            if ray in self.used_rays:
                return True
        return False

    
    def get_particles(self):
        '''Once all candidates are found, this function chooses the "best"
        matches and returns them. The reliability of the matches is considered 
        higher if
        1) they have higher number of cameras participating in the
        triangulation
        2) the RMS of distance between the crossing point and the epipolar 
        lines is smaller
        in this order. Thus, we choose the combinations of rays with highest
        number of camera participating and with the smallest RNS triangulation 
        error.
        '''
        matched_particles = []
        self.used_rays = set([])
        
        count = 0
        for k in sorted(self.candidate_dic.keys(), reverse=True):
            
            cand_k = self.candidate_dic[k]
            
            # get rid of candidates with rays that were used
            if count>0:
                cand_k = list(filter(lambda c: not(self.is_used(c)), cand_k))
            
            # triangulate all the candidate rays
            ray_crosses = [self.triangulate_rays(cand) for cand in cand_k]                
            
            # zip and sort candidates by RMS error
            dist_sorted_cands = sorted(zip(cand_k, ray_crosses),
                                       key = lambda x: x[1][2])

            
            # if the rays in the candidate have not been used, adde the
            # particle to the results list 
            for i in range(len(dist_sorted_cands)):
                used_check = self.is_used(dist_sorted_cands[i][0])                
                if not used_check:
                    p = dist_sorted_cands[i][1]
                    r_list = [(ri[0], (ri[1], self.get_eta_zeta(ri))) 
                              for ri in dist_sorted_cands[i][0]]
                    
                    new_p = [round(p[0][0], ndigits=3), 
                             round(p[0][1], ndigits=3),
                             round(p[0][2], ndigits=3),
                             r_list,
                             round(p[-1], ndigits=3)]
                    matched_particles.append(new_p)
                    self.used_rays.update(dist_sorted_cands[i][0])
            count += 1     
                 
        self.matched_particles = matched_particles
        
        
        
        
    def get_eta_zeta(self, ray):
        '''
        Retrieves the image space coordinates, eta and zeta from the 
        list self.rays.
        '''
        i = self.ray_camera_indexes[ray[0]] 
        ip1 = self.ray_camera_indexes[ray[0]+1]
        return self.rays[i:ip1][ray[1]][0], self.rays[i:ip1][ray[1]][1]

        

    def plot_ray_epipolar_lines(self, ray, zlims, ax):
        '''will plot a ray's epipolar line for a given 3D axis.'''
        import matplotlib.pyplot as plt
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        camera = self.imsys.cameras[ray[0]]
        i = self.ray_camera_indexes[ray[0]] 
        ip1 = self.ray_camera_indexes[ray[0]+1]
        eta, zeta = self.rays[i:ip1][ray[1]][:2]
        camera.plot_3D_epipolar_line(eta, zeta, zlims,
                                     ax=ax, color=colors[ray[0]])
        
        






class matching_using_time_Ray_Traversal(object):
    '''
    An implementation of a novel algorithm to improve the matching
    process using temporal information.
    '''
    
    def __init__(self, img_system, particles_dic,
                 previously_used_blobs, max_err=1e9):
        '''
        An implementation of a novel algorithm to improve the matching
        process using temporal information.
        
        inputs - 
        img_system - is an instance of the img_system object with camera 
                     objects. 
                     
        particles_dic - A dictionary with keys that are camera names, and 
                        values are lists of particle coordinates segmented in 
                        each of the cameras.
                     
        previously_used_blobs - A list. Each item in the list is a particle 
                                that was matched successfully in the previous 
                                frame.
                                
        max_err - the maximum allowable triangulation error.
        
        '''
        self.imsys = img_system
        self.pd = particles_dic
        self.prev_used_blobs = previously_used_blobs
        self.max_err = max_err
        
        # we form KDTrees for the nearest neighbour blobs search        
        #self.trees = [KDTree(self.pd[k]) for k in self.pd.keys()]
        
        
        # define an object that can take a querry and return an empty list
        class deadTree():
            
            def __init__(self):
                return None
            
            def query(self, p):
                return [None, None]
        
        
        self.trees = []
        for k in self.pd.keys():
            if len(self.pd[k])>0:
                self.trees.append(KDTree(self.pd[k]))
            else:
                self.trees.append(deadTree())
        
    
    def triangulate_candidates(self):
        '''
        Will do the triangulation of blobs that are nearest to those that 
        were successfully matched in the previous frame.
        
        1) for each successfully matched particle in the previous frame, 
           we find the nearest neighbouring blobs in the current frame
           
        2) we triangulate the blobs
        
        3) if the triangulation error is lower than the threshold max_err,
           we add the particle to a list self.matched particles.
        '''
        
        triangulated_particles = []
        for p in self.prev_used_blobs:
            
            # first, find the nearest neighbouring blobs
            p_blobs = p[3]
            nearest_blobs_num = {}
            nearest_blobs_coords = {}
            
            # for the blobs that participated in this particle
            for blb in p_blobs:
                ci,(rn, (x, y)) = blb
                cn = self.imsys.cameras[ci].name
                bn = self.trees[ci].query((x, y))[1]
                if bn is not None:
                    nearest_blobs_num[ci] = bn
                    nearest_blobs_coords[ci] = self.pd[cn][bn]
             
            # for cameras that weren't used in this particle we project
            # the particle on them and locate the nearest blob in the next 
            # frame. If the neighbour is close enough in the triangulation 
            # we will use it.
            camNum = len(self.imsys.cameras)
            usedCams = list(nearest_blobs_coords.keys())
            unusedCams = [i for i in range(camNum) if (i not in usedCams)]
            p_xyz = p[:3] 
            for ci in unusedCams:
                cn = self.imsys.cameras[ci].name
                projectionI = self.imsys.cameras[ci].projection(p_xyz)
                bn = self.trees[ci].query(projectionI)[1]
                if bn is not None:
                    nearest_blobs_num[ci] = bn
                    nearest_blobs_coords[ci] = self.pd[cn][bn]
            
                
            # second, triangulate the nearest blob neighbours:
            triangulated = self.imsys.stereo_match(nearest_blobs_coords, 1e9)
            
        
            # third, if the RMS triangulation error is low enough, add the 
            # triangulation to a list of matched particles
            if triangulated[2] < self.max_err:
                r = [(ci, 
                      (nearest_blobs_num[ci],tuple(nearest_blobs_coords[ci]))) 
                     for ci in nearest_blobs_num.keys()]
                
                p = triangulated[0]
                new_p = []
                
                new_p = [round(p[0], ndigits=3), 
                         round(p[1], ndigits=3),
                         round(p[2], ndigits=3),
                         r,
                         round(triangulated[2], ndigits=3)]
                
                triangulated_particles.append(new_p)
        
        
        # To finish off, make sure we're not using a blob more than once
        self.matched_particles = []
        key = lambda tr: tr[-1]
        triangulated_particles = sorted(triangulated_particles, key=key)
        self.used_blobs = set([])
        for p in triangulated_particles:
            test = [blb not in self.used_blobs for blb in p[3]]
            if all(test):
                self.matched_particles.append(p)
                for blob in p[3]:
                    self.used_blobs.add(blob)
        
        
        
    def return_updated_particle_dict(self):
        '''
        After finding matched particles (self.triangulate_candidates),  
        this will return an updated copy of particle_dictionary, that does
        not contain the used blobs.
        '''
        new_pd = self.pd.copy()
        for e,p in enumerate(self.matched_particles):
            p_blobs = p[3]
            for blb in p_blobs:
                ci,(rn, xy) = blb
                cn = self.imsys.cameras[ci].name
                #new_pd[cn].remove(list(xy))
                ind = new_pd[cn].index(list(xy))
                new_pd[cn][ind] = [-1,-1]
        return new_pd
                
        
        






class initiate_time_matching_Ray_Traversal(object):
    '''
    A class used in the time matching algorithm to initiate the first
    frame. This class will search for blobs that have a nearest neighbour in 
    the next frame lower than a given threshold and will first stereo-match
    only these particles.
    '''
    
    def __init__(self, img_system, particles_dic_0, particles_dic_1,
                 max_distance, RIO, voxel_size, max_err=1e9):
        '''
        input - 
        
        img_system - an instance of the imaging system class
        
        particles_dic_0 - A dictionary; keys are camera names, and values 
                         are lists of particle coordinates segmented in 
                         each of the cameras at the first frame, t=0.
                         
        particles_dic_1 - A dictionary; keys are camera names, and values 
                         are lists of particle coordinates segmented in 
                         each of the cameras at the second frame, t=0+dt.
                         
        max_distance - The maximum alowable distance between blobs to be 
                       considered neighbour. This is in image space 
                       coordinates (how many pixels blobs move in the 2D
                       images?).
        
        RIO - A nested list of 3X2 elements. The first holds the minimum and 
              maximum values of x coordinates, the second is same for y, and 
              the third for z coordinates. 
        
        voxel size - the side length of voxel cubes used in the ray traversal
                     algorithm. Given in lab coordinate scales (e.g. mm).
        
        max_err - maximum allowable RMS triangulation error.
        '''
        self.imsys = img_system
        self.pd = particles_dic_0
        self.pd1 = particles_dic_1
        self.max_dist = max_distance
        self.max_err = max_err
        self.RIO = RIO
        self.voxel_size = voxel_size
        # we form KDTrees for the nearest neighbour blobs search        
        self.trees = {}
        for k in self.pd.keys():
            self.trees[k] = KDTree(self.pd1[k])
        
    
    
    def has_neighbour(self, blob, cam):
        '''
        Return True if the given blob has a nearest neighbour and False if not.
        '''
        x,y = blob
        dist, dump = self.trees[cam].query((x, y))
        
        if dist < self.max_dist:
            return True
        
        return False
    
    
    
    def choose_blobs_with_neighbours(self):
        '''
        Will go over the blobs in particles_dic_0; if a blob has valid
        neighbours in particles_dic_1, it is added to a new
        particles_dictionary.
        '''
        
        # form the new dictionary
        self.new_pd = {}
        
        # add blobs with neighbours
        for k in self.pd.keys():
            self.new_pd[k] = []
            for blb in self.pd[k]:
                if self.has_neighbour(blb, k):
                    self.new_pd[k].append(blb)
                    
                    
    def match_blobs_with_neighbours(self):
        '''
        Once blobs that have neighbours have been found, we match them using
        the matching() class (namely, using the voxel method).
        '''
        
        # match particles using the matching object
        M = matching_Ray_Traversal(self.imsys, self.new_pd, self.RIO, self.voxel_size,
                     max_err = self.max_err)
        #return M  # <-- used for checks
        M.get_voxel_dictionary()
        M.list_candidates()
        M.get_particles()
        self.matched_particles = M.matched_particles
        

        
    def return_updated_particle_dict(self):
        '''
        After finding matched particles (self.match_blobs_with_neighbours),  
        this will return an updated copy of particle_dictionary, that does
        not contain the used blobs.
        '''
        updated_pd = self.pd.copy()
        for e,p in enumerate(self.matched_particles):
            p_blobs = p[3]
            for blb in p_blobs:
                ci,(rn, xy) = blb
                cn = self.imsys.cameras[ci].name
                #updated_pd[cn].remove(list(xy))
                ind = updated_pd[cn].index(list(xy))
                updated_pd[cn][ind] = [-1,-1]
                
        return updated_pd




    
        
        
        
