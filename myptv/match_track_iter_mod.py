# -*- coding: utf-8 -*-
"""
Created on April 27th 2022

@author: ron

Contains a class to iteratively perform matching and tracking and slowly 
depleat blobs and track more particles

"""

from numpy import loadtxt, array, savetxt
from pandas import read_csv
#from myptv.particle_matching_mod import match_blob_files
#from myptv.tracking_mod import score_optimizing_tracker
from particle_matching_mod import match_blob_files
from tracking_mod import score_optimizing_tracker




class iterative_matcher_tracker(object):
    
    
    
    def __init__(self, blob_fnames, img_system, ROI, voxel_size, 
                 max_blob_dist, mean_flow, d_max, dv_max, max_err, 
                 threshold_traj_len, frame_start=None, frame_end=None):
        
        # load the blobs
        self.blobs = []
        for fn in blob_fnames:
            #self.blobs.append(loadtxt(fn))
            self.blobs.append(array(read_csv(fn, sep='\t', header=None)))
        
        # discard blobs out of the given frame range
        if frame_start is not None:
            for i in range(len(self.blobs)):
                lst = self.blobs[i]
                self.blobs[i] = lst[lst[:,-1]>=frame_start]
        
        if frame_end is not None:
            for i in range(len(self.blobs)):
                lst = self.blobs[i]
                self.blobs[i] = lst[lst[:,-1]<=frame_end]
        
        
        # ====================================================
        #for i in range(len(self.blobs)):
        #    lst = self.blobs[i]
        #    self.blobs[i] = lst[(lst[:,0]>500)&(lst[:,1]>400)]
        # ====================================================
        
        
        self.orig_blobs = self.blobs.copy()
        
        # unpack matching and tracking parameters 
        self.imsys = img_system
        self.ROI = ROI
        self.voxel_size = voxel_size
        self.max_blob_dist = max_blob_dist
        self.max_err = max_err
        self.mean_flow = mean_flow
        self.d_max = d_max
        self.dv_max = dv_max
        self.threshold_traj_len = threshold_traj_len
        
        
        # a list that will hold the trajectories
        self.trajectories = []
        
        self.count_blobs = 0
        
        
        
    def match_and_track(self, voxel_size = None):
        '''
        performas a round of matching and then tracking the unused blobs.
        '''
        
        print('\n', 'number of available blobs:',
              sum([len(lst) for lst in self.blobs]), '\n')
        
        if voxel_size is not None:
            self.voxel_size = voxel_size
        
        # Doing the matching 
        mbf = match_blob_files([], self.imsys, self.ROI, self.voxel_size, 
                               self.max_blob_dist, blobs = self.blobs, 
                               max_err = self.max_err)
        mbf.get_particles()
        particles = mbf.return_particles()
        
        
        # Tracking the newly obtained particles
        self.sot = score_optimizing_tracker([], mean_flow = self.mean_flow, 
                                       d_max=self.d_max, dv_max=self.dv_max,
                                       particles = particles)
        self.sot.track_frames()
        
        
        # a nested list of dictionaries that hold the indexes of the blobs 
        # that were used to make trajectories. The sublist index is the camera 
        # number. The dicionaries keys are frame numbers and values are lists
        # that hold the used blob indexes.
        self.used_blobs = [{} for i in range(len(self.imsys.cameras))]
        particles_frames = set([p[-1] for p in particles])
        for dic in self.used_blobs:
            for frame in particles_frames:
                dic[frame] = []
        
        # fetching long trajs, and listing the used blobs 
        for traj in self.sot.trajs:
            if len(traj) > self.threshold_traj_len:
                self.trajectories.append(traj)
                for i in range(len(traj)):
                    blob_indexes = traj[i][4:-2]
                    frame = traj[i][-1]
                    for ci in range(len(blob_indexes)):
                        if blob_indexes[ci]==-1:
                            continue
                        
                        self.used_blobs[ci][frame].append(blob_indexes[ci])
                        self.count_blobs += 1
                        
                        # try:
                        #     self.used_blobs[ci][frame].append(blob_indexes[ci])
                        #     self.count_blobs += 1
                            
                        # except:
                        #     self.used_blobs[ci][frame] = [blob_indexes[ci]]
                        #     self.count_blobs += 1
                        
        print('')
        print('used %d blobs'%(self.count_blobs))
        print('obtained %d long trajectories in total'%len(self.trajectories))
        np = sum([len(tr) for tr in self.trajectories])
        print('with %d particles'%np, '\n')
        # updating self.blobs to remove blobs used in long trajectories
        updated_blobs_list = [[]  for i in range(len(self.imsys.cameras))]
        for ci in range(len(self.blobs)):
            for frame in particles_frames:
                fltr = lambda blb: blb[-1]==frame
                frame_blobs = list(filter(fltr, self.blobs[ci]))
                for i in range(len(frame_blobs)):
                    if i not in self.used_blobs[ci][frame]:
                        updated_blobs_list[ci].append(frame_blobs[i])
                
            # for i in range(len(self.blobs[ci])):
            #     frame = self.blobs[ci][i][-1]
            #     if i not in self.used_blobs[ci][frame]:
            #         updated_blobs_list[ci].append(self.blobs[ci][i])
                
                
        for i in range(len(updated_blobs_list)):
            self.blobs[i] = array(updated_blobs_list[i])
            
            
    def finalize(self):
        '''
        Adds the rest of "short trajectories" to the trajectory list
        '''
        for traj in self.sot.trajs:
            if len(traj) <= self.threshold_traj_len:
                self.trajectories.append(traj)



    def plot_blobs(self, camera_number=0, frames=None):
        '''
        Plots all the blobs belonging to camera with given camera_number
        and highlights the blobs that were used to generate trajectories
        '''
        import matplotlib.pyplot as plt
        
        if frames is None:
            particles_frames = self.used_blobs[0].keys()
        else:
            particles_frames = frames
        
        
        fig, ax = plt.subplots()
        for blb in self.orig_blobs[camera_number]:
            if blb[-1] in particles_frames:
                ax.plot(blb[0], blb[1], 'xk', ms=4, alpha=0.5)
        
        
        # updating the used blobs data
        self.used_blobs = [{} for i in range(len(self.imsys.cameras))]
        for dic in self.used_blobs:
            for frame in particles_frames:
                dic[frame] = []
         
        for traj in self.trajectories:
            for i in range(len(traj)):
                blob_indexes = traj[i][4:-2]
                frame = traj[i][-1]
                if frame in particles_frames:
                    for ci in range(len(blob_indexes)):
                        if blob_indexes[ci]==-1:
                            continue
                        self.used_blobs[ci][frame].append(blob_indexes[ci])
        
        
        for frame in self.used_blobs[camera_number].keys():
            fltr = lambda blb: blb[-1]==frame
            frame_blobs = list(filter(fltr, self.orig_blobs[camera_number]))
            for bi in self.used_blobs[camera_number][frame]:
                x_, y_ = frame_blobs[int(bi)][:2]
                ax.plot(x_, y_, 'rx', ms=5)
        
        




if __name__ == '__main__':
    from myptv.imaging_mod import img_system, camera
    import matplotlib.pyplot as plt
    import numpy as np
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    
    blob_fnames = ['/home/ron/working_PTV_data/1/blobs_cam1',
                   '/home/ron/working_PTV_data/1/blobs_cam2',
                   '/home/ron/working_PTV_data/1/blobs_cam3',
                   '/home/ron/working_PTV_data/1/blobs_cam4']
    
    cam1 = camera('cam1', (1280,1024))
    cam2 = camera('cam2', (1280,1024))
    cam3 = camera('cam3', (1280,1024))
    cam4 = camera('cam4', (1280,1024))
    
    for c in [cam1, cam2, cam3, cam4]: c.load('/home/ron/working_PTV_data/1')
    
    imsys = img_system([cam1, cam2, cam3, cam4])
    
    ROI = [[0.0, 70.0] ,[0.0 ,70.0], [-25.0 ,10.0]]
    #ROI = [[20.0,50.0] ,[20.0 ,50.0], [-25.0 ,10.0]]
    voxel_size = 1.5
    max_blob_dist = 1.0
    max_err=1.0
    mean_flow = 0.0 
    d_max = 0.5
    dv_max=1.0
    threshold_traj_len = 10
    frame_start = None
    frame_end = 20
    
    imt = iterative_matcher_tracker(blob_fnames, 
                                    imsys, 
                                    ROI, 
                                    voxel_size, 
                                    max_blob_dist, 
                                    mean_flow, 
                                    d_max, 
                                    dv_max,
                                    max_err,
                                    threshold_traj_len,
                                    frame_start=frame_start, 
                                    frame_end=frame_end)
    
    
    # used_blobs = []
    # for i in range(2):
    #     imt.match_and_track()
    #     used_blobs.append(imt.used_blobs)
    
    # frm_num = 2
    
    # blobs1_frm = loadtxt(blob_fnames[0])
    # blobs1_frm = blobs1_frm[blobs1_frm[:,-1]==frm_num]
    # fig, ax = plt.subplots()
    # for e,iter_ in enumerate(used_blobs[::-1]):
    #     blob_indexes = iter_[0][frm_num]
    #     for bi in blob_indexes:
    #         blb = blobs1_frm[int(bi)]
    #         ax.plot(blb[0], blb[1], 'x', color=colors[e])

            
            
    # voxel_size_list = [1.5, 2.5, 4.0, 5.0]        
    # fig, ax = plt.subplots()
    # j = 0
    # for i in range(4):
    #     imt.match_and_track(voxel_size = voxel_size_list[i])
    #     while j<len(imt.trajectories):
    #         tr = imt.trajectories[j]
    #         ax.plot(tr[:,1], tr[:,2], '-', lw=1, color=colors[i])
    #         j+=1
    #     ax.set_xlim(ROI[0][0], ROI[0][1])
    #     ax.set_ylim(ROI[1][0], ROI[1][1])
    #     ax.set_aspect('equal')
    #     fig.savefig('iter_%d.jpg'%i, dpi=300)
    
    
    voxel_size_list = [1.0, 1.5, 1.5]
    for vs in voxel_size_list[:]:
        print(vs)
        imt.match_and_track(voxel_size = vs)
    
    imt.finalize()
    
    tl = [len(tr) for tr in imt.trajectories]
    print(np.mean(tl))
    h = np.histogram(tl, bins=range(1,51))
    x_, y_ = 0.5*(h[1][1:]+h[1][:-1]), h[0]
    plt.semilogy(x_, y_, 'o-')
    
    #imt.plot_blobs(frames=[1.0,2.0])
    
    
    
    
    
    
    