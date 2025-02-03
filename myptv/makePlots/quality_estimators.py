# -*- coding: utf-8 -*-
#!/usr/bin/env python3

"""
Created on Fri May 31 15:01:48 2024

@author: ron
"""


from pandas import read_csv
from numpy import array, mean
import matplotlib.pyplot as plt
import tqdm


    
    
    
class check_matching():
    '''
    A class that containes several functions meant to help visualize 
    uncertainties and assess the quality of the analysis in the stereo 
    matching results. 
    '''
    
    def __init__(self, imsys, blob_files, particles_filename, dmax,
                 min_cam_match=3):
        '''
        A class that containes several functions meant to help visualize 
        uncertainties and assess the quality of the analysis in the stereo 
        matching results. 
        
        input -
        =======
        
        imsys - an img_system instance with calibrated cameras
        
        blob_files - a list of file names containing the resutls of the 
                     segmentation
        
        particles_filename - the name of the file containing the results of the
                             stereo matching step (the particles file)
                             
        dmax - maximum allowed triangulation uncertainty in the stereo matching
               step
               
        min_cam_match - the minimal number of cameras needed to triangulate a
                        particle in the stereo matching
        '''
        
        
        
        self.imsys = imsys
        self.dmax = dmax
        self.min_cam_match = min_cam_match
        
        # load particles and sort according to frame
        particles_data = read_csv(particles_filename, sep='\t', header=None)
        gb = particles_data.groupby(particles_data.shape[1]-1)
        self.particles = dict([(k,array(g)) for k, g in gb])
        particles_data = None
    
        # load blobs an arrange in lists according to their frame
        self.blob_dics = []
        for fn in blob_files:
            blob_data = read_csv(fn, sep='\t', header=None)
            gb = blob_data.groupby(5)
            self.blob_dics.append(dict([(k,array(g)) for k, g in gb]))
            blob_data = None
            
            
            
            
    def get_coords(self, frame, particleNum):
        '''
        Returns the coordinates dictionary of the blobs belonging to a 
        given particle.
        '''
        
        p = self.particles[frame][particleNum]
        
        cords = {}
        
        for cn in range(len(self.imsys.cameras)):
            blob_ind = int(p[cn + 3])
            if blob_ind == -1: continue
            blob = self.blob_dics[cn][frame][blob_ind]
            cords[cn] = blob[:2][::-1]
        
        return cords
            
            
            
    def rematch_particle(self, frame, particleNum):
        '''
        takes the blob data of a particle, re-matches it, and returns the
        result and the blobs.
        '''
        p = self.particles[frame][particleNum]
        cords = self.get_coords(frame, particleNum)
        
        if len(cords.keys()) < self.min_cam_match:
            return p[:3], -1, -1
        
        rm = self.imsys.stereo_match(cords, self.dmax, strict_match=True)
        
        if rm is None:
            return p[:3], -1, -1
        
        return p[:3], rm[0], rm[2]
    
    
    
    
    def get_particle_disparity(self, frame, particleNum, show_list=False):
        '''
        For a given particle, rematches it, and calculates its disparities for
        all the cameras.
        '''
        coords = self.get_coords(frame, particleNum)
        
        if len(coords.keys()) < self.min_cam_match:
            return -1
        
        p_rematch = self.imsys.stereo_match(coords, self.dmax, 
                                            strict_match=True)
        
        if p_rematch is None:
            return -1
        
        disp_lst = []
        
        for cam_ind in coords.keys():
            b = coords[cam_ind]
            proj = self.imsys.cameras[cam_ind].projection(p_rematch[0])
            disparity = sum((b-proj)**2)**0.5
            
            disp_lst.append(disparity)
        
        if show_list==True:
            print(disp_lst)
        
        return mean(disp_lst)
    
    
    
    
    def get_disparity_list(self, frames=None):
        '''
        Returns a list of the disparity between particle projection and blobs
        for all particles in the file or for a given list of frames.
        '''
        if frames is None:
            frames = list(self.particles.keys())
        
        disparity_lst = []
        for frame in frames:
            desc = 'frame: %d'%frame
            for p_ind in tqdm.tqdm(range(len(self.particles[frame])), 
                                   desc=desc):
                disparity_lst.append(self.get_particle_disparity(frame, p_ind))
        return disparity_lst
    
    
    
    
    def plot_blobs_and_particles(self, frames, camNum):
        '''
        Plots the blobs and the particle images belonging to a given list of 
        frames and camera number.
        '''
        
        fig, ax = plt.subplots()
        
        projx = [] ; projy = []
        blobx = [] ; bloby = []
        
        for frame in frames:
            for i in range(len(self.particles[frame])):
                
                coords = self.get_coords(frame, i)
                if camNum-1 in list(coords.keys()): 
                    blob = coords[camNum-1]
                    
                    p = self.particles[frame][i][:3]
                    proj = self.imsys.cameras[camNum-1].projection(p)
                    
                    projx.append(proj[0]) ; projy.append(proj[1])
                    blobx.append(blob[0]) ; bloby.append(blob[1])
                    
                    ax.plot([blob[0], proj[0]], [blob[1], proj[1]], '-k')
        
        ax.plot(projx, projy, 'or', label='particle projections')
        ax.plot(blobx, bloby, 'xb', label='segmentation results')
        ax.legend()
        
        
    def plot_disparity_histogram(self, frames=None):
        '''
        Plots a histogram for the particle disparity relative to their image
        projections.
        '''
        disp = self.get_disparity_list(frames=frames)
        fig, ax = plt.subplots()
        h = ax.hist(disp, bins='auto')
        ax.set_xlabel('disparity [pixel]')
        ax.set_ylabel('counts')
            
            
    def plot_disparity_vs_stereo_matching_err(self, frames=None):
        '''
        Generates a scatter plot of the disparity in the particle projection
        and the stereo matching error (the average distance between the
        crossing points of each ray and the mean crossing point).
        '''
        
        disp = []
        sm_err = []
        
        if frames is None:
            frames = list(self.particles.keys())
        
        for frame in frames:
            for p_ind in range(len(self.particles[frame])):
                disp.append(self.get_particle_disparity(frame, p_ind))
                sm_err.append(self.rematch_particle(frame, p_ind)[2])
            
        fig, ax = plt.subplots()
        ax.scatter(disp, sm_err, s=2)
        
        ax.set_xlabel('projection disparity [px]')
        ax.set_ylabel('stereo matching err [physical length]')
            
    
    
    
    
    
    
import tkinter as tk
from tkinter import messagebox, simpledialog

class matching_quality_GUI(check_matching):
    
    def __init__(self, imsys, blob_files, particles_filename, dmax,
                 min_cam_match=3):
        '''
        A graphical user interphase for assessing the stereo matching quality. 
        input -
        =======
        
        imsys - an img_system instance with calibrated cameras
        
        blob_files - a list of file names containing the resutls of the 
                     segmentation
        
        particles_filename - the name of the file containing the results of the
                             stereo matching step (the particles file)
                             
        dmax - maximum allowed triangulation uncertainty in the stereo matching
               step
               
        min_cam_match - the minimal number of cameras needed to triangulate a
                        particle in the stereo matching
        '''
        
        # ====================================================================
        # Setup for the functions derived from the check_matching class
        self.imsys = imsys
        self.dmax = dmax
        self.min_cam_match = min_cam_match
        
        # load particles and sort according to frame
        particles_data = read_csv(particles_filename, sep='\t', header=None)
        gb = particles_data.groupby(particles_data.shape[1]-1)
        self.particles = dict([(k,array(g)) for k, g in gb])
        particles_data = None
    
        # load blobs an arrange in lists according to their frame
        self.blob_dics = []
        for fn in blob_files:
            blob_data = read_csv(fn, sep='\t', header=None)
            gb = blob_data.groupby(5)
            self.blob_dics.append(dict([(k,array(g)) for k, g in gb]))
            blob_data = None
        
        button_width = 40
        
        
        fn = max(self.particles.keys())
        f0 = min(self.particles.keys())
        self.f_range = range(int(f0), int(fn+1))
        
        # ====================================================================
        
        
        # the GUI part 
        
        root = tk.Tk()
        root.title("Tkinter GUI Example")
        root.geometry("350x250")
        
        frame_label = tk.Label(root, 
                               text="Available frame range: %d-%d"%(f0,fn), 
                               font=("Arial", 12))
        frame_label.pack(pady=10)
        
        btn1 = tk.Button(root, text="Select Frames", 
                         command=self.select_frames,
                         width=button_width)
        btn1.pack(pady=10)
        
        btn2 = tk.Button(root, text="Plot Blobs and Particles",
                         width=button_width, 
                         command=self.btn_plot_blobs_and_particles)
        btn2.pack(pady=10)
        
        btn3 = tk.Button(root, text="Plot Disparity vs Stereo-Matching Error", 
                         command=self.btn_plot_disparity_vs_error,
                         width=button_width)
        btn3.pack(pady=10)
        
        btn4 = tk.Button(root, text="Plot Disparity Histogram", 
                         command=self.btn_plot_disparity_histogram,
                         width=button_width)
        btn4.pack(pady=10)
        
        root.mainloop()
        
    
    def select_frames(self):
        frame_range = simpledialog.askstring("Select Frames", "Enter frame range (e.g., 1-100):")
        if frame_range:
            messagebox.showinfo("Selected Frames", f"You selected frames: {frame_range}")
            
            self.f0 = int(frame_range.split('-')[0])
            self.fn = int(frame_range.split('-')[-1])
            self.f_range = list(range(self.f0, self.fn+1))

    def btn_plot_blobs_and_particles(self):
        camNum = simpledialog.askstring("Camera Number", "Choose Camera Number (e.g. '1'):")
        self.plot_blobs_and_particles(self.f_range, int(camNum))
        plt.show()
        
    def btn_plot_disparity_vs_error(self):
        self.plot_disparity_vs_stereo_matching_err(self.f_range)
        plt.show()
        
    def btn_plot_disparity_histogram(self):
        self.plot_disparity_histogram(self.f_range)
        plt.show()








# =============================================================================
# 
# if __name__ == '__main__':
# 
#         
#     folder = '../../example'
#     
#     pfn = folder + '/particles'
#     
#     from myptv.imaging_mod import camera_wrapper, img_system
#     
#     cams = [camera_wrapper('cam%d'%i, folder) for i in [1,2,3]]
#     for c in cams: c.load()
#     imsys = img_system(cams)
#     
#     bfns = [folder + '/blobs_cam%d'%i for i in [1,2,3]]
#     
#     dmax = 0.25
#     min_cam_match = 3
#     
#     cm = check_matching(imsys, bfns, pfn, dmax, min_cam_match=min_cam_match)
#     
#     
#     G = matching_quality_GUI(imsys, bfns, pfn, dmax, 
#                              min_cam_match=min_cam_match)
# =============================================================================
    

    
    
    