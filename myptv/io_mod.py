#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 08:00:17 2026

@author: owner


A module that helps dealing with saving and loading MyPTV data
as hdf5 files, where data is stored incrementaly as processing 
continues, a chunk saved for every new frame processed.

"""


import h5py
import numpy as np
import pandas as pd
import os.path




def write_to_file(fname, data, dset_type, append=True):
    '''
    A utility function used to save data in hdf5 files for MyPTV processed 
    data, representing the results of the segmentation, matching, tracking or
    smoothing steps. 
    
    - Data is saved in a file with a dataset called "dataset" 
      that containes the actual data. 
      
    - If "dataset" exists, the new data is appended to the file; otherwise
      it is generated with the given data. The idea is that every frame of the
      data is saved into the file separately.
      
    - Everytime data is added to "dataset" we store the index starting that 
      data chunk in a list "offsets". this list can later be used to querry 
      only certain segments of the data. 
      
    - For rvery data added to the file, we store a count of the samples added
      in a list "counts".
      
    - The attribute dataset.attrs['dataset type'] is a metadata, telling
      which type of data is saved here (blobs, particles, trajectories, or
      smoothed trajectories).
    
    
    input:
        
    fname (string) - the file name used
    
    data (2-dimensional numpy array) - data to be saved 
    
    dset_type (string) - Tells the type of data saved. 
    
    append (bool) - If fname exists already and append==True (default), then 
                    data is added to the existing file; if append=False, we 
                    overwrite the existing file.
    '''
    
    if type(fname)!=str:
        raise ValueError('fname must be a string')
    
    if not(fname.endswith('.hdf5')):
        fname = fname + '.hdf5'
        
    allowed_types = ['blobs', 'fiber blobs', 'fiber directions', 
                     'particles', 'trajectories', 'smoothed trajectories'] 
    if dset_type not in allowed_types:
        raise ValueError('ftypes allowed are %s'%(allowed_types))
    
    mode='a'
    if os.path.isfile(fname):
        if append==False: mode='w'
        
    
    nrows = data.shape[1]
    
    with h5py.File(fname, mode) as f:
        
        # if the dataset doesn't exist, create the dataset and metadata
        if 'dataset' not in f.keys():
            
            chunksize = max([1000, data.shape[0]*2])
            
            # 1) creat the dataset itself
            dset = f.create_dataset(
                "dataset",
                data=data,
                maxshape=(None, nrows),
                chunks=(chunksize, nrows),
                compression="gzip",
                dtype=np.float32
            )
            
            dset.attrs['dataset type'] = dset_type
            
            # 2) store the timestep offset index (e.g. [0, 105, 140,...])
            offsets = f.create_dataset(
                "offsets",
                data=[0],
                maxshape=(None,),
                dtype=np.int64
            )
            
            # 3) store the sample count per timestep (e.g. [105, 140, 136, ...])
            counts = f.create_dataset(
                "counts",
                data=[data.shape[0]],
                maxshape=(None,),
                dtype=np.int64
            )
            
            print('saved %d samples in %s'%(data.shape[0], fname))
        
        
        # if the dataset exists, update the dataset and metadata
        else:
            
            Nt = data.shape[0]
            
            # 1) append to the dataset
            dset = f['dataset']
            oldN = dset.shape[0]
            dset.resize(oldN + Nt, axis=0)
            dset[oldN:oldN+Nt] = data
            
            # 2) append the timestep offset index
            offsets = f["offsets"] 
            counts = f["counts"]
            k = offsets.shape[0]
            current_offset = offsets[k-1]
            current_count = counts[k-1]
            offsets.resize(k+1, axis=0)
            offsets[k] = current_offset + current_count
            
            # 3) append the sample counts
            counts.resize(k+1, axis=0)
            counts[k] = Nt
            
            print('added %d samples to %s'%(Nt, fname))
            
    return 






def read_from_file(fname, frame_start=None, frame_end=None):
    '''
    Reads data from a MyPTV hdf5 file of processed and returns it as a pandas
    Dataframe.
    
    Note: if frame_start = x and frame_end = y, then this returns data 
    belonging to frames with indexes from x to y
    
    input:
        
    fname (string) - name of the file to read.
    
    frame_start (int) - Index of the first frame to read. If None (defaul), 
                        it is set to the first available frame. Works like
                        list slicing.
    
    frame_end (int) - Frame index before which we end the read. If None 
                      (default), we end the reading at the end of the dataset.
                      Works like list slicing.
    '''
    
    with h5py.File(fname, "r") as f:
        dset = f['dataset']
        Nt = dset.shape[0]
        offsets = np.array(f['offsets'])
        
        if frame_start is None:
            f0=0
        else:
            f0=offsets[frame_start]
        
        if frame_end is None: # or frame_end==-1 or frame_end==Nf-1:
            fn=Nt
        elif frame_end >= len(offsets):
            fn=Nt
        else:
            fn=offsets[frame_end]
            
        if frame_start is  not None and frame_end is not None:
            if frame_start>frame_end:
                raise ValueError('frame_end should be larger than frame_start.')
            
        return pd.DataFrame(dset[f0:fn])

    


def read_file_frame_range(fname):
    '''
    Returns a set of the frame numbers within a given .hdf5 file.
    
    input:
        
    fname (string) - name of the file to read. 
    '''
    
    with h5py.File(fname, 'r') as f:
        frames = set(f['dataset'][:,-1].astype('int'))
    
    return frames



# if __name__ == '__main__':
#     fname = 'test.hdf5'
#     data = np.random.uniform(0,1,(10,3))
#     dset_type = 'smoothed trajectories'
    
#     write_to_file(fname, data, dset_type, append=True)

#     with h5py.File(fname, "r") as f:
#         print(f['dataset'].shape)
#         print(np.array(f['offsets']))
