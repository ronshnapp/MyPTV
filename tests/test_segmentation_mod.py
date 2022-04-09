# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 18:02:07 2018

@author: ron

A test function for the segmentation module

"""


from myptv.segmentation_mod import loop_segmentation



def test_segmentation():
    '''
    A test for the segmentation module, trying to segment blobs from a 
    synthetic noisy image.
    '''
    dirname = './tests/segmentation_test_files'
    
    segmentation = loop_segmentation(dirname, threshold=50, sigma=1.0)
    segmentation.segment_folder_images()
    blobs = segmentation.blobs
    test_blob_number = len(blobs)==8
    assert test_blob_number


