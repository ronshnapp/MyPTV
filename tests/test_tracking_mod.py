# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 18:02:07 2018

@author: ron

A test function for the tracking module

"""


from myptv.tracking_mod import tracker_four_frames



def test_tracking():
    '''
    A test for the tracking algorithm by using the four-frame tracking class
    over a dummy particle file with 3 particles moving over 4 frames.
    '''
    fname = './tests/tracking_test_files/test_particles'
    t4f = tracker_four_frames(fname)
    t4f.track_all_frames()
    p_list = t4f.return_connected_particles()
    
    traj_numbers = list(set([p[0] for p in p_list]))
    traj_number_test = len(traj_numbers)==3
    
    traj_lengths = [0 for i in range(len(traj_numbers))]
    for p in p_list:
        ind = int(p[0])
        traj_lengths[ind] += 1
    
    traj_length_test =  all([tl==4 for tl in traj_lengths])
    
    assert traj_length_test and traj_number_test


