#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 18:02:07 2018

@author: ron

A test function for the smoothing module

"""


from myptv.traj_smoothing_mod import smooth_trajectories
from numpy import loadtxt


def test_smoothing():
    '''
    A test for the smoothing module by smoothing three trajectories.
    '''
    fname = './tests/smoothing_test_files/trajectories'
    traj_list = loadtxt(fname)
    window = 5
    polyorder = 2
    sm = smooth_trajectories(traj_list, window, polyorder)
    sm.smooth()
    tr = sm.smoothed_trajs
    
    test_n_smoothed = len(tr)
    test_structure = all([len(tr[i])==11 for i in range(len(tr))])
    
    assert test_n_smoothed and test_structure


