#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 18:02:07 2018

@author: ron

A test function for the imaging module

"""

from myptv import imaging_mod 
from numpy import array
from math import pi


def test_imaging():
    '''
    A test for the camera and imaging system objects. We project a 
    synthetic point onto three cameras and check that the stereo mathing
    gives the same point back.
    '''
    c1 = imaging_mod.camera('1', (1000.,1000.))
    c2 = imaging_mod.camera('2', (1000.,1000.))
    c3 = imaging_mod.camera('3', (1000.,1000.))
    
    c1.O = array([400.0 , 0, 1])
    c2.O = array([0, 400.0, -1])
    c3.O = array([200.0, 400.0 ,400])
    
    c1.f = 4000
    c2.f = 4000
    c3.f = 4000
    
    c1.theta = array([0.0, -1*pi / 2.0, 0.0])
    c2.theta = array([pi / 2.0, 0., 0.])
    c3.theta = array([0.8, -0.4, 0.0])
    
    c1.calc_R()
    c2.calc_R()
    c3.calc_R()
    
    c1.xh = 1.0
    c1.yh = -1.0
    
    x = array([0.1,0.1,0.1])        
    
    imgsys = imaging_mod.img_system([c1,c2,c3])
    
    proj1 = c1.projection(x)
    proj2 = c2.projection(x)
    proj3 = c3.projection(x)

    
    coords = {0: proj1,
              1: proj2,
              2: proj3}
    
    res = imgsys.stereo_match(coords, 1e9)
    
    a, b, c = round(res[0][0], 10), round(res[0][1], 10), round(res[0][2], 10)
    assert a == 0.1 and b == 0.1, c == 0.1

