# -*- coding: utf-8 -*-
"""
Created on Fri Dec  7 18:02:07 2018

@author: ron

A test function for the calibration module

"""

from myptv.imaging_mod import camera
from myptv.calibrate_mod import calibrate 



def test_calibrate():
    '''
    A test for the calibrate class. We attempt to do a coarse calibration
    of a synthetic camera.
    '''
    cam = camera('cal_test_cam',(1280,1024), './tests/cal_test_files/cal_test_points')
    cam.load('./tests/cal_test_files')
    cal = calibrate(cam, cam.lab_points, cam.image_points)
    cal.searchCalibration()
    
    O = [1.0, 1.0, 500.0]
    O_err = sum([ (cam.O[i] - O[i])**2 for i in range(3)])**0.5
    
    theta = [0.001, 3.1415, 0.001]
    theta_err = sum([ (cam.theta[i] - theta[i])**2 for i in range(3)])**0.5
    print('mock calibration errors:', O_err, theta_err)
    assert O_err < 1.0 and theta_err < 0.1

