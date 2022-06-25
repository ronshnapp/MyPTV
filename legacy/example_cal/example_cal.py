# -*- coding: utf-8 -*-
"""
Created on Sat Nov  2 09:32:43 2019

@author: ron

example of calibrating a camera

"""



import os, sys

fpath = os.path.dirname(os.path.abspath(__file__))
#print os.path.join(str(fpath), '/..')
#sys.path.append(os.path.join(fpath, '/..'))
sys.path.append( str(fpath) + '/..' )

from utils import *
from imaging_mod import *
from calibrate_mod import *
import numpy as np
import matplotlib.pyplot as plt



#=======================
# load the image points:
#=======================

img_file = open('img_points', 'r')
a = img_file.readlines()
img_file.close()

img_coords_dic = {}

for s in a:
    s_ = s.split()
    
    if len(s_) == 3:
        n, x_, y_ = s_
        img_coords_dic[int(n)] = [float(x_), float(y_)]


#=====================
# load the lab points:
#=====================

lab_file = open('lab_coords', 'r')
a = lab_file.readlines()
lab_file.close()

lab_coords_dic = {}

for s in a:
    s_ = s.split()
    
    if len(s_) == 4:
        n, x_, y_, z_ = s_
        lab_coords_dic[int(n)] = [float(x_), float(y_), float(z_)]


lab_coords, img_coords = [], []
for k in sorted(lab_coords_dic.keys()):
    lab_coords.append(lab_coords_dic[k])
    img_coords.append(img_coords_dic[k])
    
    
    
    
#===================
# start calibration:
#===================
    
    
# initiate a camera with an initial guess:
c = camera('cam1', (3246, 2448))
c.O = np.array([1.0, 1.0, 135.0])
c.theta = np.array([3.14, 0.1, 0.1])
c.f=3000
c.calc_R()

cal = calibrate(c, lab_coords, img_coords)


print('\n')
print('D0 = %.1f'%(cal.mean_squared_err()))
print('\n')


# plot the initial guess:
fig, ax = plt.subplots()
cal.plot_proj(ax)





cal.manual_calibration()           # <-- Do some manual gross calibration
Min_res = cal.searchCalibration()  # <-- automatic solve calibration


# show the results:
print('\n')
print(c.O)
print(c.theta)
print(c.f)
print('\n')

print('D = %.1f'%(cal.mean_squared_err()))
print('\n')

cal.plot_proj()

fig, ax = plt.subplots()
ax.semilogy(cal.D_lst)
