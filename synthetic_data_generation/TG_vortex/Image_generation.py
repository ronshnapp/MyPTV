# -*- coding: utf-8 -*-
#!/usr/bin/env python3

"""
Created on Tue Jan 20 21:29:07 2026

@author: ron


This is a script used to generate synthetic PTV images for
testing purposes. We use a four camera system is an 
approximate rectangular configuration based on actual experimental
data. The cameras are calibrated using the zolof model. The
blobs in the images include a synthetic noise that is emulated
by having each blob made out of a number of image dots with
random offsets. The image are saved as 8 bit jpeg files which come 
out rather light weight. 

To generate the images, first inspect the various parameters 
that appear throughout the script, and then run it.

""" 

from tqdm import tqdm
from myptv.imaging_mod import camera_wrapper
from PIL import Image, ImageDraw
import numpy as np
from numpy.random import uniform



# ===============================
# Taylor green vortex flow field:
# ===============================


def get_u(x):
    '''
    Returns the velocity vector at a given position, 
    x = (x,y,z)
    '''
    A = 1 ; B = -1 ; C = 0
    a = np.pi*2/21 ; b = np.pi*2/21 ; c = np.pi*2/41
    x_ = np.array([a,b,c]) * x
    u = A * np.cos(x_[0]) * np.sin(x_[1]) * np.sin(x_[2])
    v = B * np.sin(x_[0]) * np.cos(x_[1]) * np.sin(x_[2])
    w = C * np.sin(x_[0]) * np.sin(x_[1]) * np.cos(x_[2])
    return np.array([u,v,w])



# ===========================================
# Constructing initial trajectory conditions:
# ===========================================


nx = 20
ny = 20
nz = 20
n_particles = nx * ny * nz
x0, xn = 5, 95
y0, yn = 5, 95
z0, zn = -50, 0

x_lst = []
v_lst = []

for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            x = np.random.uniform(x0, xn) #x0 + (i+0.5)*(xn-x0)/nx
            y = np.random.uniform(y0, yn) #y0 + (j+0.5)*(yn-y0)/ny
            z = np.random.uniform(z0, zn) #z0 + (k+0.5)*(zn-z0)/nz
            
            x = np.array([x, y, z])
            x_lst.append([x])
            v_lst.append([get_u(x)])
            
            
            

# ==========================
# Constructing trajectories:
# ==========================


dt = 0.8
steps = 50
ntraj = len(x_lst)

for j in range(ntraj):
    xi = x_lst[j][0]
    vi = v_lst[j][0]
    
    for i in range(steps):
        dx = vi * dt
        xi = xi + dx
        vi = get_u(xi)
        x_lst[j].append(xi)
        v_lst[j].append(vi)

print(len(x_lst[0]))
print(n_particles)

translations = []
for i in range(len(x_lst)):
    translations += list(np.sum(np.diff(x_lst[i], axis=0)**2, axis=1)**0.5)
print('RMS translations: ', np.mean(translations))
print('image pixels: ', 1650*1650)


# =====================
# Loading cameras:
# =====================

camNames = ['cam1', 'cam2', 'cam3', 'cam4']
cams = [camera_wrapper(cn, '.') for cn in camNames]
for cam in cams:
    cam.load()


# ======================
# Generating the images:
# ======================

resolution = 2048, 2048     # image pixel resolution
n_blobsPerParticle = 3      # number of blobs making up a particle
disparity = 0.03            # random noise of the blobs making up the particles [mm]
s = 0.8                     # sigma of a gaussian intensity point
I = 15.0                    # intensity of a fiber making point


radius = int(3*s)
x_ = range(-radius, radius+1)
X, Y = np.meshgrid(x_, x_)
blob_image = I * np.exp(-( ((X)**2 + (Y)**2) /2/s**2))


for k in tqdm(range(steps)):
    
    images = [np.zeros((resolution[1], resolution[0])) for cam in cams]
    pos_lst = []
        
    for i in range(n_particles):
        
        # particle center position
        X = x_lst[i][k]
        pos_lst.append(X)
        
        for e, cam in enumerate(cams):
            img = images[e]
            
            for j in range(n_blobsPerParticle):
                noise = np.random.normal(0, disparity, size = 3)
                cx, cy = np.array(cam.projection(X + noise)).astype('int')
                i0, i1, j0, j1 = cy-radius, cy+radius+1, cx-radius, cx+radius+1
                if (np.array([i0, i1, j0, j1]) < 0).any() or (np.array([i0, i1, j0, j1]) > 2047).any():
                    continue
                img[i0:i1, j0:j1] += blob_image
        
       
    for e, cam in enumerate(cams):
        img = images[e]
        img[img>255] = 255
        img = img.astype('int8')

        pil_im = Image.fromarray(img, mode='L')
        pil_im.save('./cam%d_images/im%03d.jpg'%(e+1,k+1))
    
    
    np.savetxt('./ground_truth/ground_truth_pos_im_%02d'%(k+1), pos_lst, delimiter='\t', fmt='%.3f')