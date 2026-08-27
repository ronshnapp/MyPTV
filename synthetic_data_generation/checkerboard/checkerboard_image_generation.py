# -*- coding: utf-8 -*-
#!/usr/bin/env python3

"""
Created on Aug 27 2026

@author: Stefano Brizzolara

This script generates synthetic images of a checkerboard target from
calibrated camera files. It is the counterpart of
TG_vortex/calibration_image_generation.py for the checkerboard calibration
procedure, and it exists so that the procedure can be tried out, and its
parameters understood, without having to run an experiment first.

The board is imaged at several positions of an imaginary translation stage,
which is the arrangement that the checkerboard procedure expects: one board,
photographed by all the cameras at the same time, at a number of known
positions along the stage.

Each square of the board is drawn by projecting its four corners through the
camera model and filling the resulting quadrilateral. Every 3D model is
therefore handled correctly, including its distortion terms, to the accuracy
with which a small square maps onto a quadrilateral.

Run it from the directory holding the camera files, for instance from
MyPTV/example:

    python ../synthetic_data_generation/checkerboard/checkerboard_image_generation.py

It writes the images into the calibration folder, and prints the parameters
that should then go into the checkerboard_calibration block of the parameters
file, including the pixel position of the board origin in each camera.

"""

import os

from numpy import array, cross
from numpy.linalg import norm

from PIL import Image, ImageDraw

from myptv.imaging_mod import camera_wrapper




# ======================================
# What to generate:
# ======================================

camera_names = ['cam1', 'cam2', 'cam3']
camera_dir = '.'

# the folder the images are written into
output_dir = './Calibration'

# the board: the number of INTERNAL corners along its two sides, and the
# physical spacing of the corners in the units of the lab coordinates
board_size = (9, 6)
square_size = 7.0

# the pose of the board at zero translation. The origin is chosen so that the
# board is centred on the lab origin, which is roughly where the measurement
# volume of the example is.
i_axis = array([1.0, 0.0, 0.0])
j_axis = array([0.0, 1.0, 0.0])

origin = array([35.0, 35.0, 0.0]) \
    - (board_size[0] - 1)/2.0*square_size*i_axis \
    - (board_size[1] - 1)/2.0*square_size*j_axis

# the positions of the translation stage, along the board normal
translations = [-8.0, -4.0, 0.0, 4.0, 8.0]

# the appearance of the images
dark, light, background = 25, 230, 128
blur = 1.2

# The polygons are drawn into an image supersample times larger in each
# direction, which is then averaged down. This matters: PIL fills polygons
# without anti-aliasing, so drawing at the final size quantizes every edge to
# a whole pixel, and the corners of the board then sit up to half a pixel away
# from where the camera model puts them. Averaging down from a larger image
# removes that, and brings the corners of the generated images to within about
# a twentieth of a pixel of the exact projection.
supersample = 4




# ======================================
# The generation:
# ======================================

def board_to_lab(u, v, translation):
    '''
    The lab coordinates of the point of the board plane at the board
    coordinates (u, v), for a board displaced by translation along its
    normal.
    '''
    e1 = i_axis/norm(i_axis)
    e2 = j_axis/norm(j_axis)
    n = cross(e1, e2)
    n = n/norm(n)
    return origin + u*e1 + v*e2 + translation*n


def render(cam, resolution, translation):
    '''
    Draws the board as seen by cam, with the board displaced by translation.
    Returns a PIL image.
    '''
    ni, nj = board_size
    s = square_size
    ss = int(supersample)

    Nx, Ny = int(resolution[0]), int(resolution[1])

    img = Image.new('L', (Nx*ss, Ny*ss), background)
    draw = ImageDraw.Draw(img)

    def project(u, v):
        '''
        The pixel position of the board point (u, v), in the coordinates of
        the supersampled image. A point at pixel eta of the final image falls
        at (eta + 0.5)*ss - 0.5 of the supersampled one, which is the mapping
        that makes averaging blocks of ss pixels put it back at eta.
        '''
        eta, zeta = cam.projection(board_to_lab(u, v, translation))
        return ((float(eta) + 0.5)*ss - 0.5, (float(zeta) + 0.5)*ss - 0.5)

    # the pattern spans one extra square beyond the internal corners on each
    # side, so that all of the internal corners are surrounded by the board
    def quad(a, b):
        '''the four projected corners of the square with indices (a, b)'''
        return [project((a + du - 1)*s, (b + dv - 1)*s)
                for du, dv in [(0, 0), (1, 0), (1, 1), (0, 1)]]

    # a light backing sheet one square larger than the pattern, as a printed
    # target has
    sheet = [project(u, v) for u, v in
             [(-2*s, -2*s), ((ni + 1)*s, -2*s),
              ((ni + 1)*s, (nj + 1)*s), (-2*s, (nj + 1)*s)]]
    draw.polygon(sheet, fill=light)

    for a in range(ni + 1):
        for b in range(nj + 1):
            if (a + b) % 2 == 0:
                draw.polygon(quad(a, b), fill=dark)

    # average the supersampled image down to the camera resolution
    if ss > 1:
        img = img.reduce(ss)

    if blur > 0:
        from PIL import ImageFilter
        img = img.filter(ImageFilter.GaussianBlur(radius=blur))

    return img


if __name__ == '__main__':

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    print('generating checkerboard images for %s\n'%str(camera_names))

    for name in camera_names:

        cam = camera_wrapper(name, camera_dir)
        cam.load()

        resolution = cam.camera.resolution

        for e, t in enumerate(translations):
            img = render(cam, resolution, t)
            fname = os.path.join(output_dir, '%s_board_%02d.tif'%(name, e))
            img.save(fname)

        # the position of the board origin in this camera, which is what the
        # origin_hint parameter needs
        X0 = board_to_lab(0.0, 0.0, translations[0])
        eta, zeta = cam.projection(X0)

        print('%s: wrote %d images; origin_hint: %.0f, %.0f'
              %(name, len(translations), eta, zeta))

    print('\nAdd a block like this to your parameters file, one per camera:\n')
    print('- checkerboard_calibration:')
    print('    camera_name: %s'%camera_names[0])
    print('    images: %s/%s_board_*.tif'%(output_dir, camera_names[0]))
    print('    board_size: %d, %d'%board_size)
    print('    square_size: %s'%square_size)
    print('    board_origin: %.4f, %.4f, %.4f'%tuple(origin))
    print('    board_i_axis: %.0f, %.0f, %.0f'%tuple(i_axis))
    print('    board_j_axis: %.0f, %.0f, %.0f'%tuple(j_axis))
    print('    translations: %s'%(', '.join(str(t) for t in translations)))
    print('    origin_hint: <the two numbers printed above for that camera>')
