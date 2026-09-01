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

It can generate images for either of the two checkerboard procedures, chosen
by the 'mode' setting below: a board on a translation stage photographed at
known positions, for checkerboard_calibration, or a board turned about freely
so that its pose is not known, for moving_board_calibration. In both cases
every camera sees the same board at the same time, which is what the
procedures expect.

Each square of the board is drawn by projecting its four corners through the
camera model and filling the resulting quadrilateral. Every 3D model is
therefore handled correctly, including its distortion terms, to the accuracy
with which a small square maps onto a quadrilateral.

Run it from the directory holding the camera files, for instance from
MyPTV/example:

    python ../synthetic_data_generation/checkerboard/checkerboard_image_generation.py

It writes the images into the calibration folder and prints the block of
parameters that should then go into the parameters file, including, for the
stage case, the pixel position of the board origin in each camera.

"""

import os

from numpy import array, cross
from numpy.linalg import norm

from PIL import Image, ImageDraw

from myptv.imaging_mod import camera_wrapper




# ======================================
# What to generate:
# ======================================

# Which of the two procedures the images are for:
#
#   'stage'   the board sits on a translation stage and is photographed at
#             known positions along it, for checkerboard_calibration
#
#   'moving'  the board is held in the volume and turned about between
#             exposures, so its pose is not known and has to be recovered,
#             for moving_board_calibration. Every camera sees the same pose
#             at the same frame, as that procedure requires.
mode = 'stage'

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

# 'stage' mode: the positions of the translation stage, along the board normal
translations = [-8.0, -4.0, 0.0, 4.0, 8.0]

# 'moving' mode: how many frames, how far the board is allowed to wander from
# the centre of the volume along each axis, and how far it may be turned.
#
# The board must be genuinely turned about, not merely slid around, because
# the intrinsics are recovered from the way the perspective changes with
# orientation. It must also travel in depth. The third component of wander is
# much the largest for that reason: the focal length trades off against how
# far away the board is, so a board that stays at one distance leaves the two
# poorly separated. Generated with the board confined to 12 mm about one
# point, the focal lengths came back 0.59 percent out and the distances
# between the cameras 0.78 percent out; with the settings below, 0.17 and
# 0.23 percent.
n_frames = 30
wander = array([30.0, 30.0, 200.0])     # mm, about the centre of the volume
max_tilt = 1.0                          # radians
volume_centre = array([35.0, 35.0, 0.0])
seed = 0

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

def stage_pose(translation):
    '''
    The pose of the board on the translation stage, at a given position along
    it. A pose is the lab position of the board's (0,0) corner together with
    the lab directions of its two sides.
    '''
    e1 = i_axis/norm(i_axis)
    e2 = j_axis/norm(j_axis)
    n = cross(e1, e2)
    n = n/norm(n)
    return (origin + translation*n, e1, e2)


def board_to_lab(u, v, pose):
    '''
    The lab coordinates of the point of the board plane at the board
    coordinates (u, v).
    '''
    P0, e1, e2 = pose
    return P0 + u*e1 + v*e2


def rotation(axis, angle):
    '''
    The rotation matrix about an axis by an angle, by Rodrigues' formula.
    '''
    from numpy import eye, sin, cos, dot
    k = array(axis, dtype=float)
    k = k/norm(k)
    K = array([[0.0, -k[2], k[1]],
               [k[2], 0.0, -k[0]],
               [-k[1], k[0], 0.0]])
    return eye(3) + sin(angle)*K + (1.0 - cos(angle))*dot(K, K)


def moving_poses(cams, resolutions):
    '''
    A set of poses for a board that is turned about freely, all of them
    visible in every camera.

    The board is turned as well as moved, because the intrinsics are recovered
    from the way the perspective changes with orientation. Poses in which the
    board would fall outside any camera's frame are rejected and drawn again,
    since moving_board_calibration can only use frames that several cameras
    saw.

    output - a list of poses, each the lab position of the board's (0,0)
             corner together with the lab directions of its two sides
    '''
    from numpy.random import default_rng

    rng = default_rng(seed)
    ni, nj = board_size
    s = square_size
    half = array([(ni - 1)/2.0*s, (nj - 1)/2.0*s, 0.0])

    # the corners of the pattern, in the board's own frame, which are what
    # must stay inside the frame
    probe = [(-s, -s), ((ni)*s, -s), ((ni)*s, (nj)*s), (-s, (nj)*s)]

    poses, tries = [], 0
    while len(poses) < n_frames and tries < 200*n_frames:
        tries += 1

        ax = rng.normal(0, 1, 3)
        R = rotation(ax, rng.uniform(0.05, max_tilt))
        e1 = R.dot(i_axis/norm(i_axis))
        e2 = R.dot(j_axis/norm(j_axis))

        centre = volume_centre + rng.uniform(-wander, wander)
        P0 = centre - (half[0]*e1 + half[1]*e2)
        pose = (P0, e1, e2)

        ok = True
        for cam, res in zip(cams, resolutions):
            for (u, v) in probe:
                eta, zeta = cam.projection(board_to_lab(u, v, pose))
                if not (0 < eta < res[0] and 0 < zeta < res[1]):
                    ok = False
                    break
            if not ok:
                break

        if ok:
            poses.append(pose)

    if len(poses) < n_frames:
        print('  only %d of the %d poses asked for fit in every camera; '
              'reduce wander or max_tilt'%(len(poses), n_frames))

    return poses


def render(cam, resolution, pose):
    '''
    Draws the board as seen by cam, in the given pose. Returns a PIL image.
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
        eta, zeta = cam.projection(board_to_lab(u, v, pose))
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

    if mode not in ('stage', 'moving'):
        raise ValueError('mode must be "stage" or "moving", not "%s"'%mode)

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    cams, resolutions = [], []
    for name in camera_names:
        cw = camera_wrapper(name, camera_dir)
        cw.load()
        cams.append(cw)
        resolutions.append(cw.camera.resolution)

    print('generating checkerboard images for %s, in "%s" mode\n'
          %(str(camera_names), mode))

    # ---------------------------------------------------------- stage
    if mode == 'stage':

        for name, cam, res in zip(camera_names, cams, resolutions):
            for e, t in enumerate(translations):
                img = render(cam, res, stage_pose(t))
                img.save(os.path.join(output_dir,
                                      '%s_board_%02d.tif'%(name, e)))

            # the position of the board origin in this camera, which is what
            # the origin_hint parameter needs
            eta, zeta = cam.projection(
                board_to_lab(0.0, 0.0, stage_pose(translations[0])))
            print('%s: wrote %d images; origin_hint: %.0f, %.0f'
                  %(name, len(translations), eta, zeta))

        print('\nAdd a block like this to your parameters file, one per '
              'camera:\n')
        print('- checkerboard_calibration:')
        print('    camera_name: %s'%camera_names[0])
        print('    images: %s/%s_board_*.tif'%(output_dir, camera_names[0]))
        print('    board_size: %d, %d'%board_size)
        print('    square_size: %s'%square_size)
        print('    board_origin: %.4f, %.4f, %.4f'%tuple(origin))
        print('    board_i_axis: %.0f, %.0f, %.0f'%tuple(i_axis))
        print('    board_j_axis: %.0f, %.0f, %.0f'%tuple(j_axis))
        print('    translations: %s'%(', '.join(str(t)
                                                for t in translations)))
        print('    origin_hint: <the two numbers printed above for that '
              'camera>')

    # ---------------------------------------------------------- moving
    else:

        poses = moving_poses(cams, resolutions)

        for name, cam, res in zip(camera_names, cams, resolutions):
            d = os.path.join(output_dir, name)
            if not os.path.isdir(d):
                os.makedirs(d)
            for e, pose in enumerate(poses):
                img = render(cam, res, pose)
                img.save(os.path.join(d, 'frame%04d.tif'%(e + 1)))
            print('%s: wrote %d images into %s'%(name, len(poses), d))

        print('\nThe same frame shows the same pose of the board in every '
              'camera,\nwhich is what moving_board_calibration needs. Add '
              'this to your\nparameters file:\n')
        print('- moving_board_calibration:')
        print('    camera_names: %s'%(', '.join(camera_names)))
        print('    images: %s/{camera}/*.tif'%output_dir)
        print('    resolution: %d, %d'%(resolutions[0][0], resolutions[0][1]))
        print('    3D_model: Tsai')
        print('    board_size: %d, %d'%board_size)
        print('    square_size: %s'%square_size)
        print('    every_nth_frame: 1')
        print('    max_frames: %d'%len(poses))
        print('    reference_camera: %s'%camera_names[0])
        print('    n_distortion: 2')
        print('    sigma: 2.0')
        print('    min_distance: 10')
        print('    min_score: 0.75')
        print('    max_iterations: 1500')
        print('    save_folder: %s'%output_dir)
        print('\nNote that the lab frame this recovers is anchored on %s, '
              'because\nnothing in the images says where the origin is. The '
              'cameras that made\nthese images are in the frame of the '
              'example, so the two will differ by\na rigid motion; distances '
              'and shapes are what should agree.'%camera_names[0])
