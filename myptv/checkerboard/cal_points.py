# -*- coding: utf-8 -*-
#!/usr/bin/env python3

"""
Created on Aug 27 2026

@author: Stefano Brizzolara


Generating MyPTV calibration points from images of a checkerboard.

This module ties together the corner detection of detection.py, the lattice
assembly of grid.py and the board geometry of board.py into the procedure
that a user actually runs. Its output is a calibration points file in the
format that MyPTV already uses, namely one row per point holding

    eta    zeta    x_lab    y_lab    z_lab

with eta and zeta the pixel coordinates of the corner. Such a file is what
myptv.utils.Cal_image_coord reads, and hence what the existing camera models
calibrate against. The checkerboard procedure is therefore a replacement for
the manual point picking stage only: it does not touch the camera models, and
a file it produces can be fed to the Tsai model, the extended Zolof model, or
any model added later, without change.


The intended experimental procedure is the one common in 3D-PTV:

    1) mount a checkerboard on a translation stage inside the measurement
       volume, oriented so that its two sides lie along known lab directions;

    2) photograph it with all the cameras simultaneously at a number of known
       stage positions, so that the corners span the depth of the volume;

    3) run this procedure once per camera, over that camera's set of images.

Because the cameras record the board at the same instant, they share one lab
frame by construction, and because the stage positions are known the depth
coordinate is known to the accuracy of the stage.

A word on the orientation hint. The lattice found in an image can be laid
over the physical board in several ways, and choosing differently in two
cameras assigns different lab coordinates to the same physical corner, which
ruins the calibration without any obvious symptom. The procedure therefore
asks for the approximate pixel position of the board origin in each camera,
which removes the ambiguity; see board.resolve_orientation. It is worth
plotting the result with plot_detection before calibrating.

"""

import os
from glob import glob

from numpy import array, zeros, isfinite

from myptv.checkerboard.detection import load_image, detect_corners
from myptv.checkerboard.grid import assemble_grid
from myptv.checkerboard.board import (resolve_orientation,
                                      board_lab_coordinates)




class checkerboard_cal_points(object):
    '''
    Finds the checkerboard corners in a set of images of one camera and
    converts them into a MyPTV calibration points file.

    Typical use:

        ch = checkerboard_cal_points(
                 './Calibration/cam1_board_*.tif',
                 board_size = (9, 6),
                 square_size = 5.0,
                 translations = [-10.0, -5.0, 0.0, 5.0, 10.0],
                 origin_hint = (212.0, 845.0))
        ch.process()
        ch.save('./Calibration/cam1_cal_points')
        ch.plot_detection(0)
    '''

    def __init__(self, images, board_size, square_size,
                 origin=(0., 0., 0.), i_axis=(1., 0., 0.),
                 j_axis=(0., 1., 0.), translations=None,
                 origin_hint=None, iaxis_hint=None,
                 sigma=2.0, min_distance=10, min_score=0.7, tol=0.35,
                 print_progress=True):
        '''
        input -

        images - either a list of paths of the images of the board, or a
                 single string holding a glob pattern such as
                 './Calibration/cam1_board_*.tif'. When a pattern is given
                 the matching files are used in sorted order, so it is worth
                 naming them with zero padded numbers.

        board_size (tuple, 2) - the number of INTERNAL corners of the board
                 along its two sides. A board of 10 by 7 squares has 9 by 6
                 internal corners, so board_size would be (9, 6). Only the
                 internal corners are used, since the corners on the outer
                 boundary of the pattern are not surrounded by the board on
                 all sides and cannot be located reliably.

        square_size (float) - the physical spacing of the corners, in the
                 length units of the lab coordinates. A pair of numbers may
                 be given if the spacing differs along the two sides.

        origin (array, 3) - the lab coordinates of the board corner (0,0),
                 for a board at zero translation.

        i_axis, j_axis (arrays, 3) - the lab space directions of the two
                 sides of the board.

        translations (list) - one entry per image, giving the displacement of
                 the board when that image was taken. Each entry is either a
                 number, taken as a distance along the board normal, or a
                 vector of three components used as it stands. If None, all
                 the images are taken to show the board at zero translation,
                 which is only meaningful for a single image.

        origin_hint (array, 2) - the approximate pixel position, in this
                 camera, of the board corner chosen as the origin. Needed to
                 fix the orientation of the board unambiguously; see the
                 module docstring.

        iaxis_hint (array, 2) - the approximate pixel position of a corner
                 lying along the first axis of the board. Needed only for a
                 square board.

        sigma (float) - the smoothing scale of the corner detector, in
                 pixels; of the order of the blur of the board edges.

        min_distance (int) - the minimum separation of two corners, in
                 pixels. It must be smaller than the spacing of the corners
                 in the image, and is also used as the size of the windows
                 in which the corners are refined.

        min_score (float) - the threshold on the corner score below which a
                 candidate is discarded.

        tol (float) - the tolerance, in units of the local corner spacing,
                 within which a corner must fall of its predicted position
                 to be added to the lattice.

        print_progress (bool) - whether to report on each image as it is
                 processed.
        '''
        if isinstance(images, str):
            self.image_files = sorted(glob(images))
            if len(self.image_files) == 0:
                raise ValueError('the pattern "%s" matched no files'%images)
        else:
            self.image_files = list(images)
            if len(self.image_files) == 0:
                raise ValueError('no images were given')

        for f in self.image_files:
            if not os.path.exists(f):
                raise ValueError('the image %s does not exist'%f)

        self.board_size = (int(board_size[0]), int(board_size[1]))
        if min(self.board_size) < 3:
            raise ValueError('board_size must be at least 3 by 3 internal '
                             'corners for the lattice to be grown reliably')

        self.square_size = square_size
        self.origin = origin
        self.i_axis = i_axis
        self.j_axis = j_axis

        if translations is None:
            translations = [0.0]*len(self.image_files)
        if len(translations) != len(self.image_files):
            raise ValueError('%d translations were given for %d images; '
                             'there must be exactly one per image'
                             %(len(translations), len(self.image_files)))
        self.translations = list(translations)

        self.origin_hint = origin_hint
        self.iaxis_hint = iaxis_hint

        self.sigma = float(sigma)
        self.min_distance = int(min_distance)
        self.min_score = float(min_score)
        self.tol = float(tol)
        self.print_progress = bool(print_progress)

        # filled in by process()
        self.rows = []
        self.results = []


    def __repr__(self):
        return ('checkerboard_cal_points instance; %d images, board %dx%d, '
                '%d points found'
                %(len(self.image_files), self.board_size[0],
                  self.board_size[1], len(self.rows)))


    def process_image(self, fname, translation, origin_hint=None,
                      iaxis_hint=None):
        '''
        Runs the detection, the lattice assembly and the lab coordinate
        assignment on a single image.

        input -
        fname (string) - path of the image
        translation - the displacement of the board for this image
        origin_hint, iaxis_hint (arrays, 2) - see the class docstring; these
                    override the values given to the constructor, which lets
                    process() carry the orientation over from one image to
                    the next.

        output - a dictionary describing what happened, with the keys:
        'file', 'ok' (bool), 'message' (string), 'n_points' (int),
        'points' (array, (N,2)) - the pixel positions of the used corners,
        'lab' (array, (N,3)) - their lab coordinates,
        'grid' (array of ints) - the oriented lattice, or None,
        'detection' (dict) - the raw output of detect_corners
        '''
        out = {'file': fname, 'ok': False, 'message': '', 'n_points': 0,
               'points': zeros((0, 2)), 'lab': zeros((0, 3)),
               'grid': None, 'detection': None}

        img = load_image(fname)

        det = detect_corners(img, sigma=self.sigma,
                             min_distance=self.min_distance,
                             min_score=self.min_score)
        out['detection'] = det

        ni, nj = self.board_size
        n_expected = ni*nj

        if len(det['points']) < 4:
            out['message'] = ('only %d corners were detected; try lowering '
                              'min_score or adjusting sigma'
                              %len(det['points']))
            return out

        idx = assemble_grid(det, tol=self.tol)

        if idx.size == 0:
            out['message'] = ('the detected corners could not be assembled '
                              'into a lattice')
            return out

        try:
            oriented = resolve_orientation(det['points'], idx,
                                           self.board_size,
                                           origin_hint=origin_hint,
                                           iaxis_hint=iaxis_hint)
        except ValueError as e:
            out['message'] = str(e)
            return out

        out['grid'] = oriented

        X = board_lab_coordinates(self.board_size, self.square_size,
                                  origin=self.origin, i_axis=self.i_axis,
                                  j_axis=self.j_axis,
                                  translation=translation)

        pts, lab = [], []
        for i in range(ni):
            for j in range(nj):
                k = oriented[i, j]
                if k < 0:          # a hole in the lattice
                    continue
                p = det['points'][k]
                if not isfinite(p).all():
                    continue
                pts.append(p)
                lab.append(X[i, j])

        out['points'] = array(pts) if len(pts) else zeros((0, 2))
        out['lab'] = array(lab) if len(lab) else zeros((0, 3))
        out['n_points'] = len(pts)
        out['ok'] = len(pts) > 0
        out['message'] = ('%d of the %d corners of the board were used'
                          %(len(pts), n_expected))

        return out


    def process(self):
        '''
        Processes every image and collects the calibration points.

        The orientation hint is carried over from one image to the next: the
        first image is oriented using the hint given to the constructor, and
        each following image is oriented using the position that the board
        origin was found at in the previous image. Since the board only
        translates between images this is both more robust than reusing a
        fixed hint and a check on consistency.

        output -
        rows (list) - the collected [eta, zeta, x, y, z] rows
        '''
        self.rows = []
        self.results = []

        hint = self.origin_hint
        ihint = self.iaxis_hint

        for fname, t in zip(self.image_files, self.translations):

            res = self.process_image(fname, t, origin_hint=hint,
                                     iaxis_hint=ihint)
            self.results.append(res)

            if self.print_progress:
                status = 'ok  ' if res['ok'] else 'FAIL'
                print('  [%s] %s: %s'%(status, os.path.basename(fname),
                                       res['message']))

            if not res['ok']:
                continue

            # carry the orientation over to the next image
            g = res['grid']
            if g is not None:
                if g[0, 0] >= 0:
                    hint = res['detection']['points'][g[0, 0]]
                if g[-1, 0] >= 0:
                    ihint = res['detection']['points'][g[-1, 0]]

            for p, X in zip(res['points'], res['lab']):
                self.rows.append([p[0], p[1], X[0], X[1], X[2]])

        n_ok = sum(1 for r in self.results if r['ok'])

        if self.print_progress:
            print('\n  %d of %d images used; %d calibration points in total'
                  %(n_ok, len(self.image_files), len(self.rows)))

        if len(self.rows) == 0:
            raise RuntimeError(
                'no calibration points were found in any of the images. The '
                'most common causes are a wrong board_size (it counts '
                'INTERNAL corners, so a board of 10x7 squares is (9,6)), a '
                'min_distance larger than the spacing of the corners in the '
                'image, or a sigma that does not match the blur of the '
                'board edges.')

        return self.rows


    def save(self, fname):
        '''
        Writes the calibration points to the disk, in the tab separated
        format that myptv.utils.Cal_image_coord reads.

        input -
        fname (string) - the path to write to, e.g. './Calibration/cam1_cal_points'
        '''
        if len(self.rows) == 0:
            raise RuntimeError('there are no calibration points to save; '
                               'call process() first')

        with open(fname, 'w') as f:
            for r in self.rows:
                f.write('%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n'
                        %(r[0], r[1], r[2], r[3], r[4]))

        if self.print_progress:
            print('  wrote %d points to %s'%(len(self.rows), fname))


    def plot_detection(self, i=0, ax=None, annotate=True):
        '''
        Plots the corners found in one image over that image, drawing the
        rows and columns of the lattice and marking the board origin. This is
        the check that the orientation of the board came out as intended,
        and it is worth looking at one such plot per camera before
        calibrating.

        input -
        i (int) - which of the images to plot
        ax - an existing matplotlib axis to draw on; a new figure is made if
             this is None
        annotate (bool) - whether to label the origin and the two axes

        output -
        ax - the axis that was drawn on
        '''
        from matplotlib.pyplot import subplots

        if len(self.results) == 0:
            raise RuntimeError('call process() before plotting')

        res = self.results[i]

        if ax is None:
            fig, ax = subplots()

        ax.imshow(load_image(res['file']), cmap='gray')
        ax.set_title(os.path.basename(res['file']))

        det = res['detection']
        if det is not None and len(det['points']):
            ax.plot(det['points'][:, 0], det['points'][:, 1], 'r.',
                    ms=3, label='detected corners')

        g = res['grid']
        if g is None:
            ax.legend(loc='upper right', fontsize=7)
            return ax

        P = det['points']

        # the rows and the columns of the lattice
        for a in range(g.shape[0]):
            ks = [k for k in g[a, :] if k >= 0]
            if len(ks) > 1:
                ax.plot(P[ks, 0], P[ks, 1], '-', lw=0.7, color='deepskyblue')
        for b in range(g.shape[1]):
            ks = [k for k in g[:, b] if k >= 0]
            if len(ks) > 1:
                ax.plot(P[ks, 0], P[ks, 1], '-', lw=0.7, color='gold')

        if annotate:
            if g[0, 0] >= 0:
                p = P[g[0, 0]]
                ax.plot(p[0], p[1], 'o', mfc='none', mec='lime', ms=12,
                        mew=2, label='board origin (0,0)')
            if g[-1, 0] >= 0:
                p = P[g[-1, 0]]
                ax.annotate('i', xy=(p[0], p[1]), color='lime',
                            fontsize=12, fontweight='bold')
            if g[0, -1] >= 0:
                p = P[g[0, -1]]
                ax.annotate('j', xy=(p[0], p[1]), color='lime',
                            fontsize=12, fontweight='bold')

        ax.legend(loc='upper right', fontsize=7)

        return ax
