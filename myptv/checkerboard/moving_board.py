# -*- coding: utf-8 -*-
#!/usr/bin/env python3

"""
Created on Aug 27 2026

@author: Stefano Brizzolara


Calibrating from a checkerboard that was moved freely.

cal_points.py handles the case where the board sits on a translation stage,
so that its position is read off the stage and the lab coordinates of every
corner are known before any image is looked at. This module handles the other
common case: the board is held in the measurement volume and turned about
between exposures, so that nothing is known about where it was. Its pose in
each frame then has to be recovered from the images, together with the
cameras.

What it needs in exchange is that the cameras are triggered together, so a
given frame shows one pose of the board seen from several directions, and
that the board really is turned about between exposures rather than merely
slid around: the intrinsics are recovered from the way the perspective
changes with orientation, and images that differ only by a translation carry
the same information over and over.

The steps are

    1) find the corners in every image                       (detection.py)
    2) order them into a lattice                             (grid.py)
    3) recover each camera's intrinsics from its homographies (pose.py)
    4) settle which way round each lattice is                 (pose.py, below)
    5) adjust cameras and board poses together against the
       reprojection error                                     (bundle.py)
    6) write out a calibration points file per camera

and the output of step 6 is an ordinary MyPTV calibration points file, so the
3D models are fitted from it afterwards exactly as they would be from hand
picked points, with the initial_calibration and final_calibration actions.


About the lab frame.

There is nothing in a set of images of a board being waved about that says
where the origin of the lab frame is or which way its axes point. The scale
is fixed, by the physical size of the squares, but the frame itself is not:
this module declares one camera to be the origin and lets everything else
follow. Distances, angles and shapes are therefore right, and the lab frame
is a convention.

For a PTV experiment that usually matters. If the coordinates have to line up
with something physical, a tank wall or a traverse, then the frame has to be
related to the apparatus afterwards, by measuring a few known points, or the
board should be placed deliberately for one exposure and that frame used as
the reference. What must not be done is to assume the frame means anything on
its own.

"""

import os
from glob import glob

from numpy import (array, zeros, eye, concatenate, linspace, median,
                   degrees)
from numpy.linalg import norm

from myptv.checkerboard.detection import load_image, detect_corners
from myptv.checkerboard.grid import assemble_grid
from myptv.checkerboard.pose import (fit_homography, homography_residual,
                                     zhang_intrinsics, pose_from_homography,
                                     relative_pose, rotation_angle,
                                     nearest_rotation, theta_from_R)
from myptv.checkerboard.bundle import multicamera_bundle




class moving_board_calibration(object):
    '''
    Recovers several cameras and the pose of a freely moved checkerboard, and
    writes MyPTV calibration points.

    Typical use:

        cal = moving_board_calibration(
                  {'cam1': './cal/cam1/*.tif',
                   'cam2': './cal/cam2/*.tif',
                   'cam3': './cal/cam3/*.tif'},
                  board_size = (8, 11),
                  square_size = 15.0)
        cal.detect()
        cal.solve()
        cal.report()
        cal.save_cal_points('./Calibration')
    '''

    def __init__(self, images, board_size, square_size,
                 sigma=2.0, min_distance=10, min_score=0.75,
                 n_dist=2, reference_camera=None, max_frames=60,
                 agreement_deg=2.0, print_progress=True):
        '''
        input -

        images (dict) - camera name -> either a list of image paths or a glob
                 pattern. The frames of different cameras are matched by
                 their position in the sorted list, so the cameras must have
                 the same number of images, in the same order, and the
                 exposures must correspond.

        board_size (tuple, 2) - the number of INTERNAL corners of the board
                 along its two sides; a board of 9 by 12 squares is (8, 11).

        square_size (float) - the physical spacing of the corners, which is
                 what fixes the scale of the whole result.

        sigma, min_distance, min_score - passed to detection.detect_corners.

        n_dist (int) - how many distortion coefficients to fit, 0, 2 or 4.

        reference_camera (string) - the camera whose pose defines the lab
                 frame. Defaults to the first.

        max_frames (int) - how many frames to use in the adjustment. All the
                 frames are detected, but a calibration does not need
                 hundreds of views of the same board and the adjustment cost
                 grows with them, so a well spread subset is taken.

        agreement_deg (float) - how closely the transformation between two
                 cameras must repeat from frame to frame for a frame's
                 labelling to be accepted.

        print_progress (bool) - whether to report as it goes.
        '''
        self.cams = list(images.keys())
        if len(self.cams) < 2:
            raise ValueError('at least two cameras are needed; a single '
                             'camera cannot fix a lab frame')

        self.files = {}
        for c, v in images.items():
            fl = sorted(glob(v)) if isinstance(v, str) else list(v)
            if len(fl) == 0:
                raise ValueError('no images for camera %s'%c)
            self.files[c] = fl

        n = set(len(v) for v in self.files.values())
        if len(n) != 1:
            raise ValueError('the cameras have different numbers of images '
                             '(%s); the frames are matched by their order, '
                             'so they must correspond'
                             %{c: len(v) for c, v in self.files.items()})

        self.board_size = (int(board_size[0]), int(board_size[1]))
        if min(self.board_size) < 3:
            raise ValueError('board_size must be at least 3 by 3 internal '
                             'corners')
        self.square_size = float(square_size)
        if self.square_size <= 0:
            raise ValueError('square_size must be positive')

        self.sigma = float(sigma)
        self.min_distance = int(min_distance)
        self.min_score = float(min_score)
        self.n_dist = int(n_dist)
        self.max_frames = int(max_frames)
        self.agreement = float(agreement_deg)
        self.ref = self.cams[0] if reference_camera is None \
            else reference_camera
        if self.ref not in self.cams:
            raise ValueError('reference_camera %s is not one of the cameras'
                             %self.ref)
        self.print_progress = bool(print_progress)

        ni, nj = self.board_size
        s = self.square_size
        B = array([[i*s, j*s, 0.0] for i in range(ni) for j in range(nj)])
        self.B = B - B.mean(axis=0)      # the board's own frame is arbitrary
        self.NPT = ni*nj

        self.grids = {}
        self.bundle = None


    def __repr__(self):
        return ('moving_board_calibration: %d cameras, %d frames, %d boards '
                'found'%(len(self.cams), len(self.files[self.ref]),
                         len(self.grids)))


    # ------------------------------------------------------------------
    def detect(self):
        '''
        Finds and orders the board in every image. Images in which it is not
        found completely are simply left out.

        output -
        grids (dict) - (camera, frame) -> (n_i, n_j, 2) corner positions
        '''
        ni, nj = self.board_size
        self.grids = {}
        self.detect_failures = {}

        for c in self.cams:
            ok = 0
            for k, fname in enumerate(self.files[c]):
                try:
                    img = load_image(fname)
                except Exception as e:
                    self.detect_failures[(c, k)] = str(e)
                    continue

                det = detect_corners(img, sigma=self.sigma,
                                     min_distance=self.min_distance,
                                     min_score=self.min_score)
                if len(det['points']) < 4:
                    self.detect_failures[(c, k)] = \
                        'only %d corners'%len(det['points'])
                    continue

                g = assemble_grid(det)
                if sorted(g.shape) != sorted(self.board_size):
                    self.detect_failures[(c, k)] = 'lattice %s'%(g.shape,)
                    continue
                if (g >= 0).sum() != self.NPT:
                    self.detect_failures[(c, k)] = \
                        'incomplete, %d of %d'%(int((g >= 0).sum()), self.NPT)
                    continue

                gg = g.T if g.shape != self.board_size else g
                self.grids[(c, k)] = det['points'][gg]
                ok += 1

            if self.print_progress:
                print('  %s: the board was found in %d of %d images'
                      %(c, ok, len(self.files[c])))

        if len(self.grids) == 0:
            raise RuntimeError(
                'the board was not found in any image. Check board_size, '
                'which counts INTERNAL corners, and min_distance, which must '
                'be smaller than the spacing of the corners in the image.')

        return self.grids


    # ------------------------------------------------------------------
    def _intrinsics(self):
        '''
        Zhang's estimate for each camera, reduced to the form the projection
        model actually uses: one focal length, no skew. Keeping the skew and
        the two separate focal lengths here and then projecting without them
        would make the starting point inconsistent with the model by however
        large they are.
        '''
        K0, f0, hres = {}, {}, {}
        BP = self.B[:, :2]

        for c in self.cams:
            Hs, rr = [], []
            for (cc, k), P in self.grids.items():
                if cc != c:
                    continue
                uv = P.reshape(-1, 2)
                H = fit_homography(BP, uv)
                Hs.append(H)
                rr.append(float((homography_residual(H, BP, uv)**2).mean())
                          ** 0.5)

            if len(Hs) < 3:
                raise RuntimeError(
                    'camera %s has only %d usable images; at least three, in '
                    'different orientations, are needed to recover its '
                    'intrinsics'%(c, len(Hs)))

            K, _ = zhang_intrinsics(Hs)
            f = 0.5*(K[0, 0] + K[1, 1])
            f0[c] = f
            K0[c] = array([[f, 0.0, K[0, 2]],
                           [0.0, f, K[1, 2]],
                           [0.0, 0.0, 1.0]])
            hres[c] = float(median(array(rr)))

        self.K0, self.f0, self.homography_residual_px = K0, f0, hres
        return K0


    def _candidates(self, c, k):
        '''
        The relabellings of the lattice that put the board in front of the
        camera with its face towards it. Two of the four survive; see the
        module docstring of pose.py.
        '''
        BP = self.B[:, :2]
        P = self.grids[(c, k)]
        out = []
        for flip in range(4):
            Q = P[::-1] if flip & 1 else P
            Q = Q[:, ::-1] if flip & 2 else Q
            uv = Q.reshape(-1, 2)
            R, t, _ = pose_from_homography(fit_homography(BP, uv), self.K0[c])
            if R[:, 2].dot(t) > 0:
                out.append((flip, R, t, uv))
        return out


    def _label(self):
        '''
        Chooses one relabelling per image so that every camera calls the same
        physical corner by the same indices.

        The transformation between two cameras does not change from frame to
        frame, and it does not care how the board is labelled provided both
        cameras are labelled the same way. So the reference camera is pinned
        arbitrarily at each frame, and the others are chosen to reproduce one
        transformation, which is found by seeing which candidate is supported
        by the most frames.
        '''
        frames = sorted(set(k for (c, k) in self.grids))
        cand = {}
        for (c, k) in self.grids:
            cand[(c, k)] = self._candidates(c, k)

        REL, support = {}, {}
        for c in self.cams:
            if c == self.ref:
                continue
            pool = [k for k in frames
                    if (self.ref, k) in cand and (c, k) in cand]
            if len(pool) < 3:
                raise RuntimeError(
                    'cameras %s and %s share only %d usable frames; they '
                    'cannot be tied together'%(self.ref, c, len(pool)))

            best = None
            for k0 in pool[:10]:
                _, R0, t0, _ = cand[(self.ref, k0)][0]
                for (_, R, t, _) in cand[(c, k0)]:
                    Rh, _th = relative_pose(R0, t0, R, t)
                    sup = 0
                    for k in pool:
                        _, Ra, ta, _ = cand[(self.ref, k)][0]
                        d = [degrees(rotation_angle(Rh.T.dot(
                                relative_pose(Ra, ta, Rb, tb)[0])))
                             for (_, Rb, tb, _) in cand[(c, k)]]
                        if d and min(d) < self.agreement:
                            sup += 1
                    if best is None or sup > best[0]:
                        best = (sup, Rh, len(pool))
            REL[c] = best[1]
            support[c] = (best[0], best[2])

        obs, pose0 = {}, {}
        for k in frames:
            if (self.ref, k) not in cand:
                continue
            _, Ra, ta, Qa = cand[(self.ref, k)][0]
            sel = {self.ref: Qa}
            for c in self.cams:
                if c == self.ref or (c, k) not in cand:
                    continue
                ds = sorted([(degrees(rotation_angle(REL[c].T.dot(
                                relative_pose(Ra, ta, Rb, tb)[0]))), Qb)
                             for (_, Rb, tb, Qb) in cand[(c, k)]],
                            key=lambda x: x[0])
                if ds[0][0] < self.agreement:
                    sel[c] = ds[0][1]
            if len(sel) < 2:
                continue
            for c, Q in sel.items():
                obs[(c, k)] = Q
            pose0[k] = (Ra, ta)

        self.labelling_support = support
        return obs, pose0


    # ------------------------------------------------------------------
    def solve(self, max_nfev=2000, loss='soft_l1', f_scale=1.0):
        '''
        Recovers the cameras and the board poses.

        input -
        max_nfev (int) - the limit on iterations of the adjustment
        loss (string) - 'soft_l1' by default, so that a single image whose
                        corners were labelled wrongly cannot drag the whole
                        solution; 'linear' for a plain least squares
        f_scale (float) - the scale of that robust loss, in pixels

        output -
        result - the scipy result object of the adjustment
        '''
        if len(self.grids) == 0:
            self.detect()

        self._intrinsics()
        obs, pose0 = self._label()

        if len(pose0) < 3:
            raise RuntimeError('only %d frames could be labelled '
                               'consistently'%len(pose0))

        frames = sorted(pose0)
        if len(frames) > self.max_frames:
            keep = set(array(frames)[linspace(
                0, len(frames)-1, self.max_frames).astype(int)])
            frames = sorted(keep)
            obs = {(c, k): v for (c, k), v in obs.items() if k in keep}
            pose0 = {k: v for k, v in pose0.items() if k in keep}

        if self.print_progress:
            print('  %d frames labelled consistently, %d used, %d '
                  'camera-frame observations'
                  %(len(pose0), len(frames), len(obs)))

        # the reference camera is the lab frame, so the board poses that came
        # out of its own homographies are already in that frame
        ext0 = {self.ref: (eye(3), zeros(3))}
        BP = self.B[:, :2]
        for c in self.cams:
            if c == self.ref:
                continue
            Rs, ts = [], []
            for k in frames:
                if (c, k) not in obs:
                    continue
                Ra, ta = pose0[k]
                Rb, tb, _ = pose_from_homography(
                    fit_homography(BP, obs[(c, k)]), self.K0[c])
                Rr, tr = relative_pose(Ra, ta, Rb, tb)
                Rs.append(Rr); ts.append(tr)
            ext0[c] = (nearest_rotation(sum(Rs)/len(Rs)),
                       sum(ts)/len(ts))

        intr0 = {c: concatenate([[self.f0[c], self.K0[c][0, 2],
                                  self.K0[c][1, 2]], zeros(self.n_dist)])
                 for c in self.cams}

        self.bundle = multicamera_bundle(self.B, obs, self.cams, intr0, ext0,
                                         pose0, n_dist=self.n_dist,
                                         fix_camera=self.ref)

        if self.print_progress:
            print('  starting the adjustment from %.3f px'
                  %self.bundle.rms())

        res = self.bundle.solve(max_nfev=max_nfev, loss=loss,
                                f_scale=f_scale)

        if self.print_progress:
            print('  finished at %.4f px after %d iterations'
                  %(self.bundle.rms(), res.nfev))

        self.result = res
        return res


    # ------------------------------------------------------------------
    def save_cal_points(self, folder, suffix='_cal_points',
                        with_view_index=True):
        '''
        Writes one calibration points file per camera, in the format that
        myptv.utils.Cal_image_coord reads, so that the 3D models can be
        fitted from them in the usual way.

        input -
        folder (string) - where to write
        suffix (string) - appended to the camera name
        with_view_index (bool) - whether to add a sixth column saying which
                    frame each point came from. Cal_image_coord reads only
                    the first five columns, so the file remains an ordinary
                    calibration points file and every existing use of it is
                    unaffected; but knowing which points share a view is what
                    allows the error to be reported per view afterwards,
                    which is how a single badly labelled image is found. See
                    myptv.makePlots.plot_calibration.

        output -
        paths (dict) - camera -> the file written
        '''
        if self.bundle is None:
            raise RuntimeError('call solve() before saving')

        lab = self.bundle.lab_coordinates()
        rows = {c: [] for c in self.cams}

        for n, (ci, k) in enumerate(self.bundle.pairs):
            c = self.cams[ci]
            uv = self.bundle.obs[n]
            X = lab[k]
            for i in range(self.NPT):
                rows[c].append((uv[i, 0], uv[i, 1], X[i, 0], X[i, 1],
                                X[i, 2], k))

        if not os.path.isdir(folder):
            os.makedirs(folder)

        paths = {}
        for c in self.cams:
            p = os.path.join(folder, c + suffix)
            with open(p, 'w') as f:
                for r in rows[c]:
                    if with_view_index:
                        f.write('%.3f\t%.3f\t%.4f\t%.4f\t%.4f\t%d\n'%r)
                    else:
                        f.write('%.3f\t%.3f\t%.4f\t%.4f\t%.4f\n'%r[:5])
            paths[c] = p
            if self.print_progress:
                print('  %s: %d points -> %s'%(c, len(rows[c]), p))

        return paths


    def camera_parameters(self):
        '''
        The recovered cameras, in the terms the Tsai model uses, which is
        what to start that model's own fit from.

        output - a dictionary, camera -> a dictionary with the keys
        'O' (the camera position in the lab frame), 'theta' (the three
        angles), 'f', 'xh', 'yh' (the principal point as an offset from the
        centre of the image, which is how camera_Tsai holds it, so the
        resolution has to be added by the caller), and 'distortion'.
        '''
        if self.bundle is None:
            raise RuntimeError('call solve() first')

        intr, ext, _ = self.bundle.unpack()
        out = {}
        for c in self.cams:
            R, t = ext[c]
            out[c] = {'O': -R.T.dot(t),
                      'theta': theta_from_R(R),
                      'f': float(intr[c][0]),
                      'principal_point': (float(intr[c][1]),
                                          float(intr[c][2])),
                      'distortion': intr[c][3:].copy()}
        return out


    def report(self):
        '''
        Prints what was recovered and the checks worth looking at, in
        particular the ones that do not depend on the choice of lab frame.
        '''
        if self.bundle is None:
            raise RuntimeError('call solve() first')

        print('\ncorners found:')
        for c in self.cams:
            n = sum(1 for (cc, _) in self.grids if cc == c)
            print('   %s: %d of %d images, homography residual %.3f px'
                  %(c, n, len(self.files[c]), self.homography_residual_px[c]))
        print('   (the homography residual is what the adjustment is trying '
              'to reach:\n    it is the error of the corners themselves, '
              'plus the lens distortion)')

        print('\nlabelling:')
        for c, (sup, tot) in self.labelling_support.items():
            print('   %s -> %s: %d of %d frames agree'
                  %(self.ref, c, sup, tot))

        print('\nreprojection: %.4f px overall'%self.bundle.rms())
        pv = self.bundle.per_view_rms()
        for c in self.cams:
            v = array([x for cc, _, x in pv if cc == c])
            if len(v):
                print('   %s: median %.4f px over %d views'
                      %(c, median(v), len(v)))
        worst = sorted(pv, key=lambda x: -x[2])[:5]
        print('   worst views: ' +
              ', '.join('%s/%d %.2f' % w for w in worst))
        print('   (a view far worse than the rest usually means its corners '
              'were\n    labelled the wrong way round; drop it and solve '
              'again)')

        pars = self.camera_parameters()
        print('\ncameras:')
        for c in self.cams:
            p = pars[c]
            print('   %s: f=%9.1f  principal point=(%.1f, %.1f)  '
                  'distortion=%s'
                  %(c, p['f'], p['principal_point'][0],
                    p['principal_point'][1],
                    ', '.join('%.4f'%x for x in p['distortion'])))
        print('   positions in the lab frame (%s is the origin):'%self.ref)
        for c in self.cams:
            print('     %s: %s'%(c, ', '.join('%9.1f'%x
                                              for x in pars[c]['O'])))

        print('   distances between them, which no choice of frame changes:')
        for i in range(len(self.cams)):
            for j in range(i+1, len(self.cams)):
                a, b = self.cams[i], self.cams[j]
                print('     %s-%s: %9.1f'
                      %(a, b, norm(pars[a]['O'] - pars[b]['O'])))

        lab = self.bundle.lab_coordinates()
        X = concatenate([lab[k] for k in lab])
        print('\nthe calibration points fill:')
        for i, ax in enumerate('xyz'):
            print('   %s from %9.1f to %9.1f'
                  %(ax, X[:, i].min(), X[:, i].max()))
        print('   Only inside this region is the calibration supported by '
              'data.\n   Make sure it covers the measurement volume.')
