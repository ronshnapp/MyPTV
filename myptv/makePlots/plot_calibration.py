# -*- coding: utf-8 -*-
#!/usr/bin/env python3

"""
Created on Aug 28 2026

@author: Stefano Brizzolara


Looking at a calibration.

A calibration is usually judged by a single number, the mean reprojection
error, and that number hides most of what can go wrong. A camera whose
average is a respectable third of a pixel may be fine everywhere except in
one corner of the frame; a set of views may be excellent apart from one whose
corners were labelled the wrong way round; and a calibration can be tight
everywhere the target went while saying nothing at all about the volume it
never reached. None of that shows up in a mean.

This module draws the things that do show it:

    per camera        how the error is distributed, and how large it is
    per view          the mean error of each image, which is how a single
                      bad image is found
    over the frame    where in the image the error lives, and which way it
                      points, which is how unmodelled distortion is
                      recognised
    in space          where the cameras are, which way they look, and where
                      the target actually went

The last two answer the questions that MATLAB's showReprojectionErrors and
showExtrinsics are usually reached for.

Everything here works from a camera and a list of (lab, image) pairs, which
is what every calibration in MyPTV produces, whether the points were picked
by hand, taken from a target on a translation stage, or recovered from a
board that was moved freely. It does not care which 3D model the camera uses:
positions come from the camera's own projection, and the field of view is
drawn from the rays it returns for the corners of its sensor, so a Tsai
camera and an extended Zolof camera are treated alike.

"""

import os

from numpy import (array, asarray, zeros, linspace, sqrt, median, mean,
                   percentile, unique, concatenate, float64, isfinite, ptp)
from numpy.linalg import norm




def read_cal_points(fname):
    '''
    Reads a MyPTV calibration points file.

    Some of them carry a sixth column saying which view, that is which image,
    each point came from. Cal_image_coord ignores it, but it is what allows
    the error to be reported per view, so it is read here when present.

    input -
    fname (string) - the path of the file

    output -
    img (array, (N,2)) - the pixel coordinates
    lab (array, (N,3)) - the lab coordinates
    view (array of ints, N) - which view each point came from, or None if the
                              file does not say
    '''
    img, lab, view = [], [], []
    with open(fname) as f:
        for ln in f:
            s = ln.split()
            if len(s) < 5:
                continue
            img.append([float(s[0]), float(s[1])])
            lab.append([float(s[2]), float(s[3]), float(s[4])])
            view.append(int(float(s[5])) if len(s) > 5 else -1)

    view = array(view, dtype=int)
    return array(img), array(lab), (None if (view < 0).any() else view)




def group_into_views(lab, tol=None):
    '''
    Splits a set of calibration points into the views they came from.

    The points of one view are the corners of the target in one position, so
    they lie on a plane and close together, while the target moves between
    views. When the file does not record which view each point came from,
    they can still be separated by that: points are put in the same view as
    the point before them unless the jump is large.

    This is a convenience for files that carry no view information. When the
    views are known, pass them in rather than relying on this.

    input -
    lab (array, (N,3)) - the lab coordinates, in the order they were written
    tol (float) - the jump that starts a new view. By default, five times the
                  median step between consecutive points.

    output -
    view (array of ints, N) - which view each point belongs to
    '''
    lab = asarray(lab, dtype=float64)
    if len(lab) < 2:
        return zeros(len(lab), dtype=int)

    step = norm(lab[1:] - lab[:-1], axis=1)
    if tol is None:
        m = median(step[step > 0]) if (step > 0).any() else 1.0
        tol = 5.0*m

    view = zeros(len(lab), dtype=int)
    v = 0
    for i in range(1, len(lab)):
        if step[i-1] > tol:
            v += 1
        view[i] = v

    return view




class calibration_report(object):
    '''
    Holds a set of cameras with the points they were calibrated against, and
    draws them.

    Typical use:

        rep = calibration_report.from_folder('./Calibration', '.',
                                             ['cam1','cam2','cam3'])
        rep.summary()
        rep.plot_error_per_camera()
        rep.plot_error_per_view()
        rep.plot_error_over_frame('cam1')
        rep.plot_extrinsics()
    '''

    def __init__(self, cameras, points, views=None):
        '''
        input -
        cameras (dict) - camera name -> any MyPTV camera, that is anything
                 with a projection(X) method and an O attribute. A
                 camera_wrapper, a camera_Tsai and a camera_extendedZolof all
                 qualify.
        points (dict) - camera name -> (img, lab), the arrays of pixel and
                 lab coordinates of that camera's calibration points.
        views (dict) - camera name -> an integer per point saying which view
                 it came from. If not given, the points are split into views
                 by group_into_views.
        '''
        self.cams = list(cameras.keys())
        self.cameras = cameras
        self.img, self.lab, self.view = {}, {}, {}

        for c in self.cams:
            if c not in points:
                raise ValueError('no calibration points for camera %s'%c)
            im, la = points[c]
            self.img[c] = asarray(im, dtype=float64)
            self.lab[c] = asarray(la, dtype=float64)
            if len(self.img[c]) != len(self.lab[c]):
                raise ValueError('camera %s has %d image points and %d lab '
                                 'points'%(c, len(self.img[c]),
                                           len(self.lab[c])))
            if views is not None and c in views:
                self.view[c] = asarray(views[c], dtype=int)
            else:
                self.view[c] = group_into_views(self.lab[c])

        self._res = {}


    # ------------------------------------------------------------------
    @classmethod
    def from_folder(cls, points_folder, camera_folder, camera_names,
                    suffix='_cal_points'):
        '''
        Builds a report from files on the disk: the camera files and the
        calibration points files that MyPTV writes.

        input -
        points_folder (string) - where the <camera><suffix> files are
        camera_folder (string) - where the camera files are
        camera_names (list) - the cameras to include
        suffix (string) - what is appended to the camera name for the points

        output - a calibration_report
        '''
        from myptv.imaging_mod import camera_wrapper

        cams, pts, views = {}, {}, {}
        for c in camera_names:
            cw = camera_wrapper(c, camera_folder)
            cw.load()
            cams[c] = cw
            p = os.path.join(points_folder, c + suffix)
            if not os.path.exists(p):
                raise ValueError('no calibration points file at %s'%p)
            im, la, vw = read_cal_points(p)
            pts[c] = (im, la)
            if vw is not None:
                views[c] = vw

        return cls(cams, pts, views if views else None)


    def __repr__(self):
        n = sum(len(self.img[c]) for c in self.cams)
        return ('calibration_report: %d cameras, %d points, %d views'
                %(len(self.cams), n, len(self.all_views())))


    # ------------------------------------------------------------------
    def residuals(self, cam):
        '''
        The reprojection error of every point of one camera: the vector, in
        pixels, from where the point was found to where the camera model puts
        it.

        input -
        cam (string) - the camera

        output -
        d (array, (N,2)) - the error vector of each point
        '''
        if cam not in self._res:
            c = self.cameras[cam]
            proj = array([c.projection(X) for X in self.lab[cam]])
            self._res[cam] = proj - self.img[cam]
        return self._res[cam]


    def errors(self, cam):
        '''
        The magnitude of the reprojection error of every point of one camera.
        '''
        return norm(self.residuals(cam), axis=1)


    def all_views(self):
        '''
        Every view index that appears, over all the cameras.
        '''
        v = [self.view[c] for c in self.cams if len(self.view[c])]
        return unique(concatenate(v)) if v else array([], dtype=int)


    def summary(self, printout=True):
        '''
        A table of the error of each camera.

        output - a list of dictionaries, one per camera
        '''
        rows = []
        for c in self.cams:
            e = self.errors(c)
            rows.append({'camera': c,
                         'points': len(e),
                         'views': len(unique(self.view[c])),
                         'mean': float(mean(e)),
                         'median': float(median(e)),
                         'p90': float(percentile(e, 90)),
                         'max': float(e.max()),
                         'rms': float(sqrt((e**2).mean()))})

        if printout:
            print('%-8s %8s %6s %8s %8s %8s %8s %8s'
                  %('camera', 'points', 'views', 'mean', 'median', 'p90',
                    'max', 'rms'))
            for r in rows:
                print('%-8s %8d %6d %8.4f %8.4f %8.4f %8.4f %8.4f'
                      %(r['camera'], r['points'], r['views'], r['mean'],
                        r['median'], r['p90'], r['max'], r['rms']))
            print('\nall in pixels.')

        return rows


    def per_view(self, cam=None, printout=False):
        '''
        The mean error of each view, which is how one bad image is found
        among many good ones.

        input -
        cam (string) - one camera, or None for all of them
        printout (bool) - whether to print the worst few

        output - a list of (camera, view, n_points, mean error) tuples
        '''
        cams = self.cams if cam is None else [cam]
        rows = []
        for c in cams:
            e = self.errors(c)
            v = self.view[c]
            for vi in unique(v):
                m = v == vi
                rows.append((c, int(vi), int(m.sum()), float(mean(e[m]))))

        if printout and rows:
            worst = sorted(rows, key=lambda r: -r[3])[:10]
            print('the worst views:')
            for c, vi, n, me in worst:
                print('   %-8s view %3d  %4d points  %.4f px'%(c, vi, n, me))

        return rows


    # ------------------------------------------------------------------
    def plot_error_per_camera(self, ax=None, bins=40):
        '''
        The distribution of the error of each camera, drawn as overlaid
        histograms with the means marked. Reading them together shows at once
        whether one camera is worse than the rest, and whether a camera has a
        long tail, which a mean would hide.
        '''
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots(figsize=(7, 4))

        hi = max(percentile(self.errors(c), 99.5) for c in self.cams)
        edges = linspace(0, hi, bins)

        for c in self.cams:
            e = self.errors(c)
            ax.hist(e, bins=edges, histtype='step', lw=1.6, label=c)
            ax.axvline(mean(e), ls=':', lw=1.0, alpha=0.6)

        ax.set_xlabel('reprojection error [px]')
        ax.set_ylabel('points')
        ax.set_title('error distribution, by camera'
                     '   (dotted lines are the means)')
        ax.legend()
        return ax


    def plot_error_per_view(self, ax=None, cam=None):
        '''
        The mean error of every view, as bars grouped by camera. This is the
        counterpart of MATLAB's showReprojectionErrors, and it is the plot to
        look at first: a view standing well above its neighbours is almost
        always one whose corners were labelled the wrong way round, and it
        should be removed and the calibration run again.
        '''
        import matplotlib.pyplot as plt

        rows = self.per_view(cam)
        cams = self.cams if cam is None else [cam]

        if ax is None:
            fig, ax = plt.subplots(figsize=(11, 4))

        allv = sorted(set(r[1] for r in rows))
        pos = {v: i for i, v in enumerate(allv)}
        w = 0.8/max(len(cams), 1)

        for k, c in enumerate(cams):
            rr = [r for r in rows if r[0] == c]
            x = [pos[r[1]] + k*w - 0.4 + w/2 for r in rr]
            y = [r[3] for r in rr]
            ax.bar(x, y, width=w, label=c)

        overall = mean([r[3] for r in rows]) if rows else 0.0
        ax.axhline(overall, color='k', ls='--', lw=1.0,
                   label='overall mean %.3f px'%overall)

        ax.set_xticks(range(len(allv)))
        ax.set_xticklabels([str(v) for v in allv], fontsize=7, rotation=90)
        ax.set_xlabel('view')
        ax.set_ylabel('mean error [px]')
        ax.set_title('mean reprojection error per view')
        ax.legend(fontsize=8)
        return ax


    def plot_error_over_frame(self, cam, ax=None, scale=None):
        '''
        Where in the image the error lives, and which way it points, drawn as
        an arrow at each calibration point.

        A cloud of arrows pointing every which way is measurement noise and
        is what a good calibration looks like. Arrows that agree with their
        neighbours, and that grow towards the edges of the frame, are
        unmodelled lens distortion. Arrows that agree over one patch only are
        usually a misplaced point.

        input -
        cam (string) - the camera
        scale (float) - how many pixels of arrow per pixel of error. Chosen
                        automatically if not given.
        '''
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 5.5))

        p = self.img[cam]
        d = self.residuals(cam)
        e = norm(d, axis=1)

        if scale is None:
            span = max(ptp(p[:, 0]), ptp(p[:, 1]))
            scale = 0.06*span/max(median(e), 1e-9)

        q = ax.quiver(p[:, 0], p[:, 1], d[:, 0]*scale, d[:, 1]*scale,
                      e, angles='xy', scale_units='xy', scale=1.0,
                      cmap='viridis', width=0.003)
        plt.colorbar(q, ax=ax, label='error [px]')

        res = self._resolution(cam)
        if res is not None:
            ax.add_patch(plt.Rectangle((0, 0), res[0], res[1], fill=False,
                                       ec='0.6', lw=1.0))
            ax.set_xlim(-0.02*res[0], 1.02*res[0])
            ax.set_ylim(1.02*res[1], -0.02*res[1])
        else:
            ax.invert_yaxis()

        ax.set_aspect('equal')
        ax.set_xlabel('eta [px]'); ax.set_ylabel('zeta [px]')
        ax.set_title('%s: error over the frame, arrows exaggerated %.0fx'
                     %(cam, scale))
        return ax


    def plot_coverage(self, ax=None):
        '''
        Where in each frame the calibration points actually are. Outside the
        region they cover, the camera model is extrapolating, and the error
        it shows on the points it was fitted to says nothing about how it
        behaves there.
        '''
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots(figsize=(7, 5))

        for c in self.cams:
            p = self.img[c]
            ax.plot(p[:, 0], p[:, 1], '.', ms=2, alpha=0.5, label=c)

        res = None
        for c in self.cams:
            res = self._resolution(c)
            if res is not None:
                break
        if res is not None:
            ax.add_patch(plt.Rectangle((0, 0), res[0], res[1], fill=False,
                                       ec='k', lw=1.2))
            ax.set_xlim(-0.03*res[0], 1.03*res[0])
            ax.set_ylim(1.03*res[1], -0.03*res[1])
        else:
            ax.invert_yaxis()

        ax.set_aspect('equal')
        ax.set_xlabel('eta [px]'); ax.set_ylabel('zeta [px]')
        ax.set_title('where the calibration points are in the frame')
        ax.legend(markerscale=4, fontsize=8)
        return ax


    # ------------------------------------------------------------------
    def _resolution(self, cam):
        '''
        The resolution of a camera, if it knows it.
        '''
        c = self.cameras[cam]
        for obj in (c, getattr(c, 'camera', None)):
            r = getattr(obj, 'resolution', None)
            if r is not None:
                return (float(r[0]), float(r[1]))
        return None


    def _frustum(self, cam, length):
        '''
        The corners of the field of view of a camera, at a given distance
        from it.

        The rays are taken from the camera itself, through its get_r, rather
        than from a rotation matrix, so this works for any 3D model: a Tsai
        camera and an extended Zolof camera both answer the question "which
        line does this pixel look along", and that is all that is needed.

        Note the word line. What get_r returns is the direction of the
        epipolar line, and its sign is not fixed by anything: the line
        through O along r is the same line as the one along -r, so
        img_system.stereo_match, which measures distances between lines,
        gives the same answer either way and has no reason to care. It does
        matter here, since a field of view drawn along the wrong sign points
        away from the scene. The Tsai model as it stands returns the sign
        pointing away from the scene. Rather than rely on that, the sign is
        settled from the data: the calibration points of this camera are in
        front of it by definition, so the ray is flipped if it disagrees.
        '''
        c = self.cameras[cam]
        res = self._resolution(cam)
        if res is None:
            return None

        O = asarray(c.O, dtype=float64)

        towards = self.lab[cam].mean(axis=0) - O
        n = norm(towards)
        if n == 0:
            return None
        towards = towards/n

        corners = [(0, 0), (res[0], 0), (res[0], res[1]), (0, res[1])]

        rays = []
        for (u, v) in corners:
            r = asarray(c.get_r(u, v), dtype=float64)
            rn = norm(r)
            if rn == 0 or not isfinite(r).all():
                return None
            rays.append(r/rn)

        # one sign for the whole frustum, taken from the middle of it
        if array(rays).mean(axis=0).dot(towards) < 0:
            rays = [-r for r in rays]

        return O, array([O + length*r for r in rays])


    def plot_extrinsics(self, ax=None, show_boards=True, frustum=None,
                        elev=22, azim=-60):
        '''
        The cameras and the target, in space: where each camera is, which way
        it looks, and where the target went. This is the counterpart of
        MATLAB's showExtrinsics.

        It is worth looking at for two reasons. It shows whether the geometry
        that came out is the geometry of the actual rig, which catches a
        camera that has been reconstructed behind the target or pointing the
        wrong way; and it shows the region the target actually visited, which
        is the region the calibration can be trusted in.

        input -
        show_boards (bool) - whether to draw the outline of the target at
                    each view
        frustum (float) - how far to draw the field of view of each camera.
                    By default, most of the way to the target.
        elev, azim (floats) - the viewing direction
        '''
        import matplotlib.pyplot as plt

        if ax is None:
            fig = plt.figure(figsize=(9, 7))
            ax = fig.add_subplot(projection='3d')

        allX = concatenate([self.lab[c] for c in self.cams])
        centre = allX.mean(axis=0)

        if frustum is None:
            d = [norm(asarray(self.cameras[c].O, dtype=float64) - centre)
                 for c in self.cams]
            frustum = 0.75*float(median(array(d)))

        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

        # the target at each view, drawn as the outline of its points
        if show_boards:
            drawn = set()
            for c in self.cams:
                v = self.view[c]
                for vi in unique(v):
                    if vi in drawn:
                        continue
                    drawn.add(vi)
                    P = self.lab[c][v == vi]
                    if len(P) < 3:
                        continue
                    ax.plot(P[:, 0], P[:, 1], P[:, 2], '.', ms=1.2,
                            color='0.55', alpha=0.5)

        # the cameras
        for k, c in enumerate(self.cams):
            col = colors[k % len(colors)]
            O = asarray(self.cameras[c].O, dtype=float64)
            ax.plot([O[0]], [O[1]], [O[2]], 'o', color=col, ms=9)
            ax.text(O[0], O[1], O[2], '  ' + c, color=col, fontsize=10)

            fr = self._frustum(c, frustum)
            if fr is None:
                continue
            O, Q = fr
            for i in range(4):
                ax.plot([O[0], Q[i, 0]], [O[1], Q[i, 1]], [O[2], Q[i, 2]],
                        '-', color=col, lw=0.8, alpha=0.8)
            for i in range(4):
                j = (i + 1) % 4
                ax.plot([Q[i, 0], Q[j, 0]], [Q[i, 1], Q[j, 1]],
                        [Q[i, 2], Q[j, 2]], '-', color=col, lw=1.2,
                        alpha=0.9)

        ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')
        ax.set_title('the cameras and the calibration target')
        ax.view_init(elev=elev, azim=azim)
        _equal_aspect_3d(ax)
        return ax




def _equal_aspect_3d(ax):
    '''
    Gives a 3D axis the same scale on all three axes, so that the geometry is
    not distorted. Without this a camera arrangement can look quite wrong.
    '''
    lims = array([ax.get_xlim3d(), ax.get_ylim3d(), ax.get_zlim3d()])
    span = (lims[:, 1] - lims[:, 0]).max()/2.0
    mid = lims.mean(axis=1)
    ax.set_xlim3d(mid[0] - span, mid[0] + span)
    ax.set_ylim3d(mid[1] - span, mid[1] + span)
    ax.set_zlim3d(mid[2] - span, mid[2] + span)
    try:
        ax.set_box_aspect((1, 1, 1))
    except Exception:
        pass
