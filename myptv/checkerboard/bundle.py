# -*- coding: utf-8 -*-
#!/usr/bin/env python3

"""
Created on Aug 27 2026

@author: Stefano Brizzolara


Adjusting several cameras and a moving board together.

The homographies of pose.py give, for each image separately, where the board
was relative to that camera. Those per image answers are noisy, and worse,
they are not consistent with one another: each camera reconstructs the board
independently, so the same physical corner ends up in slightly different
places depending on which camera is asked.

This module removes that inconsistency by solving for everything at once. The
unknowns are

    for each camera   its pose in the lab frame and its internal parameters,
    for each frame    the pose of the board,

and the quantity minimized is the distance, in pixels, between where each
corner was found and where the model puts it. The essential point is that the
board pose at a given frame is ONE set of six numbers shared by every camera
that saw it. That is what ties the cameras together: a corner is no longer
reconstructed separately by each camera but triangulated by all of them, and
cameras placed well apart constrain its depth far better than any of them
could alone.

One camera's pose is held fixed, which fixes the lab frame; without that the
whole arrangement could be rotated and translated freely without changing any
of the residuals. The scale is fixed by the physical size of the squares.


On the Jacobian.

The problem has a few hundred parameters and tens of thousands of residuals,
and it is sparse: a residual from camera c at frame f depends only on the
parameters of camera c and on the pose at frame f, and on nothing else. Given
the sparsity pattern, scipy can estimate the Jacobian by finite differences,
but it still needs of the order of fifty evaluations of the full residual to
build one, and several thousand evaluations to converge. On a real four
camera set that took several minutes and still had not converged.

So the Jacobian is written out analytically here. Every derivative is a
straightforward chain rule through

    X_lab = R_f B + t_f,   X_cam = R_c X_lab + t_c,   x = pi(X_cam),

and the only awkward piece is the derivative of a rotation matrix with
respect to its rotation vector, for which the compact expression of Gallego
and Yezzi is used,

    (1)    dR/dr_i = (r_i [r]x + [r x (I - R) e_i]x) R / |r|^2 ,

with the limit dR/dr_i = [e_i]x as |r| -> 0. With this the same problem
converges in seconds rather than not at all, which is the difference between
the method being usable and not.

The analytic Jacobian is checked against finite differences in the test
suite; get_jacobian_error is provided so that the check can be repeated on
any problem.


Reference:
Gallego, G., & Yezzi, A. (2015). A compact formula for the derivative of a
3-D rotation in exponential coordinates. Journal of Mathematical Imaging and
Vision, 51(3), 378-384.

"""

from numpy import (array, asarray, zeros, ones, eye, empty, arange, sqrt,
                   sin, cos, trace, arccos, argmax, float64, cross,
                   column_stack, concatenate, inf, abs as npabs)
from numpy.linalg import norm

from scipy.optimize import least_squares
from scipy.sparse import lil_matrix




def skew(v):
    '''
    The matrix [v]x for which [v]x w = v x w.
    '''
    return array([[0.0, -v[2], v[1]],
                  [v[2], 0.0, -v[0]],
                  [-v[1], v[0], 0.0]])




def rodrigues(r):
    '''
    The rotation matrix of a rotation vector, whose direction is the axis and
    whose length is the angle.

    input -
    r (array, 3) - the rotation vector

    output -
    R (array, (3,3))
    '''
    r = asarray(r, dtype=float64).ravel()
    th = norm(r)

    if th < 1e-12:
        return eye(3)

    K = skew(r/th)
    return eye(3) + sin(th)*K + (1.0 - cos(th))*K.dot(K)




def rvec_from_R(R):
    '''
    The rotation vector of a rotation matrix, the inverse of rodrigues.

    input -
    R (array, (3,3))

    output -
    r (array, 3)
    '''
    R = asarray(R, dtype=float64)

    c = max(-1.0, min(1.0, float((trace(R) - 1.0)/2.0)))
    th = float(arccos(c))

    if th < 1e-12:
        return zeros(3)

    if abs(3.141592653589793 - th) < 1e-6:
        # near a half turn the antisymmetric part vanishes, so the axis is
        # taken from the symmetric part instead
        A = (R + eye(3))/2.0
        i = int(argmax(A.diagonal()))
        v = A[:, i]/sqrt(max(A[i, i], 1e-30))
        return th*v/norm(v)

    w = array([R[2, 1] - R[1, 2], R[0, 2] - R[2, 0], R[1, 0] - R[0, 1]])
    return th*w/(2.0*sin(th))




def rodrigues_jacobian(r):
    '''
    The three matrices dR/dr_i of eq. (1).

    input -
    r (array, 3) - the rotation vector

    output -
    dR (list of 3 arrays, (3,3)) - the derivative of R with respect to each
                                   component of r
    '''
    r = asarray(r, dtype=float64).ravel()
    th2 = float(r.dot(r))

    if th2 < 1e-24:
        # in the limit of no rotation the derivative is just the generator
        return [skew(array([1.0, 0.0, 0.0])),
                skew(array([0.0, 1.0, 0.0])),
                skew(array([0.0, 0.0, 1.0]))]

    R = rodrigues(r)
    ImR = eye(3) - R
    rx = skew(r)

    out = []
    for i in range(3):
        e = zeros(3); e[i] = 1.0
        v = cross(r, ImR.dot(e))
        out.append((r[i]*rx + skew(v)).dot(R)/th2)

    return out




def project(X_cam, intr, n_dist):
    '''
    The pixel coordinates of points given in camera coordinates.

    The model is a pin-hole with square, unskewed pixels and an even radial
    distortion, optionally with the two tangential terms:

        x = X/Z,  y = Y/Z,  r2 = x^2 + y^2
        g = 1 + k1 r2 + k2 r2^2
        u = f (x g + 2 p1 x y + p2 (r2 + 2 x^2)) + cx
        v = f (y g + p1 (r2 + 2 y^2) + 2 p2 x y) + cy

    input -
    X_cam (array, (N,3)) - points in the camera frame
    intr (array) - [f, cx, cy] followed by n_dist distortion coefficients
    n_dist (int) - 0, 2 or 4

    output -
    uv (array, (N,2)) - pixel coordinates
    '''
    f, cx, cy = intr[0], intr[1], intr[2]
    d = intr[3:]

    x = X_cam[:, 0]/X_cam[:, 2]
    y = X_cam[:, 1]/X_cam[:, 2]
    r2 = x*x + y*y

    if n_dist == 0:
        xd, yd = x, y
    elif n_dist == 2:
        g = 1.0 + d[0]*r2 + d[1]*r2*r2
        xd, yd = x*g, y*g
    else:
        g = 1.0 + d[0]*r2 + d[1]*r2*r2
        xd = x*g + 2*d[2]*x*y + d[3]*(r2 + 2*x*x)
        yd = y*g + d[2]*(r2 + 2*y*y) + 2*d[3]*x*y

    return column_stack([f*xd + cx, f*yd + cy])




def _projection_derivatives(X_cam, intr, n_dist):
    '''
    The derivatives of the projection with respect to the point in camera
    coordinates and with respect to the intrinsic parameters.

    output -
    dXc (array, (N,2,3)) - d(u,v)/d(X,Y,Z) in the camera frame
    dI (array, (N,2,3+n_dist)) - d(u,v)/d(intrinsics)
    '''
    f = intr[0]
    d = intr[3:]
    N = len(X_cam)

    Z = X_cam[:, 2]
    x = X_cam[:, 0]/Z
    y = X_cam[:, 1]/Z
    r2 = x*x + y*y

    if n_dist == 0:
        g = ones(N)
        q = zeros(N)                      # dg/dr2
        xd, yd = x, y
        dxd_dx = ones(N); dxd_dy = zeros(N)
        dyd_dx = zeros(N); dyd_dy = ones(N)
    else:
        g = 1.0 + d[0]*r2 + d[1]*r2*r2
        q = d[0] + 2.0*d[1]*r2
        if n_dist == 2:
            xd, yd = x*g, y*g
            dxd_dx = g + 2.0*x*x*q
            dxd_dy = 2.0*x*y*q
            dyd_dx = 2.0*x*y*q
            dyd_dy = g + 2.0*y*y*q
        else:
            p1, p2 = d[2], d[3]
            xd = x*g + 2*p1*x*y + p2*(r2 + 2*x*x)
            yd = y*g + p1*(r2 + 2*y*y) + 2*p2*x*y
            dxd_dx = g + 2.0*x*x*q + 2*p1*y + p2*(2*x + 4*x)
            dxd_dy = 2.0*x*y*q + 2*p1*x + p2*(2*y)
            dyd_dx = 2.0*x*y*q + p1*(2*x) + 2*p2*y
            dyd_dy = g + 2.0*y*y*q + p1*(2*y + 4*y) + 2*p2*x

    # d(u,v)/d(x,y)
    du_dx = f*dxd_dx; du_dy = f*dxd_dy
    dv_dx = f*dyd_dx; dv_dy = f*dyd_dy

    # d(x,y)/d(X,Y,Z) = [[1/Z, 0, -x/Z], [0, 1/Z, -y/Z]]
    invZ = 1.0/Z
    dXc = empty((N, 2, 3))
    dXc[:, 0, 0] = du_dx*invZ
    dXc[:, 0, 1] = du_dy*invZ
    dXc[:, 0, 2] = -(du_dx*x + du_dy*y)*invZ
    dXc[:, 1, 0] = dv_dx*invZ
    dXc[:, 1, 1] = dv_dy*invZ
    dXc[:, 1, 2] = -(dv_dx*x + dv_dy*y)*invZ

    dI = zeros((N, 2, 3 + n_dist))
    dI[:, 0, 0] = xd            # d/df
    dI[:, 1, 0] = yd
    dI[:, 0, 1] = 1.0           # d/dcx
    dI[:, 1, 2] = 1.0           # d/dcy
    if n_dist >= 2:
        dI[:, 0, 3] = f*x*r2
        dI[:, 1, 3] = f*y*r2
        dI[:, 0, 4] = f*x*r2*r2
        dI[:, 1, 4] = f*y*r2*r2
    if n_dist >= 4:
        dI[:, 0, 5] = f*2*x*y
        dI[:, 1, 5] = f*(r2 + 2*y*y)
        dI[:, 0, 6] = f*(r2 + 2*x*x)
        dI[:, 1, 6] = f*2*x*y

    return dXc, dI




class multicamera_bundle(object):
    '''
    Holds the observations and the current estimate, and refines the estimate
    against the reprojection error.

    The observations are given as a dictionary keyed by (camera name, frame),
    each entry an (N,2) array of pixel positions of the board corners, in the
    same order as board_points. Not every camera need see every frame.
    '''

    def __init__(self, board_points, observations, camera_names,
                 intrinsics, extrinsics, poses, n_dist=2,
                 fix_camera=None, free_principal_point=True):
        '''
        input -
        board_points (array, (N,3)) - the corners in the frame of the board,
                    with the third coordinate zero
        observations (dict) - (camera, frame) -> (N,2) pixel positions
        camera_names (list) - the cameras, in a fixed order
        intrinsics (dict) - camera -> [f, cx, cy] + n_dist coefficients
        extrinsics (dict) - camera -> (R, t) taking lab to camera coordinates
        poses (dict) - frame -> (R, t) taking board to lab coordinates
        n_dist (int) - 0, 2 or 4 distortion coefficients
        fix_camera (string) - the camera whose pose defines the lab frame and
                    is therefore held fixed. Defaults to the first.
        free_principal_point (bool) - whether the principal point is refined
        '''
        self.B = asarray(board_points, dtype=float64)
        if self.B.shape[1] != 3:
            raise ValueError('board_points must be an (N,3) array')

        self.cams = list(camera_names)
        self.n_dist = int(n_dist)
        if self.n_dist not in (0, 2, 4):
            raise ValueError('n_dist must be 0, 2 or 4')

        self.fix_camera = self.cams[0] if fix_camera is None else fix_camera
        if self.fix_camera not in self.cams:
            raise ValueError('fix_camera %s is not one of the cameras'
                             %self.fix_camera)

        self.free_pp = bool(free_principal_point)

        self.frames = sorted(set(f for (_, f) in observations))
        self.pairs = [(ci, f) for f in self.frames
                      for ci, c in enumerate(self.cams)
                      if (c, f) in observations]
        if len(self.pairs) == 0:
            raise ValueError('there are no observations')

        self.obs = array([observations[(self.cams[ci], f)]
                          for (ci, f) in self.pairs], dtype=float64)

        self.NP = len(self.B)
        self.NI = 3 + self.n_dist
        self.CAMP = 6 + self.NI
        self.NC = len(self.cams)
        self.fidx = {f: j for j, f in enumerate(self.frames)}

        self.x = self._pack(intrinsics, extrinsics, poses)


    def __repr__(self):
        return ('multicamera_bundle: %d cameras, %d frames, %d observations, '
                '%d parameters'%(self.NC, len(self.frames), len(self.pairs),
                                 len(self.x)))


    # ------------------------------------------------------------------
    def _pack(self, intr, ext, poses):
        p = []
        for c in self.cams:
            R, t = ext[c]
            p += list(rvec_from_R(R)) + list(asarray(t, dtype=float64).ravel())
            p += list(asarray(intr[c], dtype=float64).ravel())
        for f in self.frames:
            R, t = poses[f]
            p += list(rvec_from_R(R)) + list(asarray(t, dtype=float64).ravel())
        return array(p, dtype=float64)


    def unpack(self, p=None):
        '''
        Returns the intrinsics, the camera poses and the board poses that a
        parameter vector stands for.
        '''
        p = self.x if p is None else p
        intr, ext, poses = {}, {}, {}

        for i, c in enumerate(self.cams):
            s = i*self.CAMP
            ext[c] = (rodrigues(p[s:s+3]), p[s+3:s+6].copy())
            intr[c] = p[s+6:s+self.CAMP].copy()

        o = self.NC*self.CAMP
        for j, f in enumerate(self.frames):
            s = o + j*6
            poses[f] = (rodrigues(p[s:s+3]), p[s+3:s+6].copy())

        return intr, ext, poses


    # ------------------------------------------------------------------
    def residuals(self, p):
        '''
        The signed differences, in pixels, between the projected and the
        observed corners, flattened.
        '''
        out = empty((len(self.pairs), self.NP, 2))

        o = self.NC*self.CAMP
        Rf, tf = {}, {}
        for j, f in enumerate(self.frames):
            s = o + j*6
            Rf[f] = rodrigues(p[s:s+3])
            tf[f] = p[s+3:s+6]

        Xlab = {f: self.B.dot(Rf[f].T) + tf[f] for f in self.frames}

        for n, (ci, f) in enumerate(self.pairs):
            s = ci*self.CAMP
            Rc = rodrigues(p[s:s+3])
            tc = p[s+3:s+6]
            intr = p[s+6:s+self.CAMP]
            Xc = Xlab[f].dot(Rc.T) + tc
            out[n] = project(Xc, intr, self.n_dist) - self.obs[n]

        return out.ravel()


    def jacobian(self, p):
        '''
        The analytic Jacobian of residuals, as a sparse matrix.
        '''
        from scipy.sparse import coo_matrix

        NP, NI, CAMP, NC = self.NP, self.NI, self.CAMP, self.NC
        nres = len(self.pairs)*NP*2
        o = NC*CAMP

        Rf, tf, dRf = {}, {}, {}
        for j, f in enumerate(self.frames):
            s = o + j*6
            rv = p[s:s+3]
            Rf[f] = rodrigues(rv)
            tf[f] = p[s+3:s+6]
            dRf[f] = rodrigues_jacobian(rv)

        Xlab = {f: self.B.dot(Rf[f].T) + tf[f] for f in self.frames}

        # A residual block depends on CAMP camera parameters and 6 pose
        # parameters, so the whole Jacobian has exactly this many non-zeros.
        # They are written into pre-allocated arrays rather than accumulated
        # in lists, which for a problem of this size is the difference
        # between the Jacobian costing more than the solve and costing
        # almost nothing.
        ncol = CAMP + 6
        nnz = len(self.pairs)*NP*2*ncol
        rows = empty(nnz, dtype=int)
        cols = empty(nnz, dtype=int)
        vals = empty(nnz, dtype=float64)

        block = 2*NP*ncol
        ru = 2*arange(NP)

        for n, (ci, f) in enumerate(self.pairs):
            s = ci*CAMP
            rv = p[s:s+3]
            Rc = rodrigues(rv)
            dRc = rodrigues_jacobian(rv)
            tc = p[s+3:s+6]
            intr = p[s+6:s+CAMP]

            Xc = Xlab[f].dot(Rc.T) + tc
            dXc, dI = _projection_derivatives(Xc, intr, self.n_dist)

            # D[:, :, k] is d(u,v)/d(parameter k), for the ncol parameters
            # this block depends on
            D = empty((NP, 2, ncol))

            for i in range(3):                       # camera rotation
                dX = Xlab[f].dot(dRc[i].T)
                D[:, 0, i] = (dXc[:, 0, :]*dX).sum(axis=1)
                D[:, 1, i] = (dXc[:, 1, :]*dX).sum(axis=1)

            D[:, :, 3:6] = dXc                       # camera translation
            D[:, :, 6:6+NI] = dI                     # intrinsics

            sf = o + self.fidx[f]*6
            for i in range(3):                       # board rotation
                dX = self.B.dot(dRf[f][i].T).dot(Rc.T)
                D[:, 0, CAMP+i] = (dXc[:, 0, :]*dX).sum(axis=1)
                D[:, 1, CAMP+i] = (dXc[:, 1, :]*dX).sum(axis=1)
            for i in range(3):                       # board translation
                D[:, 0, CAMP+3+i] = dXc[:, 0, :].dot(Rc[:, i])
                D[:, 1, CAMP+3+i] = dXc[:, 1, :].dot(Rc[:, i])

            colidx = concatenate([arange(s, s+CAMP), arange(sf, sf+6)])

            a = n*block
            rows[a:a+block] = ((n*NP*2 + ru)[:, None, None]
                               + array([0, 1])[None, :, None]).repeat(
                                   ncol, axis=2).ravel()
            cols[a:a+block] = colidx[None, None, :].repeat(
                NP, axis=0).repeat(2, axis=1).ravel()
            vals[a:a+block] = D.ravel()

        return coo_matrix((vals, (rows, cols)),
                          shape=(nres, len(p))).tocsr()


    def free_mask(self):
        '''
        A boolean array saying which parameters are free to move. The pose of
        the fixed camera is the gauge and never moves; the principal points
        are optionally held too.
        '''
        m = ones(len(self.x), dtype=bool)

        ci = self.cams.index(self.fix_camera)
        m[ci*self.CAMP: ci*self.CAMP + 6] = False

        if not self.free_pp:
            for i in range(self.NC):
                s = i*self.CAMP
                m[s+7] = False
                m[s+8] = False

        return m


    def free_idx(self):
        '''
        The indices of the parameters that are free to move.
        '''
        m = self.free_mask()
        return array([j for j in range(len(m)) if m[j]], dtype=int)


    def sparsity(self):
        '''
        The sparsity pattern of the Jacobian, for the finite difference path.
        '''
        S = lil_matrix((len(self.pairs)*self.NP*2, len(self.x)), dtype=int)
        o = self.NC*self.CAMP
        for n, (ci, f) in enumerate(self.pairs):
            r = slice(n*self.NP*2, (n+1)*self.NP*2)
            S[r, ci*self.CAMP:(ci+1)*self.CAMP] = 1
            s = o + self.fidx[f]*6
            S[r, s:s+6] = 1
        return S.tocsr()


    def get_jacobian_error(self, p=None, n_check=40, seed=0):
        '''
        Compares the analytic Jacobian with finite differences on a sample of
        its columns, and returns the largest relative discrepancy. It should
        be of the order of the finite difference truncation error, say 1e-6.
        A large value means the analytic derivative is wrong.

        input -
        p (array) - where to check; the current estimate by default
        n_check (int) - how many columns to check
        seed (int) - which columns, chosen deterministically from this

        output -
        err (float) - the largest relative difference over the columns checked
        '''
        from numpy.random import default_rng

        p = self.x.copy() if p is None else p.copy()
        J = self.jacobian(p)
        mask = self.free_mask()
        free = [j for j in range(len(p)) if mask[j]]

        rng = default_rng(seed)
        cols = rng.choice(free, size=min(n_check, len(free)), replace=False)

        r0 = self.residuals(p)
        worst = 0.0
        for j in cols:
            h = 1e-6*max(1.0, abs(p[j]))
            q = p.copy(); q[j] += h
            num = (self.residuals(q) - r0)/h
            ana = asarray(J[:, j].todense()).ravel()
            scale = max(npabs(num).max(), npabs(ana).max(), 1e-12)
            worst = max(worst, float(npabs(num - ana).max()/scale))

        return worst


    # ------------------------------------------------------------------
    def solve(self, max_nfev=200, verbose=0, loss='linear', f_scale=1.0,
              bounds=None, use_analytic=True):
        '''
        Refines the estimate.

        input -
        max_nfev (int) - the limit on residual evaluations
        verbose (int) - passed to scipy
        loss (string) - 'linear', or 'soft_l1' to be robust against a frame
                        whose corners were labelled wrongly
        f_scale (float) - the scale of the robust loss, in pixels
        bounds (tuple) - (lower, upper) arrays, or None
        use_analytic (bool) - whether to use the analytic Jacobian. Setting it
                        False falls back to finite differences with the
                        sparsity pattern, which is far slower and is here so
                        that the two can be compared.

        output -
        result - the scipy result object; the estimate is also stored on self
        '''
        # The parameters that are held fixed are removed from the problem
        # rather than merely given a zero derivative. Leaving them in makes
        # the trust region spend part of its radius on directions that cannot
        # change anything, which slows the whole solve down noticeably.
        free = self.free_idx()
        p_full = self.x.copy()

        def expand(q):
            p = p_full.copy()
            p[free] = q
            return p

        def res_free(q):
            return self.residuals(expand(q))

        if use_analytic:
            kw = {'jac': lambda q: self.jacobian(expand(q))[:, free]}
        else:
            kw = {'jac_sparsity': self.sparsity()[:, free]}

        if bounds is None:
            lo, hi = -inf, inf
        else:
            lo = asarray(bounds[0], dtype=float64)
            hi = asarray(bounds[1], dtype=float64)
            if lo.ndim:
                lo = lo[free]
            if hi.ndim:
                hi = hi[free]

        res = least_squares(res_free, self.x[free], bounds=(lo, hi),
                            x_scale='jac', loss=loss, f_scale=f_scale,
                            ftol=1e-12, xtol=1e-14, gtol=1e-12,
                            max_nfev=max_nfev, verbose=verbose, **kw)

        self.x = expand(res.x)
        return res


    def rms(self, p=None):
        '''
        The root mean square reprojection error, in pixels.
        '''
        r = self.residuals(self.x if p is None else p)
        return float(sqrt((r**2).mean()))


    def per_view_rms(self, p=None):
        '''
        The root mean square reprojection error of each (camera, frame), which
        is what to look at when hunting for a badly labelled image.

        output - a list of (camera, frame, rms) tuples
        '''
        r = self.residuals(self.x if p is None else p)
        r = r.reshape(len(self.pairs), self.NP, 2)
        out = []
        for n, (ci, f) in enumerate(self.pairs):
            out.append((self.cams[ci], f,
                        float(sqrt((r[n]**2).sum(axis=1).mean()))))
        return out


    def lab_coordinates(self):
        '''
        The lab space coordinates of every board corner at every frame, which
        is what a calibration points file needs.

        output - a dictionary, frame -> (N,3) array
        '''
        _, _, poses = self.unpack()
        return {f: self.B.dot(poses[f][0].T) + poses[f][1]
                for f in self.frames}
