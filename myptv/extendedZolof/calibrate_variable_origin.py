# -*- coding: utf-8 -*-
"""
Variable-origin extension for the Extended Zolof camera model.

The original Extended Zolof backward transformation assumes every epipolar
line for a camera passes through a single fixed point, O. This module
relaxes that assumption: it fits the origin itself as a polynomial function
of the pixel, O(x) = [C] G(x), alongside the direction, r(x) = [B] G(x),
using the same 10-term pixel basis G(x) = cam.get_xCol(x).

The approach mirrors the method used in Solof_v2_python_translation.py
(itself following Soloff et al.'s multi-camera calibration), adapted to
reuse this codebase's already-fit forward model (cam.A) and polynomial
basis functions (cam.get_XCol / cam.get_xCol), so that results stay
consistent with the rest of the Extended Zolof pipeline:

  1. For each calibration point (x_i, X_i), use the already-fit forward
     model (cam.projection, i.e. the [A] polynomial) to locate TWO points
     on the true 3D viewing ray of x_i:
       - X1: found by a spherical (theta, phi) search around X_i for the
             direction that reprojects closest to x_i.
       - X2: found by a free 3D search (initialized at X_i) for the point
             that reprojects closest to x_i.
     The local ray direction is then e_i = normalize(X2 - X1), and points
     whose reprojection residual exceeds `thresh_pix` are dropped.

  2. Build the pixel basis matrix G_i = cam.get_xCol(x_i) for the surviving
     points, and get closed-form initial guesses for C and B via
     least-squares (O(x) ~ X1, direction(x) ~ e_i).

  3. Jointly refine [C | B] with a nonlinear optimizer, minimizing the sum
     of squared distances from the TRUE calibration points X_i to the
     modeled epipolar lines l(x_i) = O(x_i) + a * direction(x_i). This is
     the same cost used in Solof_v2's dist_pt_to_line_3d_ext / indirect fit.

Use `fit_variable_origin_model()` from calibrate_extendedZolof.calibrate()
(see the `variable_origin=True` option there).
"""

import warnings

from numpy import (array, cos, sin, dot, cross, hstack, mean, median,
                    sum as npsum, abs as npabs)
from numpy.linalg import norm, lstsq
from myptv.extendedZolof.camera import (build_xCol, n_terms_for_order,
                                        order_from_n_terms)
from scipy.optimize import minimize


def _get_xCol_scaled(x, center, scale, n=10):
    '''
    Standalone version of camera.get_xCol_scaled(), used during fitting
    before px_center/px_scale are assigned to the camera object. Shares
    build_xCol() with the camera class so the two cannot drift apart.

    n - number of basis terms (3, 6, 10, 15, ...).
    '''
    x1 = (x[0] - center[0]) / scale[0]
    x2 = (x[1] - center[1]) / scale[1]
    return build_xCol(x1, x2, n)


def _get_ray_points(cam, Xi, xi):
    '''
    Given a calibration point (xi, Xi) and the camera's current forward
    model (cam.projection, i.e. cam.A), find two points, X1 and X2, that
    lie on the 3D viewing ray of xi, plus the reprojection residual (in
    pixels) at X1. Returns (X1, X2, resid).
    '''
    Xi = array(Xi, dtype=float)
    xi = array(xi, dtype=float)

    def reproj_err(X):
        return npsum((array(cam.projection(X)) - xi) ** 2)

    # --- point 1: spherical search around Xi (guarantees a unit offset) ---
    def sphere_cost(a):
        direction = array([
            cos(a[0]) * sin(a[1]),
            sin(a[0]) * sin(a[1]),
            cos(a[1])
        ])
        return reproj_err(Xi + direction)

    res_dir = minimize(sphere_cost, array([0.0, 0.0]), method='BFGS')
    a = res_dir.x
    direction = array([
        cos(a[0]) * sin(a[1]),
        sin(a[0]) * sin(a[1]),
        cos(a[1])
    ])
    X1 = Xi + direction

    # --- point 2: free 3D search, initialized at Xi ---
    res_pt = minimize(reproj_err, Xi, method='BFGS')
    X2 = res_pt.x

    resid = reproj_err(X1) ** 0.5

    return X1, X2, resid


def _dist_pt_to_line_plane(C, B, G_c, G_b, Xtrue, O0, u, v):
    '''
    Plane-constrained variant of _dist_pt_to_line_batch.

    C    : ndarray, shape (n_c, 2) - in-plane origin coefficients; column
           0 gives the plane_u offset a(x), column 1 the plane_v offset
           b(x), so O(x) = O0 + u*a(x) + v*b(x).
    B    : ndarray, shape (n_b, 3) - direction coefficients.
    G_c, G_b, Xtrue : as in _dist_pt_to_line_batch.
    O0   : ndarray (3,) - the point the plane passes through.
    u, v : ndarray (3,) each - orthonormal in-plane axes.

    Returns the point-to-line distance for every row (N,).
    '''
    ab = G_c.dot(C)                       # (N,2)
    O = O0[None, :] + ab[:, 0:1]*u[None, :] + ab[:, 1:2]*v[None, :]

    U = G_b.dot(B)                        # (N,3)
    Un = U / norm(U, axis=1, keepdims=True)

    diff = Xtrue - O
    proj = npsum(diff * Un, axis=1, keepdims=True) * Un
    perp = diff - proj

    return norm(perp, axis=1)


def _make_plane_basis(n):
    '''
    Given a plane normal n, returns two orthonormal in-plane axes (u, v)
    such that (u, v, n_hat) is a right-handed orthonormal triad.
    '''
    n_hat = array(n, dtype=float)
    n_hat = n_hat / norm(n_hat)

    # pick any axis not (nearly) parallel to n to seed the cross product
    seed = array([1.0, 0.0, 0.0])
    if abs(dot(seed, n_hat)) > 0.9:
        seed = array([0.0, 1.0, 0.0])

    u = cross(seed, n_hat)
    u = u / norm(u)
    v = cross(n_hat, u)
    v = v / norm(v)

    return u, v, n_hat


def _dist_pt_to_line_batch(C, B, G_c, G_b, Xtrue):
    '''
    C    : ndarray, shape (n_c, 3) - origin coefficients (rows = pixel
           basis terms, columns = X,Y,Z). n_c may be smaller than n_b
           when C uses a reduced polynomial order (see c_order in
           fit_variable_origin_model).
    B    : ndarray, shape (n_b, 3) - direction coefficients.
    G_c  : ndarray, shape (N, n_c) - origin pixel basis per point.
    G_b  : ndarray, shape (N, n_b) - direction pixel basis per point.
    Xtrue: ndarray, shape (N, 3) - true calibration points.

    Returns the point-to-line distance for every row (N,).
    '''
    O = G_c.dot(C)           # (N,3)
    U = G_b.dot(B)           # (N,3)
    Un = U / norm(U, axis=1, keepdims=True)

    diff = Xtrue - O
    proj = npsum(diff * Un, axis=1, keepdims=True) * Un
    perp = diff - proj

    return norm(perp, axis=1)


# number of polynomial basis terms for each supported order, using the
# same increasing-degree term ordering as get_xCol/get_xCol_scaled:
# order 1 (linear/affine) : 1, x, y                          -> 3 terms
# order 2 (quadratic)     : + x^2, y^2, xy                    -> 6 terms
# order 3 (full cubic)    : + x^2y, xy^2, x^3, y^3             -> 10 terms
_ORDER_N_TERMS = {o: n_terms_for_order(o) for o in range(1, 7)}


def fit_variable_origin_model(cam, X_list, x_list, O_init=None,
                               thresh_pix=5.0, min_points=15, c_order=3,
                               origin_model='free', freeze_C=False,
                               b_order=3):
    '''
    Fits the variable-origin epipolar model, replacing the fixed-origin
    assumption r(x) = [B]G(x), O = const.

    Two parameterizations of the per-pixel origin are supported, via
    origin_model:

      'free'  : O(x) = [C] G_c(x)          ([C] shape (n_c, 3))
                The origin floats freely in 3D. NOTE this carries a gauge
                degeneracy: sliding an origin ALONG its own ray leaves the
                epipolar line unchanged, so one of the three components is
                not identifiable from the point-to-line cost and is driven
                purely by noise - a real source of overfitting.

      'plane' : O(x) = O0 + u*a(x) + v*b(x)  ([C] shape (n_c, 2))
                The origins are constrained to a single fixed plane that
                passes through the approximate common camera center O0,
                with normal along (centroid(X) - O0) - i.e. a plane facing
                the calibration volume. (u, v) are orthonormal in-plane
                axes. This removes the degenerate direction above and
                needs 2 rather than 3 coefficient sets per basis term, so
                it is both better-conditioned and lower-variance.

    inputs -
    cam        : a camera_extendedZolof instance with an already-fit
                 forward model (cam.A must be set and up to date).
    X_list     : list/array of 3D calibration points, shape (N,3).
    x_list     : list/array of corresponding 2D image points, shape (N,2).
    O_init     : optional initial guess for a representative camera center
                 (e.g. the fixed-origin estimate from step2). Used to pick
                 a consistent sign convention for ray directions, and for
                 origin_model='plane' it is ALSO the point the plane passes
                 through (O0) and defines the plane normal. Defaults to the
                 mean of X_list, which is a poor choice for 'plane' mode -
                 pass the step2 camera-center estimate there.
    thresh_pix : calibration points whose reprojection residual (using the
                 current forward model) exceeds this are excluded from the
                 fit.
    min_points : minimum number of surviving points required to attempt
                 the fit; raises ValueError if not met.
    c_order    : polynomial order used for [C] (the ORIGIN, O(x)), as
                 opposed to [B] (the DIRECTION, r(x)), which always uses
                 the full cubic (10-term) basis. One of:
                     1 - linear/affine:  O(x) = C0 + C1*x + C2*y  (3 terms)
                     2 - quadratic (6 terms)
                     3 - full cubic (10 terms), matching [B] (default,
                         same behavior as before this option existed).
                 Since [C] predicts an unbounded lab-space origin (as
                 opposed to [B]'s bounded unit direction), it is the more
                 overfitting-prone of the two; giving it fewer free
                 parameters than [B] is a reasonable way to reduce that
                 risk when calibration points are limited or unevenly
                 distributed. Use cross_validate_origin_models() (below)
                 to check whether a given c_order actually generalizes
                 better on held-out points, rather than assuming it does.

    freeze_C   : if True, [C] is NOT fit - it is taken unchanged from
                 cam.C, and only [B] (the ray directions) is optimized.
                 c_order and origin_model are then inherited from the
                 camera as well, and so are the pixel normalization
                 (cam.px_center/px_scale) and, in 'plane' mode, the plane
                 itself (cam.O, cam.plane_u/v/n) - all of these define the
                 reference frame that [C]'s coefficients are expressed in,
                 so re-deriving any of them from a new point set would
                 silently reinterpret the frozen coefficients.

                 This is intended for refining a camera against particle
                 trajectories: those 3D target positions were themselves
                 triangulated with the current model, and the points are
                 selected for LOW triangulation error, i.e. for being most
                 consistent with the existing [C]. Refitting [C] on them
                 therefore tends to reinforce its own bias rather than
                 correct it, while [B] is far less able to exploit that
                 feedback. Requires cam.variable_origin == True.

    returns (C, B, mask, resid, px_center, px_scale, plane) -
    C     : ndarray - the fitted origin coefficients (rows = pixel basis
            terms; n_c = 3/6/10 rows depending on c_order). For
            origin_model='free' it has 3 columns (X,Y,Z); for 'plane' it
            has 2 columns (the plane_u and plane_v offsets). Layout
            matches cam.C in both cases.
    B     : ndarray, shape (n_b,3) with n_b set by b_order - the fitted
            direction coefficients,
            same layout as cam.B. NOTE: both are fit against the
            NORMALIZED pixel basis (see px_center/px_scale below) - they
            are only meaningful together with that normalization, i.e.
            via cam.get_xCol_scaled() / cam.get_ray(), not cam.get_xCol().
    mask  : boolean array, shape (N,) - which of the *input* points (in
            X_list/x_list order) survived the reprojection-residual
            rejection step and were used in the fit.
    resid : ndarray - final point-to-line distances (lab units) for the
            points used in the fit, useful as a calibration diagnostic.
    px_center, px_scale : ndarray, shape (2,) each - the pixel-coordinate
            centering/scaling used to build the basis for C and B. The
            caller should store these as cam.px_center / cam.px_scale.
    plane : None for origin_model='free'; for 'plane' a dict with keys
            'O0', 'u', 'v', 'n' (each a (3,) ndarray) defining the plane.
            The caller should store these as cam.O / cam.plane_u /
            cam.plane_v / cam.plane_n.
    origin_model : the origin model actually used. This can DIFFER from
            the requested one when freeze_C=True (the camera's own model
            wins), so callers must store THIS value on the camera rather
            than the argument they passed in - otherwise the camera would
            be labelled with a model that does not match its [C].
    '''
    if c_order not in _ORDER_N_TERMS:
        raise ValueError(f'c_order must be one of {sorted(_ORDER_N_TERMS)} '
                          f'(got {c_order!r}).')
    if origin_model not in ('free', 'plane'):
        raise ValueError("origin_model must be 'free' or 'plane' "
                          f'(got {origin_model!r}).')
    if b_order not in _ORDER_N_TERMS:
        raise ValueError(f'b_order must be one of {sorted(_ORDER_N_TERMS)} '
                          f'(got {b_order!r}).')
    n_b = _ORDER_N_TERMS[b_order]

    if freeze_C:
        # Reuse the camera's existing origin model wholesale. Everything
        # that [C]'s coefficients are expressed RELATIVE TO must be reused
        # unchanged too, or the frozen coefficients would be reinterpreted
        # against a different reference and silently produce wrong
        # origins: the pixel normalization (px_center/px_scale), and for
        # 'plane' mode also the plane's anchor point (cam.O) and axes.
        if not getattr(cam, 'variable_origin', False):
            raise ValueError(
                'freeze_C=True requires a camera that was already '
                'calibrated with the variable-origin model, but '
                f'"{getattr(cam, "name", "?")}" has variable_origin=False '
                '(so it has no [C] to freeze).')

        C_frozen = array(cam.C, dtype=float)
        n_c = C_frozen.shape[0]

        if n_c not in _ORDER_N_TERMS.values():
            raise ValueError(
                f'Camera [C] has {n_c} rows, which is not a valid '
                'polynomial basis size (3, 6 or 10).')

        cam_order = {v: k for k, v in _ORDER_N_TERMS.items()}[n_c]
        if cam_order != c_order:
            warnings.warn(
                f'freeze_C=True: using the camera\'s own C order '
                f'({cam_order}), ignoring the requested c_order={c_order}.')

        cam_model = getattr(cam, 'origin_mode', 'free')
        if cam_model != origin_model:
            warnings.warn(
                f'freeze_C=True: using the camera\'s own origin model '
                f"('{cam_model}'), ignoring the requested "
                f"origin_model='{origin_model}'.")
            origin_model = cam_model

        expected_cols = 2 if origin_model == 'plane' else 3
        if C_frozen.shape[1] != expected_cols:
            raise ValueError(
                f"Camera [C] has {C_frozen.shape[1]} columns but origin "
                f"model '{origin_model}' requires {expected_cols}. The "
                'camera file may be inconsistent.')
    else:
        n_c = _ORDER_N_TERMS[c_order]
        C_frozen = None

    X_list = array(X_list, dtype=float)
    x_list = array(x_list, dtype=float)
    N = len(X_list)

    if O_init is None:
        O_init = mean(X_list, axis=0)
    else:
        O_init = array(O_init, dtype=float)

    X1_list, e_list, X_kept, x_kept = [], [], [], []
    keep_idx = []

    for i in range(N):
        X1, X2, resid = _get_ray_points(cam, X_list[i], x_list[i])

        u = X2 - X1
        nu = norm(u)
        if nu == 0:
            continue

        if resid > thresh_pix:
            continue

        e = u / nu

        X1_list.append(X1)
        e_list.append(e)
        X_kept.append(X_list[i])
        x_kept.append(x_list[i])
        keep_idx.append(i)

    if len(keep_idx) < min_points:
        raise ValueError(
            'Too few valid rays survived reprojection-residual rejection '
            f'({len(keep_idx)} < {min_points}) to fit the variable-origin '
            'model; consider raising thresh_pix or checking the forward '
            'model (cam.A) is fit and up to date.')

    X1_arr = array(X1_list)
    e_arr = array(e_list)
    X_kept = array(X_kept)
    x_kept = array(x_kept)

    # consistent sign convention: rays should point away from O_init,
    # towards the calibration point (matters for the direction fit; the
    # point-to-line distance itself is sign-agnostic, but a consistent
    # sign keeps [B] well-behaved / interpretable).
    signs = array([
        1.0 if dot(e_arr[i], X_kept[i] - O_init) >= 0 else -1.0
        for i in range(len(keep_idx))
    ])
    e_arr = e_arr * signs[:, None]

    if freeze_C:
        # MUST reuse the camera's normalization - [C]'s coefficients are
        # expressed in these normalized pixel units.
        px_center = array(cam.px_center, dtype=float)
        px_scale = array(cam.px_scale, dtype=float)
    else:
        # --- pixel-coordinate normalization ---
        # [C] predicts unbounded lab-space coordinates (unlike [B], which
        # predicts a bounded unit direction), so fitting it against a cubic
        # polynomial in RAW pixel coordinates is poorly conditioned and prone
        # to wild extrapolation for any pixel far from the calibration points
        # (e.g. x^3 for x~1000px is ~1e9). Center/scale the pixels first so
        # every basis term stays O(1) within the calibrated footprint. These
        # constants are returned so the caller can store them on the camera
        # (cam.px_center / cam.px_scale) and use the exact same normalized
        # basis at prediction time (see camera.get_xCol_scaled / get_ray).
        px_center = mean(x_kept, axis=0)
        px_range = x_kept.max(axis=0) - x_kept.min(axis=0)
        px_scale = array([s if s > 0 else 1.0 for s in px_range / 2.0])

    # one basis long enough for both, then sliced per matrix - the terms
    # are in ascending-degree order so a prefix IS the lower-order basis.
    n_max = max(n_b, n_c)
    G_all = array([_get_xCol_scaled(xi, px_center, px_scale, n=n_max)
                   for xi in x_kept])
    G_b = G_all[:, :n_b]
    G_c = G_all[:, :n_c]   # (M, n_c) -- leading terms = the reduced-order basis,
                         # since get_xCol_scaled's terms are ordered by
                         # increasing degree (see _ORDER_N_TERMS above).

    if origin_model == 'free':
        if not freeze_C:
            # closed-form initial guesses
            C0, *_ = lstsq(G_c, X1_arr, rcond=None)   # (n_c,3): O(x) ~ X1
            B0, *_ = lstsq(G_b, e_arr, rcond=None)    # (10,3): direction(x) ~ e

            a0 = hstack([C0.ravel(), B0.ravel()])     # flat: [C terms | B terms]

        if freeze_C:
            # only [B] is free; [C] is taken from the camera as-is.
            C = C_frozen
            B0, *_ = lstsq(G_b, e_arr, rcond=None)

            def b_only_cost(b):
                B_ = b.reshape(n_b, 3)
                d = _dist_pt_to_line_batch(C, B_, G_c, G_b, X_kept)
                return npsum(d ** 2)

            res = minimize(b_only_cost, B0.ravel(), method='BFGS')
            B = res.x.reshape(n_b, 3)
        else:
            def indirect_cost(a):
                C_ = a[:n_c * 3].reshape(n_c, 3)
                B_ = a[n_c * 3:].reshape(n_b, 3)
                d = _dist_pt_to_line_batch(C_, B_, G_c, G_b, X_kept)
                return npsum(d ** 2)

            res = minimize(indirect_cost, a0, method='BFGS')

            C = res.x[:n_c * 3].reshape(n_c, 3)
            B = res.x[n_c * 3:].reshape(n_b, 3)

        final_resid = _dist_pt_to_line_batch(C, B, G_c, G_b, X_kept)
        plane = None

    else:
        if freeze_C:
            # Reuse the camera's own plane. [C]'s coefficients are offsets
            # measured from THIS anchor point along THESE axes, so
            # re-deriving the plane from the new point set would silently
            # reinterpret them.
            O0 = array(cam.O, dtype=float)
            u_ax = array(cam.plane_u, dtype=float)
            v_ax = array(cam.plane_v, dtype=float)
            n_hat = array(cam.plane_n, dtype=float)
        else:
            # ---- plane-constrained origins ----
            # (step 2) The plane passes through the approximate common
            # camera center O0 (=O_init, normally the step2 estimate), with
            # normal along the line joining O0 to the centroid of the
            # calibration points - i.e. a plane squarely facing the
            # measurement volume. Because that normal is roughly the mean
            # viewing direction, the in-plane axes (u,v) span precisely the
            # origin displacements that actually CHANGE the epipolar lines,
            # while the discarded out-of-plane direction is the
            # (unidentifiable) along-ray one.
            O0 = array(O_init, dtype=float)
            n_vec = mean(X_kept, axis=0) - O0
            if norm(n_vec) == 0:
                raise ValueError(
                    'Cannot build the origin plane: the estimated camera '
                    'center coincides with the centroid of the calibration '
                    'points. Check O_init / the step2 camera-center estimate.')
            u_ax, v_ax, n_hat = _make_plane_basis(n_vec)

        # (step 3) intersect each estimated ray with that plane, and
        # express the crossing point in the plane's own 2D coordinates.
        # These offsets (a_i, b_i) are the per-pixel deviations of the
        # true origin from the single common center O0.
        denom = e_arr.dot(n_hat)                        # (M,)
        # a ray nearly parallel to the plane has no usable (or a wildly
        # extrapolated) intersection - drop those points.
        good = npabs(denom) > 1e-6
        if good.sum() < min_points:
            raise ValueError(
                'Too few rays intersect the origin plane at a usable '
                f'angle ({int(good.sum())} < {min_points}). The plane '
                'normal may be badly estimated.')

        t = ((O0[None, :] - X1_arr).dot(n_hat)) / denom  # (M,)
        P = X1_arr + t[:, None] * e_arr                  # crossing points (M,3)

        d_plane = P - O0[None, :]
        a_off = d_plane.dot(u_ax)
        b_off = d_plane.dot(v_ax)
        ab_target = array([a_off, b_off]).T              # (M,2)

        # keep only the well-conditioned rows for fitting
        G_c_g, G_b_g = G_c[good], G_b[good]
        X_g, ab_g, e_g = X_kept[good], ab_target[good], e_arr[good]

        B0, *_ = lstsq(G_b_g, e_g, rcond=None)    # (10,3)

        if freeze_C:
            # only [B] is free; [C] is taken from the camera as-is.
            C = C_frozen

            def b_only_cost(b):
                B_ = b.reshape(n_b, 3)
                d = _dist_pt_to_line_plane(C, B_, G_c_g, G_b_g, X_g,
                                            O0, u_ax, v_ax)
                return npsum(d ** 2)

            res = minimize(b_only_cost, B0.ravel(), method='BFGS')
            B = res.x.reshape(n_b, 3)
        else:
            # (step 4) fit the in-plane deviations with the same normalized
            # polynomial basis used for the 'free' model, then jointly
            # refine C and B against the true point-to-line distance.
            C0, *_ = lstsq(G_c_g, ab_g, rcond=None)   # (n_c,2)

            a0 = hstack([C0.ravel(), B0.ravel()])

            def indirect_cost(a):
                C_ = a[:n_c * 2].reshape(n_c, 2)
                B_ = a[n_c * 2:].reshape(n_b, 3)
                d = _dist_pt_to_line_plane(C_, B_, G_c_g, G_b_g, X_g,
                                            O0, u_ax, v_ax)
                return npsum(d ** 2)

            res = minimize(indirect_cost, a0, method='BFGS')

            C = res.x[:n_c * 2].reshape(n_c, 2)
            B = res.x[n_c * 2:].reshape(n_b, 3)

        final_resid = _dist_pt_to_line_plane(C, B, G_c_g, G_b_g, X_g,
                                              O0, u_ax, v_ax)

        n_dropped = int((~good).sum())
        if n_dropped > 0:
            print(f'fit_variable_origin_model: dropped {n_dropped} points '
                  'whose rays were near-parallel to the origin plane.')

        plane = {'O0': O0, 'u': u_ax, 'v': v_ax, 'n': n_hat}

    n_rejected = N - len(keep_idx)
    if n_rejected > 0:
        print(f'fit_variable_origin_model: rejected {n_rejected}/{N} '
              f'points ({100*n_rejected/N:.1f}%) on reprojection residual.')

    mask = array([False] * N)
    mask[array(keep_idx)] = True

    return C, B, mask, final_resid, px_center, px_scale, plane, origin_model


def _mean_ray_err(ray_fn, X_arr, x_arr):
    '''Mean point-to-line distance (lab units) of ray_fn(x) -> (O, e)
    against the true 3D points X_arr, for a matching list of pixels x_arr.
    '''
    dists = []
    for Xi, xi in zip(X_arr, x_arr):
        O, e = ray_fn(xi)
        diff = array(Xi) - array(O)
        perp = diff - dot(diff, e) * e
        dists.append(norm(perp) ** 2)
    return (sum(dists) / len(dists)) ** 0.5


def cross_validate_origin_models(cam, X_list, x_list, test_fraction=0.25,
                                  thresh_pix=5.0, seed=None, c_order=3,
                                  origin_model='free', b_order=3):
    '''
    Diagnostic: splits the (already forward-model-fit, inlier) calibration
    points into a random TRAIN/TEST split, fits BOTH the original
    fixed-origin backward model and the variable-origin extension
    (with the given c_order, see fit_variable_origin_model) on TRAIN only,
    and reports the mean point-to-line ("ray") error of each on both
    TRAIN and TEST.

    This directly answers "is the variable-origin extension actually
    generalizing to pixels it wasn't fit on, or just overfitting the
    calibration points?" - if variable_test is similar to or worse than
    fixed_test (even though variable_train looks great), that's overfitting
    / extrapolation, not a genuine improvement, and you likely need more/
    better-distributed calibration points, a lower polynomial order for
    [C], or to stick with the fixed-origin model for this camera.

    Assumes cam.A is already fit (this function uses cam.projection
    internally and does not refit the forward model).

    Returns a dict:
        fixed_train, fixed_test, variable_train, variable_test - mean ray
            errors (lab units) for each model on each split.
        n_train, n_test - number of points in each split.
    '''
    from numpy.random import default_rng
    from myptv.extendedZolof.calibrate_step2_improved import \
        step2_estimate_camera_center

    rng = default_rng(seed)

    X_list = array(X_list, dtype=float)
    x_list = array(x_list, dtype=float)
    N = len(X_list)

    idx = rng.permutation(N)
    n_test = max(1, int(round(N * test_fraction)))
    test_idx = idx[:n_test]
    train_idx = idx[n_test:]

    X_train, x_train = X_list[train_idx], x_list[train_idx]
    X_test, x_test = X_list[test_idx], x_list[test_idx]

    # --- fixed-origin model, fit on TRAIN only ---
    O_train = step2_estimate_camera_center(
        cam.projection, x_train, X_train,
        thresh_pix=thresh_pix, refine=True)

    r_train = [(Xi - O_train) / norm(Xi - O_train) for Xi in X_train]
    n_b_cv = _ORDER_N_TERMS[b_order]
    Gb = array([cam.get_xCol(xi, n=n_b_cv) for xi in x_train])
    B_train, *_ = lstsq(Gb, r_train, rcond=None)

    def fixed_ray(x):
        xc = array(cam.get_xCol(x, n=n_b_cv))
        e = dot(xc, B_train)
        e = e / norm(e)
        return O_train, e

    fixed_train_err = _mean_ray_err(fixed_ray, X_train, x_train)
    fixed_test_err = _mean_ray_err(fixed_ray, X_test, x_test)

    # --- variable-origin model, fit on TRAIN only ---
    var_train_err = var_test_err = None
    try:
        C_t, B_t, _, _, px_c, px_s, plane_t, _om = fit_variable_origin_model(
            cam, X_train, x_train, O_init=O_train, thresh_pix=thresh_pix,
            min_points=max(15, int(0.5 * len(train_idx))), c_order=c_order,
            origin_model=origin_model, b_order=b_order)

        n_c = C_t.shape[0]

        def var_ray(x):
            xc = array(_get_xCol_scaled(x, px_c, px_s,
                                        n=max(B_t.shape[0], C_t.shape[0])))
            e = dot(xc[:B_t.shape[0]], B_t)
            e = e / norm(e)
            coef = dot(xc[:n_c], C_t)
            if origin_model == 'plane':
                o = (plane_t['O0'] + plane_t['u']*coef[0]
                     + plane_t['v']*coef[1])
            else:
                o = coef
            return o, e

        var_train_err = _mean_ray_err(var_ray, X_train, x_train)
        var_test_err = _mean_ray_err(var_ray, X_test, x_test)
    except ValueError as err:
        print(f'cross_validate_origin_models: variable-origin fit failed '
              f'on the train split ({err}); likely too few training '
              'points for a meaningful comparison.')

    result = {
        'fixed_train': fixed_train_err, 'fixed_test': fixed_test_err,
        'variable_train': var_train_err, 'variable_test': var_test_err,
        'n_train': len(train_idx), 'n_test': len(test_idx),
    }

    print('cross_validate_origin_models:')
    print(f'  n_train={result["n_train"]}, n_test={result["n_test"]}')
    print(f'  fixed-origin : train={fixed_train_err:.4f}  '
          f'test={fixed_test_err:.4f}')
    if var_train_err is not None:
        print(f'  variable-origin ({origin_model}, c_order={c_order}, '
              f'b_order={b_order}) : '
              f'train={var_train_err:.4f}  test={var_test_err:.4f}')
        if var_test_err > fixed_test_err:
            print('  -> variable-origin generalizes WORSE than fixed-origin '
                  'on held-out points (even if train error looked better). '
                  'This points to overfitting/extrapolation, not a genuine '
                  'improvement.')
        else:
            print('  -> variable-origin generalizes at least as well as '
                  'fixed-origin on held-out points.')

    return result
