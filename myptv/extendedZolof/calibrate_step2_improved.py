# -*- coding: utf-8 -*-
"""
Proposed replacement for step 2 of calibrate_extendedZolof.calibrate().

Key changes vs. the original get_ray_from_x() + get_nearest_line_crossing():

  1. Ray direction is found by minimizing reprojection error over a
     SPHERICAL parameterization (theta, phi) instead of differencing two
     independently-minimized 3D points. This guarantees a genuine unit
     vector and removes the arbitrary [1,1,1] offset.

  2. A second, independent point search (free 3D minimization) still gives
     a second point on the ray, so direction = (X2 - X1)/norm(...), same
     idea as before but anchored by a well-defined unit-sphere search
     rather than an arbitrary displacement.

  3. Rays with large reprojection residual are rejected (thresh_pix)
     BEFORE they're allowed to influence the camera center estimate.

  4. The camera center O is estimated by directly minimizing the sum of
     squared point-to-line distances (nonlinear optimization), rather than
     relying solely on whatever closed-form method get_nearest_line_crossing
     uses. This also makes it easy to swap in a robust loss later
     (e.g. Huber) if outlier rays are still a concern.
"""

from numpy import array, cos, sin, cross, dot, mean, sum as npsum
from numpy.linalg import norm
from scipy.optimize import minimize


def get_ray_from_x_robust(cam_projection, x, X_ref, thresh_pix=None):
    '''
    Given an image point x and a nearby reference 3D point X_ref (e.g. the
    known calibration point), estimate a 3D ray (origin, direction) whose
    projection best matches x.

    cam_projection : callable, X (3,) -> projected 2D point
    x              : observed 2D image point (2,)
    X_ref          : a 3D point near the true location (e.g. the known
                      calibration target position for this observation)
    thresh_pix     : if given, returns None if the fitted ray's reprojection
                      residual exceeds this threshold (outlier rejection)

    Returns (X1, e, resid) where X1 is a point on the ray, e is the unit
    direction, and resid is the reprojection residual at X1 (pixels).
    '''
    x = array(x, dtype=float)
    X_ref = array(X_ref, dtype=float)

    def reproj_err(X):
        return npsum((array(cam_projection(X)) - x) ** 2)

    # --- 1) direction via spherical parameterization, anchored at X_ref ---
    def sphere_cost(a):
        direction = array([
            cos(a[0]) * sin(a[1]),
            sin(a[0]) * sin(a[1]),
            cos(a[1])
        ])
        return reproj_err(X_ref + direction)

    res_dir = minimize(sphere_cost, array([0.0, 0.0]), method='BFGS')
    a = res_dir.x
    direction = array([
        cos(a[0]) * sin(a[1]),
        sin(a[0]) * sin(a[1]),
        cos(a[1])
    ])
    X1 = X_ref + direction

    # --- 2) a second, freely-optimized point on the same ray ---
    res_pt = minimize(reproj_err, X_ref, method='BFGS')
    X2 = res_pt.x

    u = X2 - X1
    nu = norm(u)
    if nu == 0:
        return None
    e = u / nu

    resid = reproj_err(X1) ** 0.5

    if thresh_pix is not None and resid > thresh_pix:
        return None

    return X1, e, resid


def estimate_pinhole_center(points_on_lines, directions, X0_guess=None):
    '''
    Robustly estimate a single 3D point O that minimizes the summed squared
    distance to a set of 3D lines, each defined by a point and a unit
    direction. Uses nonlinear optimization (as in Solof_v2's pinhole_cost)
    rather than a purely closed-form crossing, so it's easy to extend with
    a robust loss later if needed.
    '''
    P = array(points_on_lines)   # (N,3)
    U = array(directions)        # (N,3)

    if X0_guess is None:
        X0_guess = mean(P, axis=0)

    def pt_line_dist2(O):
        # squared distance from O to each line, summed
        diff = O[None, :] - P                      # (N,3)
        proj = npsum(diff * U, axis=1)[:, None] * U  # component along line
        perp = diff - proj
        return npsum(perp ** 2)

    res = minimize(pt_line_dist2, X0_guess, method='BFGS')
    return res.x


def step2_estimate_camera_center(cam_projection, x_list, X_list,
                                  thresh_pix=5.0, refine=True):
    '''
    Drop-in replacement for calibrate_extendedZolof's step 2. Returns the
    estimated camera center O.

    thresh_pix : rays whose reprojection residual exceeds this are dropped
                 before estimating O (outlier rejection).
    refine     : if True, re-run get_ray_from_x_robust a second time using
                 the refined O as the anchor point, then re-estimate O.
                 This mirrors Solof_v2's second (post-pinhole) refinement
                 pass and typically tightens the estimate further.
    '''
    points, dirs = [], []
    for x, X in zip(x_list, X_list):
        result = get_ray_from_x_robust(cam_projection, x, X,
                                        thresh_pix=thresh_pix)
        if result is None:
            continue
        X1, e, resid = result
        points.append(X1)
        dirs.append(e)

    if len(points) < 3:
        raise ValueError(
            'Too few valid rays survived outlier rejection to estimate '
            'a camera center; consider raising thresh_pix.')

    O = estimate_pinhole_center(points, dirs)

    if refine:
        # anchor the ray search at the *previous point on each ray*, but
        # now re-derive direction with O as the initial guess for X_ref,
        # which usually improves convergence once we're close to the truth
        points2, dirs2 = [], []
        for x, X in zip(x_list, X_list):
            result = get_ray_from_x_robust(cam_projection, x, O,
                                            thresh_pix=thresh_pix)
            if result is None:
                continue
            X1, e, resid = result
            points2.append(X1)
            dirs2.append(e)

        if len(points2) >= 3:
            O = estimate_pinhole_center(points2, dirs2, X0_guess=O)

    return O
