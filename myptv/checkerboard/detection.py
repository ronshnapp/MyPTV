# -*- coding: utf-8 -*-
#!/usr/bin/env python3

"""
Created on Aug 27 2026

@author: Stefano Brizzolara


Checkerboard corner detection.

This module locates the internal corners of a checkerboard (chessboard)
calibration target in an image, with sub-pixel accuracy. It is written using
only the packages that MyPTV already depends on, namely numpy, scipy and
scikit-image; in particular it does not require OpenCV.


The detection is done in four steps:


1) Candidate detection.

   An internal checkerboard corner is a saddle point of the image intensity:
   along one diagonal the intensity increases, along the other it decreases.
   Saddle points are found from the Hessian matrix of the (Gaussian smoothed)
   image,

       H = [[Ixx, Ixy], [Ixy, Iyy]] ,

   whose determinant is negative at a saddle. We therefore use

       (1)    R = -det(H) * sigma^4 = -(Ixx*Iyy - Ixy^2) * sigma^4

   as the corner response, where the sigma^4 factor makes R comparable
   between different smoothing scales. Candidates are the local maxima of R.


2) Sub-pixel refinement.

   Let c be the true corner position. Any pixel p that lies on one of the two
   checkerboard edges passing through c has an intensity gradient g(p) that is
   perpendicular to that edge, and hence perpendicular to (p - c):

       (2)    g(p) . (p - c) = 0 .

   Minimizing the sum of the squared residuals of (2) over a window of pixels
   around the corner gives a 2x2 linear system for c,

       (3)    [sum_p g g^T] c = sum_p g (g . p) ,

   which is solved directly. Pixels in flat regions have small |g| and so
   contribute little, while pixels on the two edges dominate the fit. The
   procedure is iterated a few times, re-centering the window each time.


3) Edge orientations.

   The two edge directions of a corner are recovered from a magnitude
   weighted histogram of the gradient orientations in a window around the
   corner. Since an edge direction is defined modulo pi, the histogram is
   built over [0, pi). Its two dominant modes give the two edge normals, and
   rotating them by 90 degrees gives the two edge directions. These are used
   both for scoring the corner and, later, as the local lattice directions
   when the corners are assembled into a grid (see grid.py).


4) Scoring.

   The response (1) also fires on saddle-like features that are not
   checkerboard corners. Each candidate is therefore scored by how well its
   surroundings match an ideal checkerboard corner with the edge directions
   d1 and d2 found in step 3. Writing u = (p - c).d1 and v = (p - c).d2, the
   ideal pattern alternates sign between the four quadrants,

       (4)    T(p) = sign(u) * sign(v) ,

   and the score is the absolute value of the weighted normalized
   correlation between the image and T. The absolute value is used because
   the polarity of the corner (which pair of quadrants is bright) is
   arbitrary. The score lies in [0, 1] and candidates scoring below a
   threshold are discarded.

"""

from numpy import (array, zeros, arange, meshgrid, exp, sqrt, sin, cos,
                   arctan2, pi, abs as npabs, sign, sum as npsum, argmax,
                   roll, isfinite, float64, minimum, asarray)
from numpy.linalg import solve, LinAlgError

from scipy.ndimage import gaussian_filter




def load_image(fname):
    '''
    Reads an image from the disk and returns it as a 2D float array scaled
    to the range [0, 1].

    input -
    fname (string) - path of the image file

    output -
    img (array, (Ny,Nx)) - the grayscale image
    '''
    from skimage.io import imread

    img = asarray(imread(fname), dtype=float64)

    # if the image has a colour (or alpha) axis, average the colour channels
    if img.ndim == 3:
        img = img[:, :, :3].mean(axis=2)
    elif img.ndim != 2:
        raise ValueError('cannot interpret an image with %d axes'%img.ndim)

    span = img.max() - img.min()
    if span == 0:
        raise ValueError('the image %s has no contrast'%fname)

    return (img - img.min()) / span




def saddle_response(img, sigma=2.0):
    '''
    Returns the checkerboard corner response R of eq. (1), which is positive
    at saddle points of the image intensity.

    input -
    img (array, (Ny,Nx)) - grayscale image
    sigma (float) - standard deviation, in pixels, of the Gaussian used to
                    smooth the image before differentiating it. It should be
                    of the order of the width of the blur of the
                    checkerboard edges.

    output -
    R (array, (Ny,Nx)) - the corner response
    '''
    # note on the axis order: axis 0 of the image is the row index (y) and
    # axis 1 is the column index (x), so the "order" tuples below read
    # (order in y, order in x).
    Ixx = gaussian_filter(img, sigma, order=(0, 2), mode='nearest')
    Iyy = gaussian_filter(img, sigma, order=(2, 0), mode='nearest')
    Ixy = gaussian_filter(img, sigma, order=(1, 1), mode='nearest')

    return -(Ixx*Iyy - Ixy**2) * sigma**4




def image_gradients(img, sigma=1.0):
    '''
    Returns the two components of the Gaussian smoothed image gradient.

    input -
    img (array, (Ny,Nx)) - grayscale image
    sigma (float) - standard deviation of the Gaussian, in pixels

    output -
    gx, gy (arrays, (Ny,Nx)) - derivatives with respect to x (the column
                               index) and y (the row index)
    '''
    gx = gaussian_filter(img, sigma, order=(0, 1), mode='nearest')
    gy = gaussian_filter(img, sigma, order=(1, 0), mode='nearest')
    return gx, gy




def find_candidates(img, sigma=2.0, min_distance=6, threshold_rel=0.05,
                    max_candidates=20000):
    '''
    Finds candidate checkerboard corners as the local maxima of the response
    of eq. (1).

    input -
    img (array, (Ny,Nx)) - grayscale image
    sigma (float) - smoothing scale of the response, in pixels
    min_distance (int) - minimum separation, in pixels, between two maxima.
                         It should be smaller than the spacing of the
                         checkerboard corners in the image.
    threshold_rel (float) - maxima weaker than this fraction of the strongest
                            maximum are ignored
    max_candidates (int) - a safety limit on the number of returned points

    output -
    candidates (array, (N,2)) - the (x, y) pixel coordinates of the candidates
    '''
    from skimage.feature import peak_local_max

    R = saddle_response(img, sigma=sigma)

    # peak_local_max returns (row, col) pairs, which we flip to (x, y)
    peaks = peak_local_max(R, min_distance=int(min_distance),
                           threshold_rel=threshold_rel,
                           exclude_border=int(min_distance),
                           num_peaks=int(max_candidates))

    if len(peaks) == 0:
        return zeros((0, 2))

    return array([[float(p[1]), float(p[0])] for p in peaks])




def refine_corner(gx, gy, c0, half_window=6, n_iter=4, max_shift=None):
    '''
    Refines the position of a single corner to sub-pixel accuracy by solving
    the linear system of eq. (3).

    input -
    gx, gy (arrays, (Ny,Nx)) - image gradients, as given by image_gradients
    c0 (array, 2) - the initial (x, y) estimate of the corner position
    half_window (int) - half the side of the square window of pixels used in
                        the fit. It should be somewhat smaller than the
                        spacing of the checkerboard corners in the image.
    n_iter (int) - number of refinement iterations
    max_shift (float) - if the refined position moves further than this from
                        c0 the refinement is considered to have failed and
                        None is returned. Defaults to half_window.

    output -
    c (array, 2) - the refined (x, y) position, or None if the refinement
                   failed
    '''
    Ny, Nx = gx.shape
    w = int(half_window)

    if max_shift is None:
        max_shift = float(half_window)

    x, y = float(c0[0]), float(c0[1])

    for _ in range(int(n_iter)):

        ix, iy = int(round(x)), int(round(y))

        if ix-w < 0 or iy-w < 0 or ix+w+1 > Nx or iy+w+1 > Ny:
            return None

        GX = gx[iy-w:iy+w+1, ix-w:ix+w+1]
        GY = gy[iy-w:iy+w+1, ix-w:ix+w+1]

        X, Y = meshgrid(arange(ix-w, ix+w+1, dtype=float64),
                        arange(iy-w, iy+w+1, dtype=float64))

        # a Gaussian window that suppresses the influence of distant pixels
        s = max(w/2.0, 1.0)
        W = exp(-((X - x)**2 + (Y - y)**2) / (2*s**2))

        gdotp = GX*X + GY*Y

        A = array([[npsum(W*GX*GX), npsum(W*GX*GY)],
                   [npsum(W*GX*GY), npsum(W*GY*GY)]])
        b = array([npsum(W*GX*gdotp), npsum(W*GY*gdotp)])

        try:
            c = solve(A, b)
        except LinAlgError:
            return None

        if not isfinite(c).all():
            return None

        x, y = float(c[0]), float(c[1])

        if (x - c0[0])**2 + (y - c0[1])**2 > max_shift**2:
            return None

    return array([x, y])




def corner_edge_directions(gx, gy, c, half_window=6, n_bins=36,
                           min_separation=40.0):
    '''
    Estimates the two edge directions of a checkerboard corner from a
    magnitude weighted histogram of the gradient orientations around it.

    input -
    gx, gy (arrays, (Ny,Nx)) - image gradients
    c (array, 2) - the (x, y) position of the corner
    half_window (int) - half the side of the window of pixels used
    n_bins (int) - number of histogram bins spanning [0, pi)
    min_separation (float) - minimum angle, in degrees, between the two
                             dominant modes of the histogram

    output -
    d1, d2 (arrays, 2) - unit vectors along the two edge directions, or
                         (None, None) if two modes could not be found
    '''
    Ny, Nx = gx.shape
    w = int(half_window)
    ix, iy = int(round(c[0])), int(round(c[1]))

    if ix-w < 0 or iy-w < 0 or ix+w+1 > Nx or iy+w+1 > Ny:
        return None, None

    GX = gx[iy-w:iy+w+1, ix-w:ix+w+1].ravel()
    GY = gy[iy-w:iy+w+1, ix-w:ix+w+1].ravel()

    mag = sqrt(GX**2 + GY**2)
    if mag.max() <= 0:
        return None, None

    # the orientation of an edge is defined modulo pi
    ang = arctan2(GY, GX) % pi

    bins = minimum((ang/pi * n_bins).astype(int), n_bins-1)

    hist = zeros(n_bins)
    for b, m in zip(bins, mag):
        hist[b] += m

    # smooth the histogram, wrapping around at 0 and pi
    hist = 0.5*roll(hist, -1) + hist + 0.5*roll(hist, 1)

    i1 = int(argmax(hist))

    # mask out the bins that are too close to the first mode, using the
    # circular distance modulo pi
    sep_bins = min_separation/180.0 * n_bins
    masked = hist.copy()
    for i in range(n_bins):
        d = npabs(i - i1)
        d = min(d, n_bins - d)
        if d < sep_bins:
            masked[i] = -1.0

    if masked.max() <= 0:
        return None, None

    i2 = int(argmax(masked))

    # bin centres, which are the orientations of the two edge normals
    n1 = (i1 + 0.5)/n_bins * pi
    n2 = (i2 + 0.5)/n_bins * pi

    # the edge directions are the normals rotated by 90 degrees
    d1 = array([cos(n1 + pi/2), sin(n1 + pi/2)])
    d2 = array([cos(n2 + pi/2), sin(n2 + pi/2)])

    return d1, d2




def corner_score(img, c, d1, d2, half_window=6, sigma=2.0):
    '''
    Scores how well the neighbourhood of a candidate matches an ideal
    checkerboard corner, using the template of eq. (4).

    input -
    img (array, (Ny,Nx)) - grayscale image
    c (array, 2) - the (x, y) position of the corner
    d1, d2 (arrays, 2) - the two edge directions
    half_window (int) - half the side of the window of pixels used
    sigma (float) - the blur scale of the image; pixels closer than this to
                    the corner carry little information and are down-weighted

    output -
    score (float) - a number in [0, 1]; a perfect checkerboard corner scores 1
    '''
    Ny, Nx = img.shape
    w = int(half_window)
    ix, iy = int(round(c[0])), int(round(c[1]))

    if ix-w < 0 or iy-w < 0 or ix+w+1 > Nx or iy+w+1 > Ny:
        return 0.0

    I = img[iy-w:iy+w+1, ix-w:ix+w+1]

    X, Y = meshgrid(arange(ix-w, ix+w+1, dtype=float64),
                    arange(iy-w, iy+w+1, dtype=float64))
    dx, dy = X - c[0], Y - c[1]

    u = dx*d1[0] + dy*d1[1]
    v = dx*d2[0] + dy*d2[1]
    T = sign(u) * sign(v)

    # weights: a Gaussian envelope, with the blurred core of the corner
    # suppressed
    r2 = dx**2 + dy**2
    s = max(w/2.0, 1.0)
    W = exp(-r2/(2*s**2)) * (1.0 - exp(-r2/(2*max(sigma, 0.5)**2)))

    Wsum = npsum(W)
    if Wsum <= 0:
        return 0.0

    Im = npsum(W*I)/Wsum
    Tm = npsum(W*T)/Wsum

    num = npsum(W*(I - Im)*(T - Tm))
    den = sqrt(npsum(W*(I - Im)**2) * npsum(W*(T - Tm)**2))

    if den <= 0:
        return 0.0

    return float(min(abs(num/den), 1.0))




def detect_corners(img, sigma=2.0, min_distance=6, half_window=None,
                   threshold_rel=0.05, min_score=0.7, max_candidates=20000):
    '''
    The full checkerboard corner detection pipeline: candidates are found
    from the saddle response, refined to sub-pixel accuracy, given edge
    directions, scored, and the poorly scoring ones are discarded.

    input -
    img (array, (Ny,Nx)) - grayscale image, e.g. as returned by load_image
    sigma (float) - smoothing scale, in pixels, of the corner response. This
                    is the main tuning parameter: it should be of the order
                    of the blur of the checkerboard edges.
    min_distance (int) - minimum separation, in pixels, of two corners
    half_window (int) - half the side of the windows used for the refinement,
                        the edge directions and the score. Defaults to
                        min_distance.
    threshold_rel (float) - relative threshold on the corner response
    min_score (float) - candidates whose score of eq. (4) is below this value
                        are discarded. The default of 0.7 was chosen from
                        tests on synthetic boards, in which true internal
                        corners scored above 0.85 while the spurious
                        candidates found along the outer boundary of the
                        pattern, where only two of the four quadrants belong
                        to the board, stayed below 0.55.
    max_candidates (int) - a safety limit on the number of candidates

    output - a dictionary with the keys:
    'points' (array, (N,2)) - the refined (x, y) corner positions
    'd1', 'd2' (arrays, (N,2)) - the two edge directions of each corner
    'scores' (array, N) - the score of each corner
    '''
    if half_window is None:
        half_window = int(min_distance)

    candidates = find_candidates(img, sigma=sigma, min_distance=min_distance,
                                 threshold_rel=threshold_rel,
                                 max_candidates=max_candidates)

    gx, gy = image_gradients(img, sigma=max(sigma/2.0, 0.6))

    points, D1, D2, scores = [], [], [], []

    for c0 in candidates:

        c = refine_corner(gx, gy, c0, half_window=half_window)
        if c is None:
            continue

        d1, d2 = corner_edge_directions(gx, gy, c, half_window=half_window)
        if d1 is None:
            continue

        s = corner_score(img, c, d1, d2, half_window=half_window, sigma=sigma)
        if s < min_score:
            continue

        points.append(c)
        D1.append(d1)
        D2.append(d2)
        scores.append(s)

    if len(points) == 0:
        return {'points': zeros((0, 2)), 'd1': zeros((0, 2)),
                'd2': zeros((0, 2)), 'scores': zeros(0)}

    return {'points': array(points), 'd1': array(D1), 'd2': array(D2),
            'scores': array(scores)}
