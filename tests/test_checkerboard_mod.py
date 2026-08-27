# -*- coding: utf-8 -*-
"""
Created on Aug 27 2026

@author: Stefano Brizzolara


Tests for the checkerboard calibration procedure.

The tests work on synthetic images, which makes the truth available exactly:
the images are rendered by projecting a checkerboard through a camera_Tsai
whose parameters are known, so the pixel position at which every corner
ought to be found is known to machine precision, and so are the camera
parameters that the calibration ought to recover.

Rendering is done through the homography that maps the plane of the board
onto the image. That map is exact for a pin-hole camera, which is what
camera_Tsai reduces to when its polynomial correction term E vanishes.

"""

import os
from tempfile import TemporaryDirectory

from numpy import (array, asarray, zeros, cross, floor, where, ones, hstack,
                   append, sqrt, meshgrid, arange, clip, mean, square, sin,
                   cos)
from numpy.linalg import solve, inv, norm
from numpy.random import default_rng

from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree

from myptv.TsaiModel.camera import camera_Tsai
from myptv.TsaiModel.calibrate import calibrate_Tsai
from myptv.utils import Cal_image_coord

from myptv.checkerboard.detection import detect_corners
from myptv.checkerboard.grid import assemble_grid
from myptv.checkerboard.board import (resolve_orientation,
                                      board_lab_coordinates)
from myptv.checkerboard.cal_points import checkerboard_cal_points




# the geometry shared by the tests: a camera looking along +z from a distance
# DIST at a board placed near the origin of the lab frame
RESOLUTION = (640, 512)
DISTANCE = 300.0
SQUARE = 10.0




# =============================================================================
#   helpers that render synthetic images of a checkerboard
# =============================================================================


def make_camera(name='synth_cam', O=(0., 0., -DISTANCE),
                theta=(0.06, -0.09, 0.03), f=1200.0, resolution=RESOLUTION):
    '''
    Returns a camera_Tsai with the given parameters and no polynomial
    correction, so that it is an exact pin-hole camera.
    '''
    cam = camera_Tsai(name)
    cam.O = array(O, dtype=float)
    cam.theta = array(theta, dtype=float)
    cam.f = float(f)
    cam.resolution = resolution
    cam.xh, cam.yh = 0.0, 0.0
    cam.E = zeros((3, 5))
    cam.calc_R()
    return cam


def tilted_axes(alpha=0.0, beta=0.0):
    '''
    Returns an orthonormal pair of lab space directions for the two sides of
    the board, tilted out of the plane z = 0. With alpha = beta = 0 the board
    lies in the x-y plane and faces the camera; increasing either angle
    tilts it away, so that it is seen obliquely and its image is
    foreshortened.
    '''
    e1 = array([cos(beta), 0.0, sin(beta)])
    e2 = array([sin(alpha)*sin(beta), cos(alpha), -sin(alpha)*cos(beta)])
    return e1, e2


def centered_origin(board_size, square, e1, e2):
    '''
    Returns the lab position of the board corner (0,0) that places the centre
    of the board at the origin of the lab frame, and hence on the axis of the
    camera.
    '''
    return -((board_size[0] - 1)/2.0*square*asarray(e1) +
             (board_size[1] - 1)/2.0*square*asarray(e2))


def board_view(board_size, square=SQUARE, alpha=0.0, beta=0.0, frac=0.55):
    '''
    Sets up a camera and a board pose such that the board is centred in the
    frame and its longer side spans a fraction frac of the image width. The
    focal length has to be chosen per board, since a board that is too large
    for the frame would have some of its corners cut off, and one that is too
    small would have its corners too close together to be resolved.

    output -
    cam - the camera
    origin, e1, e2 - the pose of the board, as taken by
                     board.board_lab_coordinates
    '''
    e1, e2 = tilted_axes(alpha, beta)
    origin = centered_origin(board_size, square, e1, e2)

    # the pattern spans one extra square beyond the internal corners on each
    # side, so its full extent is (max(board_size)-1+2) squares
    span = (max(board_size) + 1)*square
    f = frac*RESOLUTION[0]*DISTANCE/span

    return make_camera(f=f), origin, e1, e2


def _homography(src, dst):
    '''
    The 3x3 homography taking the four points src to the four points dst.
    '''
    A, b = [], []
    for (x, y), (X, Y) in zip(src, dst):
        A.append([x, y, 1, 0, 0, 0, -X*x, -X*y]); b.append(X)
        A.append([0, 0, 0, x, y, 1, -Y*x, -Y*y]); b.append(Y)
    h = solve(array(A, dtype=float), array(b, dtype=float))
    return append(h, 1.0).reshape(3, 3)


def render_board(cam, board_size, square, origin, i_axis, j_axis,
                 translation=0.0, blur=1.5, supersample=3, noise=0.0,
                 seed=0, dark=0.15, light=0.85, surround=1.0, background=0.5):
    '''
    Renders an image of a checkerboard as seen by cam, and returns it
    together with the exact pixel positions of the internal corners.

    input -
    cam - a camera_Tsai with E = 0
    board_size (tuple, 2) - the number of internal corners
    square (float) - the spacing of the corners in lab units
    origin, i_axis, j_axis - the pose of the board, as in
                             board.board_lab_coordinates
    translation - the displacement of the board, either a distance along its
                  normal or a vector of three components
    blur (float) - the width of the Gaussian blur applied to the image
    supersample (int) - the factor by which the image is over-sampled before
                        being averaged down, which anti-aliases the edges.
                        This is what limits how exactly the rendered image
                        matches the ideal one, and hence how closely the
                        corners can be recovered: measured against the exact
                        projection, the recovered positions have an RMS
                        error of about 0.47, 0.17, 0.07 and 0.04 pixels at
                        supersample = 1, 2, 3 and 4. The tests that measure
                        the accuracy of the detector therefore render at 4,
                        where the error is dominated by the image noise
                        rather than by the rendering.
    noise (float) - standard deviation of the Gaussian noise added
    surround (float) - the width, in squares, of the light border drawn
                       around the pattern, as a printed target has
    background (float) - the grey level outside that border

    output -
    img (array, (Ny,Nx)) - the rendered image, in [0, 1]
    truth (array, (n_i,n_j,2)) - the exact pixel position of each internal
                                 corner
    '''
    ni, nj = board_size
    s = float(square)
    ss = int(supersample)

    e1 = asarray(i_axis, dtype=float); e1 = e1/norm(e1)
    e2 = asarray(j_axis, dtype=float); e2 = e2/norm(e2)

    t = asarray(translation, dtype=float).ravel()
    if t.size == 1:
        n = cross(e1, e2)
        t = float(t[0]) * n/norm(n)

    X0 = asarray(origin, dtype=float) + t

    def board_to_lab(u, v):
        return X0 + u*e1 + v*e2

    truth = zeros((ni, nj, 2))
    for i in range(ni):
        for j in range(nj):
            truth[i, j] = cam.projection(board_to_lab(i*s, j*s))

    # the pattern spans one extra square beyond the internal corners on every
    # side, so that those corners really are internal
    u0, u1 = -s, ni*s
    v0, v1 = -s, nj*s

    quad_b = [(u0, v0), (u1, v0), (u1, v1), (u0, v1)]
    quad_i = [tuple(cam.projection(board_to_lab(u, v))) for u, v in quad_b]

    Hinv = inv(_homography(quad_b, quad_i))

    Nx, Ny = cam.resolution
    xs = (arange(Nx*ss) + 0.5)/ss - 0.5
    ys = (arange(Ny*ss) + 0.5)/ss - 0.5
    yy, xx = meshgrid(ys, xs, indexing='ij')

    p = hstack([xx.reshape(-1, 1), yy.reshape(-1, 1), ones((xx.size, 1))])
    q = p.dot(Hinv.T)
    u = q[:, 0]/q[:, 2]
    v = q[:, 1]/q[:, 2]

    on_board = (u >= u0) & (u <= u1) & (v >= v0) & (v <= v1)
    on_border = ((u >= u0 - surround*s) & (u <= u1 + surround*s) &
                 (v >= v0 - surround*s) & (v <= v1 + surround*s))

    pattern = (floor(u/s).astype(int) + floor(v/s).astype(int)) % 2

    img = where(on_board, where(pattern == 0, light, dark),
                where(on_border, light, background))
    img = img.reshape(Ny*ss, Nx*ss).reshape(Ny, ss, Nx, ss).mean(axis=(1, 3))

    img = gaussian_filter(img, blur)

    if noise > 0:
        img = img + default_rng(seed).normal(0, noise, img.shape)

    return clip(img, 0, 1), truth


def _save(img, fname):
    '''
    Writes a rendered image to the disk as a 16 bit tif.
    '''
    from skimage.io import imsave
    imsave(fname, (img*65535).astype('uint16'), check_contrast=False)




# =============================================================================
#   the tests
# =============================================================================


def test_detection_accuracy():
    '''
    Every internal corner of the board should be found, and found to well
    within a tenth of a pixel of where it is known to project.
    '''
    board_size = (9, 6)
    cam, origin, e1, e2 = board_view(board_size)

    img, truth = render_board(cam, board_size, SQUARE, origin, e1, e2,
                              noise=0.01, supersample=4)

    det = detect_corners(img, sigma=2.0, min_distance=8)

    flat = truth.reshape(-1, 2)
    d, _ = cKDTree(det['points']).query(flat)

    assert (d < 1.0).all(), 'corners missed: %d'%int((d >= 1.0).sum())

    # at supersample = 4 and this noise level the observed error is about
    # 0.055 px, most of which is the rendering rather than the detector
    rms = sqrt(mean(square(d)))
    assert rms < 0.12, 'corner position RMS error too large: %.4f px'%rms


def test_detection_rejects_pattern_boundary():
    '''
    The junctions on the outer boundary of the pattern are not surrounded by
    the board on all sides and must not be reported as corners, since they
    would extend the lattice past the declared size of the board.
    '''
    board_size = (9, 6)
    cam, origin, e1, e2 = board_view(board_size)

    img, truth = render_board(cam, board_size, SQUARE, origin, e1, e2,
                              noise=0.01)

    det = detect_corners(img, sigma=2.0, min_distance=8)

    d, _ = cKDTree(truth.reshape(-1, 2)).query(det['points'])
    assert (d < 1.0).all(), \
        '%d corners were found off the internal grid'%int((d >= 1.0).sum())


def test_grid_assembly():
    '''
    The detected corners should be assembled into a complete lattice of the
    right shape, for square and oblong boards and for boards seen at an
    angle. The lattice may come out transposed with respect to the board,
    which is the ambiguity that resolve_orientation exists to settle, so only
    the shape up to a transposition is checked here.
    '''
    for board_size in [(9, 6), (7, 7), (11, 4)]:
        for alpha, beta in [(0.0, 0.0), (0.5, -0.6)]:

            cam, origin, e1, e2 = board_view(board_size, alpha=alpha,
                                             beta=beta)

            img, truth = render_board(cam, board_size, SQUARE, origin, e1, e2,
                                      noise=0.01)

            det = detect_corners(img, sigma=2.0, min_distance=8)
            idx = assemble_grid(det)

            assert sorted(idx.shape) == sorted(board_size), \
                'board %s at tilt (%.1f,%.1f) gave a lattice of shape %s' \
                %(board_size, alpha, beta, idx.shape)

            assert (idx >= 0).sum() == board_size[0]*board_size[1], \
                'board %s at tilt (%.1f,%.1f): only %d of %d sites filled' \
                %(board_size, alpha, beta, int((idx >= 0).sum()),
                  board_size[0]*board_size[1])


def test_orientation_follows_the_hint():
    '''
    The origin hint must decide which corner of the board is taken as its
    origin, and must do so from a click that is several pixels off. Hinting
    at two different physical corners must give two different, and correct,
    answers.
    '''
    board_size = (9, 6)
    cam, origin, e1, e2 = board_view(board_size)

    img, truth = render_board(cam, board_size, SQUARE, origin, e1, e2,
                              noise=0.01)

    det = detect_corners(img, sigma=2.0, min_distance=8)
    idx = assemble_grid(det)
    P = det['points']

    # a sloppy click on the true (0,0) corner
    g = resolve_orientation(P, idx, board_size,
                            origin_hint=truth[0, 0] + array([7., -6.]))
    assert norm(P[g[0, 0]] - truth[0, 0]) < 1.0
    assert norm(P[g[-1, 0]] - truth[-1, 0]) < 1.0
    assert norm(P[g[0, -1]] - truth[0, -1]) < 1.0

    # hinting at the opposite corner of the board must relabel the lattice
    g2 = resolve_orientation(P, idx, board_size,
                             origin_hint=truth[-1, -1] + array([-5., 5.]))
    assert norm(P[g2[0, 0]] - truth[-1, -1]) < 1.0
    assert norm(P[g2[-1, 0]] - truth[0, -1]) < 1.0


def test_orientation_rejects_wrong_board_size():
    '''
    Declaring a board size that the image does not show must raise, rather
    than quietly returning a mis-indexed lattice, since a mis-indexed
    lattice attaches wrong lab coordinates to every corner.
    '''
    board_size = (9, 6)
    cam, origin, e1, e2 = board_view(board_size)

    img, truth = render_board(cam, board_size, SQUARE, origin, e1, e2,
                              noise=0.01)

    det = detect_corners(img, sigma=2.0, min_distance=8)
    idx = assemble_grid(det)

    raised = False
    try:
        resolve_orientation(det['points'], idx, (8, 5),
                            origin_hint=truth[0, 0])
    except ValueError:
        raised = True

    assert raised, 'a wrong board_size did not raise'


def test_board_lab_coordinates():
    '''
    The lab coordinates of the board corners must follow eq. (1) of board.py:
    a regular grid in the plane spanned by the two given axes, displaced
    along the normal.
    '''
    X = board_lab_coordinates((4, 3), 2.5, origin=(1., 2., 3.),
                              i_axis=(1., 0., 0.), j_axis=(0., 1., 0.),
                              translation=7.0)

    assert X.shape == (4, 3, 3)

    # the origin corner, displaced along the normal, which is +z here
    assert norm(X[0, 0] - array([1., 2., 10.])) < 1e-12

    # the spacing along the two axes
    assert norm(X[1, 0] - X[0, 0] - array([2.5, 0., 0.])) < 1e-12
    assert norm(X[0, 1] - X[0, 0] - array([0., 2.5, 0.])) < 1e-12

    # the far corner
    assert norm(X[3, 2] - array([1. + 3*2.5, 2. + 2*2.5, 10.])) < 1e-12

    # a full translation vector should be usable in place of a distance
    X2 = board_lab_coordinates((4, 3), 2.5, origin=(1., 2., 3.),
                               translation=(0., 0., 7.))
    assert norm(X2 - X) < 1e-12

    # two different spacings along the two sides
    X3 = board_lab_coordinates((3, 3), (2.0, 5.0))
    assert norm(X3[1, 0] - X3[0, 0] - array([2., 0., 0.])) < 1e-12
    assert norm(X3[0, 1] - X3[0, 0] - array([0., 5., 0.])) < 1e-12

    # a tilted board still has orthogonal, correctly spaced sides
    e1, e2 = tilted_axes(0.4, 0.5)
    X4 = board_lab_coordinates((3, 3), 2.0, i_axis=e1, j_axis=e2)
    assert abs(norm(X4[1, 0] - X4[0, 0]) - 2.0) < 1e-12
    assert abs(norm(X4[0, 1] - X4[0, 0]) - 2.0) < 1e-12
    assert abs((X4[1, 0] - X4[0, 0]).dot(X4[0, 1] - X4[0, 0])) < 1e-12


def test_end_to_end_calibration():
    '''
    The whole procedure, from images to a fitted camera.

    A board is photographed at three known positions of a translation stage
    by a camera whose parameters are known. The procedure writes a
    calibration points file, that file is read back through the same reader
    that the camera models use, and a fresh camera is fitted to it starting
    from a deliberately poor initial guess. The fit must recover the camera
    that produced the images.
    '''
    board_size = (9, 6)
    zs = [-20., 0., 20.]

    cam_true, origin, e1, e2 = board_view(board_size)

    with TemporaryDirectory() as tmp:

        files, truths = [], []
        for e, z in enumerate(zs):
            img, truth = render_board(cam_true, board_size, SQUARE, origin,
                                      e1, e2, z, noise=0.01, seed=e,
                                      supersample=4)
            fname = os.path.join(tmp, 'board_%02d.tif'%e)
            _save(img, fname)
            files.append(fname)
            truths.append(truth)

        # a sloppy click on the board origin, off by several pixels
        hint = truths[0][0, 0] + array([8., -7.])

        ch = checkerboard_cal_points(files, board_size=board_size,
                                     square_size=SQUARE, origin=origin,
                                     i_axis=e1, j_axis=e2,
                                     translations=zs, origin_hint=hint,
                                     sigma=2.0, min_distance=8,
                                     print_progress=False)
        rows = ch.process()

        n_expected = len(zs)*board_size[0]*board_size[1]
        assert len(rows) == n_expected, \
            'got %d points, expected %d'%(len(rows), n_expected)

        fname = os.path.join(tmp, 'cam_cal_points')
        ch.save(fname)

        cic = Cal_image_coord(fname)
        assert cic.N_points == len(rows)

        # the camera that made the images must project the lab coordinates of
        # the generated points back onto their pixel coordinates. This checks
        # the points themselves, independently of any fitting.
        res = [norm(cam_true.projection(array(X)) - array(x))
               for X, x in zip(cic.lab_coords, cic.image_coords)]
        assert sqrt(mean(square(res))) < 0.15, \
            'the generated points do not match the true camera: %.4f px' \
            %sqrt(mean(square(res)))

        # now fit a fresh camera, starting far from the answer
        cam = camera_Tsai('cam_fit')
        cam.resolution = cam_true.resolution
        cam.O = array([10., -12., -280.])
        cam.theta = zeros(3)
        cam.f = cam_true.f*0.96
        cam.xh, cam.yh = 0.0, 0.0
        cam.E = zeros((3, 5))
        cam.calc_R()

        cal = calibrate_Tsai(cam, cic.lab_coords, cic.image_coords)
        cal.searchCalibration(maxiter=4000, fix_f=False)

        err = cal.mean_squared_err()
        assert err < 0.5, 'calibration error too large: %.4f px'%err

        assert norm(cam.O - cam_true.O) < 5.0, \
            'camera origin off by %.3f'%norm(cam.O - cam_true.O)
        assert norm(cam.theta - cam_true.theta) < 0.01, \
            'camera angles off by %.5f'%norm(cam.theta - cam_true.theta)


def test_missing_images_raise():
    '''
    A pattern that matches nothing, or a path that does not exist, should be
    reported at once rather than after a long detection run.
    '''
    raised = False
    try:
        checkerboard_cal_points('./no_such_directory/*.tif', (9, 6), SQUARE)
    except ValueError:
        raised = True
    assert raised

    raised = False
    try:
        checkerboard_cal_points(['./no_such_file.tif'], (9, 6), SQUARE)
    except ValueError:
        raised = True
    assert raised


def test_translations_must_match_images():
    '''
    Giving a number of translations different from the number of images is a
    mistake that would attach the wrong depth to a whole image, so it must
    raise rather than pair them up silently.
    '''
    board_size = (9, 6)
    cam, origin, e1, e2 = board_view(board_size)

    with TemporaryDirectory() as tmp:

        files = []
        for e in range(2):
            img, _ = render_board(cam, board_size, SQUARE, origin, e1, e2,
                                  supersample=1)
            fname = os.path.join(tmp, 'b_%d.tif'%e)
            _save(img, fname)
            files.append(fname)

        raised = False
        try:
            checkerboard_cal_points(files, board_size, SQUARE,
                                    translations=[0.0, 1.0, 2.0])
        except ValueError:
            raised = True
        assert raised, 'a wrong number of translations did not raise'
