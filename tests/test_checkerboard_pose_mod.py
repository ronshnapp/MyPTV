# -*- coding: utf-8 -*-
"""
Created on Aug 27 2026

@author: Stefano Brizzolara


Tests for the recovery of a freely moved board: the geometry of pose.py and
the joint adjustment of bundle.py.

Everything here is checked against a synthetic arrangement whose cameras,
board poses and corner positions are all known exactly, so that what is
measured is the method and not a plausible looking number.

A word on what is checked at the end. A calibration of this kind determines
the geometry only up to a rigid motion: nothing in the images says where the
origin of the lab frame is or which way its axes point, and the procedure
simply declares one camera to be it. Comparing the recovered lab coordinates
with the true ones therefore measures the arbitrary choice of frame as much
as anything else, and on a real arrangement, where the cameras share roughly
one 'up' direction, shifting every principal point together looks so much
like translating every board that the two cannot be separated. What is
meaningful, and what is checked, is the SHAPE of the reconstruction: the
error left after the best rigid motion has been taken out, and the scale,
which the physical size of the squares fixes.

"""

from numpy import (array, asarray, zeros, eye, concatenate, cross, degrees,
                   median, diag, float64)
from numpy.linalg import norm, svd, det
from numpy.random import default_rng

from myptv.checkerboard.pose import (fit_homography, homography_residual,
                                     zhang_intrinsics,
                                     focal_from_homographies,
                                     pose_from_homography, relative_pose,
                                     rotation_angle, nearest_rotation,
                                     theta_from_R, R_from_theta)
from myptv.checkerboard.bundle import (skew, rodrigues, rvec_from_R,
                                       rodrigues_jacobian, project,
                                       multicamera_bundle)


# the arrangement used by most of the tests: a modest board seen by four
# cameras from well separated directions
BOARD = (6, 8)
SQUARE = 15.0
W, H = 2000, 1400
CAMS = ['c0', 'c1', 'c2', 'c3']


def board_points(board=BOARD, square=SQUARE):
    '''
    The corners of the board, in the board's own frame, centred on it.
    '''
    B = array([[i*square, j*square, 0.0] for i in range(board[0])
               for j in range(board[1])])
    return B - B.mean(axis=0)


def look_at(O, target=None, up=None):
    '''
    A camera at O looking at target; returns (R, t) taking lab to camera
    coordinates.
    '''
    O = asarray(O, dtype=float64)
    target = zeros(3) if target is None else asarray(target, dtype=float64)
    up = array([0.0, 0.0, 1.0]) if up is None else asarray(up, dtype=float64)

    z = target - O; z = z/norm(z)
    x = cross(up, z); x = x/norm(x)
    y = cross(z, x)
    R = array([x, y, z])
    return R, -R.dot(O)


def synthetic_setup(n_frames=14, noise=0.1, seed=0, n_dist=2, f=5000.0,
                    spread=90.0):
    '''
    Builds a four camera arrangement, a set of board poses, and the corner
    positions they produce, with noise.

    output - a dictionary holding the truth and the observations
    '''
    rng = default_rng(seed)
    B = board_points()

    intr = {'c0': array([f,       W/2 - 18, H/2 + 12, -0.04, 0.25]),
            'c1': array([f - 60,  W/2 + 14, H/2 - 9,  -0.03, 0.18]),
            'c2': array([f + 70,  W/2 - 6,  H/2 + 4,  -0.05, 0.30]),
            'c3': array([f - 25,  W/2 + 22, H/2 + 17, -0.02, 0.12])}
    intr = {c: intr[c][:3 + n_dist] for c in CAMS}

    Os = {'c0': array([0.0, -900.0, 120.0]),
          'c1': array([-260.0, -880.0, 120.0]),
          'c2': array([700.0, -640.0, 190.0]),
          'c3': array([-700.0, -660.0, 170.0])}
    ext = {c: look_at(Os[c]) for c in CAMS}

    # The board is held facing the cameras, which sit towards -y. Its own
    # third axis, e1 x e2, therefore points along +y, away from them. This
    # matters for more than realism: the test that separates the index flips
    # asks which way that axis points, so a board left lying in the x-y
    # plane, edge-on to cameras that look horizontally, would make the
    # question meaningless.
    R0 = rodrigues(array([-1.5707963267948966, 0.0, 0.0]))

    poses, obs = {}, {}
    for k in range(n_frames):
        ax = rng.normal(0, 1, 3); ax = ax/norm(ax)
        Rf = rodrigues(ax*rng.uniform(0.05, 0.6)).dot(R0)
        tf = rng.uniform(-spread, spread, 3)
        poses[k] = (Rf, tf)
        Xl = B.dot(Rf.T) + tf
        for c in CAMS:
            Rc, tc = ext[c]
            uv = project(Xl.dot(Rc.T) + tc, intr[c], n_dist)
            if (uv[:, 0].min() < 10 or uv[:, 0].max() > W - 10 or
                    uv[:, 1].min() < 10 or uv[:, 1].max() > H - 10):
                continue
            obs[(c, k)] = uv + rng.normal(0, noise, uv.shape)

    return dict(B=B, intr=intr, ext=ext, Os=Os, poses=poses, obs=obs,
                n_dist=n_dist, noise=noise)


def align_rigid(P, Q):
    '''
    The rotation and translation taking P as close as possible to Q, and the
    scale that would do best if it were allowed. Used to separate the
    arbitrary choice of lab frame from a genuine distortion of the shape.
    '''
    P, Q = asarray(P, dtype=float64), asarray(Q, dtype=float64)
    mp, mq = P.mean(axis=0), Q.mean(axis=0)
    A = (P - mp).T.dot(Q - mq)/len(P)
    U, S, Vt = svd(A)
    D = eye(3)
    if det(U)*det(Vt) < 0:
        D[2, 2] = -1.0
    R = Vt.T.dot(D).dot(U.T)
    s = float((S*diag(D)).sum()*len(P)/((P - mp)**2).sum())
    t = mq - R.dot(mp)
    return R, t, s




# =============================================================================
#   rotations
# =============================================================================


def test_rodrigues_round_trip():
    '''
    Turning a rotation vector into a matrix and back must return it.
    '''
    rng = default_rng(1)
    worst = 0.0
    for _ in range(200):
        r = rng.normal(0, 1, 3)
        r = r/norm(r)*rng.uniform(0.0, 3.0)
        R = rodrigues(r)
        worst = max(worst, norm(rodrigues(rvec_from_R(R)) - R))
    assert worst < 1e-10, 'worst round trip error %.3e'%worst


def test_rodrigues_jacobian():
    '''
    The analytic derivative of a rotation with respect to its rotation
    vector must agree with finite differences, including near the identity
    where the formula has a removable singularity.
    '''
    rng = default_rng(2)
    worst = 0.0
    for k in range(120):
        r = rng.normal(0, 1, 3)
        if k % 8 == 0:
            r = r*1e-10
        dR = rodrigues_jacobian(r)
        for i in range(3):
            h = 1e-7
            q = r.copy(); q[i] += h
            num = (rodrigues(q) - rodrigues(r))/h
            worst = max(worst, float(abs(num - dR[i]).max()))
    assert worst < 1e-5, 'worst derivative error %.3e'%worst


def test_skew():
    '''
    The skew matrix must reproduce the cross product.
    '''
    rng = default_rng(3)
    for _ in range(20):
        a, b = rng.normal(0, 1, 3), rng.normal(0, 1, 3)
        assert norm(skew(a).dot(b) - cross(a, b)) < 1e-12


def test_theta_conversion_matches_the_Tsai_model():
    '''
    The three angles of the Tsai model must be recovered from a rotation
    matrix and rebuild it, and must agree with the model's own calc_R.
    '''
    from myptv.TsaiModel.camera import camera_Tsai

    rng = default_rng(4)
    worst = 0.0
    for _ in range(200):
        th = rng.uniform(-1.4, 1.4, 3)
        R = R_from_theta(th)
        worst = max(worst, norm(R_from_theta(theta_from_R(R)) - R))
    assert worst < 1e-10, 'worst round trip %.3e'%worst

    cam = camera_Tsai('t')
    cam.theta = array([0.31, -0.22, 0.47])
    cam.calc_R()
    assert norm(R_from_theta(cam.theta) - cam.R) < 1e-12




# =============================================================================
#   the geometry of a single view
# =============================================================================


def test_homography_is_exact_for_a_plane():
    '''
    A plane seen by a pin-hole camera maps to the image by a homography
    exactly, so with no distortion and no noise the residual must vanish.
    '''
    s = synthetic_setup(n_frames=6, noise=0.0, n_dist=0, seed=5)
    B = s['B']
    for (c, k), uv in list(s['obs'].items())[:6]:
        H = fit_homography(B[:, :2], uv)
        r = homography_residual(H, B[:, :2], uv)
        assert r.max() < 1e-6, 'residual %.3e on a noiseless plane'%r.max()


def test_zhang_intrinsics_recovers_the_camera():
    '''
    The intrinsics must come back from homographies of the board in several
    orientations.
    '''
    s = synthetic_setup(n_frames=25, noise=0.0, n_dist=0, seed=6)
    B = s['B']
    for c in CAMS:
        Hs = [fit_homography(B[:, :2], uv) for (cc, k), uv in s['obs'].items()
              if cc == c]
        if len(Hs) < 6:
            continue
        K, cond = zhang_intrinsics(Hs)
        ftrue, cxtrue, cytrue = s['intr'][c][:3]
        assert abs(K[0, 0] - ftrue)/ftrue < 0.01, \
            '%s: alpha %.1f vs %.1f'%(c, K[0, 0], ftrue)
        assert abs(K[1, 1] - ftrue)/ftrue < 0.01, \
            '%s: beta %.1f vs %.1f'%(c, K[1, 1], ftrue)
        assert abs(K[0, 2] - cxtrue) < 0.02*W, \
            '%s: u0 %.1f vs %.1f'%(c, K[0, 2], cxtrue)
        assert abs(K[1, 2] - cytrue) < 0.02*H, \
            '%s: v0 %.1f vs %.1f'%(c, K[1, 2], cytrue)


def test_zhang_needs_several_orientations():
    '''
    Two views cannot determine five intrinsic parameters, and the function
    must say so rather than return something arbitrary.
    '''
    s = synthetic_setup(n_frames=4, noise=0.0, n_dist=0, seed=7)
    B = s['B']
    Hs = [fit_homography(B[:, :2], uv)
          for (cc, k), uv in list(s['obs'].items())[:2]]
    raised = False
    try:
        zhang_intrinsics(Hs)
    except ValueError:
        raised = True
    assert raised, 'two homographies did not raise'


def test_focal_from_homographies():
    '''
    With the principal point given, the focal length alone must come back.
    '''
    s = synthetic_setup(n_frames=25, noise=0.0, n_dist=0, seed=8)
    B = s['B']
    for c in CAMS:
        Hs = [fit_homography(B[:, :2], uv) for (cc, k), uv in s['obs'].items()
              if cc == c]
        if len(Hs) < 6:
            continue
        ftrue, cx, cy = s['intr'][c][:3]
        K, f = focal_from_homographies(Hs, (cx, cy))
        assert abs(f - ftrue)/ftrue < 0.005, \
            '%s: f %.1f vs %.1f'%(c, f, ftrue)


def test_pose_from_homography():
    '''
    Given the intrinsics, the pose of the board must come back from its
    homography.
    '''
    s = synthetic_setup(n_frames=12, noise=0.0, n_dist=0, seed=9)
    B = s['B']
    worst_ang, worst_t = 0.0, 0.0
    for (c, k), uv in s['obs'].items():
        ftrue, cx, cy = s['intr'][c][:3]
        K = array([[ftrue, 0, cx], [0, ftrue, cy], [0, 0, 1.0]])
        H = fit_homography(B[:, :2], uv)
        R, t, ok = pose_from_homography(H, K)

        Rc, tc = s['ext'][c]
        Rf, tf = s['poses'][k]
        R_true = Rc.dot(Rf)
        t_true = Rc.dot(tf) + tc

        worst_ang = max(worst_ang, degrees(rotation_angle(R.T.dot(R_true))))
        worst_t = max(worst_t, norm(t - t_true))

    assert worst_ang < 0.5, 'worst rotation error %.3f deg'%worst_ang
    assert worst_t < 1.0, 'worst translation error %.3f mm'%worst_t


def test_handedness_separates_the_index_flips():
    '''
    The lattice of grid.py is indexed arbitrarily, and for a board that is
    turned about there is no pixel position to hint at. Two of the four
    flips describe the board seen from behind and must be rejected by
    requiring its front face to point at the camera, which is the condition
    dot(r1 x r2, t) > 0. Exactly two candidates must survive, and the correct
    one must be among them.
    '''
    s = synthetic_setup(n_frames=10, noise=0.0, n_dist=0, seed=10)
    B = s['B']

    checked = 0
    for (c, k), uv in s['obs'].items():
        ftrue, cx, cy = s['intr'][c][:3]
        K = array([[ftrue, 0, cx], [0, ftrue, cy], [0, 0, 1.0]])

        Rc, tc = s['ext'][c]
        Rf, tf = s['poses'][k]
        R_true = Rc.dot(Rf)

        P = uv.reshape(BOARD[0], BOARD[1], 2)
        kept, correct = 0, 0
        for flip in range(4):
            Q = P[::-1] if flip & 1 else P
            Q = Q[:, ::-1] if flip & 2 else Q
            H = fit_homography(B[:, :2], Q.reshape(-1, 2))
            R, t, _ = pose_from_homography(H, K)
            if R[:, 2].dot(t) > 0:
                kept += 1
                if degrees(rotation_angle(R.T.dot(R_true))) < 1.0:
                    correct += 1

        assert kept == 2, '%d flips survived, expected 2'%kept
        assert correct == 1, 'the correct flip was not among those kept'
        checked += 1

    assert checked > 20, 'only %d views were checked'%checked


def test_relative_pose_is_the_same_at_every_frame():
    '''
    The cameras do not move, so the transformation between two of them,
    worked out from the board seen in both, must come out the same whichever
    frame it is taken from.
    '''
    s = synthetic_setup(n_frames=12, noise=0.0, n_dist=0, seed=11)
    B = s['B']

    def pose(c, k):
        ftrue, cx, cy = s['intr'][c][:3]
        K = array([[ftrue, 0, cx], [0, ftrue, cy], [0, 0, 1.0]])
        return pose_from_homography(fit_homography(B[:, :2], s['obs'][(c, k)]),
                                    K)[:2]

    for c in CAMS[1:]:
        Rs = []
        for k in s['poses']:
            if ('c0', k) not in s['obs'] or (c, k) not in s['obs']:
                continue
            R0, t0 = pose('c0', k)
            R1, t1 = pose(c, k)
            Rs.append(relative_pose(R0, t0, R1, t1)[0])
        if len(Rs) < 4:
            continue
        Rm = nearest_rotation(sum(Rs)/len(Rs))
        angs = [degrees(rotation_angle(Rm.T.dot(R))) for R in Rs]
        assert max(angs) < 0.5, \
            'c0 -> %s varies by up to %.3f deg between frames'%(c, max(angs))




# =============================================================================
#   the joint adjustment
# =============================================================================


def _build(s, seed=99, n_dist=None):
    '''
    A bundle starting from a deliberately perturbed version of the truth.
    '''
    rng = default_rng(seed)
    n_dist = s['n_dist'] if n_dist is None else n_dist

    intr0, ext0, poses0 = {}, {}, {}
    for c in CAMS:
        intr0[c] = concatenate([[s['intr'][c][0]*1.01, W/2, H/2],
                                zeros(n_dist)])
        R, t = s['ext'][c]
        ext0[c] = (rodrigues(rng.normal(0, 0.01, 3)).dot(R),
                   t + rng.normal(0, 8.0, 3))
    ext0['c0'] = s['ext']['c0']
    for k, (R, t) in s['poses'].items():
        poses0[k] = (rodrigues(rng.normal(0, 0.02, 3)).dot(R),
                     t + rng.normal(0, 6.0, 3))

    return multicamera_bundle(s['B'], s['obs'], CAMS, intr0, ext0, poses0,
                              n_dist=n_dist, fix_camera='c0')


def test_analytic_jacobian_matches_finite_differences():
    '''
    The Jacobian is written out by hand, so it is checked against finite
    differences. A wrong Jacobian does not raise; it quietly makes the solve
    fail or converge somewhere wrong, which is exactly the sort of thing a
    test has to catch.
    '''
    s = synthetic_setup(n_frames=8, noise=0.05, seed=12)
    ba = _build(s)
    err = ba.get_jacobian_error(n_check=45, seed=1)
    assert err < 1e-4, 'largest relative discrepancy %.3e'%err


def test_fixed_camera_does_not_move():
    '''
    The pose of the camera that defines the lab frame is the gauge and must
    come back untouched, otherwise the frame has drifted and the result means
    something different from what it claims.
    '''
    s = synthetic_setup(n_frames=8, noise=0.1, seed=13)
    ba = _build(s)
    before = ba.x.copy()
    ba.solve(max_nfev=40)
    ci = CAMS.index('c0')
    moved = abs(ba.x[ci*ba.CAMP:ci*ba.CAMP + 6]
                - before[ci*ba.CAMP:ci*ba.CAMP + 6]).max()
    assert moved == 0.0, 'the fixed camera moved by %.3e'%moved


def test_bundle_recovers_the_geometry():
    '''
    The whole point: starting from a perturbed estimate, the adjustment must
    reproduce the arrangement that made the images.

    The shape of the reconstruction is what is checked, for the reason set
    out at the top of this file, together with the scale and the distances
    between the cameras, both of which are independent of the lab frame.
    '''
    s = synthetic_setup(n_frames=16, noise=0.1, seed=14)
    ba = _build(s)

    assert ba.rms() > 5.0, 'the starting point was supposed to be perturbed'

    ba.solve(max_nfev=400)

    # the fit itself should reach the noise
    assert ba.rms() < 2.0*s['noise'], \
        'reprojection %.4f px against %.2f px of noise'%(ba.rms(), s['noise'])

    intr, ext, poses = ba.unpack()

    # the focal lengths
    for c in CAMS:
        rel = abs(intr[c][0] - s['intr'][c][0])/s['intr'][c][0]
        assert rel < 0.02, '%s: focal length off by %.3f%%'%(c, 100*rel)

    # the distances between the cameras, which no choice of frame can change
    Or = {c: -ext[c][0].T.dot(ext[c][1]) for c in CAMS}
    for i in range(len(CAMS)):
        for j in range(i+1, len(CAMS)):
            a, b = CAMS[i], CAMS[j]
            got = norm(Or[a] - Or[b])
            want = norm(s['Os'][a] - s['Os'][b])
            assert abs(got - want) < 0.02*want, \
                '%s-%s: %.1f mm instead of %.1f mm'%(a, b, got, want)

    # and the shape of the reconstructed corners
    ks = sorted(poses)
    Xr = concatenate([s['B'].dot(poses[k][0].T) + poses[k][1] for k in ks])
    Xt = concatenate([s['B'].dot(s['poses'][k][0].T) + s['poses'][k][1]
                      for k in ks])
    R, t, sc = align_rigid(Xr, Xt)
    resid = norm(Xr.dot(R.T) + t - Xt, axis=1)

    assert median(resid) < 0.5, \
        'the shape is distorted by %.4f mm'%median(resid)
    assert abs(sc - 1.0) < 5e-3, 'the scale is off by %.2e'%(sc - 1.0)


def test_lab_coordinates_are_a_rigid_board():
    '''
    Whatever else it gets wrong, the adjustment cannot deform the board: the
    corners come from one rigid object, so their spacing must be exactly what
    was put in.
    '''
    s = synthetic_setup(n_frames=8, noise=0.1, seed=15)
    ba = _build(s)
    ba.solve(max_nfev=60)

    lab = ba.lab_coordinates()
    for k, X in lab.items():
        G = X.reshape(BOARD[0], BOARD[1], 3)
        d1 = norm(G[1:, :] - G[:-1, :], axis=2)
        d2 = norm(G[:, 1:] - G[:, :-1], axis=2)
        assert abs(d1 - SQUARE).max() < 1e-9
        assert abs(d2 - SQUARE).max() < 1e-9


def test_bundle_tolerates_a_missing_camera():
    '''
    A camera that did not see the board in some frames must simply
    contribute nothing there, not break the bookkeeping.
    '''
    s = synthetic_setup(n_frames=12, noise=0.1, seed=16)
    dropped = [k for k in list(s['poses'])[:4]]
    for k in dropped:
        s['obs'].pop(('c2', k), None)
        s['obs'].pop(('c3', k), None)

    ba = _build(s)
    n_before = len(ba.pairs)
    ba.solve(max_nfev=200)
    assert len(ba.pairs) == n_before
    assert ba.rms() < 2.0*s['noise'], \
        'reprojection %.4f px with cameras missing'%ba.rms()


def test_per_view_rms_finds_a_bad_view():
    '''
    A view whose corners were labelled wrongly should stand out clearly in
    the per view errors, which is how such a view is meant to be found.
    '''
    s = synthetic_setup(n_frames=12, noise=0.1, seed=17)
    victim = sorted(k for (c, k) in s['obs'] if c == 'c1')[3]
    P = s['obs'][('c1', victim)].reshape(BOARD[0], BOARD[1], 2)
    s['obs'][('c1', victim)] = P[::-1, ::-1].reshape(-1, 2)   # a half turn

    ba = _build(s)
    ba.solve(max_nfev=200, loss='soft_l1', f_scale=1.0)

    pv = {(c, k): v for c, k, v in ba.per_view_rms()}
    bad = pv[('c1', victim)]
    others = sorted(v for kk, v in pv.items() if kk != ('c1', victim))
    assert bad > 10*median(array(others)), \
        'the mislabelled view scored %.3f against a median of %.3f' \
        %(bad, median(array(others)))
