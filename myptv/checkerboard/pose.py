# -*- coding: utf-8 -*-
#!/usr/bin/env python3

"""
Created on Aug 27 2026

@author: Stefano Brizzolara


Recovering the pose of a checkerboard that was moved freely.

The procedure in cal_points.py takes the pose of the board from the
experiment: the board sits on a translation stage, and the stage readings say
where it is. That is the tidiest way to calibrate, but it is not the only way
people photograph a target. It is at least as common to hold the board in the
volume and turn it about between exposures, in which case the pose of the
board in each image is not known and has to be recovered from the images
themselves before any lab coordinates can be written down. This module does
that recovery. What it needs in exchange is that the cameras are triggered
together, so that a given frame shows one board pose seen from several
directions.

The route is the classical one of Zhang (2000), extended over several
cameras.


1) A homography per image.

   The board is flat, so for a pin-hole camera the map from board
   coordinates to pixels is a homography H, which is fitted to the detected
   corners. Lens distortion is not a homography, and shows up as a residual;
   it is left for the camera model to absorb later.

2) The intrinsics of each camera.

   Writing H = [h1 h2 h3] and H = s K [r1 r2 t], the fact that r1 and r2 are
   orthonormal gives two linear constraints per image on the entries of

       (1)    B = K^-T K^-1 ,

   namely h1' B h2 = 0 and h1' B h1 = h2' B h2. Three images of the board in
   different orientations therefore over-determine the six entries of the
   symmetric B, which is solved for in the least-squares sense and factored
   to give K. This is why the board has to be turned about between exposures:
   images that differ only by a translation give the same constraints over
   and over, and B stays singular.

3) The pose of the board in each camera.

   With K known, s = 1/|K^-1 h1| and

       (2)    r1 = s K^-1 h1,  r2 = s K^-1 h2,  t = s K^-1 h3,

   with r3 = r1 x r2. The r1 and r2 that come out of (2) are not exactly
   orthonormal, because H was fitted to noisy corners, so the nearest proper
   rotation is taken.

4) One lab frame for all the cameras.

   Step 3 gives the pose of the board in the frame of each camera
   separately. Since the cameras fire together, the board pose at frame f is
   one physical thing, and the transformation between two cameras follows:

       (3)    T(c1 -> c2) = T(c2,f) T(c1,f)^-1 ,

   which must come out the same for every f, because the cameras do not
   move. That it does is a strong check that the whole chain is right, and
   this module reports the spread over frames so it can be looked at. The
   lab frame is then fixed by declaring the board at one chosen frame to
   define it, which makes every camera's extrinsics and every board pose
   expressible in one frame.

   Note what this does and does not give. The geometry of the cameras
   relative to each other and to the board is recovered, and the scale comes
   from the physical size of the squares, so lengths are right. The
   orientation and origin of the lab frame, however, are whatever the board
   happened to be doing in the reference frame. If lab coordinates have to
   line up with something physical, a tank wall or a traverse, then either
   the board should be placed deliberately for the reference exposure, or the
   frame has to be related to the apparatus afterwards.

5) The remaining index ambiguity.

   The lattice out of grid.py is indexed arbitrarily, and for a board that
   is turned about there is no fixed pixel position to hint at. The
   ambiguity is instead resolved from the geometry, in three stages. The
   shape of the board rules out the transposition unless it is square. Of
   the four flips that remain, two describe the board seen from behind; they
   are recognisable because the pose that comes out of (2) puts the board
   behind the camera or gives an improper rotation, and are discarded. The
   last two differ by a half turn about the board normal, and are separated
   by requiring the pose to vary smoothly along the sequence, a half turn
   between neighbouring frames being far larger than any real motion.


Reference:
Zhang, Z. (2000). A flexible new technique for camera calibration. IEEE
Transactions on Pattern Analysis and Machine Intelligence, 22(11), 1330-1334.

"""

from numpy import (array, asarray, zeros, ones, eye, hstack, vstack, sqrt,
                   cross, dot, arcsin, arctan2, float64, isfinite, trace)
from numpy.linalg import svd, inv, norm, det




def fit_homography(board_pts, img_pts):
    '''
    Fits the homography taking planar board coordinates to pixels, in the
    least-squares sense, using the normalized direct linear transformation.

    The normalization matters: the entries of the design matrix otherwise mix
    quantities of order one with quantities of order 10^6, and the smallest
    singular value is then set by rounding rather than by the data.

    input -
    board_pts (array, (N,2)) - coordinates of the corners on the board, in
                               physical units
    img_pts (array, (N,2)) - the corresponding pixel coordinates

    output -
    H (array, (3,3)) - the homography, normalized so that H[2,2] = 1
    '''
    X = asarray(board_pts, dtype=float64)
    U = asarray(img_pts, dtype=float64)

    if X.shape[0] != U.shape[0] or X.shape[0] < 4:
        raise ValueError('at least four corresponding points are needed, got '
                         '%d and %d'%(X.shape[0], U.shape[0]))

    def normalizer(P):
        '''the similarity taking P to zero mean and mean distance sqrt(2)'''
        m = P.mean(axis=0)
        d = sqrt(((P - m)**2).sum(axis=1)).mean()
        s = sqrt(2.0)/d if d > 0 else 1.0
        T = array([[s, 0, -s*m[0]],
                   [0, s, -s*m[1]],
                   [0, 0, 1.0]])
        return T

    Tx, Tu = normalizer(X), normalizer(U)

    Xn = (hstack([X, ones((len(X), 1))]).dot(Tx.T))[:, :2]
    Un = (hstack([U, ones((len(U), 1))]).dot(Tu.T))[:, :2]

    N = len(Xn)
    A = zeros((2*N, 9))
    A[0::2, 0] = Xn[:, 0]; A[0::2, 1] = Xn[:, 1]; A[0::2, 2] = 1.0
    A[0::2, 6] = -Un[:, 0]*Xn[:, 0]
    A[0::2, 7] = -Un[:, 0]*Xn[:, 1]
    A[0::2, 8] = -Un[:, 0]
    A[1::2, 3] = Xn[:, 0]; A[1::2, 4] = Xn[:, 1]; A[1::2, 5] = 1.0
    A[1::2, 6] = -Un[:, 1]*Xn[:, 0]
    A[1::2, 7] = -Un[:, 1]*Xn[:, 1]
    A[1::2, 8] = -Un[:, 1]

    _, _, Vt = svd(A)
    Hn = Vt[-1].reshape(3, 3)

    H = inv(Tu).dot(Hn).dot(Tx)

    if H[2, 2] == 0:
        raise ValueError('the homography came out singular')

    return H/H[2, 2]




def homography_residual(H, board_pts, img_pts):
    '''
    The distances, in pixels, between the corners and where the homography
    puts them. This is the quantity that certifies a set of detected corners
    without any knowledge of where the board was: a plane seen by a pin-hole
    camera maps to the image by an exact homography, so what is left over is
    the error in locating the corners plus the lens distortion.

    input -
    H (array, (3,3)) - the homography
    board_pts (array, (N,2)) - board coordinates
    img_pts (array, (N,2)) - pixel coordinates

    output -
    r (array, N) - the distance for each point
    '''
    X = hstack([asarray(board_pts, dtype=float64), ones((len(board_pts), 1))])
    P = X.dot(H.T)
    P = P[:, :2]/P[:, 2:3]
    return sqrt(((P - asarray(img_pts, dtype=float64))**2).sum(axis=1))




def _v_ij(H, i, j):
    '''
    The row vector v such that v . b = h_i' B h_j, where b holds the six
    distinct entries of the symmetric B of eq. (1) in the order
    B11, B12, B22, B13, B23, B33.
    '''
    hi, hj = H[:, i], H[:, j]
    return array([hi[0]*hj[0],
                  hi[0]*hj[1] + hi[1]*hj[0],
                  hi[1]*hj[1],
                  hi[0]*hj[2] + hi[2]*hj[0],
                  hi[1]*hj[2] + hi[2]*hj[1],
                  hi[2]*hj[2]])




def zhang_intrinsics(homographies):
    '''
    Recovers the intrinsic matrix of a camera from homographies of a planar
    target in several orientations, as in step 2 of the module docstring.

    input -
    homographies (list of arrays, (3,3)) - at least three homographies, from
                board coordinates to pixels, of the board in DIFFERENT
                orientations. Orientations that differ only by a translation
                carry no new information and will leave the problem singular.

    output -
    K (array, (3,3)) - the intrinsic matrix, [[a, g, u0], [0, b, v0],
                       [0, 0, 1]]
    cond (float) - the ratio of the second smallest to the smallest singular
                   value of the constraint matrix. It measures how well the
                   orientations of the board determine the answer; a value
                   close to one means they do not, and the result should not
                   be trusted. Well spread orientations give tens or more.
    '''
    Hs = [asarray(H, dtype=float64) for H in homographies]

    if len(Hs) < 3:
        raise ValueError('at least three homographies are needed to recover '
                         'the intrinsics, got %d'%len(Hs))

    rows = []
    for H in Hs:
        rows.append(_v_ij(H, 0, 1))
        rows.append(_v_ij(H, 0, 0) - _v_ij(H, 1, 1))

    V = vstack(rows)

    _, s, Vt = svd(V)
    b = Vt[-1]

    cond = float(s[-2]/s[-1]) if s[-1] > 0 else float('inf')

    B11, B12, B22, B13, B23, B33 = b

    # factor B = K^-T K^-1; see Zhang (2000), appendix B
    d = B11*B22 - B12**2
    if d == 0 or B11 == 0:
        raise ValueError('the constraints on the intrinsics came out '
                         'degenerate; the board was probably not turned '
                         'about enough between exposures')

    v0 = (B12*B13 - B11*B23)/d
    lam = B33 - (B13**2 + v0*(B12*B13 - B11*B23))/B11

    if lam/B11 <= 0 or lam*B11/d <= 0:
        raise ValueError('the recovered intrinsics are not physical; the '
                         'board was probably not turned about enough between '
                         'exposures, or the corner indexing is inconsistent')

    alpha = sqrt(lam/B11)
    beta = sqrt(lam*B11/d)
    gamma = -B12*alpha**2*beta/lam
    u0 = gamma*v0/beta - B13*alpha**2/lam

    K = array([[alpha, gamma, u0],
               [0.0, beta, v0],
               [0.0, 0.0, 1.0]])

    return K, cond




def focal_from_homographies(homographies, principal_point):
    '''
    Recovers only the focal length, with the principal point given and the
    pixels assumed square and unskewed, from the same constraints that
    zhang_intrinsics uses.

    This exists because the full five parameter estimate is often not usable
    in practice. The two constraints per image determine the focal length far
    better than they determine the principal point, and when the field of
    view is narrow and the target stays in one region of the frame the two are
    so strongly correlated that the unconstrained solution puts the principal
    point hundreds of pixels outside anything physical and absorbs the error
    into a large skew. On a real three camera set the unconstrained estimate
    gave principal points more than 600 pixels off centre and skews of up to
    32, while the focal lengths were sensible. Fixing what is badly determined
    and solving for what is well determined gives a much better starting point
    for a subsequent refinement, which can then free the principal point again
    against the full reprojection error rather than against these two
    algebraic constraints.

    With the principal point (cx, cy) fixed and w = 1/f^2, the matrix of
    eq. (1) is

        B = [[w, 0, -cx w], [0, w, -cy w],
             [-cx w, -cy w, (cx^2 + cy^2) w + 1]]

    so each of the two constraints is linear in the single unknown w, and all
    the images are combined by least squares.

    input -
    homographies (list of arrays, (3,3)) - homographies from board
                coordinates to pixels, of the board in different orientations
    principal_point (tuple, 2) - the pixel coordinates to fix the principal
                point at, normally the centre of the image

    output -
    K (array, (3,3)) - the intrinsic matrix
    f (float) - the focal length, in pixels
    '''
    cx, cy = float(principal_point[0]), float(principal_point[1])

    rows, rhs = [], []

    for H in homographies:
        H = asarray(H, dtype=float64)

        def pq(i, j):
            hi, hj = H[:, i], H[:, j]
            p = (hi[0]*hj[0] + hi[1]*hj[1]
                 - cx*(hi[0]*hj[2] + hi[2]*hj[0])
                 - cy*(hi[1]*hj[2] + hi[2]*hj[1])
                 + (cx**2 + cy**2)*hi[2]*hj[2])
            q = hi[2]*hj[2]
            return p, q

        p12, q12 = pq(0, 1)
        p11, q11 = pq(0, 0)
        p22, q22 = pq(1, 1)

        rows.append(p12);         rhs.append(-q12)
        rows.append(p11 - p22);   rhs.append(-(q11 - q22))

    A = array(rows, dtype=float64).reshape(-1, 1)
    b = array(rhs, dtype=float64)

    denom = float((A[:, 0]**2).sum())
    if denom == 0:
        raise ValueError('the constraints on the focal length came out '
                         'degenerate')

    w = float((A[:, 0]*b).sum()/denom)

    if w <= 0:
        raise ValueError('the recovered focal length is not physical (w=%g); '
                         'the board was probably not turned about enough '
                         'between exposures, or the corner indexing is '
                         'inconsistent'%w)

    f = 1.0/sqrt(w)

    K = array([[f, 0.0, cx],
               [0.0, f, cy],
               [0.0, 0.0, 1.0]])

    return K, f




def nearest_rotation(M):
    '''
    The proper rotation closest to M in the Frobenius sense, obtained from
    the singular value decomposition. Needed because the two columns that
    come out of eq. (2) are only approximately orthonormal.

    input -
    M (array, (3,3))

    output -
    R (array, (3,3)) - a matrix with R' R = I and det(R) = +1
    '''
    U, _, Vt = svd(asarray(M, dtype=float64))
    R = U.dot(Vt)

    if det(R) < 0:
        U[:, -1] *= -1.0
        R = U.dot(Vt)

    return R




def pose_from_homography(H, K):
    '''
    The pose of the board in the frame of the camera, from the homography and
    the intrinsics, as in step 3 of the module docstring.

    The sign of the scale is fixed by requiring the board to be in front of
    the camera. If it is not, the homography describes the board seen from
    behind, which happens when the corner indexing has been flipped, and the
    returned pose is marked as such.

    input -
    H (array, (3,3)) - homography from board coordinates to pixels
    K (array, (3,3)) - the intrinsic matrix

    output -
    R (array, (3,3)) - rotation taking board coordinates to camera coordinates
    t (array, 3) - the origin of the board in camera coordinates
    ok (bool) - False if the board came out behind the camera, in which case
                the pose is meaningless and the indexing should be flipped
    '''
    Kinv = inv(asarray(K, dtype=float64))
    H = asarray(H, dtype=float64)

    h1, h2, h3 = H[:, 0], H[:, 1], H[:, 2]

    a1, a2, a3 = Kinv.dot(h1), Kinv.dot(h2), Kinv.dot(h3)

    n = norm(a1)
    if n == 0 or not isfinite(n):
        return eye(3), zeros(3), False

    s = 1.0/n

    r1, r2, t = s*a1, s*a2, s*a3

    # the board must be in front of the camera
    ok = bool(t[2] > 0)
    if not ok:
        r1, r2, t = -r1, -r2, -t

    r3 = cross(r1, r2)

    R = nearest_rotation(vstack([r1, r2, r3]).T)

    return R, t, ok




def relative_pose(R1, t1, R2, t2):
    '''
    The rigid transformation taking the frame of camera 1 to that of camera 2,
    given the pose of the same board in both, as in eq. (3).

    input -
    R1, t1 - pose of the board in camera 1
    R2, t2 - pose of the board in camera 2

    output -
    R (array, (3,3)), t (array, 3) - the transformation, x2 = R x1 + t
    '''
    R = dot(R2, R1.T)
    t = asarray(t2, dtype=float64) - dot(R, asarray(t1, dtype=float64))
    return R, t




def rotation_angle(R):
    '''
    The angle, in radians, of the rotation R. Used to compare two rotations
    by the angle of R1' R2, which is the natural distance between them.

    input -
    R (array, (3,3))

    output -
    angle (float) - in [0, pi]
    '''
    c = (trace(asarray(R, dtype=float64)) - 1.0)/2.0
    c = max(-1.0, min(1.0, float(c)))
    from math import acos
    return acos(c)




def theta_from_R(R):
    '''
    The three angles of the Tsai camera model that correspond to a rotation
    matrix, that is, the angles for which

        R = Rx(theta1) Ry(theta2) Rz(theta3)

    with the convention used by myptv.TsaiModel.camera.camera_Tsai.calc_R.
    Writing out that product gives

        R[0,2] = sin(theta2)
        R[0,0] = cos(theta2) cos(theta3),  R[0,1] = -cos(theta2) sin(theta3)
        R[1,2] = -sin(theta1) cos(theta2), R[2,2] = cos(theta1) cos(theta2)

    from which the angles follow. The decomposition is not unique when
    cos(theta2) vanishes, that is when the second angle is a quarter turn;
    the branch with theta3 = 0 is returned there.

    input -
    R (array, (3,3))

    output -
    theta (array, 3)
    '''
    R = asarray(R, dtype=float64)

    s2 = max(-1.0, min(1.0, float(R[0, 2])))
    theta2 = arcsin(s2)

    c2 = sqrt(max(0.0, 1.0 - s2**2))

    if c2 < 1e-8:
        # degenerate: only the sum or difference of the other two angles is
        # determined, so one of them is set to zero
        theta3 = 0.0
        theta1 = arctan2(R[1, 0], R[1, 1])
    else:
        theta3 = arctan2(-R[0, 1], R[0, 0])
        theta1 = arctan2(-R[1, 2], R[2, 2])

    return array([float(theta1), float(theta2), float(theta3)])




def R_from_theta(theta):
    '''
    The rotation matrix of the three Tsai angles, the inverse of
    theta_from_R. Kept here so that the two can be checked against each
    other.

    input -
    theta (array, 3)

    output -
    R (array, (3,3))
    '''
    from math import sin, cos
    t1, t2, t3 = [float(t) for t in theta]

    Rx = array([[1, 0, 0],
                [0, cos(t1), -sin(t1)],
                [0, sin(t1), cos(t1)]])
    Ry = array([[cos(t2), 0, sin(t2)],
                [0, 1, 0],
                [-sin(t2), 0, cos(t2)]])
    Rz = array([[cos(t3), -sin(t3), 0],
                [sin(t3), cos(t3), 0],
                [0, 0, 1]])

    return dot(dot(Rx, Ry), Rz)
