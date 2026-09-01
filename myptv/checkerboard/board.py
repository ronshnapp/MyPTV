# -*- coding: utf-8 -*-
#!/usr/bin/env python3

"""
Created on Aug 27 2026

@author: Stefano Brizzolara


The geometry of a checkerboard calibration target.

Once the corners found in an image have been assembled into a lattice (see
grid.py), each of them is identified by a pair of indices (a, b) that runs
over the lattice as it happens to be oriented in that particular image. To
turn those into lab space coordinates two further things are needed, and this
module supplies both.


1) Resolving the orientation of the lattice.

   The indices (a, b) coming out of the growth are arbitrary: the site that
   ended up being called (0,0) is whichever corner the growth was seeded
   from, and the direction called "a" is whichever edge direction the
   detector happened to return first. A lattice can therefore be laid over
   the physical board in eight different ways, obtained by combining a
   transposition with a flip of either index.

   This ambiguity is not harmless. In a multi camera experiment, mapping the
   lattice differently in two cameras assigns different lab coordinates to
   the same physical corner, and the resulting calibration is silently and
   completely wrong. It must therefore be resolved explicitly rather than by
   a convention that happens to depend on the viewing direction.

   The ambiguity is resolved with a hint: the approximate pixel position, in
   that camera, of the corner that the user has chosen to be the origin of
   the board. Of the eight mappings, the ones whose shape does not match the
   declared board size are discarded, and of those that remain the one whose
   origin falls closest to the hint is chosen. For a board with an unequal
   number of corners along its two sides this fixes the mapping completely.
   For a square board two mappings survive, differing in which side of the
   board is called the first axis, and a second hint pointing along that
   axis is needed.

   When no hint is given a deterministic fallback is used, and a warning is
   issued, because that fallback is only safe for a single camera.


2) Assigning lab coordinates.

   The board is flat, so the corner with indices (i, j) lies at

       (1)    X(i, j) = X0 + i * s * e1 + j * s * e2 + t

   where X0 is the lab position of the board origin, s is the physical
   spacing of the corners, e1 and e2 are unit vectors along the two sides of
   the board in lab space, and t is the translation of the board for the
   image in question.

   The intended use is the one common in 3D-PTV, in which the board is
   mounted on a translation stage and photographed by all cameras
   simultaneously at a number of known positions along the stage. Then e1,
   e2 and X0 are fixed by how the board is mounted, and t is a known
   displacement along the stage, which the caller may give either as a
   distance along the board normal or as a full three component vector. The
   several positions of the board together provide the spread in depth that
   a calibration needs, and because the cameras see the board at the same
   instant they share one lab frame by construction.

   Note that this module does not attempt to infer the pose of the board
   from the images. The pose is taken from the experiment, where it is known
   to the accuracy of the translation stage, which is far better than what
   fitting it to the images would give.

"""

from warnings import warn

from numpy import asarray, cross, float64, zeros, sqrt




def _normalize(v, name):
    '''
    Returns v as a unit vector, raising an informative error if it is null.
    '''
    v = asarray(v, dtype=float64).ravel()

    if v.shape != (3,):
        raise ValueError('%s must have three components, got %d'
                         %(name, v.size))

    n = sqrt((v**2).sum())
    if n == 0:
        raise ValueError('%s must not be the null vector'%name)

    return v/n




def lattice_variants(idx):
    '''
    Returns the eight ways of laying a grown lattice over the board, obtained
    by combining a transposition with a flip of either index.

    input -
    idx (array of ints, (n1,n2)) - the lattice, as returned by
                                  grid.grid_to_array

    output - a list of (variant, key) pairs, in which variant is the
             re-indexed lattice and key is the (transpose, flip1, flip2)
             tuple that produced it
    '''
    variants = []

    for transpose in [False, True]:
        A = idx.T if transpose else idx

        for flip1 in [False, True]:
            for flip2 in [False, True]:
                B = A[::-1, :] if flip1 else A
                B = B[:, ::-1] if flip2 else B
                variants.append((B, (transpose, flip1, flip2)))

    return variants




def resolve_orientation(points, idx, board_size, origin_hint=None,
                        iaxis_hint=None):
    '''
    Chooses which of the eight mappings of the lattice onto the board is the
    right one, as described in the module docstring.

    input -
    points (array, (N,2)) - the (x, y) corner positions
    idx (array of ints, (n1,n2)) - the lattice, as returned by
                                   grid.grid_to_array
    board_size (tuple, 2) - the number of internal corners of the board along
                            its two sides, (n_i, n_j)
    origin_hint (array, 2) - the approximate pixel position of the corner
                             chosen as the board origin. Strongly recommended;
                             required in practice for a multi camera
                             calibration.
    iaxis_hint (array, 2) - the approximate pixel position of a corner lying
                            along the first axis of the board, away from the
                            origin. Needed only for a square board, where the
                            origin alone leaves two possibilities.

    output -
    idx_oriented (array of ints, (n_i,n_j)) - the lattice re-indexed so that
                            site (i, j) is the board corner (i, j)
    '''
    ni, nj = int(board_size[0]), int(board_size[1])

    candidates = [(v, k) for v, k in lattice_variants(idx)
                  if v.shape == (ni, nj)]

    if len(candidates) == 0:
        raise ValueError(
            'the lattice found in the image has shape %s, which does not '
            'match the declared board_size of %s (nor its transpose). Either '
            'the board_size parameter is wrong, or the detection found only '
            'part of the board, or it grew past it.'
            %(str(idx.shape), str((ni, nj))))

    # with no hint we cannot do better than a convention, which is only safe
    # when a single camera is involved
    if origin_hint is None:
        warn('no origin_hint was given, so the orientation of the board is '
             'being fixed by a convention that depends on the direction from '
             'which the board is viewed. This is NOT safe for a multi camera '
             'calibration, where it can assign different lab coordinates to '
             'the same physical corner in different cameras. Give '
             'origin_hint to remove the ambiguity.', RuntimeWarning)

        # the corner nearest the origin of the image, breaking ties towards
        # the first axis running along the image x axis
        def fallback_cost(v):
            p00 = points[v[0, 0]]
            p10 = points[v[-1, 0]]
            return (p00[0]**2 + p00[1]**2,
                    -abs(p10[0] - p00[0]))

        best = min(candidates, key=lambda c: fallback_cost(c[0]))
        return best[0]

    origin_hint = asarray(origin_hint, dtype=float64).ravel()[:2]

    def cost(v):
        c = ((points[v[0, 0]] - origin_hint)**2).sum()
        if iaxis_hint is not None:
            h = asarray(iaxis_hint, dtype=float64).ravel()[:2]
            c += ((points[v[-1, 0]] - h)**2).sum()
        return c

    ordered = sorted(candidates, key=lambda c: cost(c[0]))
    best, best_cost = ordered[0][0], cost(ordered[0][0])

    # for a square board the origin hint alone leaves two mappings, which
    # differ in which side of the board is taken as the first axis
    if ni == nj and iaxis_hint is None and len(ordered) > 1:
        if cost(ordered[1][0]) <= 4.0*best_cost:
            warn('the board is square and no iaxis_hint was given, so which '
                 'side of the board is taken as its first axis is ambiguous. '
                 'Give iaxis_hint to remove the ambiguity.', RuntimeWarning)

    return best




def board_lab_coordinates(board_size, square_size, origin=(0., 0., 0.),
                          i_axis=(1., 0., 0.), j_axis=(0., 1., 0.),
                          translation=0.0):
    '''
    Returns the lab space coordinates of every internal corner of the board,
    using eq. (1).

    input -
    board_size (tuple, 2) - the number of internal corners along the two
                            sides of the board, (n_i, n_j)
    square_size (float) - the physical spacing of the corners, in the length
                          units used for the lab coordinates. If the two
                          sides of the board have different spacings, a pair
                          of values may be given.
    origin (array, 3) - the lab coordinates of the board corner (0, 0)
    i_axis (array, 3) - the lab space direction of the first side of the
                        board; normalized internally
    j_axis (array, 3) - the lab space direction of the second side
    translation (float or array of 3) - the displacement of the board for
                        this image. A single number is taken as a distance
                        along the board normal, i_axis x j_axis; a three
                        component vector is used as it stands.

    output -
    X (array, (n_i,n_j,3)) - X[i,j] holds the lab coordinates of the board
                             corner (i, j)
    '''
    ni, nj = int(board_size[0]), int(board_size[1])

    e1 = _normalize(i_axis, 'i_axis')
    e2 = _normalize(j_axis, 'j_axis')

    if abs(float(e1.dot(e2))) > 1e-6:
        warn('i_axis and j_axis are not perpendicular (their normalized dot '
             'product is %.3g). They are used as given, but a checkerboard '
             'normally has perpendicular sides.'%float(e1.dot(e2)),
             RuntimeWarning)

    s = asarray(square_size, dtype=float64).ravel()
    if s.size == 1:
        s1 = s2 = float(s[0])
    elif s.size == 2:
        s1, s2 = float(s[0]), float(s[1])
    else:
        raise ValueError('square_size must be one or two numbers')

    if s1 <= 0 or s2 <= 0:
        raise ValueError('square_size must be positive')

    X0 = asarray(origin, dtype=float64).ravel()
    if X0.shape != (3,):
        raise ValueError('origin must have three components')

    t = asarray(translation, dtype=float64).ravel()
    if t.size == 1:
        normal = cross(e1, e2)
        n = sqrt((normal**2).sum())
        if n == 0:
            raise ValueError('i_axis and j_axis are parallel, so the board '
                             'normal is undefined and a scalar translation '
                             'cannot be interpreted')
        t = float(t[0]) * normal/n
    elif t.size != 3:
        raise ValueError('translation must be a number or a vector of three '
                         'components')

    X = zeros((ni, nj, 3))
    for i in range(ni):
        for j in range(nj):
            X[i, j] = X0 + i*s1*e1 + j*s2*e2 + t

    return X
