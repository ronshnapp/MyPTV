# -*- coding: utf-8 -*-
"""
Created on Nov 14 2026

Tests for the matching of calibration blobs to target points by growing from
a few marked points (cal_point_matching_mod).

The tests use synthetic data, so the truth is known: a target is projected
through a camera that has a strong radial distortion and a wavy perturbation
on top of it (so that no low order model describes it), a few of the dots are
dropped, noise is added to the rest, and spurious blobs are thrown in.
"""

from numpy import array, cos, sin, sqrt, vstack, zeros
from numpy.linalg import norm
from numpy.random import default_rng

from myptv.cal_point_matching_mod import grow_calibration_matches



def distorted_camera(X, k1=-0.4, f=2600., res=(1600, 1200),
                     angles=(0.45, -0.35, 0.1), distance=420.):
    '''projects the lab points X (N,3) on a strongly distorted image'''
    a, b, c = angles
    Rx = array([[1, 0, 0], [0, cos(a), -sin(a)], [0, sin(a), cos(a)]])
    Ry = array([[cos(b), 0, sin(b)], [0, 1, 0], [-sin(b), 0, cos(b)]])
    Rz = array([[cos(c), -sin(c), 0], [sin(c), cos(c), 0], [0, 0, 1]])
    Xc = (X - X.mean(axis=0)).dot((Rz.dot(Ry).dot(Rx)).T)
    Xc = Xc + array([0, 0, distance])
    u = Xc[:, :2] / Xc[:, 2:3]
    r2 = (u**2).sum(axis=1)[:, None] * f**2 / (res[0]/2)**2
    u = u * (1 + k1*r2)
    u = u + 0.004*sin(18*u[:, ::-1])
    return u*f + array(res)/2



def stepped_target(nx=24, ny=18, d=4.0, levels=(0., 6., 12.), width=32.):
    '''a dot target made of flat strips that are at different depths'''
    return array([(i*d, j*d, levels[int(i*d//width) % len(levels)])
                  for i in range(nx) for j in range(ny)], dtype=float)



def make_blobs(targets, seed=0, drop=0.05, noise=0.3, n_spurious=20,
               res=(1600, 1200)):
    '''
    Returns the blobs, and a dictionary with the true blob position of each
    target point that is seen in the image.
    '''
    rng = default_rng(seed)
    x = distorted_camera(targets, res=res)
    seen = ((x[:, 0] > 20) & (x[:, 0] < res[0]-20) &
            (x[:, 1] > 20) & (x[:, 1] < res[1]-20) &
            (rng.uniform(size=len(targets)) > drop))
    xb = x[seen] + rng.normal(0, noise, size=(seen.sum(), 2))
    spurious = rng.uniform([0, 0], res, size=(n_spurious, 2))
    truth = {tuple(t): p for t, p in zip(targets[seen], xb)}
    return vstack([xb, spurious]), truth



def marked(truth, lab_points):
    '''the rows a user would have produced by marking the given points'''
    return [list(truth[tuple(X)]) + list(X) for X in lab_points]



def count(matcher, truth):
    '''the number of right and of wrong pairs'''
    right, wrong = 0, 0
    for p in matcher.point_pairs:
        key = tuple(p[2:])
        if key in truth and norm(truth[key] - p[:2]) < 1e-6:
            right += 1
        else:
            wrong += 1
    return right, wrong



def first_seen(truth, candidates):
    '''the first of the candidate lab points that is seen in the image'''
    for X in candidates:
        if tuple(X) in truth:
            return tuple(X)
    raise ValueError('none of the candidates is seen')



def test_flat_target_from_three_marked_points():
    '''Three neighbouring marked points are enough for a flat target.'''
    targets = stepped_target(levels=(0.,))
    blobs, truth = make_blobs(targets)
    seeds = [(40., 32., 0.), (44., 32., 0.), (40., 36., 0.)]
    seeds = [first_seen(truth, [s]) for s in seeds]

    m = grow_calibration_matches(blobs, targets, marked(truth, seeds),
                                 verbose=False)
    m.pair_points()
    right, wrong = count(m, truth)

    assert wrong == 0
    assert right >= 0.99 * len(truth)
    assert len(m.conflicts) == 0



def test_flat_part_is_not_left_without_a_marked_point():
    '''
    With the marked points all on one level of a stepped target, the other
    levels are left unpaired and are not guessed.
    '''
    targets = stepped_target()
    blobs, truth = make_blobs(targets)
    seeds = [(8., 32., 0.), (12., 32., 0.), (8., 36., 0.)]

    m = grow_calibration_matches(blobs, targets, marked(truth, seeds),
                                 verbose=False)
    m.pair_points()
    right, wrong = count(m, truth)

    assert wrong == 0
    assert all([p[4] == 0.0 for p in m.point_pairs])
    assert len(m.unmatched_blobs) > 0



def test_stepped_target_with_one_marked_point_per_extra_part():
    '''Three marked points on one strip and a single one on each other.'''
    targets = stepped_target()
    blobs, truth = make_blobs(targets)
    seeds = [(8., 32., 0.), (12., 32., 0.), (8., 36., 0.)]
    for x0, z in [(44., 6.), (76., 12.)]:
        seeds.append(first_seen(truth, [(x0, y, z) for y in
                                        (32., 36., 28., 40.)]))

    m = grow_calibration_matches(blobs, targets, marked(truth, seeds),
                                 verbose=False)
    m.pair_points()
    right, wrong = count(m, truth)

    assert wrong == 0
    assert right >= 0.99 * len(truth)
    assert len(m.conflicts) == 0



def test_irregular_target():
    '''The target does not have to be a regular grid.'''
    rng = default_rng(3)
    targets = stepped_target(levels=(0.,))
    targets[:, :2] += rng.uniform(-1.2, 1.2, size=(len(targets), 2))
    blobs, truth = make_blobs(targets)

    # three points near the middle that are seen and not on one line
    centre = targets.mean(axis=0)
    order = sorted(truth.keys(), key=lambda X: norm(array(X) - centre))
    seeds = [order[0]]
    for X in order[1:]:
        v = [array(s) - array(seeds[0]) for s in seeds[1:]] + \
            [array(X) - array(seeds[0])]
        if len(v) == 2:
            area = abs(v[0][0]*v[1][1] - v[0][1]*v[1][0])
            if area < 0.3 * norm(v[0]) * norm(v[1]):
                continue
        seeds.append(X)
        if len(seeds) == 3:
            break

    m = grow_calibration_matches(blobs, targets, marked(truth, seeds),
                                 verbose=False)
    m.pair_points()
    right, wrong = count(m, truth)

    assert wrong == 0
    assert right >= 0.98 * len(truth)



def test_marked_point_that_is_not_in_the_target_file_is_refused():
    targets = stepped_target(levels=(0.,))
    blobs, truth = make_blobs(targets)
    seeds = marked(truth, [(40., 32., 0.), (44., 32., 0.), (40., 36., 0.)])

    # a typing error: 63 instead of 36
    seeds[2][3] = 63.0

    m = grow_calibration_matches(blobs, targets, seeds, verbose=False)
    m.pair_points()

    assert len(m.seed_targets) == 2
    assert any(['not in the target file' in msg for msg in m.messages])



def test_marked_point_with_coordinates_of_its_neighbour_is_reported():
    '''
    A marked point that was given the lab coordinates of the dot next to
    it has to be reported, as it cannot be seen in the errors of the result.
    '''
    targets = stepped_target(levels=(0.,))
    blobs, truth = make_blobs(targets)
    seeds = marked(truth, [(40., 32., 0.), (44., 32., 0.), (40., 36., 0.),
                           (60., 48., 0.)])

    # the fourth point is given the coordinates of its neighbour
    seeds[3][2] = 64.0

    m = grow_calibration_matches(blobs, targets, seeds, verbose=False)
    m.pair_points()

    assert 4 in [c[0] for c in m.conflicts]



def test_results_can_be_saved_and_read_as_calibration_points(tmp_path):
    from myptv.utils import Cal_image_coord

    targets = stepped_target(levels=(0.,))
    blobs, truth = make_blobs(targets)
    seeds = marked(truth, [(40., 32., 0.), (44., 32., 0.), (40., 36., 0.)])

    m = grow_calibration_matches(blobs, targets, seeds, verbose=False)
    m.pair_points(verify=False)
    fname = str(tmp_path / 'cam_cal_points')
    m.save_results(fname)

    cic = Cal_image_coord(fname)
    assert cic.N_points == len(m.point_pairs)
    assert len(cic.lab_coords[0]) == 3 and len(cic.image_coords[0]) == 2
