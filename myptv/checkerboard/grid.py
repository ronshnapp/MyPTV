# -*- coding: utf-8 -*-
#!/usr/bin/env python3

"""
Created on Aug 27 2026

@author: Stefano Brizzolara


Assembling detected checkerboard corners into an ordered grid.

The detection step (see detection.py) returns an unordered cloud of corner
positions. To be of any use for a calibration, each corner has to be given a
pair of integer lattice indices (i, j) identifying which corner of the
physical board it is. This module performs that assignment.

The difficulty is that the mapping from the board to the image is not known
in advance: it involves perspective, so the corners in the image are not a
regular grid, and it may involve lens distortion, so the rows and columns are
not even straight. The approach taken here is therefore a local one, in the
spirit of the region growing method of Geiger et al. (2012): rather than
fitting a global model, the grid is grown outwards one corner at a time, and
each new corner is predicted from its already placed neighbours.


The growth works as follows.

1) A seed corner is chosen and given the indices (0, 0). Its four lattice
   neighbours are found by looking, along each of the two edge directions
   returned by the detector, for the nearest corner lying within a narrow
   angular cone. These are given the indices (1,0), (-1,0), (0,1) and (0,-1).

2) The grid is then grown repeatedly. At each pass, every empty lattice site
   adjacent to an occupied one is considered, its position in the image is
   predicted from the occupied sites around it, and the nearest detected
   corner that has not yet been used is assigned to it, provided it lies
   within a tolerance of the prediction.

   Two kinds of prediction are used, both of which are exact for an affine
   map and accurate for a smooth one over the span of a single cell:

       collinear    p(i+1,j) = 2 p(i,j) - p(i-1,j)
       parallelogram    p(i,j) = p(i-1,j) + p(i,j-1) - p(i-1,j-1)

   All the available predictions are averaged, which makes the growth robust
   to a single badly placed neighbour.

3) The pass is repeated until no further corner can be added. The occupied
   sites are then packed into a rectangular array of lattice indices.

Because the growth is local, it follows the curvature of a distorted image
and stops of its own accord at the edge of the board.

Reference:
Geiger, A., Moosmann, F., Car, O., & Schuster, B. (2012). Automatic camera
and range sensor calibration using a single shot. IEEE International
Conference on Robotics and Automation, 3936-3943.

"""

from numpy import (array, zeros, full, sqrt, median, dot, argsort, cos,
                   inf, isfinite)

from scipy.spatial import cKDTree




def estimate_spacing(points):
    '''
    Estimates the typical spacing, in pixels, between neighbouring corners of
    the board, as the median distance from a corner to its nearest neighbour.
    For a checkerboard the nearest neighbour of a corner is one of its
    lattice neighbours, so this is a robust estimate even when part of the
    cloud is spurious.

    input -
    points (array, (N,2)) - the (x, y) corner positions

    output -
    spacing (float) - the estimated lattice spacing in pixels
    '''
    if len(points) < 2:
        raise ValueError('at least two corners are needed to estimate the '
                         'lattice spacing')

    tree = cKDTree(points)
    d, _ = tree.query(points, k=2)

    return float(median(d[:, 1]))




def find_neighbor_in_direction(points, tree, i, direction, spacing,
                              used=None, max_angle=0.45, r_min=0.5,
                              r_max=1.7):
    '''
    Finds the lattice neighbour of corner i lying in a given direction.

    Candidates are the corners at a distance between r_min*spacing and
    r_max*spacing from corner i whose bearing differs from the given
    direction by less than max_angle. Of these, the nearest is returned.

    input -
    points (array, (N,2)) - the corner positions
    tree (cKDTree) - a tree built on points
    i (int) - index of the corner whose neighbour is sought
    direction (array, 2) - a unit vector giving the direction to search in
    spacing (float) - the estimated lattice spacing, in pixels
    used (set) - indices that are already assigned and may not be reused
    max_angle (float) - half the opening angle of the search cone, in radians
    r_min, r_max (floats) - the range of accepted distances, in units of
                            spacing

    output -
    j (int) - the index of the neighbour, or None if there is none
    '''
    if used is None:
        used = set()

    cos_max = cos(max_angle)

    candidates = tree.query_ball_point(points[i], r_max*spacing)

    best, best_d = None, inf

    for j in candidates:
        if j == i or j in used:
            continue

        w = points[j] - points[i]
        d = sqrt(w[0]**2 + w[1]**2)

        if d < r_min*spacing or d > r_max*spacing:
            continue

        if dot(w/d, direction) < cos_max:
            continue

        if d < best_d:
            best, best_d = j, d

    return best




def _predictions(grid, points, i, j):
    '''
    Returns the predicted image position of the empty lattice site (i, j),
    together with a local length scale, based on the occupied sites around
    it. See the module docstring for the two prediction rules used.

    input -
    grid (dict) - maps (i,j) tuples to indices into points
    points (array, (N,2)) - the corner positions
    i, j (ints) - the lattice site to predict

    output -
    p (array, 2) - the predicted (x, y) position, or None if the site cannot
                   be predicted
    scale (float) - a local estimate of the lattice spacing in pixels, or
                    None
    '''
    preds, scales = [], []

    # collinear extrapolation along the four directions
    for di, dj in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
        a = (i + di, j + dj)
        b = (i + 2*di, j + 2*dj)
        if a in grid and b in grid:
            pa, pb = points[grid[a]], points[grid[b]]
            preds.append(2*pa - pb)
            scales.append(sqrt(((pa - pb)**2).sum()))

    # parallelogram completion over the four diagonal quadrants
    for di, dj in [(-1, -1), (-1, 1), (1, -1), (1, 1)]:
        a = (i + di, j)
        b = (i, j + dj)
        c = (i + di, j + dj)
        if a in grid and b in grid and c in grid:
            pa, pb, pc = points[grid[a]], points[grid[b]], points[grid[c]]
            preds.append(pa + pb - pc)
            scales.append(0.5*(sqrt(((pa - pc)**2).sum()) +
                               sqrt(((pb - pc)**2).sum())))

    if len(preds) == 0:
        return None, None

    p = array(preds).mean(axis=0)
    scale = float(median(array(scales)))

    if not isfinite(p).all() or not isfinite(scale) or scale <= 0:
        return None, None

    return p, scale




def grow_grid(points, d1, d2, seed, spacing=None, tol=0.35, max_sites=100000):
    '''
    Grows a lattice of corners outwards from a seed corner, as described in
    the module docstring.

    input -
    points (array, (N,2)) - the corner positions
    d1, d2 (arrays, (N,2)) - the two edge directions of each corner
    seed (int) - index of the corner to start from, which is given the
                 indices (0, 0)
    spacing (float) - the lattice spacing in pixels; estimated from points if
                      not given
    tol (float) - a corner is accepted at a lattice site if it lies within
                  tol times the local spacing of the predicted position
    max_sites (int) - a safety limit on the size of the grid

    output -
    grid (dict) - maps (i,j) tuples to indices into points. Empty if the seed
                  could not be given four neighbours.
    '''
    if spacing is None:
        spacing = estimate_spacing(points)

    tree = cKDTree(points)

    grid = {(0, 0): int(seed)}
    used = {int(seed)}

    # step 1: the four neighbours of the seed, along +-d1 and +-d2
    e1, e2 = d1[seed], d2[seed]

    for ij, direction in [((1, 0), e1), ((-1, 0), -e1),
                          ((0, 1), e2), ((0, -1), -e2)]:

        j = find_neighbor_in_direction(points, tree, seed, direction,
                                       spacing, used=used)
        if j is not None:
            grid[ij] = int(j)
            used.add(int(j))

    # the seed needs at least one neighbour along each axis for the growth to
    # have somewhere to go
    have_1 = (1, 0) in grid or (-1, 0) in grid
    have_2 = (0, 1) in grid or (0, -1) in grid

    if not (have_1 and have_2):
        return {}

    # step 2: grow until nothing more can be added
    while len(grid) < max_sites:

        # the empty sites next to an occupied one
        frontier = set()
        for (i, j) in grid:
            for di, dj in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
                if (i+di, j+dj) not in grid:
                    frontier.add((i+di, j+dj))

        if len(frontier) == 0:
            break

        # try the best predicted sites first, so that well constrained sites
        # claim their corner before poorly constrained ones do
        proposals = []
        for (i, j) in frontier:
            p, scale = _predictions(grid, points, i, j)
            if p is None:
                continue
            proposals.append((i, j, p, scale))

        if len(proposals) == 0:
            break

        added = 0
        for (i, j, p, scale) in proposals:

            if (i, j) in grid:
                continue

            d, k = tree.query(p, k=1)

            if k >= len(points) or int(k) in used:
                continue

            if d > tol*scale:
                continue

            grid[(i, j)] = int(k)
            used.add(int(k))
            added += 1

        if added == 0:
            break

    return grid




def grid_to_array(grid):
    '''
    Packs the dictionary returned by grow_grid into a rectangular array of
    indices, with -1 marking the lattice sites at which no corner was found.

    input -
    grid (dict) - maps (i,j) tuples to indices into the corner array

    output -
    idx (array of ints, (n1,n2)) - idx[a,b] is the index of the corner at
                                   lattice site (i_min+a, j_min+b), or -1
    '''
    if len(grid) == 0:
        return zeros((0, 0), dtype=int) - 1

    keys = list(grid.keys())
    i_min = min(k[0] for k in keys)
    i_max = max(k[0] for k in keys)
    j_min = min(k[1] for k in keys)
    j_max = max(k[1] for k in keys)

    idx = full((i_max - i_min + 1, j_max - j_min + 1), -1, dtype=int)

    for (i, j), k in grid.items():
        idx[i - i_min, j - j_min] = k

    return idx




def assemble_grid(detection, spacing=None, tol=0.35, n_seeds=12):
    '''
    Assembles the corners found by detection.detect_corners into an ordered
    grid. Several seeds are tried, in order of decreasing corner score, and
    the largest resulting grid is kept; this makes the result insensitive to
    a poor choice of seed.

    input -
    detection (dict) - the dictionary returned by detection.detect_corners
    spacing (float) - the lattice spacing in pixels; estimated if not given
    tol (float) - the snapping tolerance, in units of the local spacing
    n_seeds (int) - how many seeds to try

    output -
    idx (array of ints, (n1,n2)) - the lattice, as returned by grid_to_array
    '''
    points = detection['points']
    d1, d2 = detection['d1'], detection['d2']
    scores = detection['scores']

    if len(points) < 4:
        return zeros((0, 0), dtype=int) - 1

    if spacing is None:
        spacing = estimate_spacing(points)

    # the strongest corners make the most reliable seeds
    order = argsort(scores)[::-1][:int(n_seeds)]

    best = {}
    for seed in order:
        grid = grow_grid(points, d1, d2, int(seed), spacing=spacing, tol=tol)
        if len(grid) > len(best):
            best = grid

    return grid_to_array(best)
