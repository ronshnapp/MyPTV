# -*- coding: utf-8 -*-
"""
Created on Nov 14 2026

Matching segmented calibration blobs to the points of a target file by
growing outwards from a few hand marked points.


The problem.

To calibrate a camera, each dot that was segmented in the calibration image
has to be told which point of the target file it is. The traditional way of
doing so in MyPTV (utils.match_calibration_blobs_and_points) is to fit a
rough camera model to 10-15 hand marked points, project the whole target file
through it, and pair every blob with the nearest projection. That only works
if the rough model is accurate to better than half a dot spacing over the
whole image, which a low order model fitted to a dozen points cannot deliver
when the optics are complicated (strong lens distortion, refraction through
windows and so on). The user is then left to mark more and more points.


The idea.

Identifying the dots does not require a model that is accurate everywhere.
However complicated the optics are, the mapping from the target to the image
is smooth, and so over the distance between neighbouring dots it is very
nearly affine. So, instead of projecting the whole target at once, we
start from the few marked points and work outwards:

    1) Local growth. An unmatched target point that sits next to points
       that were already matched is projected on the image with an affine
       map fitted to those matched neighbours only. If a free blob lies
       within a small fraction of the local dot spacing from the predicted
       position, the two are paired. The newly paired points then serve
       as neighbours to the points beyond them, and so on, until nothing
       can be added. The affine map is only ever extrapolated by about one
       dot spacing, which is why its simplicity does not matter.

    2) Cleaning. Once the growth stalls, pairs that do not agree with
       their surroundings are removed: each matched point is predicted from
       its neighbours as in step 1, and the pair is dropped if its blob is
       not where the neighbours say it should be.

    3) Global refit. A polynomial model is then fitted to everything that
       was matched so far, with an order that is raised as more points
       become available. This model has two roles: pairs that disagree with
       it are thrown away, and target points that were skipped by the local
       growth are picked up through it, provided that they are surrounded
       by matched points, so the model is not extrapolated. If steps 2 or 3
       changed anything, we go back to step 1.

Nothing here assumes a particular target shape. The target points do not need
to form a regular grid, to have a uniform spacing, or to lie on planes; the
only assumptions are that neighbouring target points are also neighbours in
the image, and that the blobs are separated from each other by more than the
error of a one-spacing affine extrapolation.


What the user has to mark.

The local affine map can only predict the image of a target point from
neighbours that surround it in the directions that matter:

    - to start growing over a flat part of a target, mark 3 points on it
      that are close to each other and not on a single line. (Points that
      are spread over the target, as one would mark for an initial
      calibration, work too as long as the image is not distorted much,
      but neighbouring points do not have this limitation);

    - a flat part that is displaced out of the plane of the points matched so
      far (another level of a stepped target, or another position of a
      translated plate) cannot be predicted from them, because nothing that
      was matched says how the image changes in the third direction. One
      marked point on each such part is usually enough: the in-plane
      behaviour is then borrowed from the already matched parts for the
      first steps. If nothing grows from it, mark two more next to it. Parts
      that are close enough to two matched parts at different depths (e.g.
      the ridges of a stepped target) are reached without marking anything
      on them, since the line of sight is then known.

Marking more points than this minimum does no harm.

A warning is in place regarding periodic targets. If a marked point is given
the lab coordinates of the dot next to the one that was clicked, the growth
from it gives a perfectly self consistent, and wrong, result. This cannot be
detected from the data alone, since a regular grid looks the same when it is
shifted by one period. Two safeguards are used: marked points whose lab
coordinates are not found in the target file are refused, and marked points
that disagree with the final global model are reported.
"""

from numpy import array, asarray, zeros, sqrt, median
from numpy import loadtxt, savetxt, hstack, abs as npabs
from numpy.linalg import norm, svd, lstsq, pinv
from scipy.spatial import cKDTree




def load_blobs_xy(blob_fname):
    '''
    Reads a MyPTV blob file and returns an (N,2) array of the blob
    centers as (x, y) image coordinates. (In blob files, the first column is
    the y coordinate and the second one is x.)
    '''
    blobs = loadtxt(blob_fname, ndmin=2)
    return blobs[:, 1::-1].copy()




def polynomial_terms(order):
    '''
    Returns the list of exponents (i,j,k) of all the monomials
    X^i Y^j Z^k with i+j+k <= order.
    '''
    terms = []
    for n in range(order+1):
        for i in range(n, -1, -1):
            for j in range(n-i, -1, -1):
                terms.append((i, j, n-i-j))
    return terms




class polynomial_projection(object):
    '''
    A plain polynomial map from lab space to image space, fitted by least
    squares. It is used here as a model-independent stand in for a camera:
    all that the matching needs is a smooth 3D -> 2D function.
    '''

    def __init__(self, X, x, order):
        '''
        X - (N,3) array, lab space points
        x - (N,2) array, their image space coordinates
        order - int, the polynomial order
        '''
        X = asarray(X, dtype=float)
        x = asarray(x, dtype=float)
        self.order = order
        self.terms = polynomial_terms(order)

        # the coordinates are normalized so the problem is well conditioned
        self.X0 = X.mean(axis=0)
        scale = (X - self.X0).std(axis=0)
        scale[scale < 1e-12 * max(scale.max(), 1e-300)] = 1.0
        self.scale = scale

        # lstsq returns the minimum norm solution, so directions that the
        # points do not constrain (e.g. Z^3 for a target with 3 levels) get
        # no contribution instead of blowing up.
        A = self._design(X)
        self.coefs = lstsq(A, x, rcond=1e-10)[0]
        self._pinv = pinv(A, rcond=1e-10)

    def _design(self, X):
        Xn = (asarray(X, dtype=float) - self.X0) / self.scale
        cols = [Xn[:, 0]**i * Xn[:, 1]**j * Xn[:, 2]**k
                for i, j, k in self.terms]
        return array(cols).T

    def __call__(self, X):
        '''returns the (N,2) image coordinates of the (N,3) lab points'''
        X = asarray(X, dtype=float)
        if X.ndim == 1:
            return self._design(X[None, :]).dot(self.coefs)[0]
        return self._design(X).dot(self.coefs)

    def leverage(self, X):
        '''
        Returns, for each of the (N,3) lab points, the factor by which the
        noise of the fitted image points is amplified in the prediction made
        at that point (in terms of the variance). It is of the order of
        n_terms/n_points, and so much smaller than 1, at places that are
        surrounded by fitted points, and it shoots up where the polynomial
        is being extrapolated.
        '''
        a = self._design(asarray(X, dtype=float).reshape(-1, 3))
        return (a.dot(self._pinv)**2).sum(axis=1)




def _fmt(X):
    '''a short string for a point's coordinates'''
    return '(' + ', '.join(['%g' % (v + 0.0) for v in X]) + ')'




class grow_calibration_matches(object):
    '''
    Pairs segmented calibration blobs with target file points, starting from
    a few hand marked points and growing outwards. See the module docstring
    for the description of the method.

    To do the pairing, use
    grow_calibration_matches.pair_points()

    To save the results, use
    grow_calibration_matches.save_results(fname)

    The attribute point_pairs and the method save_results are the same as in
    utils.match_calibration_blobs_and_points, so the two can be used
    interchangeably.
    '''

    def __init__(self, blobs_xy, targets, marked_points, tol=0.3,
                 n_neighbors=8, max_reach=3.0, global_tol=0.25,
                 max_order=3, verbose=True):
        '''
        input -

        blobs_xy - (N,2) array of the (x,y) image coordinates of the segmented
                   blobs, or a string with the path of a blob file.

        targets - (M,3) array with the lab space coordinates of the target
                  points, or a string with the path of a target file.

        marked_points - list of hand marked points, each being
                        [x_image, y_image, x_lab, y_lab, z_lab], or a string
                        with the path of a file holding them (namely, the
                        "manualPoints" file of the initial calibration GUI).

        tol - float. In the local growth, a blob is paired with a target
              point if its distance from the predicted position is less than
              tol times the distance between the blob and its nearest
              neighbour. Should be well under 0.5.

        n_neighbors - int, the number of matched neighbours used to fit each
                      local affine map.

        max_reach - float. A target point is only predicted from matched
                    points that are closer to it than max_reach times its
                    distance to its nearest neighbour in the target file.

        global_tol - float, same as tol but for the points paired through the
                     global polynomial model.

        max_order - int, the highest order of the global polynomial (1-3).
        '''
        if isinstance(blobs_xy, str):
            blobs_xy = load_blobs_xy(blobs_xy)
        if isinstance(targets, str):
            targets = loadtxt(targets, ndmin=2)
        if isinstance(marked_points, str):
            marked_points = loadtxt(marked_points, ndmin=2)

        self.blobs = asarray(blobs_xy, dtype=float).reshape(-1, 2)
        self.targets = asarray(targets, dtype=float).reshape(-1, 3)
        self.marked = asarray(marked_points, dtype=float).reshape(-1, 5)

        self.tol = tol
        self.k = n_neighbors
        self.k_wide = 8 * n_neighbors
        self.near_reach = 1.5
        self.max_out_of_span = 0.1
        self.max_amplification = 5.0
        self.max_leverage = 1.0
        self.max_support = 6.0
        self.max_reach = max_reach
        self.global_tol = global_tol
        self.max_order = max_order
        self.verbose = verbose

        if len(self.blobs) < 3 or len(self.targets) < 3:
            raise ValueError('Too few blobs or target points to match.')

        # the local spacing of the blobs and of the target points, taken as
        # the distance of each one from its nearest neighbour
        self.blob_tree = cKDTree(self.blobs)
        self.blob_spacing = self.blob_tree.query(self.blobs, k=2)[0][:, 1]
        self.target_tree = cKDTree(self.targets)
        self.target_spacing = self.target_tree.query(self.targets, k=2)[0][:, 1]

        self.messages = []
        self.conflicts = []
        self.unverified = []
        self.point_pairs = []
        self.model = None
        self.model_rms = None

        self._reset()


    def _say(self, msg):
        self.messages.append(msg)
        if self.verbose:
            print(msg)


    def _reset(self):
        # t2b[i] = index of the blob paired with target i (-1 if not paired)
        self.t2b = zeros(len(self.targets), dtype=int) - 1
        self.b2t = zeros(len(self.blobs), dtype=int) - 1
        # how each pair was found: 0-marked, 1-local growth, 2-global model,
        # 3-local growth using a wide neighbourhood
        self.origin = zeros(len(self.targets), dtype=int) - 1
        self.rejected = set()
        self.seed_targets = []
        self.seed_numbers = []


    def _pair(self, ti, bi, origin):
        self.t2b[ti] = bi
        self.b2t[bi] = ti
        self.origin[ti] = origin


    def _unpair(self, ti):
        bi = self.t2b[ti]
        self.rejected.add((int(ti), int(bi)))
        self.b2t[bi] = -1
        self.t2b[ti] = -1
        self.origin[ti] = -1



    def set_marked_points(self, leave_out=None):
        '''
        Ties each hand marked point to a blob and to a target file point.
        Points that cannot be tied unambiguously are left out with a message,
        since a wrong starting point spoils everything that grows from it.
        '''
        self._reset()

        for n, p in enumerate(self.marked):
            if leave_out is not None and n+1 == leave_out:
                continue
            x_img, X_lab = p[:2], p[2:]
            label = 'Marked point %d, lab: %s,' % (n+1, _fmt(X_lab))

            # the lab coordinates must be those of a target file point
            d, ti = self.target_tree.query(X_lab)
            if d > 0.05 * self.target_spacing[ti]:
                self._say('%s is not in the target file (the nearest target '
                          'point is %s); it is ignored. Check for a typing '
                          'error.' % (label, _fmt(self.targets[ti])))
                continue

            # the image coordinates must fall on a segmented blob
            d, bi = self.blob_tree.query(x_img)
            if d > 0.5 * self.blob_spacing[bi]:
                self._say('%s was not marked on a segmented blob; it is '
                          'ignored.' % label)
                continue

            if self.t2b[ti] == bi:
                continue

            if self.t2b[ti] != -1 or self.b2t[bi] != -1:
                self._say('%s conflicts with another marked point (same blob '
                          'or same lab coordinates); it is ignored.' % label)
                continue

            self._pair(ti, bi, 0)
            self.seed_targets.append(int(ti))
            self.seed_numbers.append(n+1)

        return len(self.seed_targets)



    def _local_prediction(self, ti, nbrs):
        '''
        Predicts the image position of target point ti with an affine map
        fitted to the matched target points nbrs (indices, sorted by their
        distance from ti). The map is made to pass exactly through the
        nearest of them. Returns None if the neighbours do not determine the
        map in the direction that leads to ti. Otherwise, returns the
        predicted position and a flag that is True if the prediction leans
        on far away neighbours (see _step_coefficients).
        '''
        a = nbrs[0]
        if len(nbrs) < 2:
            return None

        dXt = self.targets[ti] - self.targets[a]
        xa = self.blobs[self.t2b[a]]
        dX = self.targets[nbrs[1:]] - self.targets[a]
        dx = self.blobs[self.t2b[nbrs[1:]]] - xa

        # neighbours nearer to ti count more, yet far ones are never
        # neglected altogether, as they may be the only source for the
        # derivatives in some direction.
        dist = norm(self.targets[nbrs[1:]] - self.targets[ti], axis=1)
        w = 1.0 / (1.0 + (dist / max(dist.min(), 1e-300))**2)

        res = self._step_coefficients(dX, w, dXt, dist,
                                      self.target_spacing[ti])
        if res is None:
            return None
        c, far = res
        return xa + c.dot(dx), far



    def _step_coefficients(self, dX, w, dXt, dist, spacing):
        '''
        Solves the weighted least squares problem behind the local
        predictions. Given the steps dX from the anchor to its neighbours and
        their weights w, it returns the coefficients c such that the image of
        the step dXt is predicted as sum_i c_i dx_i, with dx_i the image
        steps to the neighbours. dist holds the distances of the neighbours
        from the point being predicted and spacing is the local dot spacing.
        Returns None if the step dXt is not determined by the neighbours, if
        the prediction is too sensitive to the noise in their positions, or
        if it is a long step that leans on neighbours that are far away.
        Otherwise, returns c and a flag that is True if the prediction leans
        on far away neighbours.
        '''
        U, s, Vt = svd(dX * w[:, None], full_matrices=False)
        if s[0] == 0:
            return None
        r = int((s > 1e-4 * s[0]).sum())
        U, s, Vt = U[:, :r], s[:r], Vt[:r]

        # The step from the anchor to ti has to be within the subspace that
        # the neighbours span; otherwise the prediction is not determined
        # (e.g. a point above a plane on which all the neighbours lie). A
        # small part of the step is allowed to be out of it, to allow for
        # targets that are slightly curved.
        if (norm(dXt - Vt.T.dot(Vt.dot(dXt))) >
                self.max_out_of_span * spacing):
            return None

        # The prediction is a linear combination of the positions of the
        # neighbouring blobs. If the coefficients are large the prediction
        # amplifies the noise of the blob positions and cannot be trusted.
        # That happens when the step relies on a direction in which the
        # neighbours hardly extend. For reference, the plain collinear
        # extrapolation, x3 = 2*x2 - x1, has an amplification of sqrt(5).
        c = (dXt.dot(Vt.T) / s).dot(U.T) * w
        if sqrt((c**2).sum() + (1.0 - c.sum())**2) > self.max_amplification:
            return None

        # An affine map is good over a few dot spacings and not much more.
        # So, the neighbours that the prediction actually leans on, i.e.
        # those with the large coefficients, should be near the point that
        # is predicted. Otherwise, with few matched points around (as it is
        # right at the start), the map would be set by far away points and
        # still look perfectly determined.
        #
        # There is one case in which leaning on far away points is
        # tolerable, and that is a step to an immediate neighbour of the
        # anchor. To be sent to a wrong blob there the map has to be off by
        # about 100%, whereas the image scale and orientation hardly ever
        # change that much across a target. This is what lets the growth
        # start from marked points that are spread over the target, or from
        # a single point marked on a part of it. The caller is told about it
        # through the flag, and tightens the tolerance.
        far = (npabs(c).dot(dist) >
               self.max_support * spacing * max(npabs(c).sum(), 1e-300))
        if far and norm(dXt) > self.near_reach * spacing:
            return None

        return c, far



    def _wide_prediction(self, ti, nbrs):
        '''
        Predicts the image position of target point ti in the case that its
        nearest matched neighbours are not enough, e.g. when they all lie on
        a plane and ti is off it, or when there is just one of them. A wider
        set of matched points, nbrs, is then used.

        Fitting a plain 3D affine map to a wide neighbourhood is not good
        enough: how the image changes along the directions that only far
        points tell about would be measured where the far points are, and
        with a distorted image this is not how it changes where ti is. But,
        whatever the lens does, it does it to the image, and so the lines of
        sight in lab space are hardly affected and they vary slowly over the
        target. Therefore, this is done in three steps:

            1) a 3D affine map is fitted to the wide neighbourhood, and the
               direction that it maps to zero is taken as the local line of
               sight, d. A step along d does not move the image.

            2) all the steps from the anchor are projected along d on the
               plane perpendicular to it, which makes the problem two
               dimensional whatever the shape of the target is.

            3) a 2D affine map is fitted to the projected steps with weights
               that favor the points nearest to ti, and it is used only if
               the points it leans on are indeed near ti.

        As a safeguard, the result has to agree with that of the plain 3D
        affine map of step 1. Where the two disagree the extrapolation is
        not to be trusted, and it is better to leave the point alone than to
        pair it with a wrong blob.

        Returns the predicted position, the length of the step along the
        line of sight in units of the dot spacing, and the flag telling if
        the prediction leans on far away neighbours; or None.
        '''
        a = nbrs[0]
        if len(nbrs) < 4:
            return None

        dXt = self.targets[ti] - self.targets[a]
        xa = self.blobs[self.t2b[a]]
        dX = self.targets[nbrs[1:]] - self.targets[a]
        dx = self.blobs[self.t2b[nbrs[1:]]] - xa
        dist = norm(self.targets[nbrs[1:]] - self.targets[ti], axis=1)
        w = 1.0 / (1.0 + (dist / max(dist.min(), 1e-300))**2)
        spacing = self.target_spacing[ti]

        # 1) the 3D affine map, dx = dX.dot(J), and the line of sight
        U, s, Vt = svd(dX * w[:, None], full_matrices=False)
        if s[0] == 0 or s[-1] < 1e-4 * s[0]:
            return None
        J = (Vt.T / s).dot(U.T.dot(dx * w[:, None]))        # (3,2)
        x_affine = xa + dXt.dot(J)

        Uj, sj, Vjt = svd(J, full_matrices=True)
        if sj[-1] < 1e-6 * sj[0]:
            return None
        E = Uj[:, :2]             # spans the plane perpendicular to d

        # 2+3) the 2D affine map between the projected steps and the image
        w2 = w**2
        # (the line of sight may be taken from far away, but how the image
        # changes across it may not, as this is what the lens distorts; the
        # check for that is within _step_coefficients.)
        res = self._step_coefficients(dX.dot(E), w2, dXt.dot(E), dist,
                                      spacing)
        if res is None:
            return None
        c, far = res
        if far and norm(dXt) > self.near_reach * spacing:
            return None
        x_pred = xa + c.dot(dx)

        # the safeguard
        bi = self.blob_tree.query(x_pred)[1]
        if norm(x_pred - x_affine) > 0.5 * self.tol * self.blob_spacing[bi]:
            return None

        # how far along the line of sight the step goes, in dot spacings
        depth = abs(dXt.dot(Uj[:, 2])) / spacing

        return x_pred, depth, far



    def grow_local(self):
        '''
        The local growth (step 1 of the module docstring). Returns the number
        of pairs that were added.
        '''
        n_added = 0

        while True:
            matched = (self.t2b != -1).nonzero()[0]
            free = (self.t2b == -1).nonzero()[0]
            if len(matched) < 2 or len(free) == 0:
                break

            # the candidates are the free target points adjacent to matched
            tree = cKDTree(self.targets[matched])
            k = min(self.k, len(matched))
            dist, ind = tree.query(self.targets[free], k=k)
            if k == 1:
                dist, ind = dist[:, None], ind[:, None]

            # the distance of each candidate from the nearest matched point,
            # in units of the local spacing of the target points
            reach = dist[:, 0] / self.target_spacing[free]

            proposals = []
            for row in (reach <= self.max_reach).nonzero()[0]:
                ti = free[row]
                # The proposals are ranked by how safe they are:
                #   0 - an immediate neighbour of matched points, predicted
                #       from its nearest neighbours;
                #   1 - a longer step predicted likewise (over a missing dot
                #       or a gap);
                #   2 - a step that needed a wide neighbourhood but stays at
                #       the depth of the point it starts from (the first
                #       steps from a point marked alone on a part of the
                #       target);
                #   3 - a step to a different depth (to a part of the target
                #       that has no marked point on it).
                # The tolerance is halved for rank 3, and for steps that
                # lean on far away points, which are given rank 2 at least.
                x_pred, wide, far, tol = None, False, False, self.tol
                rank = 0 if reach[row] <= self.near_reach else 1

                res = self._local_prediction(ti, matched[ind[row]])
                if res is not None:
                    x_pred, far = res

                elif len(matched) > k:
                    # the nearest neighbours do not tell how the image changes
                    # in the direction of this point (e.g. they are all on
                    # one plane, and the point is off it). Matched points
                    # further away might.
                    k_wide = min(self.k_wide, len(matched))
                    ind_wide = tree.query(self.targets[ti], k=k_wide)[1]
                    res = self._wide_prediction(ti, matched[ind_wide])
                    if res is not None:
                        x_pred, depth, far = res
                        wide = True
                        rank = 2
                        if depth > 0.5:
                            rank, tol = 3, 0.5 * self.tol

                if x_pred is None:
                    continue

                if far:
                    rank, tol = max(rank, 2), 0.5 * self.tol

                d, bi = self.blob_tree.query(x_pred)
                if self.b2t[bi] != -1 or (int(ti), int(bi)) in self.rejected:
                    continue

                err = d / self.blob_spacing[bi]
                if err < tol:
                    proposals.append((err, rank, ti, bi, wide))

            # Safe steps go first: in each pass only the proposals of the
            # safest rank that is present are accepted. So, the riskier steps
            # get their turn only once the safer ones have run out, by which
            # time there are more matched points around to predict them from
            # and fewer free blobs to confuse them with.
            if len(proposals) > 0:
                best_rank = min([p[1] for p in proposals])
                proposals = [p for p in proposals if p[1] == best_rank]

            # accept the proposals, best first, each blob being used once
            n_pass = 0
            for err, rank, ti, bi, wide in sorted(proposals):
                if self.b2t[bi] == -1 and self.t2b[ti] == -1:
                    self._pair(ti, bi, 3 if wide else 1)
                    n_pass += 1

            n_added += n_pass
            if n_pass == 0:
                break

        return n_added



    def remove_local_outliers(self):
        '''
        Removes pairs that do not agree with their surroundings: each matched
        target point is predicted from its matched neighbours, as in the
        local growth but with the point itself left out, and the pair is
        removed if the blob is far from the prediction, both relative to how
        well this works for the other pairs and relative to the dot spacing.
        This catches blobs that are not dots of the target but happen to be
        near where a dot should be, and dots whose center came out badly in
        the segmentation. The marked points are not removed. Returns the
        number of pairs removed.
        '''
        matched = (self.t2b != -1).nonzero()[0]
        if len(matched) < 8:
            return 0

        tree = cKDTree(self.targets[matched])
        k = min(self.k + 1, len(matched))
        ind = tree.query(self.targets[matched], k=k)[1]

        res = zeros(len(matched)) - 1.0
        for row, ti in enumerate(matched):
            pred = self._local_prediction(ti, matched[ind[row, 1:]])
            if pred is not None and not pred[1]:
                res[row] = norm(pred[0] - self.blobs[self.t2b[ti]])

        checked = res >= 0
        if checked.sum() < 8:
            return 0

        sigma = 1.4826 * median(res[checked])
        spacing = self.blob_spacing[self.t2b[matched]]
        bad = checked & (res > 5.0 * sigma) & (res > 0.1 * spacing)

        n_removed = 0
        for ti in matched[bad]:
            if self.origin[ti] != 0:
                self._unpair(ti)
                n_removed += 1
        return n_removed



    def fit_global_model(self):
        '''
        Fits the global polynomial to the present pairs. The order is set by
        the number of pairs that are available. Returns False if there are too
        few of them or if they all lie on a single plane, in which case
        a 3D model cannot be had.
        '''
        self.model, self.model_low, self.model_rms = None, None, None
        matched = (self.t2b != -1).nonzero()[0]
        X = self.targets[matched]
        x = self.blobs[self.t2b[matched]]

        if len(matched) < 12:
            return False

        # coplanar points cannot say how the image varies off their plane
        s = svd(X - X.mean(axis=0), compute_uv=False)
        if s[-1] < 1e-6 * s[0]:
            return False

        order = 1
        for o in range(2, self.max_order+1):
            if len(matched) >= 4 * len(polynomial_terms(o)):
                order = o

        self.model = polynomial_projection(X, x, order)
        self.model_low = None
        if order > 1:
            self.model_low = polynomial_projection(X, x, order-1)
        self.model_residuals = norm(self.model(X) - x, axis=1)
        self.model_matched = matched
        self.model_rms = sqrt((self.model_residuals**2).mean())
        return True



    def refine_global(self):
        '''
        The global refit (step 3 of the module docstring). Returns the
        number of pairs that were removed and the number that were added.
        '''
        if not self.fit_global_model():
            return 0, 0

        # (a) remove the pairs that disagree with the model. The marked
        # points are never removed here; they are reported at the end.
        res, matched = self.model_residuals, self.model_matched
        spacing = self.blob_spacing[self.t2b[matched]]
        sigma = 1.4826 * median(res)
        bad = (res > 4.0 * sigma) & (res > 0.5 * spacing)
        n_removed = 0
        for ti in matched[bad]:
            if self.origin[ti] != 0:
                self._unpair(ti)
                n_removed += 1

        if n_removed > 0:
            if not self.fit_global_model():
                return n_removed, 0
            res, matched = self.model_residuals, self.model_matched
            spacing = self.blob_spacing[self.t2b[matched]]

        # (b) pick up free target points through the model, but only if the
        # model describes the matched points well in the first place.
        #
        # A polynomial that is extrapolated may be off by a whole dot spacing
        # while fitting the matched points perfectly, and pairing a point
        # with the neighbour of its blob is far worse than leaving it
        # unpaired. Two conditions are therefore set. First, the point has to
        # be among the matched points and not beyond them, which is measured
        # by its leverage (see polynomial_projection.leverage). Second, the
        # model of one order lower has to send it to the very same blob.
        # Parts of the target that are beyond the matched points are thus
        # left for the user to mark a point on.
        n_added = 0
        if (self.model_low is not None and
                self.model_rms < 0.5 * self.global_tol * median(spacing)):
            free = (self.t2b == -1).nonzero()[0]
            proposals = []
            if len(free) > 0:
                d, bi = self.blob_tree.query(self.model(self.targets[free]))
                err = d / self.blob_spacing[bi]
                d_low, bi_low = self.blob_tree.query(
                                        self.model_low(self.targets[free]))
                agree = (bi_low == bi) & (d_low < 0.5 * self.blob_spacing[bi])
                inside = (self.model.leverage(self.targets[free]) <
                          self.max_leverage)
                for i in ((err < self.global_tol) & agree & inside).nonzero()[0]:
                    if self.b2t[bi[i]] != -1:
                        continue
                    if (int(free[i]), int(bi[i])) in self.rejected:
                        continue
                    proposals.append((err[i], free[i], bi[i]))

            for err, ti, bi in sorted(proposals):
                if self.b2t[bi] == -1 and self.t2b[ti] == -1:
                    self._pair(ti, bi, 2)
                    n_added += 1

        return n_removed, n_added



    def _run(self, max_iterations=25):
        '''alternates the local growth and the global refit until settled'''
        for it in range(max_iterations):
            self.grow_local()
            n_removed_local = self.remove_local_outliers()
            n_removed, n_global = self.refine_global()
            if n_removed_local == 0 and n_removed == 0 and n_global == 0:
                break



    def pair_points(self, verify=True):
        '''
        Does the actual pairing of blobs and target points. The paired points
        are stored in the attribute self.point_pairs, as rows of
        [x_image, y_image, x_lab, y_lab, z_lab].

        If verify is True, the marked points are afterwards checked against
        each other (see verify_marked_points).
        '''
        self.messages = []
        self.conflicts = []
        n_seeds = self.set_marked_points()

        if n_seeds < 2:
            self._say('At least 2 valid marked points are needed to start '
                      'the matching, and 3 that are not on one line are '
                      'recommended (%d found).' % n_seeds)
            self._collect()
            return self.point_pairs

        self._run()
        self._collect()
        self._report()

        if verify and n_seeds > 2:
            self.verify_marked_points()

        return self.point_pairs



    def verify_marked_points(self):
        '''
        Checks the marked points against each other. A marked point with
        wrong lab coordinates (e.g. those of the neighbouring dot) starts a
        growth that is self consistent, so it cannot be spotted by looking at
        the errors of the result. What can be done is to leave each
        marked point out in turn, repeat the matching with the others, and
        see which blob they assign to the target point of the one left out.
        If it is not the blob that the user had marked, the marked points
        contradict each other. Then, either the point that was left out
        or one of those from which the growth reached it is wrong.

        Points that the others cannot reach (e.g. the only marked point on a
        level of the target) can obviously not be checked like this.

        The results are stored in self.conflicts as a list of
        (marked point number, lab coordinates, marked (x,y), implied (x,y)),
        and in self.unverified as a list of marked point numbers.
        '''
        state = (self.t2b.copy(), self.b2t.copy(), self.origin.copy(),
                 set(self.rejected), list(self.seed_targets),
                 list(self.seed_numbers), self.messages, self.verbose)
        seeds = list(zip(self.seed_numbers, self.seed_targets,
                         self.t2b[self.seed_targets]))

        self.verbose = False
        conflicts, unverified = [], []
        for number, ti, bi in seeds:
            self.messages = []
            self.set_marked_points(leave_out=number)
            if len(self.seed_targets) >= 2:
                self._run()
            if self.t2b[ti] == -1 and self.b2t[bi] == -1:
                unverified.append(number)
            elif self.t2b[ti] != bi:
                implied = None
                if self.t2b[ti] != -1:
                    implied = tuple(self.blobs[self.t2b[ti]])
                conflicts.append((number, tuple(self.targets[ti]),
                                  tuple(self.blobs[bi]), implied))

        # back to the results obtained with all the marked points
        (self.t2b, self.b2t, self.origin, self.rejected, self.seed_targets,
         self.seed_numbers, self.messages, self.verbose) = state
        self.fit_global_model()
        self.conflicts, self.unverified = conflicts, unverified

        if len(conflicts) == 0:
            self._say('The marked points are consistent with each other.')
        else:
            self._say('\nWARNING: the marked points do not agree with each '
                      'other, so some of the pairs may be WRONG. Check the '
                      'lab coordinates of the following marked points:')
            for number, X, x, implied in conflicts:
                msg = ('  marked point %d, lab: %s, image: %s'
                       % (number, _fmt(X), _fmt(x)))
                if implied is not None:
                    msg += ('; the other marked points place it at %s'
                            % _fmt(implied))
                else:
                    msg += ('; the other marked points assign this blob to '
                            'a different target point')
                self._say(msg)
            self._say('(A single wrong point makes the points that it '
                      'contradicts be listed along with it. If they are all '
                      'correct, mark another point next to each of them and '
                      'match again.)')

        if len(unverified) > 0:
            self._say('Marked points that were not cross checked, since the '
                      'matching cannot reach them without them: %s.'
                      % ', '.join(['%d' % n for n in unverified]))

        return conflicts



    def _collect(self):
        matched = (self.t2b != -1).nonzero()[0]
        self.point_pairs = [hstack([self.blobs[self.t2b[ti]],
                                    self.targets[ti]]) for ti in matched]
        self.unmatched_targets = self.targets[self.t2b == -1]
        self.unmatched_blobs = self.blobs[self.b2t == -1]



    def _report(self):
        '''Summarizes the results and points out what is suspicious.'''
        n = int((self.t2b != -1).sum())
        self._say('\nMatched %d out of %d target points (%d blobs); %d by '
                  'local growth and %d through the global model.'
                  % (n, len(self.targets), len(self.blobs),
                     int(((self.origin == 1) | (self.origin == 3)).sum()),
                     int((self.origin == 2).sum())))

        if self.fit_global_model():
            self._say('Global polynomial model (order %d): RMS error of '
                      '%.2f pixels.' % (self.model.order, self.model_rms))

            # marked points that disagree with everything else
            res = dict(zip(self.model_matched, self.model_residuals))
            sigma = 1.4826 * median(self.model_residuals)
            for ti in self.seed_targets:
                bi = self.t2b[ti]
                if (res[ti] > 4.0 * sigma and
                        res[ti] > 0.5 * self.blob_spacing[bi]):
                    self._say('Warning: the marked point with lab coordinates '
                              '%s is %.1f pixels off the model fitted to all '
                              'the matched points. Check it.'
                              % (_fmt(self.targets[ti]), res[ti]))
        else:
            self._say('Note: the matched points are too few or all lie on one '
                      'plane, so the result could not be cross checked with '
                      'a global model. If the target has other levels, mark '
                      'one point on each of them.')

        # marked points from which nothing grew
        if n > len(self.seed_targets):
            tree = cKDTree(self.targets[self.t2b != -1])
            for ti in self.seed_targets:
                d = tree.query(self.targets[ti], k=2)[0][1]
                if d > self.max_reach * self.target_spacing[ti]:
                    self._say('Warning: nothing was matched around the '
                              'marked point with lab coordinates %s. Either '
                              'its coordinates are wrong, or the image is '
                              'too distorted to grow from it alone; in that '
                              'case mark two more dots next to it (not on '
                              'one line with it).'
                              % _fmt(self.targets[ti]))

        n_free = len(self.unmatched_targets)
        if len(self.unmatched_blobs) == 0:
            self._say('All the blobs were paired.')
        elif n_free > 0:
            self._say('%d target points and %d blobs were left unpaired. If '
                      'these blobs are dots of the target, mark one of them '
                      'and match again.' % (n_free, len(self.unmatched_blobs)))



    def save_results(self, fname):
        '''
        Saves the pairs of blobs and target points in a given file name.
        '''
        fmt = ['%.2f', '%.2f', '%.2f', '%.2f', '%.2f']
        savetxt(fname, self.point_pairs, fmt=fmt, delimiter='\t')
