# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 12:15:21 2019

@author: ron


calibration module for Extended Zolof camera instances. This is used to 
obtain the [A], [B], and O Extended Zolof model parameters, and optionally
the [C] parameters of the variable-origin extension (see
calibrate_variable_origin.py).

"""

from numpy import array, dot
from numpy.linalg import lstsq, norm
from scipy.optimize import minimize
from numpy import median, abs as npabs
import warnings

from myptv.extendedZolof.camera import camera_extendedZolof
from myptv.utils import line, get_nearest_line_crossing

from pandas import read_csv



class calibrate_extendedZolof(camera_extendedZolof):
    '''
    This object is used to calibrate cameras against a given list
    of lab and camera point coordinates. 
    '''
    
    
    def __init__(self, camera, x_list, X_list, quadratic=False):
        '''
        Given a list of 2D points, x=(x,y), and a list of 3D point X=(X,Y,Z), 
        we assume that given a point X, we can compute x by a polynomial 
        of degree 3, as - 
        
        x = A0 + A1*X + A2*Y + A3*Z +
            A4*X^2 + A5*Y^2 + A6*Z^2 + A7*XY + A8*YZ + A9*ZX + A108XYZ
            A11*XY^2 + A12*XZ^2 + A13*YX^2 + A14*YZ^2 + A15*ZX^2 + A16*ZY^2 +
            A17*X^3 + A18*Y^3
            
        if quadratic==True, then only the quadratic terms are used. This is 
        used in the initial calibration.
        '''
        self.cam = camera
        self.A = [[0.0 for i in range(19)] for j in [0,1]]
        self.B = [[0.0 for i in range(10)] for j in [0, 1, 2]]
        self.C = [[0.0 for i in range(10)] for j in [0, 1, 2]]
        self.x_list = x_list
        self.X_list = X_list
        self.quadratic = quadratic
        
        
        
    
    def calibrate(self, variable_origin=False, thresh_pix=5.0, c_order=3,
                  origin_model='free'):
        '''
        Given a list of points, x and X, this function attempts to determine
        the A, B (and O) coefficients, or A, B, C (and O) coefficients if
        variable_origin==True.
        
        variable_origin - if False (default), fits the original model where
                          all epipolar lines pass through a single fixed
                          camera center, self.O. If True, additionally fits
                          [C] coefficients so that the origin of each
                          epipolar line is itself a polynomial function of
                          the pixel, O(x) = [C] G(x) (see
                          calibrate_variable_origin.py). self.O is still
                          computed either way: with variable_origin==True it
                          is used only as an initial estimate / sign-
                          convention anchor for the ray fit.
        
        thresh_pix       - reprojection-residual threshold (pixels) used for
                          outlier rejection in both the camera-center
                          estimate (step 2) and, if variable_origin==True,
                          the per-pixel origin/direction fit.

        c_order           - only used if variable_origin==True. Polynomial
                          order for [C] (1=linear/affine, 2=quadratic,
                          3=full cubic, default, matching [B]). Fewer terms
                          means less flexibility - and less capacity to
                          overfit - for the per-pixel origin specifically.
                          See calibrate_variable_origin.fit_variable_origin_model
                          and .cross_validate_origin_models for more detail
                          and a way to check generalization on held-out
                          points.

        origin_model      - only used if variable_origin==True. Either
                          'free' (default; O(x) floats freely in 3D) or
                          'plane' (O(x) is constrained to a fixed plane
                          through the approximate common camera center,
                          facing the calibration volume). 'plane' uses 2
                          rather than 3 coefficient sets per basis term
                          and removes the along-ray gauge degeneracy of
                          the 'free' model, so it is typically better
                          conditioned and less prone to overfitting.
        '''
        # 1) finding the A coefficients - 
        #if self.quadratic==False:
        #    XColumns = [self.cam.get_XCol(Xi) for Xi in self.X_list]
        
        #elif self.quadratic==True:
        #    XColumns = []
        #    for Xi in self.X_list:
        #        Xcol_i = self.cam.get_XCol(Xi)
                # (Here, -9 is for quadratic, and -15 is linear)
        #        for i in range(-9,0): Xcol_i[i] = 0  
        #        XColumns.append(Xcol_i)
        
        #res = lstsq(XColumns, self.x_list, rcond=None)
        #self.A = res[0]
        
        self.A, A_mask = fit_A_robust(
            self.cam, self.X_list, self.x_list,
            quadratic=self.quadratic,
            sigma_thresh=3.5
            )
        
        self.A_mask = A_mask  # keep for diagnostics/plotting
        self.x_inlayers = array(self.x_list)[A_mask]
        self.X_inlayers = array(self.X_list)[A_mask]

        # IMPORTANT: update the camera's forward model *now*, not at the
        # end. Steps 2 (camera-center estimate) and the variable-origin fit
        # below both call self.cam.projection(), which uses self.cam.A
        # internally - if we don't set it here they would use a stale A
        # (whatever the camera object held before this calibration call).
        self.cam.A = self.A
        
        # 2) finding the best camera center -
        #line_list = []
        #for i in range(0, len(self.X_list)):
        #    O, e = self.get_ray_from_x(self.x_list[i], X0=self.X_list[i])
        #    line_list.append(line(O, e)) 
        #self.O = get_nearest_line_crossing(line_list)
        
        from myptv.extendedZolof.calibrate_step2_improved import step2_estimate_camera_center
        self.O = step2_estimate_camera_center(
            self.cam.projection, #self.x_list, self.X_list,
            self.x_inlayers, self.X_inlayers,
            thresh_pix=thresh_pix,
            refine=True
        )
        self.cam.O = self.O

        if not variable_origin:
            # 3) finding the unit vector for each X (fixed-origin model) -
            r_list = []
            for Xi in self.X_inlayers: #self.X_list:
                r = (Xi - self.O)/norm(Xi - self.O)
                r_list.append(r)

            # 4) finding the B coefficients -
            xColumns = [self.cam.get_xCol(xi) for xi in self.x_inlayers] #self.x_list]
            res = lstsq(xColumns, r_list, rcond=None)
            self.B = res[0]

            self.cam.B = self.B
            self.cam.variable_origin = False

        else:
            # 3'+4') fitting a per-pixel origin O(x)=[C]G(x) together with
            # the direction B(x)=[B]G(x), instead of assuming a single
            # fixed O for every epipolar line.
            from myptv.extendedZolof.calibrate_variable_origin import \
                fit_variable_origin_model

            (self.C, self.B, self.ray_mask, self.ray_resid,
             self.px_center, self.px_scale, self.plane) = \
                fit_variable_origin_model(
                    self.cam, self.X_inlayers, self.x_inlayers,
                    O_init=self.O,
                    thresh_pix=thresh_pix,
                    c_order=c_order,
                    origin_model=origin_model
                )

            self.cam.C = self.C
            self.cam.B = self.B
            self.cam.px_center = self.px_center
            self.cam.px_scale = self.px_scale
            self.cam.origin_mode = origin_model

            if self.plane is not None:
                # the plane passes through O0; keep cam.O consistent with
                # it, since get_ray() builds O(x) relative to cam.O.
                self.O = self.plane['O0']
                self.cam.O = self.O
                self.cam.plane_u = self.plane['u']
                self.cam.plane_v = self.plane['v']
                self.cam.plane_n = self.plane['n']

            self.cam.variable_origin = True
        
        
            
    def get_ray_from_x(self, x, X0=None):
        '''
        Given a point in 2D image space, this function returns a line in 3D
        that passes through this point. The line is represented with six 
        parameters: one point in 3D, O, and one unit vector in 3D, e.
        '''
        
        func = lambda X: sum((array(self.projection(X)) - array(x))**2)
        
        if X0 is None:
            X0 = array([0,0,0])
        
        X02 = array(X0) + array([1,1,1])
            
        O = minimize(func, X0).x
        dX = minimize(func, X02).x
        e = (O-dX)/sum((O-dX)**2)**0.5
        
        return O, e
        
    
    
    def mean_squared_err(self):
        '''
        Calculates and returns the mean square of the deviations in camera
        space (forward, 3D -> 2D projection error).
        '''
        errorsSquard = []
        
        try:
            x_lst = array(self.x_inlayers)
        except:
            x_lst = array(self.x_list)
        
        try:
            X_lst = array(self.X_inlayers)
        except:
            X_lst = array(self.X_list)
                          

        for i in range(len(X_lst)):
            xProj = dot(self.get_XCol(X_lst[i]), self.cam.A)
            errorsSquard.append( norm(array(xProj)-array(x_lst[i]))**2 )
        
        return (sum(errorsSquard)/len(errorsSquard))**0.5


    def mean_ray_err(self):
        '''
        Calculates and returns the mean point-to-line distance (lab units)
        between the calibration points and their fitted epipolar lines
        (backward, 2D -> 3D ray error).

        For variable_origin fits this uses the residuals already computed
        during fit_variable_origin_model (self.ray_resid). For fixed-origin
        fits it recomputes the distance from each inlier point to its ray
        using self.cam.get_ray().
        '''
        if getattr(self.cam, 'variable_origin', False) and \
                hasattr(self, 'ray_resid'):
            resid = array(self.ray_resid)
            return (sum(resid**2)/len(resid))**0.5

        dists = []
        for i in range(len(self.X_inlayers)):
            Xi = array(self.X_inlayers[i])
            xi = self.x_inlayers[i]
            O, e = self.cam.get_ray(xi[0], xi[1])
            diff = Xi - O
            perp = diff - dot(diff, e)*e
            dists.append(norm(perp)**2)

        return (sum(dists)/len(dists))**0.5
        
        
        
    def plot_err_distribution(self, ax = None):
        import matplotlib.pyplot as plt
        from numpy import sum as npsum
        
        if ax == None:
            fig, ax = plt.subplots()
        
        #imc = array(self.x_list)
        #z_lst = array([self.cam.projection(x) for x in self.X_list])
        try:
            imc = array(self.x_inlayers)
        except:
            imc = array(self.x_list)
        
        try:
            z_lst = array([self.cam.projection(x) for x in self.X_inlayers])
        except:
            z_lst = array([self.cam.projection(x) for x in self.X_list])
        
        err = npsum((imc-z_lst)**2, axis=1)**0.5
        
        h = ax.hist( err, bins='auto')
        ax.set_xlabel('Camera projection err [px]')
        ax.set_ylabel('Counts')
        
        
        
        
    def plot_proj(self, ax = None):
        import matplotlib.pyplot as plt
        
        if ax == None:
            fig, ax = plt.subplots()
        
        try:
            imc = array(self.x_inlayers)
        except:
            imc = array(self.x_list)
        ax.plot(imc[:,0], imc[:,1], 'ob')
        for i in range(imc.shape[0]):
            ax.text(imc[i,0], imc[i,1], '%d'%i, color = 'b')
        
        try:
            z_lst = array([self.cam.projection(x) for x in self.X_inlayers])
        except:
            z_lst = array([self.cam.projection(x) for x in self.X_list])
        ax.plot( z_lst[:,0], z_lst[:,1], 'xr' )
        for i in range(z_lst.shape[0]):
            ax.text(z_lst[i,0], z_lst[i,1], '%d'%i, color = 'r')
            
        ax.set_aspect('equal')
        
        
        
        
        
        





 
 
def fit_A_robust(cam, X_list, x_list, quadratic=False,
                  sigma_thresh=3.0, max_iterations=3,
                  max_reject_fraction=0.25):
    '''
    Robustly fits the A (3D -> 2D) polynomial coefficients with iterative
    outlier rejection.
 
    sigma_thresh        : points with residual > median + sigma_thresh*MAD
                           are rejected. ~3.0 is a reasonable starting
                           point; lower it to reject more aggressively.
    max_iterations       : safety cap on the number of refit iterations.
    max_reject_fraction  : if more than this fraction of points would be
                           rejected, stop early and warn -- usually means
                           thresh is too strict or something else is wrong
                           (e.g. bad initial correspondences), not just
                           "a few outliers".
 
    Returns (A, mask) where mask is a boolean array marking which points
    were kept in the final fit.
    '''
    # --- build the design matrix, same as the original code ---
    if quadratic == False:
        XColumns = [cam.get_XCol(Xi) for Xi in X_list]
    else:
        XColumns = []
        for Xi in X_list:
            Xcol_i = cam.get_XCol(Xi)
            # (-9 is for quadratic, -15 is linear)
            for i in range(-9, 0):
                Xcol_i[i] = 0
            XColumns.append(Xcol_i)
 
    XColumns = array(XColumns)
    xTarget = array(x_list)
    N = len(xTarget)
 
    mask = array([True] * N)
    A = None
 
    for iteration in range(max_iterations):
 
        res = lstsq(XColumns[mask], xTarget[mask], rcond=None)
        A = res[0]
 
        # residuals for ALL points, so rejected ones can be reconsidered
        residuals = norm(XColumns.dot(A) - xTarget, axis=1)
 
        med = median(residuals[mask])
        mad = median(npabs(residuals[mask] - med)) * 1.4826  # -> ~std
 
        if mad == 0:
            # residuals are (numerically) identical; nothing more to reject
            break
 
        new_mask = residuals < (med + sigma_thresh * mad)
 
        if new_mask.sum() < (1 - max_reject_fraction) * N:
            warnings.warn(
                'Outlier rejection would remove more than '
                f'{max_reject_fraction*100:.0f}% of points; stopping '
                'early and keeping the previous iteration\'s fit. '
                'Consider raising sigma_thresh or checking your '
                'correspondences.')
            break
 
        if (new_mask == mask).all():
            break  # converged: no change in kept/rejected points
 
        mask = new_mask
 
    n_rejected = N - mask.sum()
    if n_rejected > 0:
        print(f'fit_A_robust: rejected {n_rejected}/{N} points '
              f'({100*n_rejected/N:.1f}%) as outliers.')
 
    return A, mask




        
        
        
        
        
        

class calibrate_with_particles_EZ(object):
    '''
    A class used to refine the calibration using particles data. In short,
    after the primary clibration is done, matching and tracking can be used
    to obtain trajectories from the experimental data. Here, we can leverage
    the trajectories obtained to minimize further the calibration error.
    The assumption is that longer trajectories are considered more reliable 
    as compared to shorter trajectories, so we use only "long" trajectories 
    in this process.
    '''
    
    def __init__(self, traj_filename, camera, cam_number, blobs_fname, 
                 min_traj_len = 15, max_point_number = 1000):
        '''
        input -
        
        traj_filename - the name of the file containing the trajectories from
                        which the calibration points are taken.
        
        camera - an instance of the camera we wish to try and re-calibrate
        
        cam_number - int, >= 1; the number (index) of the camera to be
                     calibrated. For example, if this is 1 then it takes
                     the blobs from the 4th column of the trajectory file.
        
        blobs_fname - The name of the file that contains the segmented 
                      particles' data.
        
        min_traj_len - only trajectories longer then this number will be used
                       in the calibration
                       
        max_point_number - the maximum number of points that shall be taken
                           to re-calibrate the camera. Note that too many 
                           points might lead to long calculation times.
        '''
        
        self.traj_fname = traj_filename
        self.camera = camera
        self.cam_number = cam_number
        self.blobs_fname = blobs_fname
        self.min_traj_len = min_traj_len
        self.max_point_number = max_point_number
        
        # gathering points from trajectory file and matching them to blobs
        self.fetch_points()
        
        
    def fetch_points(self):
        '''
        This will fetch the calibration points from the trajectory file
        '''
        print('fetching data from trajectories and blobs file...')
        
        # load the trajectories and sort them in a dictionary according to id
        fltr = lambda k,g: k!=-1 and len(g)>=self.min_traj_len
        tr_data = read_csv(self.traj_fname, delimiter='\t', header=None)
        self.trajs = dict([(k,array(g)) for k, g in tr_data.groupby(by=0) 
                                                                if fltr(k,g)])

        # load blobs an arrange in lists according to their frame
        blob_data = read_csv(self.blobs_fname, delimiter='\t', header=None)
        self.blobs = dict([(k,array(g)) for k, g in blob_data.groupby(by=5)])

        # the index of the blobs column in the trajectories file
        ind = self.cam_number + 3
        
        skip = int(self.min_traj_len/4)
        
        # extract the data from blobs and trajs dictionaries
        all_valid_points = []
        for k in self.trajs.keys():
            if len(self.trajs[k])>=self.min_traj_len:
                for p in self.trajs[k][::skip]:     # <--- taking only once some 
                    blob_num = int(p[ind])     #      data points in each traj
                    if blob_num==-1:continue
                    frame = p[-1]
                    err = p[-2]
                    blob_coords = self.blobs[frame][blob_num][:2]
                    
                    point = p[1:4]
                    
                    all_valid_points.append( (err,(point, blob_coords[::-1])) )
        
        
        if len(all_valid_points)==0:
            msg1 = 'No trajectories were found. min_traj_len may be too high.'
            raise ValueError(msg1)
        
        
        # sort the points according to their triangulation error
        all_valid_points = sorted(all_valid_points, key=lambda x: x[0])
        
        # get the best N points into a final list
        self.cal_points = [p[1] 
                           for p in all_valid_points[:self.max_point_number]]


    def get_particle_disparity(self):
        '''
        Returns a list with the discrepancies between the segmented blob 
        positions and the projections of their 3D positions onto the camera.
        '''
        cam = self.camera
        disparities = [cam.projection(p[0]) - p[1] for p in self.cal_points]
        return disparities
        

    def get_calibrate_instance(self):
        
        # initiating a calibrate object using this data
        self.cal = calibrate_extendedZolof(self.camera.camera, 
                                           [p[1] for p in self.cal_points], 
                                           [p[0] for p in self.cal_points])
        return self.cal
        
        
