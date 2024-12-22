# -*- coding: utf-8 -*-
"""

@author: agambino

An implementation of the trajectory segments stitching algorithm, for fibers. We are 
looking to connect together trajectories that, for some reason, were broken 
during the tracking process, and to connect them together. Once the trajectories are
stitched, the fiber orientations are connected with a polynomial fitting and then 
smoothed.

The algorithm is defined by Xu (2008): Haitao Xu, Tracking Lagrangian 
trajectories in position–velocity space, Meas. Sci. Technol. 19 (2008)
075105 (10pp).
"""

from numpy import array, gradient, dot, savetxt, linalg, vstack, lexsort, unique, where, zeros, diff, square, sum, zeros_like, isnan #agambino
from numpy import append as NPappend
from myptv.utils import fit_polynomial
from myptv.traj_smoothing_mod import smooth_traj_poly
from itertools import groupby


class traj_ori_stitching(object):    
    '''
    Finds a list of trajectory pairs to connect to each
    other based on the algorithm of Xu(2008). For now we are not
    using trajectory acceleration information for the stitching.
    In addition, we also interpolate the "missing" samples at the connected
    interval by using a cubic polynomial that is fitted with the last and 
    first two data points of the connected trajectories.
    '''
    
    def __init__(self, traj_list, Ts, dm, polyorder, window):
        '''
        inputs -
        traj_list - list, holding the trajectories. The algorithm relies
                    on the trajectories' velocity, so we input here
                    the results of smoothed trajectories. Use the same
                    format as the one saved in 
                    traj_smooth_mod.smooth_trajectories.save_results().
        Ts - float, the maximum allowable time separation allowed 
             for stitching
        dm - maximum distance (in r,v*dt space) for connection
        '''
        self.traj_list = traj_list
        self.Ts = Ts
        self.dm = dm
        self.polyorder = polyorder
        self.window = window

        # ===============================================================
        # for now we will not be using acceleration data in the stitching
        # but the framework is prepared through the parameter wa. For now
        # it is zero, so acceleration is not used.
        self.wa = 0.0
        # ===============================================================

    def get_traj(self, i):
        '''
        given a trajectory number, i, this will return the trajectory data, 
        with samples sorted according to time.
        
        input -
        i - int, trajectory number
        '''
        key = lambda s: s[-1]
        #tr = array(sorted(self.traj_list[self.traj_list[:,0] == i], key=key))
        tr = array(sorted(self.traj_list[self.index_id_hash[i],:], key=key))
        return tr


    def calc_d(self, tr_i, tr_j):
        '''
        Taking 2 trajectories with trajectory numbers i and j. Trajectory i
        ends at time t_ie and trajectory j starts at time t_js; it is 
        necessary that t_ie < t_js.
        
        We project the position and velocity of particle i onto time t_js 
        assuming constant velocity: 
            x_i(t_js) = x_i(t_ie) + u_i(t_ie) * (t_js - t_ie)
            u_i(t_js) = u_i(t_ie)  (+ possible future acceleration term)
            
        Then, we calculate the distance between particle j and the projection
        of particle i at the start time of particle j
            dij = [ (x_i(t_js) - x_j(t_js))^2 + (u_i(t_js) - u_j(t_js))^2 ]^0.5
        
        * if a trajectory has more than one sample, we estimate its velocity
          based on simple forward differences. 
          
        ** For now, we do not attempt to connect single samples. This is 
           something that could be added in the future.
        
        
        return - if dij<=dm and t_js - t_ie < self.Ts: returns dij. 
                 Else, returns -1.
        '''
        # note the end and starting times
        t_ie = tr_i[-1,-1]
        t_js = tr_j[0,-1]
        
        # note their time separation:
        dt = t_js - t_ie
        if dt > self.Ts or dt<0:
            return -1
        
        # taking the position velocity and acceleration of trajectory i
        x_ie = tr_i[-1,1:4]
        v_ie = tr_i[-1,4:7]
        a_ie = tr_i[-1,7:10]
        
        # projecting particle i to time t_js
        x_i_js = x_ie + dt * v_ie
        v_i_js = v_ie + self.wa * dt * a_ie
        
        # calculating the distance between trajectories at t_js
        dx = sum((x_i_js - tr_j[0,1:4])**2)**0.5
        dv = sum((v_i_js - tr_j[0,4:7])**2)**0.5
        
        # calculate d_ij; if it's small enough return it, else return -1.
        dij = (dx**2 + dv**2)**0.5
        
        if dij<=self.dm: return dij
        else: return -1

    
    
    def calc_dij(self):
        '''
        Calculates dij for all possible trajectories and lists them.
        '''
        traj_ids = list(set(self.traj_list[:,0]))
        if -1 in traj_ids:
            traj_ids.remove(-1)
        
        # make a dictionary whose keys are trajectory ids and values are
        # the trajectories themselves
# =============================================================================
#         Legacy code:
#
#         traj_dic = {}
#         for id_ in traj_ids:
#             tr = self.get_traj(id_)
#             v = tr[:,4:7]
#             if (v==0).all():
#                 x = tr[:,1:4]
#                 v = gradient(x, axis=0)
#                 a = gradient(v, axis=0)
#                 tr[:,4:7] = v
#                 tr[:,7:10] = a
#                 
#             traj_dic[id_] = tr 
# =============================================================================
        
        traj_dic = {}
        trajs = [array(sorted(list(g), key=lambda a: a[-1])) for k,g in 
                 groupby(self.traj_list, key=lambda x: x[0]) if k!=-1]
        for tr in trajs:
            id_ = tr[0,0]
            v = tr[:,4:7]
            #print(len(tr))
            if (v==0).all() and len(tr)>1: #agambino
                x = tr[:,1:4]
                v = gradient(x, axis=0)
                a = gradient(v, axis=0)
                tr[:,4:7] = v
                tr[:,7:10] = a
            traj_dic[id_] = tr 
        
        # make a dictionary whose keys are frame numbers and values are lists
        # of the trajectories ids that start at these times
        t0_dic = {}
        for id_ in traj_ids:
            t0_i = traj_dic[id_][0,-1]
            try:
                t0_dic[t0_i].append(id_)
            except:
                t0_dic[t0_i] = [id_]
            
        
        # make the list of d_ij
        self.dij_list = []
        count, N = 0, len(traj_ids)
        for id_i in traj_ids:
            print('', end='\r')
            print('calculating d_ij: %d / %d'%(count,N), end='\r') 
            
            traj_i = traj_dic[id_i]
            ts_i = traj_i[-1,-1]
            t0_to_search = [ts_i + i for i in range(1, 1 + self.Ts)]
            
            for t0 in t0_to_search:
                try:
                    for id_j in t0_dic[t0]:
                        traj_j = traj_dic[id_j]
                        dij = self.calc_d(traj_i, traj_j)
                        if dij != -1:
                            dt = traj_j[0,-1] - traj_i[-1,-1]
                            self.dij_list.append((id_i, id_j, dt, dij))
                except:
                    continue
                
            count += 1
        print('') 

        
    
    
    def find_best_stitch_candidates(self):
        '''
        After generating the dij list, if there are more than one possible
        connections to be made for any single trajectory, this function 
        discards the connections with the higher dij values.
        '''
        sort_d = lambda s: s[-1]
        dij_sorted = sorted(self.dij_list, key = sort_d)
        
        new_dij_list = []
        i_added = []
        j_added = []
        for dij in dij_sorted:
            i, j = dij[0], dij[1]
            
            if i in i_added or j in j_added:
                continue
            
            else:
                new_dij_list.append(dij)
                i_added.append(i)
                j_added.append(j)
                
        self.dij_list = new_dij_list
                
        
    def connect_traj(self):
        '''
        After finding the best candidates for stitching, this connects the
        trajectories. The connection is made by changing the trajectory number
        of the earlier trajectory (i) to that of the later trajectory (j).
        
        We also interpolate the "missing" samples at the connected interval
        by using a cubic polynomial that is fitted with the last and first 
        two data points of the connected trajectories.
        '''
        connected_traj_i = []
        connected_traj_j = []
        
        count = 1
        for i,j,dt,dij in self.dij_list:
            print('', end='\r')
            print(' %d / %d'%(count, len(self.dij_list)), end='\r')
            
            tr_i = self.get_traj(i)
            tr_j = self.get_traj(j)
            
            #if j was previously connected:
            while len(tr_j)==0:
                j = connected_traj_j[connected_traj_i.index(j)]
                tr_j = self.get_traj(j)
            
            # find the times that need to be filled in
            t_ie = tr_i[-1,-1]
            
            # get the trajectory positions needed for the polynomial fitting
            x_i = tr_i[-2:,1:4]   # <-- last two data points of i
            x_j = tr_j[:2, 1:4]   # <-- first two data points of j
            p_i = tr_i[-2:,10:13]   #agambino
            p_j = tr_j[:2, 10:13]   #agambino

            
            tm_fitting = [-1.0, 0.0, dt, dt+1.0]  # <-- time for fitting
            #print('tm_fitting: ', tm_fitting)
            
            # the time needed to interpolate
            tm_interp = [float(i) for i in range(1,int(dt))]
            #print('ecco: ', len(tm_interp))
            if len(tm_interp)>1: # if len(tm_interp)=1, orientation stitching becomes bad, so we avoid it.
                # in each direction fit a polynomial and fill missing samples
                interp_samples_j = [[j,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,tm_+t_ie] #agambino
                                  for tm_ in tm_interp]
                interp_p_plus_j = [[0,0,0,tm_+t_ie] #agambino
                                  for tm_ in tm_interp]
                interp_p_minus_j = [[0,0,0,tm_+t_ie] #agambino
                                  for tm_ in tm_interp]

                for k in range(3):
                    #polyfit positions
                    x_ = list(x_i[:,k]) + list(x_j[:,k])
                    poly_coeffs_x = fit_polynomial(tm_fitting, x_, 3)

                    #polyfit orientations 
                    p_plus = list(p_i[:,k]) + list(p_j[:,k]) #agambino      
                    poly_coeffs_p_plus = fit_polynomial(tm_fitting, p_plus, 3) #agambino
                    p_minus = list(p_i[:,k]) + list(-p_j[:,k]) #agambino      
                    poly_coeffs_p_minus = fit_polynomial(tm_fitting, p_minus, 3) #agambino
                    #print('ecco p_plus: ', p_plus)
                    #print('ecco p_minus: ', p_minus)
                    #print('k: ', k)
                    #print('ecco len(tm_interp): ', len(tm_interp))

                    for e, tm_ in enumerate(tm_interp):
                        # interpolate position
                        tm_vect = [tm_**3, tm_**2, tm_, 1.0]
                        x_interp = dot(poly_coeffs_x, tm_vect)
                        p_plus_interp = dot(poly_coeffs_p_plus, tm_vect) #agambino
                        p_minus_interp = dot(poly_coeffs_p_minus, tm_vect) #agambino
                        interp_samples_j[e][1+k] = x_interp
                        interp_p_plus_j[e][k] = p_plus_interp #agambino
                        interp_p_minus_j[e][k] = p_minus_interp #agambino

                # to solve the sign ambiguity, we fit the polynomial over two set of orientations from different segments
                p_incr_plus = sum(sum(square(diff(array(interp_p_plus_j)[:, :3], axis=0)), axis=1))
                p_incr_minus = sum(sum(square(diff(array(interp_p_minus_j)[:, :3], axis=0)), axis=1))
                #print('ecco p_incr_plus: ', p_incr_plus)
                #print('ecco p_incr_minus: ', p_incr_minus)
                #print('ecco interp_p_plus_j: ', interp_p_plus_j)
                #print('ecco interp_p_minus_j: ', interp_p_minus_j)

                if p_incr_plus >= p_incr_minus:
                    for a in range(len(interp_samples_j)):  #agambino: normalization of new interpolated orientations
                        interp_samples_j[a][10:13] = interp_p_minus_j[a][0:3]
                        vnorm = linalg.norm(interp_samples_j[a][10:13])
                        interp_samples_j[a][10:13] = interp_samples_j[a][10:13] / vnorm
                elif p_incr_plus < p_incr_minus:
                    for a in range(len(interp_samples_j)):  #agambino: normalization of new interpolated orientations
                        interp_samples_j[a][10:13] = interp_p_plus_j[a][0:3]
                        vnorm = linalg.norm(interp_samples_j[a][10:13])
                        interp_samples_j[a][10:13] = interp_samples_j[a][10:13] / vnorm

            else:
                continue
            '''
            if len(tm_interp)>0:
                # in each direction fit a polynomial and fill missing samples
                interp_samples_j = [[j,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,tm_+t_ie] #agambino
                                  for tm_ in tm_interp]
                for k in range(3):
                    #polyfit positions
                    x_ = list(x_i[:,k]) + list(x_j[:,k])
                    poly_coeffs_x = fit_polynomial(tm_fitting, x_, 3)

                    #polyfit orientations 
                    p_ = list(p_i[:,k]) + list(p_j[:,k]) #agambino      
                    poly_coeffs_p = fit_polynomial(tm_fitting, p_, 3) #agambino


                    for e, tm_ in enumerate(tm_interp):
                        # interpolate position
                        tm_vect = [tm_**3, tm_**2, tm_, 1.0]
                        x_interp = dot(poly_coeffs_x, tm_vect)
                        p_interp = dot(poly_coeffs_p, tm_vect) #agambino
                        interp_samples_j[e][1+k] = x_interp
                        interp_samples_j[e][10+k] = p_interp #agambino
                
                
                for a in range(len(interp_samples_j)): #agambino: normalization of new interpolated orientations
                    vnorm = linalg.norm(interp_samples_j[a][10:13])
                    interp_samples_j[a][10:13] = interp_samples_j[a][10:13] / vnorm
                    #print('ecco interp samples: ', interp_samples_j[a][10:13])
                '''
            # add the interpolated samples to the traj_list and update
            # index hash table
            Ntl = len(self.traj_list)
            Ninterp = len(interp_samples_j)
            self.traj_list = NPappend(self.traj_list, 
                                          interp_samples_j, axis=0)
            self.index_id_hash[j] += list(range(Ntl, Ninterp + Ntl))
            
            # change the traj_number of traj i to j
            for k in range(len(self.traj_list)):
                if self.traj_list[k][0] == i:
                    self.traj_list[k][0] = j
            i_indexes = self.index_id_hash[i]
            self.index_id_hash[j] += i_indexes
            self.index_id_hash[i] = []
            
            connected_traj_i.append(i)
            connected_traj_j.append(j)
            
            count += 1

            
        # finish by calculating the velocity and acceleration of the 
        # the stitched trajectory by using the smoothing function
# =============================================================================
#         for id_ in list(set(self.traj_list[:,0])):
#             
#             # (we add the single samples at the end)
#             if id_==-1:
#                 continue
#             
#             traj = self.get_traj(id_)
#             
# =============================================================================

        self.new_traj_list = []
        trajs = [array(sorted(list(g), key=lambda a: a[-1])) for k,g in 
                 groupby(self.traj_list, key=lambda x: x[0]) if k!=-1]
        trajs = vstack(trajs)
        trajs = trajs[lexsort((trajs[:, 19], trajs[:, 0]))]
        
        trajs[:, 10:13] = self.traj_sign_switch(trajs[:, 0], trajs[:, 10:13])
        #self.new_traj_list = array(trajs)
        
        
        # smooth everything
        spacing = int(self.window/2)+1
        indunique = unique(trajs[:, 0])
        for index in indunique:
            traj = trajs[trajs[:, 0] == index, :]
            #mask = isnan(traj[:, 10:14]).any(axis=1)
            #traj = traj[~mask]            
            if len(traj)<=int(self.window):
                continue
            new_x, new_v, new_acc = smooth_traj_poly(traj[:,1:4].T, self.window, self.polyorder, 1) #agambino
            new_p, new_pdot, new_pdotdot = smooth_traj_poly(traj[:,10:13].T, self.window, self.polyorder, 1) #agambino
            new_x = array(new_x).T
            new_v = array(new_v).T
            new_acc = array(new_acc).T
            new_p = array(new_p).T
            new_pdot = array(new_pdot).T
            new_pdotdot = array(new_pdotdot).T
            #new_v[:spacing,:] = 0
            #new_v[-spacing:,:] = 0
            #new_acc[:spacing,:] = 0
            #new_acc[-spacing:,:] = 0
            #new_pdot[:spacing,:] = 0
            #new_pdot[-spacing:,:] = 0
            #new_pdotdot[:spacing,:] = 0
            #new_pdotdot[-spacing:,:] = 0
            #traj[:, 1:4] = new_x
            #traj[:, 4:7] = new_v
            #traj[:, 7:10] = new_acc
            #traj[:, 10:13] = new_p
            #traj[:, 13:16] = new_pdot
            #traj[:, 16:19] = new_pdotdot

            mask = (traj[spacing:-spacing, 4:7] == 0).all(axis=1)

            traj[spacing:-spacing][mask, 1:4] = new_x[spacing:-spacing][mask,:]
            traj[spacing:-spacing][mask, 4:7] = new_v[spacing:-spacing][mask,:]
            traj[spacing:-spacing][mask, 7:10] = new_acc[spacing:-spacing][mask,:]
            traj[spacing:-spacing][mask, 10:13] = new_p[spacing:-spacing][mask,:]
            traj[spacing:-spacing][mask, 13:16] = new_pdot[spacing:-spacing][mask,:]
            traj[spacing:-spacing][mask, 16:19] = new_pdotdot[spacing:-spacing][mask,:]
            trajs[trajs[:, 0] == index, :] = traj    

        self.new_traj_list = array(trajs)    
        
    def stitch_trajectories(self):
        '''
        This performs all the steps of the trajectory stitching process.
        '''
        whr = self.traj_list[:,0] != -1
        traj_ids = list(set(self.traj_list[whr,0]))
        ntr = len(traj_ids)
        print('starting at %d trajectories'%(ntr))
        nsmp = len(self.traj_list[whr])*1.0
        print('with %.1f samples per trajectory on average'%(nsmp/ntr),'\n')
        
        
        # making a hash table of id -> index list of trajectories
        self.index_id_hash = {}
        for i in range(len(self.traj_list)):
            id_ = self.traj_list[i][0]
            try:
                self.index_id_hash[id_].append(i)
            except:
                self.index_id_hash[id_] = [i]
        

        print('searching for candidates to connect:')
        self.calc_dij()
        print('fetching the best matches...')
        self.find_best_stitch_candidates()
        

        N = len(self.dij_list)
        print('found %d connections to be made,'%N)
        print('connecting...')
        ntraj0 = len(self.traj_list)
        self.connect_traj()
        print('finished connecting trajectories', '\n')
        
        print('interpolated %d new samples'%(len(self.traj_list) - ntraj0))
        traj_ids = list(set(self.traj_list[:,0]))
        ntr = len(traj_ids)
        print('finished with %d trajectories'%(ntr))
        whr = self.traj_list[:,0] != -1
        nsmp = len(self.traj_list[whr])*1.0
        print('at %.1f samples per trajectory on average'%(nsmp/ntr),'\n')    


            
    def traj_sign_switch(self, indd, ori): #agambino
        '''
        Function that restores orientation signs over a trajectory
        '''
        ori_straight = zeros_like(ori)
        indd_unique = unique(indd)
        for index in indd_unique:
            traj_indd = where(indd == index)[0]
            traj_ori = ori[traj_indd,:]
            traj_ori_straight = zeros_like(traj_ori)
            traj_ori_straight[0, :] = traj_ori[0, :]
            for i in range(1, traj_ori.shape[0]):
                dot_product = dot(traj_ori_straight[i-1,:], traj_ori[i,:])
            
                if dot_product < 0:
                    traj_ori_straight[i,:] = -traj_ori[i,:]
                else:
                    traj_ori_straight[i,:] = traj_ori[i,:]
            ori_straight[traj_indd,:] = traj_ori_straight
        return ori_straight

    def save_results(self, fname):
        '''
        saves the results on the hard drive as a text file with a given file
        name.
        '''
        fmt = ['%d', '%.3f', '%.3f', '%.3f', '%.3f', '%.3f', '%.3f', '%.3f', 
               '%.3f', '%.3f', '%.3f']
        savetxt(fname, self.new_traj_list, fmt=fmt, delimiter='\t')
            
            