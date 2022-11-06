# -*- coding: utf-8 -*-
"""
Created on Sat March 19 2022

@author: ron


Implementation of a method for smoothing trajectories.
The method uses polynomial fitting after Luethi et al 2005.

"""

from myptv.utils import fit_polynomial
from numpy import dot, array, savetxt






class smooth_trajectories(object):
    '''
    A class used to smooth trajectories in a list of trajectories. 
    During smoothing, we also calculate the velocity and acceleration
    of the trajectories.
    The input trajectory list structure is the same as the files produced by
    the classes in tracking_mod.py.
    
    Note - only trajectories whose length is larger than the window size
    will be smoothed and saved. Shorter trajectories are saved with zero
    velocity and accelerations.
    '''
    
    def __init__(self, traj_list, window, polyorder, repetitions=1):
        '''
        Input -
        traj_list - a nested list of NX5 elements. N is the number of samples,
                    then the following indices are trajectory number, 
                    x position, y position, z position, frame number.
        window - odd position integer, window size used in the polynomial
                 fitting
        polyorder - integer, degree of the polynomial used in the smoothing.
        repetitions - the number of times that the smoothing should be 
                      repeated for each trajectory.
        '''
        
        self.traj_list = traj_list
        self.window = window
        self.polyorder = polyorder
        self.repetitions = repetitions
        
        
    def smooth(self):
        '''
        Performs the smoothing and returns the results. 
        '''
        
        # organizing trajectories in a dictionary:
        traj_dic = {}
        zero_length_trajs = []
        for i in range(len(self.traj_list)):
            tr = self.traj_list[i]
            
            # for unconnected samples, put zero velocity and acceleration:
            if tr[0] == -1: 
                new_tr = [tr[0], tr[1], tr[2], tr[3], 
                          0.0, 0.0, 0.0, 0.0, 0.0, 0.0, tr[-1]]
                zero_length_trajs.append(new_tr)
            
            # from the samples, make a trajectory dictionary
            if tr[0] in traj_dic.keys():
                traj_dic[tr[0]].append(tr)
            else:
                traj_dic[tr[0]] = [tr]
        
        
        short_trajs = []
        smoothed_traj_list = []
        N = len(traj_dic.keys())
        count = 0
        for tr_num in traj_dic.keys():
            print('', end='\r')
            print(' progress: %.1f%%'%(count/N*100), end='\r') 
            
            # for too short trajectories mark zero velocity and acceleration:
            if len(traj_dic[tr_num]) < self.window:
                for i in range(len(traj_dic[tr_num])):
                    tr = traj_dic[tr_num][i]
                    new_tr = [tr[0], tr[1], tr[2], tr[3], 
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, tr[-1]]
                    short_trajs.append(new_tr)
                continue
            
            # sort samples according to time:
            traj = sorted(traj_dic[tr_num], key=lambda s: s[-1])
            
            # smoothing
            p, v, a = smooth_traj_poly(array(traj).T[1:4,:], 
                                             self.window, 
                                             self.polyorder,
                                             repetitions=self.repetitions)
            
            # setting a new trajectories
            new_traj = []
            for i in range(len(traj_dic[tr_num])):
                new_traj.append([])
                new_traj[-1].append(traj[i][0])
                new_traj[-1].append(p[0][i])
                new_traj[-1].append(p[1][i])
                new_traj[-1].append(p[2][i])
                new_traj[-1].append(v[0][i])
                new_traj[-1].append(v[1][i])
                new_traj[-1].append(v[2][i])
                new_traj[-1].append(a[0][i])
                new_traj[-1].append(a[1][i])
                new_traj[-1].append(a[2][i])
                new_traj[-1].append(traj[i][-1])
                
            smoothed_traj_list += new_traj
            count+=1
            
        print('')
        print('smoothed samples: %d'%(len(smoothed_traj_list)))
        print('too short to smooth: %d'%(len(short_trajs)))
        print('single samples: %d'%(len(zero_length_trajs)))
        
        smoothed_traj_list += short_trajs    
        smoothed_traj_list += zero_length_trajs
        self.smoothed_trajs = smoothed_traj_list
        
        
    def save_results(self, fname):
        '''
        Will save the smoothed trajectories in a text file.
        '''
        fmt = ['%d', '%.3f', '%.3f', '%.3f', '%.3f', '%.3f', '%.3f', '%.3f', 
               '%.3f', '%.3f', '%.3f']
        savetxt(fname, self.smoothed_trajs, fmt=fmt, delimiter='\t')








def smooth_traj_poly(traj, window, polyorder, repetitions=1):
    '''
    will smooth the particle position using a moving polynomial, where
    the velocities and accelerations are calculated by
    differentiating the polynomial coefficients analytically.
    
    input - 
    traj - a nested list (or numpy array), with shapes (3,N), N being
           the number of position samples.
    windows - integer or tuple with length 3. The window sliding size
              used for smoothing.
    polyorder - integer, the order of the polynomial used for the fitting.
    repetitions - integer, the number of times the smoothing should be 
                  repeated for the trajectory 
    
    
    returns-
    new_pos - the smoothed particle's position
    new_vel - the particle velocity
    new_acc - the particle acceleration
    '''
    
    if type(window) == int:
        window = (window, window, window)
    
    if type(polyorder) == int:
        polyorder = (polyorder, polyorder, polyorder)
    
    if type(window) != tuple: 
        raise TypeError('window must be either integer or tuple')
    
    if type(polyorder) != tuple: 
        raise TypeError('polyorder must be either a number or tuple')
    
    test = [w%2 != 1 for w in window]
    if sum(test) != 0:
        raise ValueError('windows must be a positive odd integer.')
    
    test = [polyorder[i] > window[i] for i in [0,1,2]]
    if sum(test) != 0:
        raise ValueError('polyorder cannot be larger than window.')
    
    
    N = len(traj[0])
    
    test = [w > N for w in window]
    if sum(test) != 0:
        raise ValueError('window cannot be larger than the trajectory length.')
    
    
    for j in range(repetitions):
        
        new_pos = [[0.0 for i in range(N)] for j in range(3)] 
        new_vel = [[0.0 for i in range(N)] for j in range(3)] 
        new_acc = [[0.0 for i in range(N)] for j in range(3)]
        
        sequence = zip(range(3), window, polyorder)
        for e,win,po in sequence:
            time_ = range(win)
            ev_point = [float(win/2)**i for i in range(po + 1)][::-1]
            
            Deriv_mat = []
            for i in range(po+1):
                a = [0.0 for j in range(po+1)]
                if i!=0:
                    a[i-1] = po - (i-1)
                Deriv_mat.append(a)
            
            
            for i in range(N):
                if i < win/2:            # get the portion of size widow for
                    p_ = traj[e][:win]   # fitting the polynomial
                elif N - (1+i) < win/2:
                    p_ = traj[e][-1*win:]
                else:
                    p_ = traj[e][i-int(win/2) : i+int(win/2)+1]
                
                C = fit_polynomial(time_, p_, po)      
                C1 = dot(Deriv_mat, C)
                C2 = dot(Deriv_mat, C1)
                
                
                if i < win/2: 
                    ev_point = [float(i%(win/2))**k for k in range(po + 1)[::-1]]
                elif N - (1+i) < win/2:
                    ev_point = [float(win + i - N)**k for k in range(po + 1)[::-1]]
                else:
                    ev_point = [float(win/2-0.5)**k for k in range(po + 1)[::-1]]
                
                new_pos[e][i] = dot(ev_point, C)
                new_vel[e][i] = dot(ev_point, C1)
                new_acc[e][i] = dot(ev_point, C2)
            
        traj = array(new_pos)

    return new_pos, new_vel, new_acc

