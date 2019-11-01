#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 13:33:57 2018

@author: ron


Utility code to use for the MyPTV package.

"""



from numpy import dot, array
from numpy.linalg import inv, norm


def line_dist(O1, r1, O2, r2):
    '''
    2 lines are defined as (O1 + a r1) and  (O2 + b r2), where O are origins,
    r are direction vectors, and a,b are variables (in n dimensions). 
    This utility calculates the minimal distance between these 2 lines.
    
    input - 
    O1,O2,r1,r2 (arrays, n) - line parameters
    
    output - 
    dist (float) -the minimum distance between the lines
    x (array, n)- the point that is nearest to the two points crossing
    '''
    
    # find the a,b that minimize the distance:
    A = array([[dot(r1,r1), -dot(r1,r2)],
               [dot(r1,r2), -dot(r2,r2)]])
    
    B = array([[dot(r1,O2-O1)],[dot(r2,O2-O1)]])
    
    try:
        a,b = dot(inv(A), B)
    except:
        a,b = 0.0, 0.0
    
    # use the a,b to calc the minimum distance:
    l1,l2 = O1 + a*r1 , O2 + b*r2
    dist = norm(l1 - l2)
    x = (l1+l2)*0.5
    
    return dist, x



def point_line_dist(O,r,P):
    '''
    for a line (O + a r) and a point P, this returns the distance between
    the line and the point.
    '''
    a = dot(r, P - O) / dot(r, r)
    l = O + a*r
    d = norm(l-P)
    return d


