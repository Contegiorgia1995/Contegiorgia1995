# -*- coding: utf-8 -*-
"""
Created on Sat Aug  6 14:06:12 2022

@author: Giorgia
"""
from types import SimpleNamespace
import time
import itertools as it
import numpy as np
from scipy import optimize
from scipy import stats
from __future__ import print_function
from numpy import *
from statistics import NormalDist



options.AdultRisk       = 'Y';



def discretize_nonstationary_ar_rouwen(rho, sigma_z_l,sigma_z_t, n):
 #Rouwen method to discretize AR(1), non-stationary, process w symmetric innovations
 # Based on WP Fella, Gallipoli and Pan (2015)
 # It matches persistence of the original process exactly
    
    q = 1/2 * (1+rho*sigma_z_l/sigma_z_t)
    nu = sqrt(n-1) * sigma_z_t

    P = [[q, 1-q],[1-q, q]]

    for i in range (2,n):
        P = np.dot(q,[np.append(np.append(np.reshape(P,[i,i]), np.zeros([i,1]),axis = 1), np.zeros([1,i+1]),  axis = 0)]) + np.dot((1-q),[np.append(np.append(np.zeros([i,1]), np.reshape(P,[i,i]), axis = 1),np.zeros([1,i+1]), axis = 0)]) + np.dot((1-q),[np.append(np.zeros([1,i+1]), np.append(np.reshape(P,[i,i]), np.zeros([i,1]), axis = 1), axis = 0)]) + np.dot(q,[np.append(np.zeros([1,i+1]), np.append(np.zeros([i,1]), np.reshape(P,[i,i]), axis = 1), axis = 0)])
        P = np.reshape(P, [i+1,i+1])
        P[1:i,:] = P[1:i,:]/2

    zgrid = linspace(-nu, nu, n).T
    
    return   P, zgrid

