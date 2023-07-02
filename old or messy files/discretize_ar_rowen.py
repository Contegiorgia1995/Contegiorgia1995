# -*- coding: utf-8 -*-
"""
Created on Sat Aug  6 13:59:13 2022

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
import pyimsl
from statistics import NormalDist
import par_income_process
from par_income_process import par_income_process

model = SimpleNamespace() #model
model.par = SimpleNamespace() #param
model.sol = SimpleNamespace() 
model.options = SimpleNamespace()
options = model.options
par = model.par
par_est = model.par

options.AdultRisk       = 'Y'; #from start

def discretize_ar_rouwen(mu_eps,rho, sigma_eps, n):
    
# Rouwen method to discretize AR(1) process w symmetric innovations
# It matches persistence of the original process exactly

#     mu_eps = 0
    
    q = (rho+1)/2
    nu = sqrt(n-1) * sigma_eps/(1-rho**2)**0.5

    P = [[q, 1-q],[1-q, q]]
    
    for i in range (2, n):
        P = np.dot(q,[np.append(np.append(np.reshape(P,[i,i]), np.zeros([i,1]),axis = 1), np.zeros([1,i+1]),  axis = 0)]) + np.dot((1-q),[np.append(np.append(np.zeros([i,1]), np.reshape(P,[i,i]), axis = 1),np.zeros([1,i+1]), axis = 0)]) + np.dot((1-q),[np.append(np.zeros([1,i+1]), np.append(np.reshape(P,[i,i]), np.zeros([i,1]), axis = 1), axis = 0)]) + np.dot(q,[np.append(np.zeros([1,i+1]), np.append(np.zeros([i,1]), np.reshape(P,[i,i]), axis = 1), axis = 0)])
        P = np.reshape(P, [i+1,i+1])
        P[1:i,:] = P[1:i,:]/2


    zgrid = np.linspace(mu_eps/(1-rho)-nu, mu_eps/(1-rho)+nu, n).T
    


