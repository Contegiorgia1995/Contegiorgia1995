
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

model = SimpleNamespace() #model
model.par = SimpleNamespace() #param
model.sol = SimpleNamespace() 
model.options = SimpleNamespace()
options = model.options
par = model.par
par_est = model.par

educ  = [0,1,2]
educ = educ[-1]

w                   = par.w = 1
AGE_PROF            = inc.age_prof
FE                  = par.inc_fe
FE_pos              = inc_fe_pos
#fe = FE_pos #### added myself to fit with prob_ret
EDUC                = par.educ
    
def ret_rep(par,fe,educ,options):
# keyboard
# Pension replacement rate: 
# Source: Krueger, Ludwig 2015
# Parameters
    
    
    
## 1. Average lifetime income (AIME)
    age_ret             = range(par.Je2_pos,par.Jr_pos-1)
    totinc              = w * np.sum(exp(np.reshape(np.repeat([FE[educ][:len(fe)][:]],len(age_ret),axis = 0),[len(age_ret),len(fe)]) + np.reshape(np.repeat([AGE_PROF[educ][0][age_ret]],len(fe)),[len(age_ret),len(fe)])), axis = 0)##### get rid of rrays in first element of the sum

    prob_educ           = par.dist_educ
    den                 = 0
    h0_pos              = math.ceil(par.N_fe/2) # To find individual that starts with h0 = average
    
    for ie in range (0,len(par.educ)):
        den             =  den + prob_educ[ie] * w *np.sum(exp(FE[educ][h0_pos-1].T + AGE_PROF[ie][0][age_ret]))
    y                   = totinc/den

## 2. Marginal replacement rates
    if options.ExPolRetirement == 'Y':
        tau1 = options.PolRetMult[0]
        tau2 = options.PolRetMult[1]
        tau3 = options.PolRetMult[2]
    elif options.ExPolRetirement == 'N':
        tau1   = 0.9
        tau2   = 0.32
        tau3   = 0.15    

##3. Bendpoints
    b1                  = 0.24
    b2                  = 1.35
    b3                  = 1.99

## 4. Replacement rate
    rep_rate            = (tau1*y)*(y<= b1) + (tau1*b1 + tau2*(y-b1))*(y>b1)*( y<= b2) + (tau1*b1 + tau2*(b2-b1) + tau3*(y- b2))*(y>b2)*( y<= b3) + (tau1*b1 + tau2*(b2-b1) + tau3*(b3-b2))*(y>b3);
                
## Replacement Benefits:
    avg_inc             = den/len(age_ret)
    rep                 = avg_inc * rep_rate
    return rep
    #         fprintf('educ = %i, fe = %i, y=%3.3f,ret_rep = %3.3f , rep = %3.3f, avg inc = %3.3f \n',educ,fe,y,rep_rate,rep,totinc/length(age_ret));
