# -*- coding: utf-8 -*-
"""
Created on Sat Aug 13 16:49:46 2022

@author: Giorgia
"""

from types import SimpleNamespace
import time
import itertools as it
import numpy as np
from scipy import optimize
#from scipy.interpolate import RegularGridInterpolator
from scipy.io import loadmat
import os
import scipy.interpolate as spi

def prob_workT(par,options,Vp,Cp,j_pos):
    r_sav     = par.r_sav
    r_debt    = par.r_debt
    beta      = par.beta
    gammac    = par.gammac
    w         = par.w
    Lambda    = par.Lambda
    
    S = par.grids[0][j_pos]## matlab uses {} because par.grids is a cell array which I do not manage to convert in python
    EDUC      = par.educ
    AGE_PROF  = par.inc.age_prof
    INNO_pos  = par.inc.inno_pos
    FE_pos    = par.inc.fe_pos
    
    Sp = np.zeros([len(INNO_pos),len(EDUC),len(S),len(FE_pos)])
    C = np.zeros([len(INNO_pos),len(EDUC),len(S),len(FE_pos)])
    V = np.zeros([len(INNO_pos),len(EDUC),len(S),len(FE_pos)])
    boundgrid = np.zeros([len(INNO_pos),len(FE_pos),len(EDUC)])
    
    for educ in range (0,len(EDUC)):
        S       = par.grids[educ][j_pos]
        Spgrid  = par.grids[educ][j_pos+1]
        
        FE        = par.inc.fe[educ]
        INNO      = par.inc.z_val[educ][j_pos,:]
        
        r         = (r_sav * (Spgrid>=0) + r_debt * (Spgrid<0)).T
        rs        = (r_sav * (S>=0) + r_debt * (S<0))
        
        for ife in range (0,len(FE_pos)):
            splVp = spi.interp1d(Spgrid,Vp[educ,:,ife],axis=0, fill_value="extrapolate")
            # People don't work next period so only disutility of consumption. 
            Cpp = Cp[educ,:,ife]
            ucp = beta*(1+r)*Cpp**(-gammac)
            
            for inno  in range (0,len(INNO_pos)):
                dispinc = (1-Lambda) * exp(AGE_PROF[educ][0][j_pos]) * exp(FE[ife]) * exp(INNO[inno]) #AGE_PROF[educ][0][j_pos], i THINK THE [0] IS NEEDED ONLY BECAUSE IT IS A 1 DIMENSION DICTIONARY(FLEXIBLE) 
    
                feasible = (dispinc + (1+rs)*S - Spgrid[0]>0)
                not_feasible = (dispinc + (1+rs)*S - Spgrid[0]<=0)
    
                C[inno,educ,feasible,ife],Sp[inno,educ,feasible,ife],boundgrid[inno,ife,educ] = EGM(par,ucp,dispinc,Spgrid,S[feasible])
    
                vp = splVp(Sp[inno,educ,feasible,ife])
                V[inno,educ,feasible,ife] = (C[inno,educ,feasible,ife]**(1-gammac))/(1-gammac) + beta*vp
                C[inno,educ,not_feasible,ife] = 0
                V[inno,educ,not_feasible,ife] = (-(10**(5/gammac))**(1-gammac))/(1-gammac)
                Sp[inno,educ,not_feasible,ife]= Spgrid[0]

    
    errors = sum(boundgrid[:])/np.size(C[:])
    if options.timer_on == 'Y':
            if sum(boundgrid[:]) >= len(EDUC)*len(FE_pos)*len(INNO_pos):
                 print(f"j : {par.age[j_pos-1]}, Share of errors (increase grid) = {errors}")
##################################################################################################            
    return V, C, Sp
    
#########################THE PROBLEM IS THAT i AM MIZING TYPES IN A LIST
#compare with matlab calling C[32][4][2] or any other combination of [inno][educ], first row [0] is 00000 for C and S, -2 for V. if we comment out " C[inno,educ,not_feasible,ife] = 0", everything is fine
#it works great except C[32][:][:][0],
#C(not_feasible,ife,educ,inno)