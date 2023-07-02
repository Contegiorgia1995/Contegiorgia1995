# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 16:56:07 2022

@author: Giorgia
"""
from types import SimpleNamespace
import time
import itertools as it
import numpy as np
from scipy import optimize
from scipy.interpolate import RegularGridInterpolator
import scipy.interpolate as spi


##all coming from file parameters



def prob_ret(par,options,Vp,Cp,j_pos):
                ##########    r_sav = par.r_sav
    r_sav = par.r_sav
    r_debt = par.r_debt

    beta = par.beta
    gammac = par.gammac
    EDUC = par.educ
    
    par.N_fe = 5
    N_fe         = par.N_fe # from parameters par.N_fe = 5
    par.inc.fe_pos   = np.linspace(0,N_fe-1,N_fe) ## again from parameters
    inc.fe_pos = par.inc.fe_pos
    FE_pos = inc.fe_pos #form par_income _process we have inc_fe_pos = np.arange(0,N_fe,N_fe) (still need N_fe)
    
##########################indexing using j_pos are right only if we use j_pos = j in solve_model!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if j_pos == par.Jd_pos-2: #consume all, next period is dead
        S = par.grids[0][j_pos]## matlab uses {} because par.grids is a cell array which I do not manage to convert in python
        Sp = np.zeros([len(EDUC),len(S),len(FE_pos)])
        C = np.zeros([len(EDUC),len(S),len(FE_pos)])
        V = np.zeros([len(EDUC),len(S),len(FE_pos)])
        
        for educ in range(0,len(EDUC)):
            S = par.grids[educ][j_pos]
            C[educ,:,:,] = ((1+r_sav)*np.repeat([S],len(FE_pos),axis =0) + np.reshape(np.repeat(ret_rep(par,FE_pos,educ,options),len(S)),[len(FE_pos),len(S)])).T
            V[educ,:,:,]= (C[educ,:,:,]**(1-gammac))/(1-gammac)
    else:# solve Euler equation
        S = par.grids[0][j_pos-1]## matlab uses {} because par.grids is a cell array which I do not manage to convert in python
        Sp = np.zeros([len(EDUC),len(S),len(FE_pos)])
        C = np.zeros([len(EDUC),len(S),len(FE_pos)])
        V = np.zeros([len(EDUC),len(S),len(FE_pos)])

        boundgrid = np.zeros([len(FE_pos),len(EDUC)])

        for educ in range (0,len(EDUC)):
            S       = par.grids[educ][j_pos]
            Spgrid  = par.grids[educ][j_pos+1]
        
            for ife in range (0,len(FE_pos)):
                ##############################################
                #splVp = RegularGridInterpolator((Spgrid,),Vp[educ][ife][:]) 
                splVp = spi.interp1d(Spgrid,Vp[educ,:,ife],axis=0, fill_value="extrapolate")
                #splVp = scipy.interpolate.interp1d(Spgrid,Vp[educ][ife][:])
                ##############################################
                Cpp   = Cp[educ,:,ife] #actually Cp, not C

                ucp     = beta*(1+r_sav)*Cpp**(-gammac)

                dispinc = ret_rep(par,FE_pos,educ,options)
                dispinc = dispinc[ife]

                C[educ,:,ife],Sp[educ,:,ife],boundgrid[ife,educ] = EGM(par,ucp,dispinc,Spgrid,S)

                V[educ,:,ife] = C[educ,:,ife]**(1-gammac)/(1-gammac) + beta*splVp(Sp[educ,:,ife])
                
                errors = sum(boundgrid[:])/np.size(C[:])
    
###################################################################    
        if options.timer_on == 'Y':
             if sum(boundgrid[:]) >= len(EDUC)*len(FE_pos):
                 print(f"j : {par.age[j_pos-1]}, Share of errors (increase grid) = {errors}")
##################################################################################################            
    return V, C, Sp
