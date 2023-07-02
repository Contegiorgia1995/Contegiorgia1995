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
import some_stream


import ret_rep
from ret_rep import ret_rep

import EGM
from EGM import EGM

model = SimpleNamespace() #model
model.par = SimpleNamespace() #param
model.sol = SimpleNamespace() 
model.options = SimpleNamespace()
options = model.options
par = model.par
par_est = model.par

##all coming from file parameters
par = model.par
par.time_period     = 2;
prime        = 3/100; # Smets and Wouters AER 2007
par.r        = (1+prime)**par.time_period-1
par.r_sav     = (1+prime)**par.time_period-1
r_sav      = par.r_sav            # Roys, Seshadri - 2014 XXX CONSIDER ADDING TAXES
#par_est.int_iota = 0.1; # XXX Add to estimation? XXX
par.int_iota = 0.1 ###### my modeification from par_est which is part of a more complicated function
par.r_debt   =  ((1+prime+par.int_iota)**par.time_period - 1)
par.gammac = 0.5 #from Matlab file parameters

par.educ = np.array([1,2,3])

par.age             = np.linspace(0,80,41)
par.Jd              = 80     
                               
par.Jd_pos          = list(par.age).index(par.Jd)+1 #np.where(par.age == par.Jd)


#par.grids = np.zeros((3,len(par.age))) #from 8.4 Savings in the file parameters


def prob_ret(par,options,Vp,Cp,j_pos):
                ##########    r_sav = par.r_sav
    beta = par.beta
    gammac = par.gammac
    EDUC = par.educ
    
    par.N_fe = 5
    N_fe         = par.N_fe # from parameters par.N_fe = 5
    par.inc.fe_pos   = np.linspace(0,N_fe-1,N_fe) ## again from parameters
    inc.fe_pos = par.inc.fe_pos
    FE_pos = inc.fe_pos #form par_income _process we have inc_fe_pos = np.arange(0,N_fe,N_fe) (still need N_fe)
    


    if j_pos == par.Jd_pos-1: #consume all, next period is dead
        S = par.grids[0][j_pos-1]## matlab uses {} because par.grids is a cell array which I do not manage to convert in python
        Sp = np.zeros([len(EDUC),len(FE_pos),len(S)])
        C = np.zeros([len(EDUC),len(FE_pos),len(S)])
        V = np.zeros([len(EDUC),len(FE_pos),len(S)])
        
        for educ in range(0,len(EDUC)):
            S = par.grids[educ][j_pos-1]
            C[educ,:,:,] = (1+r_sav)*np.repeat([S],len(FE_pos),axis =0) + np.reshape(np.repeat(ret_rep(par,FE_pos,educ,options),len(S)),[len(FE_pos),len(S)]) ## numbers slightly wrong

    else:# solve Euler equation
        S = par.grids[0][j_pos-1]## matlab uses {} because par.grids is a cell array which I do not manage to convert in python
        Sp = np.zeros([len(EDUC),len(FE_pos),len(S)])
        C = np.zeros([len(EDUC),len(FE_pos),len(S)])
        V = np.zeros([len(EDUC),len(FE_pos),len(S)])

        boundgrid = np.zeros([len(FE_pos),len(EDUC)])

        for educ in range (0,len(EDUC)):
            S       = par.grids[educ][j_pos-1]
            Spgrid  = par.grids[educ][j_pos]
        
            for ife in range (0,len(FE_pos)):
                ##############################################
                splVp = RegularGridInterpolator((Spgrid,),Vp[educ][ife][:]) 
                ##############################################
                Cpp   = Cp[educ][ife][:]

                ucp     = beta*(1+r_sav)* Cpp**(-gammac)

                dispinc = ret_rep(par,FE_pos,educ,options)
                dispinc = dispinc[ife]

                C[educ][ife][:],Sp[educ][ife][:],boundgrid[ife][educ] = EGM(par,ucp,dispinc,Spgrid,S)

                V[educ,ife,:] = C[educ,ife,:]**(1-gammac)/(1-gammac) + beta*splVp(Sp[educ,ife,:])
                
    errors = sum(boundgrid[:])/np.size(C[:])
    
###################################################################    
     if options.timer_on == 'Y':
         if sum(boundgrid[:]) >= len(EDUC)*len(FE_pos):
             print(f"j : {par.age[j_pos-1]}, Share of errors (increase grid) = {errors}")
##################################################################################################            
    return V, C, Sp
name = "Fred"
>>> f"He said his name is {name}."