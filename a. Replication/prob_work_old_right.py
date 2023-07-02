# -*- coding: utf-8 -*-
"""
Created on Sun Aug 14 15:48:14 2022

@author: Giorgia
"""

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
from scipy.interpolate import RegularGridInterpolator
from scipy.io import loadmat
import os
import scipy.interpolate as spi
from scipy.interpolate import griddata
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import interpn
from scipy.interpolate import interp2d
from scipy.interpolate import NearestNDInterpolator

def prob_work_old(par,options,Vp,Cp,j_pos):
    r_sav     = par.r_sav
    r_debt    = par.r_debt
    beta      = par.beta
    gammac    = par.gammac
    w         = par.w
    Lambda    = par.Lambda
    #parameters(par_est, options)
    
    S = par.grids[0][j_pos]## matlab uses {} because par.grids is a cell array which I do not manage to convert in python
    EDUC      = par.educ
    AGE_PROF  = par.inc.age_prof
    INNO_pos  = par.inc.inno_pos
    FE_pos    = par.inc.fe_pos
    
    Sp = np.zeros([len(INNO_pos),len(EDUC),len(S),len(FE_pos)])
    C = np.zeros([len(INNO_pos),len(EDUC),len(S),len(FE_pos)])
    V = np.zeros([len(INNO_pos),len(EDUC),len(S),len(FE_pos)])
    boundgrid = np.zeros([len(INNO_pos),len(FE_pos),len(EDUC)])
    INNO_posT = np.reshape(INNO_pos,[1,5])
    inno3 = np.repeat(INNO_posT,len(S),axis=0)

    
    for educ in range (0,len(EDUC)):
        S       = par.grids[educ][j_pos]
        FE        = par.inc.fe[educ]
        INNO      = par.inc.z_val[educ][j_pos,:]
        Spgrid  = par.grids[educ][j_pos+1]
        r         = (r_sav * (Spgrid>=0) + r_debt * (Spgrid<0)).T
        rs        = (r_sav * (S>=0) + r_debt * (S<0))
        
        INNOp_prob = par.inc.z_prob[educ][:,j_pos+1,:].T # Note I need j_pos + 1 here!
        """it would be good to understand why I need to transpose par.inc.z_prob all the times, as I had to transpose P when I filled the matrix in "par_income_process"     
        """
        Spgrid     = par.grids[educ][j_pos+1]
        
        FE        = par.inc.fe[educ]
        INNO      = par.inc.z_val[educ][j_pos,:]
        
        r         = (r_sav * (Spgrid>=0) + r_debt * (Spgrid<0)).T
        rs        = (r_sav * (S>=0) + r_debt * (S<0))
        
        INNO2,S2  = np.meshgrid(INNO_pos,Spgrid)

    
        for ife in range (0,len(FE_pos)):
            ##points in grid given by INNO2
            ##values in Vp[]
            # splVp = RegularGridInterpolator(S2,Vp[:,educ,:,ife].T,method = 'linear')#########may need to invert S2 and INNO2 ????
            # splVp = LinearNDInterpolator((S2,INNO2), Vp[:,educ,:,ife].T)
            ######
            splVp = GlobalSpline2D(Spgrid, INNO_pos,Vp[:,educ,:,ife])
            
            # for i in range(0,len(S)):
            #     if Vp[inno,educ,i,ife] != 0:
            #         Vp[inno,educ,i,ife] = Vp[inno,educ,i,ife] 
            #         splVp = interp2d(Spgrid, INNO_pos,Vp[:,educ,:,ife],  'linear' )
            #     else:
            #         Vp[inno,educ,i,ife] = (-(10**(5/gammac))**(1-gammac))/(1-gammac)
            #         splVp = interp2d(Spgrid, INNO_pos,Vp[:,educ,:,ife],  'linear' )
            
            #splVp = NearestNDInterpolator(S2,Vp[:,educ,:,ife].T,  'linear' )
            #matplotlib.pyplot.contourf(S2, INNO2, Vp[:,educ,:,ife].T, cmap = 'jet')
            #splVp = interpn(INNO2,Vp[:,educ,:,ife],S2, )
            ## People don't work next period so only disutility of consumption. 
            for inno  in range (0,len(INNO_pos)):
                
                innop_prob = INNOp_prob[inno,:]
                Cpp = Cp[:,educ,:,ife].T
                if min(np.ravel(Cpp.any())<0 and abs(min(np.ravel(Cpp.any()))))>1e-6:
                    Cpp = max(Cpp,0)
                ucp = beta*np.reshape(((1+r)),[80,])*np.dot(Cpp**(-gammac), np.reshape(innop_prob,[5,]))
            
            
                dispinc = (1-Lambda) * exp(AGE_PROF[educ][0][j_pos]) * exp(FE[ife]) * exp(INNO[inno]) #AGE_PROF[educ][0][j_pos], i THINK THE [0] IS NEEDED ONLY BECAUSE IT IS A 1 DIMENSION DICTIONARY(FLEXIBLE) 
    
                feasible = (dispinc + (1+rs)*S - Spgrid[0]>0)
                not_feasible = (dispinc + (1+rs)*S - Spgrid[0]<=0)
                

                C[inno,educ,feasible,ife],Sp[inno,educ,feasible,ife],boundgrid[inno,ife,educ] = EGM(par,ucp,dispinc,Spgrid,S[feasible])
                
                #c_final, sp_final, boundgrid = EGM_old_worker(par,ucp,dispinc,Spgrid,S[feasible])
                Sp3 = np.repeat(np.reshape(Sp[inno,educ,:,ife],[80,1]),len(INNO_pos),axis = 1)

                vp = np.dot(splVp(Sp3[:,0],inno3[0,:]).T,innop_prob)
                #vp = splVp(Sp3)
                
                V[inno,educ,feasible,ife] = C[inno,educ,feasible,ife]**(1-gammac)/(1-gammac) + beta*vp[feasible]
                C[inno,educ,not_feasible,ife] = 0
                V[inno,educ,not_feasible,ife] = (-(10**(5/gammac))**(1-gammac))/(1-gammac)
                Sp[inno,educ,not_feasible,ife] = Spgrid[0]
    
    errors = sum(boundgrid[:])/size(C[:])
    if options.timer_on == 'Y':
        if sum(boundgrid[:]) >= len(EDUC)*len(FE_pos)*len(INNO_pos):
            print(f"j : {par.age[j_pos-1]}, Share of errors (increase grid) = {errors}")
##################################################################################################            
    return V, C, Sp
    
#########################THE PROBLEM IS THAT i AM MIZING TYPES IN A LIST
#compare with matlab calling C[32][4][2] or any other combination of [inno][educ], first row [0] is 00000 for C and S, -2 for V. if we comment out " C[inno,educ,not_feasible,ife] = 0", everything is fine
#it works great except C[32][:][:][0],
#C(not_feasible,ife,educ,inno)