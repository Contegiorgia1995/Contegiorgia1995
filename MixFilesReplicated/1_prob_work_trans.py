# -*- coding: utf-8 -*-
"""
Created on Sat Oct 22 14:25:49 2022

@author: Giorgia
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 21:34:08 2022

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
import math 

def prob_work_trans(par,options,Vp,Cp,Vc0,j_pos):
    
    r_sav     = par.r_sav
    r_debt    = par.r_debt
    beta      = par.beta
    gammac    = par.gammac
    w         = par.w
    Lambda    = par.Lambda
    lambdan   = par.lambdan
    gamman    = par.gamman
    
    S         = par.grids[0][j_pos]
    EDUC      = par.educ
    
    AGE_PROF  = par.inc.age_prof
    INNO_pos  = par.inc.inno_pos
    INNO_posT = np.reshape(INNO_pos,[1,5])
    FE_pos    = par.inc.fe_pos
    
    N         = np.arange(0,4,1)
    PHI_pos   = np.arange(0,len(par.PHI))
    
    Sp = np.zeros([len(INNO_pos),len(EDUC),len(N),len(PHI_pos),len(S),len(FE_pos)])
    Sp[:] = np.nan
    C  = np.zeros([len(INNO_pos),len(EDUC),len(N),len(PHI_pos),len(S),len(FE_pos)])
    C[:] = np.nan
    V  = np.zeros([len(INNO_pos),len(EDUC),len(N),len(PHI_pos),len(S),len(FE_pos)])
    V[:] = np.nan
    Ck = np.zeros([len(INNO_pos),len(EDUC),len(N),len(PHI_pos),len(S),len(FE_pos)])
    Ck[:] = np.nan
    Tp = np.zeros([len(INNO_pos),len(EDUC),len(N),len(PHI_pos),len(S),len(FE_pos)])
    Tp[:] = np.nan
    boundgrid = np.zeros([len(INNO_pos),len(N),len(PHI_pos),len(FE_pos),len(EDUC)])
    
    ## Without children and phi = 1 (phi is transfer to each child)
    if options.Fertility == "Endo": # Only need to solve case without children in benchmark case where we allow them to choose endogenous fertility
        i_n   = 0
        iphi =0
        ig   = 0 #Group does not matter in value function tomorrow if no children
        
        for educ in range (0,len(EDUC)):
            
            S          = par.grids[educ][j_pos]
            Spgrid  = par.grids[educ][j_pos+1]
            
            FE        = par.inc.fe[educ]
            INNO      = par.inc.z_val[educ][j_pos,:]
            
            
            r         = (r_sav * (Spgrid>=0) + r_debt * (Spgrid<0)).T
            rs        = (r_sav * (S>=0) + r_debt * (S<0))
            
            inno3 = np.repeat(INNO_posT,len(S),axis=0)
            INNOp_prob = par.inc.z_prob[educ][:,j_pos+1,:].T # Note I need j_pos + 1 here!

            INNO2,S2  = np.meshgrid(INNO_pos,Spgrid)
            
            for ife in range (0,len(FE_pos)):

                splVp = GlobalSpline2D(Spgrid, INNO_pos,Vp[:,educ,:,ife] )
           
                for inno  in range (0,len(INNO_pos)):
                
                    innop_prob = INNOp_prob[inno,:]
                    Cpp = Cp[:,educ,:,ife].T
                    
                    if min(np.ravel(Cpp.any())<0 and abs(min(np.ravel(Cpp.any()))))>1e-6:
                        Cpp = max(Cpp,0)
                    ucp = beta*np.reshape(((1+r)),[80,])*np.dot(Cpp**(-gammac), np.reshape(innop_prob,[5,]))
            
                    
                    h   = exp(AGE_PROF[educ][0][j_pos]) * exp(FE[ife])* exp(INNO[inno])
                    dispinc = w*h*(1-Lambda)
                    
                    feasible = (dispinc + (1+rs)*S - Spgrid[0]>0)
                    #not_feasible = (1-feasible)
                    not_feasible = (dispinc + (1+rs)*S - Spgrid[0]<=0)

                    C[inno,educ,i_n,iphi,feasible,ife],Sp[inno,educ,i_n,iphi,feasible,ife],boundgrid[inno,i_n,iphi,ife,educ] = EGM(par,ucp,dispinc,Spgrid,S[feasible])
                    Ck[inno,educ,i_n,iphi,feasible,ife] = np.zeros(len(C[inno,educ,i_n,iphi,feasible,ife]))                    
                    
                    Sp3 = np.repeat(np.reshape(Sp[inno,educ,i_n,iphi,feasible,ife],[80,1]),len(INNO_pos),axis = 1)
                    
                    vp = np.dot(splVp(Sp3[:,0],inno3[0,:]).T,innop_prob)
                    # vp = np.dot((splVp(Sp3[:,0],inno3[0,:]).T),innop_prob)
                    # for j in range(0,len(vp)):
                    #     if math.isnan(vp[j]):
                    #         vp[j] = -inf
                    # vp = np.sort(vp)
                    # for j in range(0,len(vp)):
                    #     if math.isinf(vp[j]):
                    #         vp[j] = nan
                    
                    V[inno,educ,i_n,iphi,feasible,ife] = C[inno,educ,i_n,iphi,feasible,ife]**(1-gammac)/(1-gammac) + beta*vp[feasible]
                    
                    C[inno,educ,i_n,iphi,not_feasible,ife] = 0;
                    V[inno,educ,i_n,iphi,not_feasible,ife] = -(10**(5/gammac))**(1-gammac)/(1-gammac)
                    Sp[inno,educ,i_n,iphi,not_feasible,ife]= Spgrid[0]
                    
                    Tp[inno,educ,i_n,iphi,:,ife] = np.zeros(len(C[inno,educ,i_n,iphi,:,ife]))
        
        # Fill cell for all PHI;
        for iphi in range (1,len(PHI_pos)):
            for educ in range (0,len(EDUC)):
                for ife in range (0,len(FE_pos)):
                    for inno in range (0,len(FE_pos)):
                        C[inno,educ,i_n,iphi,:,ife]  = C[inno,educ,i_n,0,:,ife]
                        Ck[inno,educ,i_n,iphi,:,ife] = Ck[inno,educ,i_n,0,:,ife]
                        Sp[inno,educ,i_n,iphi,:,ife] = Sp[inno,educ,i_n,0,:,ife]
                        V[inno,educ,i_n,iphi,:,ife]  = V[inno,educ,i_n,0,:,ife]
                        Tp[inno,educ,i_n,iphi,:,ife] = Tp[inno,educ,i_n,0,:,ife]
                        boundgrid[inno,i_n,0,ife,educ] = boundgrid[inno,i_n,0,ife,educ]

#############
#PROBLEMS: 
# vp is worng as consequence of wrong Vp which is wrong from wrong V wrong itself because splVp is not correct and gives a shitty aproximantion

## Case with Children
    if options.Fertility == 'Endo': # Only need to solve case without children in benchmark case where we allow them to choose endogenous fertility
            in_1    = 1
    elif options.Fertility == 'Exo':
            in_1    = 0
            
    for educ in range (0,len(EDUC)):
    
        INNOp_prob  = par.inc.z_prob[educ][:,j_pos+1,:].T  # Note I need j_pos + 1 here!
        
        S          = par.grids[educ][j_pos]
        Spgrid  = par.grids[educ][j_pos+1]
    
        r         = (r_sav * (Spgrid>=0) + r_debt * (Spgrid<0)).T
        rs        = (r_sav * (S>=0) + r_debt * (S<0))
        
        inno3 = np.repeat(INNO_posT,len(S),axis=0)
        
        FE        = par.inc.fe[educ]
        INNO      = par.inc.z_val[educ][j_pos,:]
        
        INNO2,S2  = np.meshgrid(INNO_pos,Spgrid)
        
        for i_n in range(in_1,len(N)):
        #for i_n in range(in_1,in_1+1):
            #i_n = 1
            n             = par.N[i_n]##### careful here that par.N is [0,1,2,3] also in Matlab
            for iphi in range (0,len(PHI_pos)):
                phi        = par.PHI[iphi]
                for ife in range (0,len(FE_pos)):
                    splVp = GlobalSpline2D(Spgrid, INNO_pos,Vp[:,educ,:,ife])
                    for inno in range (0,len(INNO_pos)):
                        innop_prob      = INNOp_prob[inno,:]
                        h   = exp(AGE_PROF[educ][0][j_pos]) * exp(FE[ife])* exp(INNO[inno])*(1-Lambda)
                        
                        Cpp = Cp[:,educ,:,ife].T
                        if min(np.ravel(Cpp.any())<0 and abs(min(np.ravel(Cpp.any()))))>1e-6:
                            Cpp = max(Cpp,0)
                        ucp = beta*np.reshape(((1+r)),[80,])*np.dot(Cpp**(-gammac), np.reshape(innop_prob,[5,]))
                        
                        
                        labor_inc = w*h*(1-Lambda)
                        
                        n_final = n* (par.fam_size/2)# 3 as in matlab
                        ChildCost_0 = ChildCost(par,labor_inc,n_final,0,options)
                        ChildCost_1 = ChildCost(par,labor_inc,n_final,1,options)
                        ChildCost_opt = min(ChildCost_0,ChildCost_1)
                        Tp[inno,educ,i_n,iphi,:,ife]   = (ChildCost_1>ChildCost_0)
                        
                        dispinc = labor_inc- ChildCost_opt -  n_final  *phi
                        feasible = (dispinc + (1+rs)*S - Spgrid[0]>0)
                        not_feasible = (dispinc + (1+rs)*S - Spgrid[0]<=0)
                    
         ########################### Altruism
                        hp              = exp(AGE_PROF[educ][0][j_pos]) * exp(FE[ife])* exp(INNO[inno])
                        Vc0_aux     = np.reshape(Vc0[:,iphi,:],[par.N_fe*len(par.psy_val_hs)])
                        
                        
                        Vc0_prob_h0 = PG(hp)
                        Vc0_prob_psy= par.psy_prob[educ,:].T
                        Vc0_prob    = gridmake(Vc0_prob_psy,Vc0_prob_h0).T## 500 x2
                        Vc0_prob    = Vc0_prob[:,1] * Vc0_prob[:,0]
                        
                        Gn = beta * lambdan * ( n_final )**(gamman) * np.dot(Vc0_aux.T, Vc0_prob)
                        
                        if options.Ck == 'Yes':
                                C[inno,educ,i_n,iphi,feasible,ife],Sp[inno,educ,i_n,iphi,feasible,ife],boundgrid[inno,i_n,iphi,ife,educ] = EGM_withchild(par,ucp,dispinc,Spgrid,S[feasible],n_final) 

                                f_n     = (lambdan/n_final**(1-gamman))**(1/gammac)
                                Ck[inno,educ,i_n,iphi,feasible,ife]     = f_n * C[inno,educ,i_n,iphi,feasible,ife]
                                                          
                                Sp3 = np.repeat(np.reshape(Sp[inno,educ,i_n,iphi,:,ife],[len(feasible),1]),len(INNO_pos),axis = 1)   
                                
                                vp = np.dot((splVp(Sp3[:,0],inno3[0,:]).T),innop_prob)
                                # for j in range(0,len(vp)):
                                #     if math.isnan(vp[j]):
                                #         vp[j] = -inf
                                # vp = np.sort(vp)
                                # for j in range(0,len(vp)):
                                #     if math.isinf(vp[j]):
                                #         vp[j] = nan

                                  
                    
   
                                V[inno,educ,i_n,iphi,feasible,ife] = (C[inno,educ,i_n,iphi,feasible,ife]**(1-gammac))/(1-gammac)*(C[inno,educ,i_n,iphi,feasible,ife]>=0) + lambdan* n_final**(gamman) * (Ck[inno,educ,i_n,iphi,feasible,ife]**(1-gammac))/(1-gammac)*(C[inno,educ,i_n,iphi,feasible,ife]>=0) -((10**(5/gammac))**(1-gammac))/(1-gammac)*(C[inno,educ,i_n,iphi,feasible,ife]<0) + beta*vp[feasible]*(C[inno,educ,i_n,iphi,feasible,ife]>=0) + Gn*(C[inno,educ,i_n,iphi,feasible,ife]>=0)
                                # Fix bad extrapolation:
                                Sp[inno,educ,i_n,iphi,feasible,ife]  = Sp[inno,educ,i_n,iphi,feasible,ife]*(C[inno,educ,i_n,iphi,feasible,ife] >=0) + Spgrid[0]*(C[inno,educ,i_n,iphi,feasible,ife]<0)
                                C[inno,educ,i_n,iphi,feasible,ife]  = C[inno,educ,i_n,iphi,feasible,ife]*(C[inno,educ,i_n,iphi,feasible,ife] >= 0) + 0*(C[inno,educ,i_n,iphi,feasible,ife] <0)
                                Ck[inno,educ,i_n,iphi,feasible,ife] = Ck[inno,educ,i_n,iphi,feasible,ife]*(Ck[inno,educ,i_n,iphi,feasible,ife]>= 0) + 0*(Ck[inno,educ,i_n,iphi,feasible,ife]<0)

                                C[inno,educ,i_n,iphi,not_feasible,ife] = 0
                                Ck[inno,educ,i_n,iphi,not_feasible,ife] = 0
                                V[inno,educ,i_n,iphi,not_feasible,ife] = -(10**(5/gammac))**(1-gammac)/(1-gammac)
                                Sp[inno,educ,i_n,iphi,not_feasible,ife]= Spgrid[0]
                                
                        elif options.Ck == 'No':
                                C[inno,educ,i_n,iphi,feasible,ife],Sp[inno,educ,i_n,iphi,feasible,ife],boundgrid[inno,i_n,iphi,ife,educ] = EGM(par,ucp,dispinc,Spgrid,S[feasible])
                                Ck[inno,educ,i_n,iphi,feasible,ife] =np.zeros(len(C[inno,educ,i_n,iphi,feasible,ife]))
                                Sp3         = np.repeat(np.reshape(Sp[inno,educ,i_n,iphi,:,ife],[len(feasible),1]),len(INNO_pos),axis = 1)    
                                vp          = np.dot(splVp(Sp3[:,0],inno3[0,:]).T,innop_prob)### wrong interpolation
                                vp = np.dot((splVp(Sp3[:,0],inno3[0,:]).T),innop_prob)
                                # for j in range(0,len(vp)):
                                #     if math.isnan(vp[j]):
                                #         vp[j] = -inf
                                # vp = np.sort(vp)
                                # for j in range(0,len(vp)):
                                #     if math.isinf(vp[j]):
                                #         vp[j] = nan
                                        
                                        
                                V[inno,educ,i_n,iphi,feasible,ife] = C[inno,educ,i_n,iphi,feasible,ife]**(1-gammac)/(1-gammac)*(C[inno,educ,i_n,iphi,feasible,ife]>=0) -(10**(5/gammac))**(1-gammac)/(1-gammac)*(C[inno,educ,i_n,iphi,feasible,ife]<0) + beta*vp[feasible]*(C[inno,educ,i_n,iphi,feasible,ife]>=0) + Gn*(C[inno,educ,i_n,iphi,feasible,ife]>=0)
                                
                                # Fix bad extrapolation:
                                Sp[inno,educ,i_n,iphi,feasible,ife] = Sp[inno,educ,i_n,iphi,feasible,ife]*(C[inno,educ,i_n,iphi,feasible,ife] >=0) + Spgrid[0]*(C[inno,educ,i_n,iphi,feasible,ife]<0)
                                C[inno,educ,i_n,iphi,feasible,ife] = C[inno,educ,i_n,iphi,feasible,ife]*(C[inno,educ,i_n,iphi,feasible,ife] >=0) + 0*(C[inno,educ,i_n,iphi,feasible,ife]<0)
                                C[inno,educ,i_n,iphi,not_feasible,ife] = 0
                                Ck[inno,educ,i_n,iphi,not_feasible,ife] = 0
                                V[inno,educ,i_n,iphi,not_feasible,ife] = -(10**(5/gammac))**(1-gammac)/(1-gammac)
                                Sp[inno,educ,i_n,iphi,not_feasible,ife]= Spgrid[0]


##Search Grid:
    V_2    = np.zeros([len(INNO_pos),len(EDUC),len(N),len(S),len(FE_pos)])
    V_2[:] = np.nan
    C_2    = np.zeros([len(INNO_pos),len(EDUC),len(N),len(S),len(FE_pos)])
    C_2[:] = np.nan
    Ck_2   = np.zeros([len(INNO_pos),len(EDUC),len(N),len(S),len(FE_pos)])
    Ck_2[:] = np.nan
    Sp_2   = np.zeros([len(INNO_pos),len(EDUC),len(N),len(S),len(FE_pos)])
    Sp_2[:] = np.nan
    PHIp_2 = np.zeros([len(INNO_pos),len(EDUC),len(N),len(S),len(FE_pos)])
    PHIp_2[:] = np.nan
    Tp_2   = np.zeros([len(INNO_pos),len(EDUC),len(N),len(S),len(FE_pos)])
    Tp_2[:] = np.nan
    
    posP = np.zeros([len(INNO_pos),len(EDUC),len(N),len(S),len(FE_pos)])
    for i_s in range (0,len(S)):
        
        for ife in range (0,len(FE_pos)):
                
            for inno in range (0,len(INNO_pos)):
                for educ in range (0,len(EDUC)):
                    for i_n in range (0,len(N)):
                        Vaux                           = np.squeeze(V[inno,educ,i_n,:,i_s,ife])     
                        posauxP                        = np.argmax(Vaux)### otherwise it would be zero and in Matlab it is 1
                        posP[inno,educ,i_n,i_s,ife]    = int64(posauxP)
                        V_2[inno,educ,i_n,i_s,ife]     = V[inno,educ,i_n,int64(posP[inno,educ,i_n,i_s,ife]),i_s,ife]################understand!!!!!!!!!!!!!!!!!!!!!!!!!!
                        C_2[inno,educ,i_n,i_s,ife]     = C[inno,educ,i_n,int64(posP[inno,educ,i_n,i_s,ife]),i_s,ife]
                        Ck_2[inno,educ,i_n,i_s,ife]    = Ck[inno,educ,i_n,int64(posP[inno,educ,i_n,i_s,ife]),i_s,ife]
                        Sp_2[inno,educ,i_n,i_s,ife]    = Sp[inno,educ,i_n,int64(posP[inno,educ,i_n,i_s,ife]),i_s,ife]
                        Tp_2[inno,educ,i_n,i_s,ife]    = Tp[inno,educ,i_n,int64(posP[inno,educ,i_n,i_s,ife]),i_s,ife]
                        PHIp_2[inno,educ,i_n,i_s,ife]  = posP[inno,educ,i_n,i_s,ife]
                        

    
    errors = sum(boundgrid[:])/len(C[:])
    if options.timer_on == 'Y':
            if sum(boundgrid[:]) >= len(EDUC)*len(FE_pos)*len(INNO_pos)*len(PHI_pos)*len(N):
                printf('j : {par.age(j_pos)}, Share of errors (increase grid) ={errors} \n')

    return V_2, C_2, Ck_2, Sp_2, PHIp_2, Tp_2
####
#Problems:
    
#in search Grid we use Vaux = np.squeeze(V[inno,educ,i_n,:,i_s,ife]) which determines the grid. Since my V[] is output of a bad approximation, numbers are messed up: all V_2, SP_2, Ck_2, PHIp_2

#########I am using C instead of V in Vaux but need to change that