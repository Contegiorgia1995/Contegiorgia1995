# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 14:10:17 2022

@author: Giorgia
"""
import numpy as np
import scipy.interpolate as spi
from scipy.interpolate import griddata
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import interpn
from scipy.interpolate import interp2d

def prob_fertility(par,options,Vp,Cp,j_pos):
    
# keyboard
# Solve household problem at fertility period: choice on number of children
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
    
    AGE_PROF   = par.inc.age_prof
    INNO_pos   = par.inc.inno_pos
    INNO_posT = np.reshape(INNO_pos,[1,5])

    FE_pos     = par.inc.fe_pos
    
    N          = np.arange(0,4,1)
    
    Sp = np.zeros([len(INNO_pos),len(EDUC),len(N),len(S),len(FE_pos)])
    Sp[:] = np.nan
    C  = np.zeros([len(INNO_pos),len(EDUC),len(N),len(S),len(FE_pos)])
    C[:] = np.nan
    V  = np.zeros([len(INNO_pos),len(EDUC),len(N),len(S),len(FE_pos)])
    V[:] = np.nan
    Ck = np.zeros([len(INNO_pos),len(EDUC),len(N),len(S),len(FE_pos)])
    Ck[:] = np.nan
    Tp = np.zeros([len(INNO_pos),len(EDUC),len(N),len(S),len(FE_pos)])
    Tp[:] = np.nan
    boundgrid = np.zeros([len(INNO_pos),len(N),len(FE_pos),len(EDUC)])
    
    
    ## Without children and phi = 1
    if options.Fertility == 'Endo':
    # Only need to solve case without children in benchmark case where we allow them to choose endogenous fertility
            i_n   = 0
            
            for educ in range (0,len(EDUC)):
                
                S          = par.grids[educ][j_pos]
                Spgrid  = par.grids[educ][j_pos+1]
                
                FE        = par.inc.fe[educ]
                INNO      = par.inc.z_val[educ][j_pos,:]
                
                inno3 = np.repeat(INNO_posT,len(S),axis=0)
                INNOp_prob = par.inc.z_prob[educ][:,j_pos+1,:].T # Note I need j_pos + 1 here!
                
                r         = (r_sav * (Spgrid>=0) + r_debt * (Spgrid<0))
                rs        = (r_sav * (S>=0) + r_debt * (S<0))
                
                INNO2,S2  = np.meshgrid(INNO_pos,Spgrid)
                
                
                for ife in range (0,len(FE_pos)):
                    splVp = GlobalSpline2D(Spgrid, INNO_pos,Vp[:,educ,i_n,:,ife] )
                    for inno  in range (0,len(INNO_pos)):
                    
                        innop_prob = INNOp_prob[inno,:]
                        
                        h   = exp(AGE_PROF[educ][0][j_pos]) * exp(FE[ife])* exp(INNO[inno])
                        
                        Cpp = Cp[:,educ,i_n,:,ife].T
                        
                        if min(np.ravel(Cpp.any())<0 and abs(min(np.ravel(Cpp.any()))))>1e-6:
                            Cpp = max(Cpp,0)
                        ucp = beta*np.reshape(((1+r)),[80,])*np.dot(Cpp**(-gammac), np.reshape(innop_prob,[5,]))
                        
                        dispinc = w*h*(1-Lambda) + par.fert_trans[i_n]
                        
                        feasible = (dispinc + (1+rs)*S - Spgrid[0]>0)
                        
                        not_feasible = (dispinc + (1+rs)*S - Spgrid[0]<=0)
                        
                        C[inno,educ,i_n,feasible,ife],Sp[inno,educ,i_n,feasible,ife],boundgrid[inno,i_n,ife,educ] = EGM(par,ucp.T,dispinc,Spgrid,S[feasible])
                        Ck[inno,educ,i_n,feasible,ife] = np.zeros(len(C[inno,educ,i_n,feasible,ife]))
                        
                        Sp3 = np.repeat(np.reshape(Sp[inno,educ,i_n,feasible,ife],[80,1]),len(INNO_pos),axis = 1)
                        
                        vp = np.dot((splVp(Sp3[:,0],inno3[0,:]).T),innop_prob)
                        # for j in range(0,len(vp)):
                        #     if math.isnan(vp[j]):
                        #         vp[j] = -inf
                        # vp = np.sort(vp)
                        # for j in range(0,len(vp)):
                        #     if math.isinf(vp[j]):
                        #         vp[j] = nan

                        
                        V[inno,educ,i_n,feasible,ife] = C[inno,educ,i_n,feasible,ife]**(1-gammac)/(1-gammac) + beta*vp[feasible]
                        
                        C[inno,educ,i_n,not_feasible,ife] = 0;
                        V[inno,educ,i_n,not_feasible,ife] = -(10**(5/gammac))**(1-gammac)/(1-gammac)
                        Sp[inno,educ,i_n,not_feasible,ife]= Spgrid[0]

                        Tp[inno,educ,i_n,:,ife] = np.zeros(len(C[inno,educ,i_n,:,ife]))
            
    
    ## Case with Children
    if options.Fertility == 'Endo': #Only need to solve case without children in benchmark case where we allow them to choose endogenous fertility
            in_1    = 1
    elif options.Fertility == 'Exo':
            in_1    = 0

    
    for educ in range (0,len(EDUC)):
        INNOp_prob = par.inc.z_prob[educ][:,j_pos+1,:].T
        S          = par.grids[educ][j_pos]
        Spgrid  = par.grids[educ][j_pos+1]
    
        r         = (r_sav * (Spgrid>=0) + r_debt * (Spgrid<0)).T
        rs        = (r_sav * (S>=0) + r_debt * (S<0))
        
        inno3 = np.repeat(INNO_posT,len(S),axis=0)
        FE        = par.inc.fe[educ]
        INNO      = par.inc.z_val[educ][j_pos,:]
        INNO2,S2  = np.meshgrid(INNO_pos,Spgrid)

        for i_n in range(in_1,len(N)):
            n             = par.N[i_n]
            for ife in range (0,len(FE_pos)):
                splVp = GlobalSpline2D(S2[:,0], INNO2[0,:],Vp[:,educ,i_n,:,ife])
                for inno in range (0,len(INNO_pos)):
                    innop_prob      = INNOp_prob[inno,:]
                    
                    h   = exp(AGE_PROF[educ][0][j_pos]) * exp(FE[ife])* exp(INNO[inno])
                    
                    Cpp = Cp[:,educ,i_n,:,ife].T
                    if min(np.ravel(Cpp.any())<0 and abs(min(np.ravel(Cpp.any()))))>1e-6:
                        Cpp = max(Cpp,0)
                    ucp = beta*np.reshape(((1+r)),[80,1])*np.dot(Cpp**(-gammac), np.reshape(innop_prob,[5,1]))
                    
                    labor_inc = w*h*(1-Lambda)
                    
                    n_final = n* (par.fam_size/2)

                    ChildCost_0 = ChildCost(par,labor_inc,n_final,0,options)
                    ChildCost_1 = ChildCost(par,labor_inc,n_final,1,options)
                    ChildCost_opt = min(ChildCost_0,ChildCost_1)
                    Tp[inno,educ,i_n,:,ife]   = (ChildCost_1>ChildCost_0)
                    
                    dispinc = labor_inc- ChildCost_opt + par.fert_trans[i_n]
                    feasible = (dispinc + (1+rs)*S - Spgrid[0]>0)
                    not_feasible = (dispinc + (1+rs)*S - Spgrid[0]<=0)
                    
                    if options.Ck == 'Yes':
                        
                        C[inno,educ,i_n,feasible,ife],Sp[inno,educ,i_n,feasible,ife],boundgrid[inno,i_n,ife,educ] = GEGM_withchild(par,ucp.T,dispinc,Spgrid,S[feasible],splVp,innop_prob,n_final) 

                        f_n     = (lambdan/n_final**(1-gamman))**(1/gammac)
                        Ck[inno,educ,i_n,feasible,ife]     = f_n * C[inno,educ,i_n,feasible,ife]
                        
                        Sp3 = np.repeat(np.reshape(Sp[inno,educ,i_n,:,ife],[len(feasible),1]),len(INNO_pos),axis = 1)                    
                        vp = np.dot(splVp(Sp3[:,0],inno3[0,:]).T,innop_prob)
                        
                        #vp = np.dot((splVp(Sp3[:,0],inno3[0,:]).T),innop_prob)
                        # for j in range(0,len(vp)):
                        #     if math.isnan(vp[j]):
                        #         vp[j] = -inf
                        # vp = np.sort(vp)
                        # for j in range(0,len(vp)):
                        #     if math.isinf(vp[j]):
                        #         vp[j] = nan
                                
                                
                        V[inno,educ,i_n,feasible,ife] = (C[inno,educ,i_n,feasible,ife]**(1-gammac))/(1-gammac)*(C[inno,educ,i_n,feasible,ife]>=0) + lambdan* n_final**(gamman) * (Ck[inno,educ,i_n,feasible,ife]**(1-gammac))/(1-gammac)*(C[inno,educ,i_n,feasible,ife]>=0) -((10**(5/gammac))**(1-gammac))/(1-gammac)*(C[inno,educ,i_n,feasible,ife]<0) + beta*vp[feasible]
                        
                        # Fix bad extrapolation:
                        Sp[inno,educ,i_n,feasible,ife]  = Sp[inno,educ,i_n,feasible,ife]*(C[inno,educ,i_n,feasible,ife] >=0) + Spgrid[0]*(C[inno,educ,i_n,feasible,ife]<0)
                        C[inno,educ,i_n,feasible,ife]  = C[inno,educ,i_n,feasible,ife]*(C[inno,educ,i_n,feasible,ife] >= 0) + 0*(C[inno,educ,i_n,feasible,ife] <0)
                        Ck[inno,educ,i_n,feasible,ife] = Ck[inno,educ,i_n,feasible,ife]*(Ck[inno,educ,i_n,feasible,ife]>= 0) + 0*(Ck[inno,educ,i_n,feasible,ife]<0)
                        
                        C[inno,educ,i_n,not_feasible,ife] = 0
                        Ck[inno,educ,i_n,not_feasible,ife] = 0
                        V[inno,educ,i_n,not_feasible,ife] = -(10**(5/gammac))**(1-gammac)/(1-gammac)
                        Sp[inno,educ,i_n,not_feasible,ife]= Spgrid[0]
                            
                    if options.Ck == 'No':
                        C[inno,educ,i_n,feasible,ife],Sp[inno,educ,i_n,feasible,ife],boundgrid[inno,i_n,ife,educ]= GEGM(par,ucp.T,dispinc,Spgrid,S[feasible],splVp,innop_prob)
                        Ck[inno,educ,i_n,feasible,ife] =np.zeros(len(C[inno,educ,i_n,iphi,feasible,ife]))
                        Ck[:] = np.nan   
                        Sp3 = np.repeat(np.reshape(Sp[inno,educ,i_n,feasible,ife],[80,1]),len(INNO_pos),axis = 1)
                        
                        vp = np.dot((splVp(Sp3[:,0],inno3[0,:]).T),innop_prob)
                        # for j in range(0,len(vp)):
                        #     if math.isnan(vp[j]):
                        #         vp[j] = -inf
                        # vp = np.sort(vp)
                        # for j in range(0,len(vp)):
                        #     if math.isinf(vp[j]):
                        #         vp[j] = nan

                        
                        V[inno,educ,i_n,feasible,ife] = C[inno,educ,i_n,feasible,ife]**(1-gammac)/(1-gammac) + beta*vp[feasible]
                        
                        # Fix bad extrapolation:
                        Sp[inno,educ,i_n,feasible,ife]  = Sp[inno,educ,i_n,feasible,ife]*(C[inno,educ,i_n,feasible,ife] >=0) + Spgrid[0]*(C[inno,educ,i_n,feasible,ife]<0)
                        C[inno,educ,i_n,feasible,ife]  = C[inno,educ,i_n,feasible,ife]*(C[inno,educ,i_n,feasible,ife] >= 0) + 0*(C[inno,educ,i_n,feasible,ife] <0)
                        
                        C[inno,educ,i_n,not_feasible,ife] = 0
                        V[inno,educ,i_n,not_feasible,ife] = -(10**(5/gammac))**(1-gammac)/(1-gammac)
                        Sp[inno,educ,i_n,not_feasible,ife]= Spgrid[0]
    
    
    
    ## Search Grid: only with endogenous fertility
    
    V_2    = np.zeros([len(INNO_pos),len(EDUC),len(S),len(FE_pos)])
    V_2[:] = np.nan
    C_2    = np.zeros([len(INNO_pos),len(EDUC),len(S),len(FE_pos)])
    C_2[:] = np.nan
    Ck_2   = np.zeros([len(INNO_pos),len(EDUC),len(S),len(FE_pos)])
    Ck_2[:] = np.nan
    Sp_2   = np.zeros([len(INNO_pos),len(EDUC),len(S),len(FE_pos)])
    Sp_2[:] = np.nan
    Np_2 = np.zeros([len(INNO_pos),len(EDUC),len(S),len(FE_pos)])
    Np_2[:] = np.nan
    Tp_2   = np.zeros([len(INNO_pos),len(EDUC),len(S),len(FE_pos)])
    Tp_2[:] = np.nan
    
    posN = np.zeros([len(INNO_pos),len(EDUC),len(S),len(FE_pos)])
    
    if options.Fertility == 'Endo': #Only need to solve case without children in benchmark case where we allow them to choose endogenous fertility
            for i_s in range (0,len(S)):
                for ife in range (0,len(FE_pos)):    
                    for inno in range (0,len(INNO_pos)):
                        for educ in range (0,len(EDUC)):
                            for i_n in range (0,len(N)):
                                Vaux                           = np.squeeze(V[inno,educ,:,i_s,ife])     
                                # Max wrt N:
                                posauxN                        = np.argmax(Vaux)### otherwise it would be zero and in Matlab it is 1
                                posN[inno,educ,i_s,ife]        = int64(posauxN)
                                V_2[inno,educ,i_s,ife]     = V[inno,educ,int64(posN[inno,educ,i_s,ife]),i_s,ife]################understand!!!!!!!!!!!!!!!!!!!!!!!!!!
                                C_2[inno,educ,i_s,ife]     = C[inno,educ,int64(posN[inno,educ,i_s,ife]),i_s,ife]
                                Ck_2[inno,educ,i_s,ife]    = Ck[inno,educ,int64(posN[inno,educ,i_s,ife]),i_s,ife]
                                Sp_2[inno,educ,i_s,ife]    = Sp[inno,educ,int64(posN[inno,educ,i_s,ife]),i_s,ife]
                                Tp_2[inno,educ,i_s,ife]    = Tp[inno,educ,int64(posN[inno,educ,i_s,ife]),i_s,ife]
                                Np_2[inno,educ,i_s,ife]    = posN[inno,educ,i_s,ife]

    
    elif options.Fertility == 'Exo':
            for i_s in range (0,len(S)):
                for ife in range (0,len(FE_pos)):    
                    for inno in range (0,len(INNO_pos)):
                        for educ in range (0,len(EDUC)):
                            V_2[inno,educ,i_s,ife]  = V[inno,educ,0,i_s,ife]
                            C_2[inno,educ,i_s,ife]  = C[inno,educ,0,i_s,ife]
                            Ck_2[inno,educ,i_s,ife] = Ck[inno,educ,0,i_s,ife]
                            Sp_2[inno,educ,i_s,ife] = Sp[inno,educ,0,i_s,ife]
                            Tp_2[inno,educ,i_s,ife] = Tp[inno,educ,0,i_s,ife]
                            Np_2[inno,educ,i_s,ife] = 0
    
    
    errors = sum(boundgrid[:])/len(C[:])
    if options.timer_on == 'Y':
        if sum(boundgrid[:]) >= len(EDUC)*len(FE_pos)*len(INNO_pos)*len(N):
            print('j : {par.age(j_pos)}, Share of errors (increase grid) ={errors} \n')

    return V_2,C_2,Ck_2,Sp_2,Np_2,Tp_2
