# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 13:13:07 2022

@author: Giorgia
"""
import numpy as np
import scipy.interpolate as spi

def prob_educ_hs_drop(par,options,Vp,Cp):
    # Solve household problem when it is HS dropout, for age 12 and 16
    # HS dropout:
    educ      = 0
    PSY    = par.psy_val_hs
    r_sav     = par.r_sav
    r_debt    = par.r_debt
    beta      = par.beta
    gammac    = par.gammac
    w         = par.w
    Lambda    = par.Lambda
    
    V       = {}
        
    for i in range (0, 1):
        a = [0]
        V[i]= [a]*3
   
    C       = {}
    # for i in range(0,3):
    #     V[i] = [[]]
        
    for i in range (0, 1):
        a = [0]
        C[i]= [a]*3
        
    Sp       = {}
    # for i in range(0,3):
    #     V[i] = [[]]
        
    for i in range (0, 1):
        a = [0]
        Sp[i]= [a]*3
    
    FE_pos     = par.inc.fe_pos
    AGE_PROF   = par.inc.age_prof
    INNO_pos   = par.inc.inno_pos
    INNO_posT = np.reshape(INNO_pos,[1,5])

    FE         = par.inc.fe[educ]
    
    Vpaux     = np.squeeze(Vp[:,educ,:,:])
    Cpaux     = np.squeeze(Cp[:,educ,:,:])
    
    for j_pos in range (par.Je2_pos,par.Je1_pos-2,-1):
        INNOp_prob = par.inc.z_prob[educ][:,j_pos+1,:].T
        S          = par.grids[educ][j_pos] # This grids do not change by education (at this age)
        inno3      = np.repeat(INNO_posT,len(S),axis=0)
        
        V2         = np.zeros([len(INNO_pos),len(S),len(FE_pos)])
        C2         = np.zeros([len(INNO_pos),len(S),len(FE_pos)])
        Sp2        = np.zeros([len(INNO_pos),len(S),len(FE_pos)])
        boundgrid  = np.zeros([len(INNO_pos),len(FE_pos)])
        
        Spgrid     = par.grids[educ][j_pos+1]
        r         = (r_sav * (Spgrid>=0) + r_debt * (Spgrid<0)).T
        rs        = (r_sav * (S>=0) + r_debt * (S<0))

        INNO2,S2  = np.meshgrid(INNO_pos,Spgrid) 
        
        INNO      = par.inc.z_val[educ][j_pos,:]
        
        for ife in range(0,len(FE_pos)):
            splVp = GlobalSpline2D(S2[:,0], INNO2[0,:], np.squeeze(Vpaux[:,:,ife])) ####bCAREFUL WITH INPUTS HERE!!!!!
            for inno  in range (0,len(INNO_pos)):
            
                innop_prob = INNOp_prob[inno,:]
                
                h   = exp(AGE_PROF[educ][0][j_pos]) * exp(FE[ife])* exp(INNO[inno])
                
                Cpp           = np.squeeze(Cpaux[:,:,ife]).T
                if min(np.ravel(Cpp.any())<0 and abs(min(np.ravel(Cpp.any()))))>1e-6:
                    Cpp  = max(Cpp,0);

                
                
                #ucp = beta*np.reshape(((1+r)),[len(r),])*np.dot(Cpp**(-gammac), np.reshape(innop_prob,[5,]))
                
                ucp = beta*np.reshape(((1+r)),[len(r),1])*np.dot(Cpp**(-gammac), np.reshape(innop_prob,[size(innop_prob),1]))
                #ucp = np.reshape(ucp,[size(ucp),1])  
                dispinc = w*h*(1-Lambda) + par.init_trans[educ][j_pos]
                
                feasible = (dispinc + (1+rs)*S - Spgrid[0]>0)
                
                not_feasible = (dispinc + (1+rs)*S - Spgrid[0]<=0)
                                
                C2[inno,feasible,ife],Sp2[inno,feasible,ife],boundgrid[inno,ife] = GEGM(par,ucp.T,dispinc,Spgrid,S[feasible],splVp,innop_prob)
                
                C2[inno,not_feasible,ife] = 0
                Sp2[inno,not_feasible,ife] = -Spgrid[0]
                
                Sp3 = np.repeat(np.reshape(Sp2[inno,feasible,ife],[len(Sp2[inno,feasible,ife]),1]),len(INNO_pos),axis = 1)
                vp = np.dot(splVp(Sp3[:,0],inno3[0,:]).T,innop_prob)
                
                V2[inno,feasible,ife]=C2[inno,feasible,ife]**(1-gammac)/(1-gammac) + beta*vp[feasible]
                V2[inno,not_feasible,ife] = -(10**(5/gammac))**(1-gammac)/(1-gammac)

        
        errors = sum(boundgrid[:])/len(C2[:])
        if options.timer_on == 'Y':
                if sum(boundgrid[:]) >= len(FE_pos)*len(INNO_pos):
                    print('j : {par.age(j_pos)}, Share of errors (increase grid) ={errors} \n')

        
        V[0][j_pos - (par.Je1_pos-1)] = V2
        C[0][j_pos - (par.Je1_pos-1)] = C2
        Sp[0][j_pos - (par.Je1_pos-1)] = Sp2
        
        Vpaux = copy.deepcopy(V2)
        Cpaux = copy.deepcopy(C2)

    
    # Expectation about innovation
    V0         = np.zeros([len(S),len(FE_pos)])
    exp_prob   = np.squeeze(par.inc.z_prob[educ][:,j_pos,0])
    for i_s in range(0,len(S)):
        for fe in range(0,len(FE_pos)):
            vp        = np.squeeze(V2[:,i_s,fe])
            V0[i_s,fe] = np.dot(np.reshape(vp,[1,len(vp)]),exp_prob)

    return V0,V,C,Sp 
