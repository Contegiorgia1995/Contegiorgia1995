# -*- coding: utf-8 -*-
"""
Created on Sat Sep 17 12:57:07 2022

@author: Giorgia
"""
import numpy.matlib
def prob_educ_hs_grad(par,options,Vp,Cp):
    
# Solve household problem when it is HS Graduate, for age 12 and 16

# HS graduate:
    educ      = 1
    PSY       = par.psy_val_hs
    r_sav     = par.r_sav
    r_debt    = par.r_debt
    beta      = par.beta
    gammac    = par.gammac
    w         = par.w
    Lambda    = par.Lambda
    
    # for i in range (0, 1):
    #     a = [0]
    #     V[i]= [a]*3
    
    V       = {}  
        
    for i in range (0, 1):
        a = [0]
        V[i]= [a]*3
        V[0][0] = [a]*len(PSY)
########################################CAN'T FIGURE OUT    
    C       = {}  
        
    for i in range (0, 1):
        a = [0]
        C[i]= [a]*3
        C[0][0] = [a]*len(PSY)
        
    # C =  {}
    # for i in range (1, 3):
    #     a = [0]
    #     C[0] = [a] * 100
    #     C[i]= [a]*3
#############################################
        
    Sp       = {} 
        
    for i in range (0, 1):
        a = [0]
        Sp[i]= [a]*3
        Sp[0][0] = [a]*len(PSY)
        
    # for i in range (0, 1):
    #     a = [0]
    #     Sp[i]= [a]*3
        
    FE_pos     = par.inc.fe_pos
    AGE_PROF   = par.inc.age_prof
    INNO_pos   = par.inc.inno_pos
    INNO_posT = np.reshape(INNO_pos,[1,5])
    
    FE         = par.inc.fe[educ]
    
    Vpaux     = np.squeeze(Vp[:,educ,:,:])
    Cpaux     = np.squeeze(Cp[:,educ,:,:])
    
    
    for j_pos in range (par.Je2_pos,par.Je1_pos-1,-1):
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
            #splVp = spi.interp2d(S2, INNO2, np.squeeze(Vpaux[:,:,ife]) )
            splVp = GlobalSpline2D(S2[:,0], INNO2[0,:], np.squeeze(Vpaux[:,:,ife]))
            for inno  in range (0,len(INNO_pos)):
            
                innop_prob = INNOp_prob[inno,:]
                
                h   = exp(AGE_PROF[educ][0][j_pos]) * exp(FE[ife])* exp(INNO[inno])
                
                Cpp           = np.squeeze(Cpaux[:,:,ife]).T
                if min(np.ravel(Cpp.any())<0 and abs(min(np.ravel(Cpp.any()))))>1e-6:
                    Cpp  = max(Cpp,0)

                ucp = beta*np.reshape(((1+r)),[len(r),1])*np.dot(Cpp**(-gammac), np.reshape(innop_prob,[5,1]))
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
        
        Vpaux = V2
        Cpaux = C2

   

    ## Age 16 (age of HS)
    # Expectation about innovation
    exp_prob   = np.squeeze(par.inc.z_prob[educ][:,j_pos,0])
    
    j_pos = par.Je1_pos-1
        
    S          = par.grids[educ][j_pos] # This grids do not change by education (at this age)
    inno3      = np.repeat(INNO_posT,len(S),axis=0)
    
    V2         = np.zeros([len(S),len(FE_pos)])
    C2         = np.zeros([len(S),len(FE_pos)])
    Sp2        = np.zeros([len(S),len(FE_pos)])
    boundgrid  = np.zeros([len(FE_pos)])
    
    Spgrid     = par.grids[educ][j_pos+1]
    r         = (r_sav * (Spgrid>=0) + r_debt * (Spgrid<0)).T
    rs        = (r_sav * (S>=0) + r_debt * (S<0))

    INNO2,S2  = np.meshgrid(INNO_pos,Spgrid) 
    
    INNO      = par.inc.z_val[educ][j_pos,:]
    
    
    for ife in range(0,len(FE_pos)):
        splVp = GlobalSpline2D(S2[:,0], INNO2[0,:], np.squeeze(Vpaux[:,:,ife]))
            
        h   = exp(AGE_PROF[educ][0][j_pos]) * exp(FE[ife])* exp(INNO[inno])
        
        Cpp           = np.squeeze(Cpaux[:,:,ife]).T
        if min(np.ravel(Cpp.any())<0 and abs(min(np.ravel(Cpp.any()))))>1e-6:
            Cpp  = max(Cpp,0)
        
        ucp = beta*np.reshape(((1+r)),[len(r),1])*np.dot(Cpp**(-gammac), np.reshape(exp_prob,[5,1]))
        #dispinc = -par.pe1 + par.init_trans[educ][j_pos]
        dispinc = w*h*(1-Lambda) + par.init_trans[educ][j_pos]
        feasible = (dispinc + (1+rs)*S - Spgrid[0]>0)
        
        not_feasible = (dispinc + (1+rs)*S - Spgrid[0]<=0)
        
        C2[feasible,ife],Sp2[feasible,ife],boundgrid[ife] = GEGM(par,ucp.T,dispinc,Spgrid,S[feasible],splVp,exp_prob)
        
        C2[not_feasible,ife] = 0
        Sp2[not_feasible,ife] = 0
        
        Sp3 = np.repeat(np.reshape(Sp2[:,ife],[len(Sp2[:,ife]),1]),len(INNO_pos),axis = 1)
        vp = np.dot(splVp(Sp3[:,0],inno3[0,:]).T,exp_prob)
        
        V2[feasible,ife]=C2[feasible,ife]**(1-gammac)/(1-gammac) + beta*vp[feasible]
        V2[not_feasible,ife] = -(10**(5/gammac))**(1-gammac)/(1-gammac)
    
    if j_pos >= par.Je1_pos: #only case j_pos = 10,9
        V[0][j_pos - (par.Je1_pos-1)] = V2
        C[0][j_pos - (par.Je1_pos-1)] = C2
        Sp[0][j_pos - (par.Je1_pos-1)] = Sp2
    else: #case j_pos = 8
                                    
        for i in range(0, len(PSY) ):    
            V[0][j_pos - (par.Je1_pos-1)][i] = V2 - np.reshape(np.repeat(np.repeat(PSY[i],size(V2,0)),size(V2,1)),[size(V2,0),size(V2,1)])
            C[0][j_pos - (par.Je1_pos-1)][i] =  C2
            Sp[0][j_pos - (par.Je1_pos-1)][i] =  Sp2

    
    Vpaux = V2
    Cpaux = C2
    
    return V,C, Sp