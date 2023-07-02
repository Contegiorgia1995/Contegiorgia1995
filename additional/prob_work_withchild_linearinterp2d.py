# -*- coding: utf-8 -*-
"""
Created on Sat Sep 24 16:09:29 2022

@author: Giorgia
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 17:45:09 2022

@author: Giorgia
"""

import numpy as np
from numba import njit, boolean, int32, double, void
@njit(void(double[:],double[:],double[:,:],double[:],double[:],double[:]),fastmath=True)
def interp_2d_vec(grid1,grid2,value,xi1,xi2,yi):
    """ 2d interpolation for vector of points
        
    Args:

        grid1 (numpy.ndarray): 1d grid
        grid2 (numpy.ndarray): 1d grid
        value (numpy.ndarray): value array (2d)
        xi1 (numpy.ndarray): input vector
        xi2 (numpy.ndarray): input vector
        yi (numpy.ndarray): output vector

    """

    for i in range(yi.size):
        yi[i] = interp_2d(grid1,grid2,value,xi1[i],xi2[i])
        
        
def prob_work_with_child(par,options,Vp,Cp,j_pos):
# Solve household problem when j = Jc+1: children consume at home
# keyboard
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
    
    N         = np.arange(0,4,1)
    
    Sp = np.zeros([len(INNO_pos),len(EDUC),len(N),len(S),len(FE_pos)])
    C  = np.zeros([len(INNO_pos),len(EDUC),len(N),len(S),len(FE_pos)])
    V  = np.zeros([len(INNO_pos),len(EDUC),len(N),len(S),len(FE_pos)])
    Ck = np.zeros([len(INNO_pos),len(EDUC),len(N),len(S),len(FE_pos)])
    Tp = np.zeros([len(INNO_pos),len(EDUC),len(N),len(S),len(FE_pos)])
    boundgrid = np.zeros([len(INNO_pos),len(N),len(FE_pos),len(EDUC)])
    
    ## Without children and phi = 1
    if options.Fertility == 'Endo': # Only need to solve case without children in benchmark case where we allow them to choose endogenous fertility
            i_n   = 0
            
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
                   # interp_2d_vec(grid1,grid2,value,Sp3[:,0],inno3[0,:],yi)
                    #splVp = interp2d(Spgrid, INNO_pos,Vp[:,educ,i_n,:,ife] )
                    
                    
                    for inno  in range (0,len(INNO_pos)):
                        
                   
                        innop_prob = INNOp_prob[inno,:]
                        
                        
                        h   = exp(AGE_PROF[educ][0][j_pos]) * exp(FE[ife])* exp(INNO[inno])
                        
                        Cpp = Cp[:,educ,i_n,:,ife].T
                        
                        if min(np.ravel(Cpp.any())<0 and abs(min(np.ravel(Cpp.any()))))>1e-6:
                            Cpp = max(Cpp,0)
                        ucp = beta*np.reshape(((1+r)),[80,])*np.dot(Cpp**(-gammac), np.reshape(innop_prob,[5,]))

                        dispinc = w*h*(1-Lambda)
                        feasible = ((dispinc + (1+rs)*S - Spgrid[0])>0)
                        
                        not_feasible = (dispinc + (1+rs)*S - Spgrid[0]<=0)
                        
                        C[inno,educ,i_n,feasible,ife],Sp[inno,educ,i_n,feasible,ife],boundgrid[inno,i_n,ife,educ] = EGM(par,ucp.T,dispinc,Spgrid,S[feasible])
                        Ck[inno,educ,i_n,feasible,ife] = np.zeros(len(C[inno,educ,i_n,feasible,ife]))
                        
                        #splVp = spi.interp1d(Spgrid,Vp[inno,educ,i_n,:,ife],axis=0, fill_value="extrapolate")
                        
                        Sp3 = np.repeat(np.reshape(Sp[inno,educ,i_n,feasible,ife],[80,1]),len(INNO_pos),axis = 1)
                        
                        splVp =  interp_2d_vec(Spgrid, INNO_pos, Vp[:,educ,i_n,:,ife],Sp3[:,0],inno3[0,:],yi)
                        
                        
                        #
                        vp = splVp(Sp[inno,educ,i_n,feasible,ife])*innop_prob
                        #vp = np.dot(splVp(Sp3[:,0],inno3[0,:]).T,innop_prob)

                        for j in range(0,len(vp)):
                            if math.isnan(vp[j]):
                                vp[j] = -inf
                        vp = np.sort(vp)
                        for j in range(0,len(vp)):
                            if math.isinf(vp[j]):
                                vp[j] = nan
                        
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
            n             = par.N[i_n]
            for ife in range (0,len(FE_pos)):
                splVp = interp2d(Spgrid, INNO_pos,Vp[:,educ,i_n,:,ife])
                for inno in range (0,len(INNO_pos)):
                    innop_prob      = INNOp_prob[inno,:]
                    h   = exp(AGE_PROF[educ][0][j_pos]) * exp(FE[ife])* exp(INNO[inno])
                    
                    Cpp = np.squeeze(Cp[:,educ,i_n,:,ife].T)
                    if min(np.ravel(Cpp.any())<0 and abs(min(np.ravel(Cpp.any()))))>1e-6:
                        Cpp = max(Cpp,0)
                    ucp = beta*np.reshape(((1+r)),[80,1])*np.dot(Cpp**(-gammac), np.reshape(innop_prob,[5,1]))
                    
                    labor_inc = w*h*(1-Lambda)
                    
                    n_final = n* (par.fam_size/2)
                    
                    ChildCost_0 = ChildCost(par,labor_inc,n_final,0,options)
                    ChildCost_1 = ChildCost(par,labor_inc,n_final,1,options)
                    ChildCost_opt = min(ChildCost_0,ChildCost_1)
                    Tp[inno,educ,i_n,:,ife]   = (ChildCost_1>ChildCost_0)
                    
                    dispinc = labor_inc- ChildCost_opt
                    feasible = (dispinc + (1+rs)*S - Spgrid[0]>0)
                    not_feasible = (dispinc + (1+rs)*S - Spgrid[0]<=0)
                    
                    
                    if options.Ck == 'Yes':
                        
                        C[inno,educ,i_n,feasible,ife],Sp[inno,educ,i_n,feasible,ife],boundgrid[inno,i_n,ife,educ] = GEGM_withchild(par,ucp.T,dispinc,Spgrid,S[feasible],splVp,innop_prob,n_final) 
                        #                        EGM_withchild(par,ucp,dispinc,Spgrid,S,n_final); %%%% To do: Generalized EGM with child
                        f_n     = (lambdan/n_final**(1-gamman))**(1/gammac)
                        Ck[inno,educ,i_n,feasible,ife]     = f_n * C[inno,educ,i_n,feasible,ife]
                        
                        
                        Sp3 = np.repeat(np.reshape(Sp[inno,educ,i_n,:,ife],[len(Spgrid),1]),len(INNO_pos),axis = 1)   #######wrong!! but Sp[] ok                 
                        vp = np.dot(splVp(Sp3[:,0],inno3[0,:]).T,innop_prob)
                        #vp = np.dot((splVp(Sp3[:,0],inno3[0,:]).T),innop_prob)
                        # for j in range(0,len(vp)):
                        #     if math.isnan(vp[j]):
                        #         vp[j] = -inf
                                
                        # for j in range(0,len(vp)):
                        #     if vp[j] == -inf:
                        #         vp[j] = np.sort(vp[j])
                                
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


                        
                    elif options.Ck == 'No':
                        C[inno,educ,i_n,feasible,ife],Sp[inno,educ,i_n,feasible,ife],boundgrid[inno,i_n,ife,educ]= GEGM(par,ucp.T,dispinc,Spgrid,S[feasible],splVp,innop_prob)
                        
                        Ck[inno,educ,i_n,feasible,ife] = np.zeros(len(C[inno,educ,i_n,feasible,ife]))
                        
                        Sp3         = np.repeat(np.reshape(Sp[inno,educ,i_n,:,ife],[len(feasible),1]),len(INNO_pos),axis = 1)    
                        vp          = np.dot(splVp(Sp3[:,0],inno3[0,:]).T,innop_prob)### wrong interpolation
                        
                        for j in range(0,len(vp)):
                            if math.isnan(vp[j]):
                                vp[j] = -inf
                        vp = np.sort(vp)
                        for j in range(0,len(vp)):
                            if math.isinf(vp[j]):
                                vp[j] = nan
                        
                        
                        # for j in range(0,len(C[inno,educ,i_n,feasible,ife])):
                        #     if C[inno,educ,i_n,feasible,ife][j]>=0:
                        #         V[inno,educ,i_n,feasible,ife][j] = C[inno,educ,i_n,feasible,ife][j]**(1-gammac)/(1-gammac)  + beta*vp[feasible][j]
                        #     else:
                        #         V[inno,educ,i_n,feasible,ife][j] = -(10**(5/gammac))**(1-gammac)/(1-gammac) + beta*vp[feasible][j]

                        V[inno,educ,i_n,feasible,ife] = C[inno,educ,i_n,feasible,ife]**(1-gammac)/(1-gammac)*(C[inno,educ,i_n,feasible,ife]>=0) -(10**(5/gammac))**(1-gammac)/(1-gammac)*(C[inno,educ,i_n,feasible,ife]<0) + beta*vp[feasible]

                        
                        C[inno,educ,i_n,not_feasible,ife] = 0
                        Ck[inno,educ,i_n,not_feasible,ife] = 0
                        V[inno,educ,i_n,not_feasible,ife] = -(10**(5/gammac))**(1-gammac)/(1-gammac)
                        Sp[inno,educ,i_n,not_feasible,ife]= Spgrid[0]
                        


    errors = sum(boundgrid[:])/len(C[:])
    if options.timer_on == 'Y':
        if sum(boundgrid[:]) >= len(EDUC)*len(FE_pos)*len(INNO_pos)*len(N):
            print('j : {par.age(j_pos)}, Share of errors (increase grid) ={errors} \n')
    return V,C,Ck,Sp,Tp

####To Do:
    #outcome of V[] is messed up for small or negative numbers in first rows. OFTEN IT SKIPS NUMBERS!!
