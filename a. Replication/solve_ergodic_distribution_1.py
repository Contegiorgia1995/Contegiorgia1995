# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 13:49:56 2022

@author: Giorgia
"""#https://github.com/aprsa/ndpolator/blob/main/ndpolator/ndpolator.py
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import interpn
import types, copy
def solve_ergodic_distribution1(par,pol,options):
# Interpolate policy functions into dense grid

    ##%% Distribute policy functions
    S           = copy.deepcopy(pol.S)
    Se          = copy.deepcopy(pol.Se)
    C           = copy.deepcopy(pol.C)
    Ce          = copy.deepcopy(pol.Ce)
    Ck          = copy.deepcopy(pol.Ck)
    Np          = copy.deepcopy(pol.Np)
    PHIp        = copy.deepcopy(pol.PHIp)
    Tp          = copy.deepcopy(pol.Tp)
    tau0        = copy.deepcopy(pol.tau0)
    
    #del pol
    
    
   ## Interpolate policy functions in a dense grid
    # Parameters
    EDUC                            = par.educ
    N                               = par.N
    PSY                             = par.psy_val_hs
    FE_pos                          = par.inc.fe_pos
    FE_pos                          = np.array(FE_pos)
    INNO_pos                        = par.inc.inno_pos
    INNO_pos                        = np.array(INNO_pos)
    
    # Arrays for new policies
    S_dense           = {}
    for i in range (0, par.Jd_pos):
        S_dense[i] = [[]]
    
    

    Se_dense                  = {}
    for i in range(0,3):
         Se_dense[i] = [[],[],[]]
   
    C_dense           = {}
    for i in range (0, par.Jd_pos):
        C_dense[i] = [[]]
   
    Ck_dense           = {}
    for i in range (0, par.Jd_pos):
        Ck_dense[i] = [[]]
    
    Ce_dense                  = {}
    for i in range(0,3):
         Ce_dense[i] = [[],[],[]]   
    # for i in range (0, 3):
    #     a = [0]
    #     Ce_dense[i]= [a]*3
    
    Tp_dense           = {}
    for i in range (0, par.Jd_pos):
        Tp_dense[i] = [[]]
    
    ## 1. Create dense grid
    curv                            = 3
    
    par.grids_dense                 = {}
    for i in range (0, len(EDUC)):
        a = [0]
        par.grids_dense[i]= [a]*par.Jd_pos ##################may need to transpose
    
    for educ in range (0, len(EDUC)):
        par.grids_dense[educ][par.Je1_pos-1] = par.grids[educ][par.Je1_pos-1]

    n_2_s                           = 100
    
    
    for j_pos in range (par.Je1_pos,par.Jd_pos):
        for educ in range (0, len(EDUC)):
            S0                          = par.grids[educ][j_pos]
            s_min                       = abs(np.min(S0))
            s_max                       = abs(np.max(S0))
            par.grids_dense[educ][j_pos] = np.linspace(-s_min**(1/curv),s_max**(1/curv),n_2_s)**curv

    
    ## 2. Interpolate policy function of savings for education
    educ  = 0
    for t in range(0,3):
        S_dense_aux      = par.grids_dense[educ][par.Je1_pos+t-1]
        S_orig           = par.grids[educ][par.Je1_pos+t-1]
        Sp_max           = np.max(par.grids_dense[educ][par.Je1_pos+t])
    
    
        S_pol            = Se[educ][t]
        C_pol            = Ce[educ][t]
        if np.max(S_pol[:]) > Sp_max:### probably will need to se for i in S_pol, S_pol[i] > Sp_max
            j_pos        = par.Je1_pos+t-1
            print(f'Attention: age = {par.age[j_pos]}, educ = {educ}, extrapolation in ergodic distribution, max extrap = {log(np.max(S_pol[:])/Sp_max)}')
            
            S_pol        = np.where(S_pol>Sp_max, Sp_max, S_pol)
     
        Sp_min           = np.min(par.grids_dense[educ][par.Je1_pos+t])
        if np.min(S_pol[:]) < Sp_min:##### or any??--- no , min and max should be ok
            j_pos        = par.Je1_pos+t-1
            print(f'Attention: age ={par.age(j_pos)}, educ = {educ}, extrapolation in ergodic distribution, min extrap = {log(np.min(S_pol[:])/Sp_min)}')
            S_pol        = np.where(S_pol<Sp_min, Sp_min, S_pol)##########################check!!
            
        pol = np.zeros([len(INNO_pos),len(S_dense_aux),len(FE_pos)])
        for i in range (0, len(INNO_pos)):
            for j in range(0, len(S_dense_aux)):
                for k in range (0, len(FE_pos)):
                    pol[i,j,k] = interp_3d(INNO_pos,S_orig,FE_pos,S_pol,INNO_pos[i],S_dense_aux[j],FE_pos[k])
                    
        Se_dense[educ][t] = np.reshape(pol,[len(INNO_pos),len(S_dense_aux),len(FE_pos)])   ### probably no need to reshape         

        
        pol = np.zeros([len(INNO_pos),len(S_dense_aux),len(FE_pos)])
        for i in range (0, len(INNO_pos)):
            for j in range(0, len(S_dense_aux)):
                for k in range (0, len(FE_pos)):
                    pol[i,j,k] = interp_3d(INNO_pos,S_orig,FE_pos,C_pol,INNO_pos[i],S_dense_aux[j],FE_pos[k])
        Ce_dense[educ][t] = np.reshape(pol,[len(INNO_pos),len(S_dense_aux),len(FE_pos)])  
    
        
    
    educ  = 1
    for t in range(0,3):
        S_dense_aux      = par.grids_dense[educ][par.Je1_pos+t-1]
        S_orig           = par.grids[educ][par.Je1_pos+t-1]
        Sp_max           = np.max(par.grids_dense[educ][par.Je1_pos+t])
        
        S_pol            = Se[educ][t]
        C_pol            = Ce[educ][t]
        if np.max(np.ravel(S_pol[:])) > Sp_max:
            j_pos        = par.Je1_pos+t-1
            print(f'Attention: age = {par.age[j_pos]}, educ = {educ}, extrapolation in ergodic distribution, max extrap = {log(np.max(S_pol[:])/Sp_max)}')
            S_pol        = np.where(S_pol>Sp_max, Sp_max, S_pol)
        
        Sp_min           = min(par.grids_dense[educ][par.Je1_pos+t])
        if min(np.ravel(S_pol[:])) < Sp_min:
            j_pos        = par.Je1_pos+t-1;
            print(f'Attention: age ={par.age[j_pos]}, educ = {educ}, extrapolation in ergodic distribution, min extrap = {log(np.min(S_pol[:])/Sp_min)}')
            S_pol        = np.where(S_pol<Sp_min, Sp_min, S_pol)################until here
        
        if t <= 0:
            psy_orig         = PSY
            S1,FE1,PSY1    = np.meshgrid(S_orig,FE_pos,psy_orig)
            S2,FE2,PSY2   = np.meshgrid(S_dense_aux,FE_pos,psy_orig)
            
            pol = np.zeros([len(psy_orig),len(S_dense_aux),len(FE_pos)])
            for i in range (0, len(psy_orig)):
                for j in range(0, len(S_dense_aux)):
                    for k in range (0, len(FE_pos)):
                        pol[i,j,k] = interp_3d(psy_orig,S_orig,FE_pos,S_pol,psy_orig[i],S_dense_aux[j],FE_pos[k])
                        
            Se_dense[educ][t] = pol#np.reshape(pol,[len(INNO_pos),len(S_dense_aux),len(FE_pos)])
            #Se_dense[educ][t][:,:,:]    = reshape(pol(S2,FE2,PSY2), length(S_dense_aux),length(FE_pos),length(psy_orig))
            
            
            pol = np.zeros([len(psy_orig),len(S_dense_aux),len(FE_pos)])
            for i in range (0, len(psy_orig)):
                for j in range(0, len(S_dense_aux)):
                    for k in range (0, len(FE_pos)):
                        pol[i,j,k] = interp_3d(psy_orig,S_orig,FE_pos,C_pol,psy_orig[i],S_dense_aux[j],FE_pos[k])
            
            Ce_dense[educ][t] = pol
            
        else:
            S1,FE1,INNO1   = np.meshgrid(S_orig,FE_pos,INNO_pos)
            S2,FE2,INNO2   = np.meshgrid(S_dense_aux,FE_pos,INNO_pos)
            
            pol = np.zeros([len(INNO_pos),len(S_dense_aux),len(FE_pos)])
            for i in range (0, len(INNO_pos)):
                for j in range(0, len(S_dense_aux)):
                    for k in range (0, len(FE_pos)):
                        pol[i,j,k] = interp_3d(INNO_pos,S_orig,FE_pos,S_pol,INNO_pos[i],S_dense_aux[j],FE_pos[k])
            Se_dense[educ][t] = pol
            
            pol = np.zeros([len(INNO_pos),len(S_dense_aux),len(FE_pos)])
            for i in range (0, len(INNO_pos)):
                for j in range(0, len(S_dense_aux)):
                    for k in range (0, len(FE_pos)):
                        pol[i,j,k] = interp_3d(INNO_pos,S_orig,FE_pos,C_pol,INNO_pos[i],S_dense_aux[j],FE_pos[k])
            Ce_dense[educ][t] = pol

        
    educ  = 2
    for t in range(0,3):
        S_dense_aux      = par.grids_dense[educ][par.Je1_pos+t-1]
        S_orig           = par.grids[educ][par.Je1_pos+t-1]
        Sp_max           = np.max(par.grids_dense[educ][par.Je1_pos+t])
        
        psy_orig             = PSY
        
        S_pol            = Se[educ][t]
        C_pol            = Ce[educ][t]
        
        if np.max(np.ravel(S_pol[:])) > Sp_max:
            j_pos        = par.Je1_pos+t-1
            print(f'Attention: age = {par.age[j_pos]}, educ = {educ}, extrapolation in ergodic distribution, max extrap = {log(np.max(S_pol[:])/Sp_max)}')
            S_pol        = np.where(S_pol>Sp_max, Sp_max, Sp_pol)
        
        Sp_min           = np.min(par.grids_dense[educ][par.Je1_pos+t])
        if np.min(np.ravel(S_pol[:])) < Sp_min:
            j_pos        = par.Je1_pos+t-1;
            print(f'Attention: age ={par.age[j_pos]}, educ = {educ}, extrapolation in ergodic distribution, min extrap = {log(np.min(S_pol[:])/Sp_min)}')
            S_pol = np.where(S_pol<Sp_min, Sp_min, S_pol)
    
        if t <= 0:
            psy_orig         = PSY
            S1,FE1,PSY1    = np.meshgrid(S_orig,FE_pos,psy_orig)
            S2,FE2,PSY2   = np.meshgrid(S_dense_aux,FE_pos,psy_orig)
            
            pol = np.zeros([len(psy_orig),len(S_dense_aux),len(FE_pos)])
            for i in range (0, len(psy_orig)):
                for j in range(0, len(S_dense_aux)):
                    for k in range (0, len(FE_pos)):
                        pol[i,j,k] = interp_3d(psy_orig,S_orig,FE_pos,S_pol,psy_orig[i],S_dense_aux[j],FE_pos[k])
                        
            Se_dense[educ][t] = pol#np.reshape(pol,[len(INNO_pos),len(S_dense_aux),len(FE_pos)])
            #Se_dense[educ][t][:,:,:]    = reshape(pol(S2,FE2,PSY2), length(S_dense_aux),length(FE_pos),length(psy_orig))
            
            
            pol = np.zeros([len(psy_orig),len(S_dense_aux),len(FE_pos)])
            for i in range (0, len(psy_orig)):
                for j in range(0, len(S_dense_aux)):
                    for k in range (0, len(FE_pos)):
                        pol[i,j,k] = interp_3d(psy_orig,S_orig,FE_pos,C_pol,psy_orig[i],S_dense_aux[j],FE_pos[k])
            
            Ce_dense[educ][t] = pol
            
            
        else:
            # S1,FE1   = np.meshgrid(S_orig,FE_pos)
            # S2,FE2   = np.meshgrid(S_dense_aux)
            pol = np.zeros([len(S_dense_aux),len(FE_pos)])
            
            for j in range(0, len(S_dense_aux)):
                for k in range (0, len(FE_pos)):
                    pol[j,k] = interp_2d(S_orig,FE_pos,S_pol,S_dense_aux[j],FE_pos[k])
            
            Se_dense[educ][t]= pol
            
            
            pol = np.zeros([len(S_dense_aux),len(FE_pos)])
            for j in range(0, len(S_dense_aux)):
                for k in range (0, len(FE_pos)):
                    pol[j,k] = interp_2d(S_orig,FE_pos,C_pol,S_dense_aux[j],FE_pos[k])
            
            Ce_dense[educ][t] = pol
        
        #del Se
    
   ## 3. Interpolate policy function of savings for working young
    for j_pos in range (par.Je2_pos+1,par.Jc_pos-1):
        S_dense[j_pos]        = np.zeros([len(INNO_pos),len(EDUC),len(par.grids_dense[0][j_pos]),len(FE_pos)])
        S_dense[j_pos][:] = np.nan
        C_dense[j_pos]        = np.zeros([len(INNO_pos),len(EDUC),len(par.grids_dense[0][j_pos]),len(FE_pos)])
        C_dense[j_pos][:] = np.nan

        for educ in range (0,len(EDUC)):
            
            S_dense_aux      = par.grids_dense[educ][j_pos]
            S_orig           = par.grids[educ][j_pos]#### to transpose
            
            Sp_max           = np.max(par.grids_dense[educ][j_pos+1])
            S_pol            = np.squeeze(S[j_pos][:,educ,:,:])
            C_pol            = np.squeeze(C[j_pos][:,educ,:,:])
            
            if np.max(np.ravel(S_pol[:])) > Sp_max:
                print(f'Attention: age = {par.age[j_pos]}, educ = {educ}, extrapolation in ergodic distribution, max extrap = {log(np.max(S_pol[:])/Sp_max)}')
                S_pol        = np.where(S_pol>Sp_max, Sp_max, S_pol)
    
            
            Sp_min           = np.min(par.grids_dense[educ][j_pos+1])
            if np.min(np.ravel(S_pol[:])) > Sp_min:
                print(f'Attention: age = {par.age[j_pos]}, educ = {educ}, extrapolation in ergodic distribution, min extrap = {log(np.min(S_pol[:])/Sp_min)}')
                S_pol = np.where(S_pol<Sp_min, Sp_min, S_pol)
    
            
            # S1,FE1,INNO1   = np.meshgrid(S_orig,FE_pos,INNO_pos)
            # S2,FE2,INNO2   = np.meshgrid(S_dense_aux,FE_pos,INNO_pos)
            
            pol = np.zeros([len(INNO_pos),len(S_dense_aux),len(FE_pos)])
            for i in range (0, len(INNO_pos)):
                for j in range(0, len(S_dense_aux)):
                    for k in range (0, len(FE_pos)):
                        pol[i,j,k] = interp_3d(INNO_pos,S_orig,FE_pos,S_pol,INNO_pos[i],S_dense_aux[j],FE_pos[k])
                        
                        
            S_dense[j_pos][:,educ,:,:]  = pol
            
            pol = np.zeros([len(INNO_pos),len(S_dense_aux),len(FE_pos)])
            for i in range (0, len(INNO_pos)):
                for j in range(0, len(S_dense_aux)):
                    for k in range (0, len(FE_pos)):
                        pol[i,j,k] = interp_3d(INNO_pos,S_orig,FE_pos,C_pol,INNO_pos[i],S_dense_aux[j],FE_pos[k])
            
            C_dense[j_pos][:,educ,:,:]  = pol
        

    
    ## 4. Interpolate policy function of savings for fertility
    j_pos                   = par.Jc_pos-1
    S_dense[j_pos]          = np.zeros([len(INNO_pos),len(EDUC),len(par.grids_dense[0][j_pos]),len(FE_pos)])
    S_dense[j_pos][:]       = np.nan
    C_dense[j_pos]          = np.zeros([len(INNO_pos),len(EDUC),len(par.grids_dense[0][j_pos]),len(FE_pos)])
    C_dense[j_pos][:]       = np.nan
    Ck_dense[j_pos]         = np.zeros([len(INNO_pos),len(EDUC),len(par.grids_dense[0][j_pos]),len(FE_pos)])
    Ck_dense[j_pos][:]      = np.nan
    Np_dense                = np.zeros([len(INNO_pos),len(EDUC),len(par.grids_dense[0][j_pos]),len(FE_pos)])
    Np_dense[:]             = np.nan
    Tp_dense[j_pos]         = np.zeros([len(INNO_pos),len(EDUC),len(par.grids_dense[0][j_pos]),len(FE_pos)])
    Tp_dense[j_pos][:]        = np.nan
    
    for educ in range(0,len(EDUC)):      
        S_dense_aux         = par.grids_dense[educ][j_pos]
        S_orig              = par.grids[educ][j_pos]
        Sp_max              = np.max(par.grids_dense[educ][j_pos+1])
    
        S_pol               = np.squeeze(S[j_pos][:,educ,:,:])
        C_pol               = np.squeeze(C[j_pos][:,educ,:,:])
        Ck_pol              = np.squeeze(Ck[j_pos][:,educ,:,:])
        
        if np.max(np.ravel(S_pol[:])) > Sp_max:
            print(f'Attention: age = {par.age[j_pos]}, educ = {educ}, extrapolation in ergodic distribution, min extrap = {log(np.min(S_pol[:])/Sp_min)}')
            S_pol           = np.where(S_pol>Sp_max, Sp_max, S_pol)
        
            Sp_min          = np.min(par.grids_dense[educ][j_pos+1])
        
        if np.min(S_pol[:]) > Sp_min:
            print(f'Attention: age = {par.age[j_pos]}, educ = {educ}, extrapolation in ergodic distribution, min extrap = {log(np.min(S_pol[:])/Sp_min)}')
            S_pol           = np.where(S_pol<Sp_min, Sp_min, S_pol)
            
        pol = np.zeros([len(INNO_pos),len(S_dense_aux),len(FE_pos)])
        for i in range (0, len(INNO_pos)):
            for j in range(0, len(S_dense_aux)):
                for k in range (0, len(FE_pos)):
                    pol[i,j,k] = interp_3d(INNO_pos,S_orig,FE_pos,S_pol,INNO_pos[i],S_dense_aux[j],FE_pos[k])
                    
                    
        S_dense[j_pos][:,educ,:,:]  = pol
        
        pol = np.zeros([len(INNO_pos),len(S_dense_aux),len(FE_pos)])
        for i in range (0, len(INNO_pos)):
            for j in range(0, len(S_dense_aux)):
                for k in range (0, len(FE_pos)):
                    pol[i,j,k] = interp_3d(INNO_pos,S_orig,FE_pos,C_pol,INNO_pos[i],S_dense_aux[j],FE_pos[k])
        
        C_dense[j_pos][:,educ,:,:]  = pol
        
        
        pol = np.zeros([len(INNO_pos),len(S_dense_aux),len(FE_pos)])
        for i in range (0, len(INNO_pos)):
            for j in range(0, len(S_dense_aux)):
                for k in range (0, len(FE_pos)):
                    pol[i,j,k] = interp_3d(INNO_pos,S_orig,FE_pos,Ck_pol,INNO_pos[i],S_dense_aux[j],FE_pos[k])
        
        Ck_dense[j_pos][:,educ,:,:]  = pol
        
        Np_pol                      = np.squeeze(Np[:,educ,:,:])
        pol = np.zeros([len(INNO_pos),len(S_dense_aux),len(FE_pos)])
        for i in range (0, len(INNO_pos)):
            for j in range(0, len(S_dense_aux)):
                for k in range (0, len(FE_pos)):
                    points = (INNO_pos,S_orig,FE_pos)
                    point = np.array([INNO_pos[i],S_dense_aux[j],FE_pos[k]])
                    pol[i,j,k] = int(interpn(points,Np_pol,point,method = "nearest"))###################################careful to check if results are the same, here I use nearest neighboorhood extrapolation, and int() to round numbers, but matlba may want to have a nearest neightborhood interpolation
        
        
        Np_dense[:,educ,:,:]        = pol
        
        Tp_pol                      = np.squeeze(Tp[j_pos][:,educ,:,:])
        pol = np.zeros([len(INNO_pos),len(S_dense_aux),len(FE_pos)])
        for i in range (0, len(INNO_pos)):
            for j in range(0, len(S_dense_aux)):
                for k in range (0, len(FE_pos)):
                    points = (INNO_pos,S_orig,FE_pos)
                    point = np.array([INNO_pos[i],S_dense_aux[j],FE_pos[k]])
                    pol[i,j,k] = (interpn(points,Tp_pol,point,method = "nearest"))###################################careful to check if results are the same, here I use nearest neighboorhood extrapolation, and int() to round numbers, but matlba may want to have a nearest neightborhood interpolation
                    
        Tp_dense[j_pos][:,educ,:,:]  = pol
        
        ## 5. Interpolate policy function of savings for work with child
    for j_pos  in range  (par.Jc_pos,par.Jc_pos+par.Je1_pos-3):
        
        S_dense[j_pos]          = np.zeros([len(INNO_pos),len(EDUC),len(N),len(par.grids_dense[0][j_pos]),len(FE_pos)])
        S_dense[j_pos][:]       = np.nan
        C_dense[j_pos]          = np.zeros([len(INNO_pos),len(EDUC),len(N),len(par.grids_dense[0][j_pos]),len(FE_pos)])
        C_dense[j_pos][:]       = np.nan
        Ck_dense[j_pos]         = np.zeros([len(INNO_pos),len(EDUC),len(N),len(par.grids_dense[0][j_pos]),len(FE_pos)])
        Ck_dense[j_pos][:]      = np.nan
        Tp_dense[j_pos]         = np.zeros([len(INNO_pos),len(EDUC),len(N),len(par.grids_dense[0][j_pos]),len(FE_pos)])
        Tp_dense[j_pos][:]        = np.nan
        
        for educ in range(0,len(EDUC)):
            S_dense_aux       = par.grids_dense[educ][j_pos]
            S_orig            = par.grids[educ][j_pos]
            Sp_max            = np.max(par.grids_dense[educ][j_pos+1])
                
            S_pol             = np.squeeze(S[j_pos][:,educ,:,:,:])
            C_pol             = np.squeeze(C[j_pos][:,educ,:,:,:])
            Ck_pol            = np.squeeze(Ck[j_pos][:,educ,:,:,:])
            Tp_pol            = np.squeeze(Tp[j_pos][:,educ,:,:,:])
            
            if np.max(ravel(S_pol[:])) > Sp_max:
                print(f'Attention: age = {par.age[j_pos]}, educ = {educ}, extrapolation in ergodic distribution, min extrap = {log(np.max(S_pol[:])/Sp_max)}')
                S_pol           = np.where(S_pol>Sp_max, Sp_max, S_pol)
            
                Sp_min           = np.min(par.grids_dense[educ][j_pos+1])
            if np.min(np.ravel(S_pol[:])) > Sp_min:
                print(f'Attention: age = {par.age[j_pos]}, educ = {educ}, extrapolation in ergodic distribution, min extrap = {log(np.min(S_pol[:])/Sp_min)}')
                S_pol           = np.where(S_pol<Sp_min, Sp_min, S_pol)
        
            if options.Fertility == 'Endo':
                    
                pol = np.zeros([len(INNO_pos),len(N),len(S_dense_aux),len(FE_pos)])
                for i in range (0, len(INNO_pos)):
                    for n in range (0,len(N)):
                        for j in range(0, len(S_dense_aux)):
                            for k in range (0, len(FE_pos)):
                                pol[i,n,j,k] = interp_4d(INNO_pos,N,S_orig,FE_pos,S_pol,INNO_pos[i],N[n],S_dense_aux[j],FE_pos[k])
                            
                            
                S_dense[j_pos][:,educ,:,:,:]  = pol
                
                pol = np.zeros([len(INNO_pos),len(N),len(S_dense_aux),len(FE_pos)])
                for i in range (0, len(INNO_pos)):
                    for n in range (0,len(N)):
                        for j in range(0, len(S_dense_aux)):
                            for k in range (0, len(FE_pos)):
                                pol[i,n,j,k] = interp_4d(INNO_pos,N,S_orig,FE_pos,C_pol,INNO_pos[i],N[n],S_dense_aux[j],FE_pos[k])
                            
                C_dense[j_pos][:,educ,:,:,:]   = pol
                
                pol = np.zeros([len(INNO_pos),len(N),len(S_dense_aux),len(FE_pos)])
                for i in range (0, len(INNO_pos)):
                    for n in range (0,len(N)):
                        for j in range(0, len(S_dense_aux)):
                            for k in range (0, len(FE_pos)):
                                pol[i,n,j,k] = interp_4d(INNO_pos,N,S_orig,FE_pos,Ck_pol,INNO_pos[i],N[n],S_dense_aux[j],FE_pos[k])

                Ck_dense[j_pos][:,educ,:,:,:]  = pol
                
                pol = np.zeros([len(INNO_pos),len(N),len(S_dense_aux),len(FE_pos)])
                for i in range (0, len(INNO_pos)):
                    for n in range (0,len(N)):
                        for j in range(0, len(S_dense_aux)):
                            for k in range (0, len(FE_pos)):
                                points = (INNO_pos,N,S_orig,FE_pos)
                                point = np.array([INNO_pos[i],N[n],S_dense_aux[j],FE_pos[k]])
                                pol[i,n,j,k] = (interpn(points,Tp_pol,point,method = "nearest"))
                             
                Tp_dense[j_pos][:,educ,:,:,:]  = pol
             
            elif options.Fertility == 'Exo':  ##########################################################################     
                pol = np.zeros([len(INNO_pos),len(S_dense_aux),len(FE_pos)])
                for i in range (0, len(INNO_pos)):
                    for j in range(0, len(S_dense_aux)):
                        for k in range (0, len(FE_pos)):
                            pol[i,j,k] = interp_3d(INNO_pos,S_orig,FE_pos,S_pol,INNO_pos[i],S_dense_aux[j],FE_pos[k])
                S_dense[j_pos][:,educ,0,:,:]  = pol
                
                pol = np.zeros([len(INNO_pos),len(S_dense_aux),len(FE_pos)])
                for i in range (0, len(INNO_pos)):
                    for j in range(0, len(S_dense_aux)):
                        for k in range (0, len(FE_pos)):
                            pol[i,j,k] = interp_3d(INNO_pos,S_orig,FE_pos,C_pol,INNO_pos[i],S_dense_aux[j],FE_pos[k])
                C_dense[j_pos][:,educ,0,:,:] = pol
                
                pol = np.zeros([len(INNO_pos),len(S_dense_aux),len(FE_pos)])
                for i in range (0, len(INNO_pos)):
                    for j in range(0, len(S_dense_aux)):
                        for k in range (0, len(FE_pos)):
                            pol[i,j,k] = interp_3d(INNO_pos,S_orig,FE_pos,Ck_pol,INNO_pos[i],S_dense_aux[j],FE_pos[k])
                Ck_dense[j_pos][:,educ,0,:,:]   = pol
                
                Tp_pol                      = np.squeeze(Tp[j_pos][:,educ,:,:])
                
                pol = np.zeros([len(INNO_pos),len(S_dense_aux),len(FE_pos)])
                for i in range (0, len(INNO_pos)):
                    for j in range(0, len(S_dense_aux)):
                        for k in range (0, len(FE_pos)):
                            points = (INNO_pos,S_orig,FE_pos)
                            point = np.array([INNO_pos[i],S_dense_aux[j],FE_pos[k]])
                            pol[i,j,k] = (interpn(points,Tp_pol,point,method = "nearest"))
                            
                Tp_dense[j_pos][:,educ,0,:,:]   = reshape(pol(S2,FE2,INNO2),length(S_dense_aux),length(FE_pos),1,length(INNO_pos))
                    
        
            ## 6. Interpolate policy function of savings for work + transfers to child
            j_pos                       = par.Jc_pos+par.Je1_pos-3
            
            S_dense[j_pos]          = np.zeros([len(INNO_pos),len(EDUC),len(N),len(par.grids_dense[0][j_pos]),len(FE_pos)])
            S_dense[j_pos][:]       = np.nan
            C_dense[j_pos]          = np.zeros([len(INNO_pos),len(EDUC),len(N),len(par.grids_dense[0][j_pos]),len(FE_pos)])
            C_dense[j_pos][:]       = np.nan
            Ck_dense[j_pos]         = np.zeros([len(INNO_pos),len(EDUC),len(N),len(par.grids_dense[0][j_pos]),len(FE_pos)])
            Ck_dense[j_pos][:]      = np.nan
            Tp_dense[j_pos]         = np.zeros([len(INNO_pos),len(EDUC),len(N),len(par.grids_dense[0][j_pos]),len(FE_pos)])
            Tp_dense[j_pos][:]      = np.nan
            PHIp_dense              = np.zeros([len(INNO_pos),len(EDUC),len(N),len(par.grids_dense[0][j_pos]),len(FE_pos)])
            PHIp_dense[:]           = np.nan
            
            for educ in range (0,len(EDUC)):
                S_dense_aux     = par.grids_dense[educ][j_pos]
                S_orig          = par.grids[educ][j_pos]
                Sp_max          = np.max(par.grids_dense[educ][j_pos+1])
                
            
                S_pol           = np.squeeze(S[j_pos][:,educ,:,:,:])
                C_pol           = np.squeeze(C[j_pos][:,educ,:,:,:])
                Ck_pol          = np.squeeze(Ck[j_pos][:,educ,:,:,:])
                Tp_pol          = np.squeeze(Tp[j_pos][:,educ,:,:,:])
                PHIp_pol        = np.squeeze(PHIp[:,educ,:,:,:])
                
                #     aux               = reshape(1*(H_dense_aux<= CUTOFFS(educ,1)) + ...
                #         2.*(H_dense_aux> CUTOFFS(educ,1)).*(H_dense_aux<= CUTOFFS(educ,2)) +...
                #         3.*(H_dense_aux> CUTOFFS(educ,2)),1,length(H_dense_aux));
                #     Grs_dense{j_pos}(:,:,educ,:)  = repmat(aux,length(par.grids_dense{1,j_pos}),1,1,length(N));

    
                if np.max(ravel(S_pol[:])) > Sp_max:
                    print(f'Attention: age = {par.age[j_pos]}, educ = {educ}, extrapolation in ergodic distribution, min extrap = {log(np.max(S_pol[:])/Sp_max)}')
                    S_pol           = np.where(S_pol>Sp_max, Sp_max, S_pol)
                
                    Sp_min           = np.min(par.grids_dense[educ][j_pos+1])
                if np.min(np.ravel(S_pol[:])) > Sp_min:
                    print(f'Attention: age = {par.age[j_pos]}, educ = {educ}, extrapolation in ergodic distribution, min extrap = {log(np.min(S_pol[:])/Sp_min)}')
                    S_pol           = np.where(S_pol<Sp_min, Sp_min, S_pol)
    

                if options.Fertility == 'Endo':
                    
                    pol = np.zeros([len(INNO_pos),len(N),len(S_dense_aux),len(FE_pos)])
                    for i in range (0, len(INNO_pos)):
                        for n in range (0,len(N)):
                            for j in range(0, len(S_dense_aux)):
                                for k in range (0, len(FE_pos)):
                                    pol[i,n,j,k] = interp_4d(INNO_pos,N,S_orig,FE_pos,S_pol,INNO_pos[i],N[n],S_dense_aux[j],FE_pos[k])
                    S_dense[j_pos][:,educ,:,:,:]  = pol
                    
                    pol = np.zeros([len(INNO_pos),len(N),len(S_dense_aux),len(FE_pos)])
                    for i in range (0, len(INNO_pos)):
                        for n in range (0,len(N)):
                            for j in range(0, len(S_dense_aux)):
                                for k in range (0, len(FE_pos)):
                                    pol[i,n,j,k] = interp_4d(INNO_pos,N,S_orig,FE_pos,C_pol,INNO_pos[i],N[n],S_dense_aux[j],FE_pos[k])
                                
                    C_dense[j_pos][:,educ,:,:,:]   = pol
                    
                    pol = np.zeros([len(INNO_pos),len(N),len(S_dense_aux),len(FE_pos)])
                    for i in range (0, len(INNO_pos)):
                        for n in range (0,len(N)):
                            for j in range(0, len(S_dense_aux)):
                                for k in range (0, len(FE_pos)):
                                    pol[i,n,j,k] = interp_4d(INNO_pos,N,S_orig,FE_pos,Ck_pol,INNO_pos[i],N[n],S_dense_aux[j],FE_pos[k])

                    Ck_dense[j_pos][:,educ,:,:,:]  = pol
                    
                    pol = np.zeros([len(INNO_pos),len(N),len(S_dense_aux),len(FE_pos)])
                    for i in range (0, len(INNO_pos)):
                        for n in range (0,len(N)):
                            for j in range(0, len(S_dense_aux)):
                                for k in range (0, len(FE_pos)):
                                    points = (INNO_pos,N,S_orig,FE_pos)
                                    point = np.array([INNO_pos[i],N[n],S_dense_aux[j],FE_pos[k]])
                                    pol[i,n,j,k] = (interpn(points,Tp_pol,point,method = "nearest"))
                                 
                    Tp_dense[j_pos][:,educ,:,:,:]  = pol
                    
                    pol = np.zeros([len(INNO_pos),len(N),len(S_dense_aux),len(FE_pos)])
                    for i in range (0, len(INNO_pos)):
                        for n in range (0,len(N)):
                            for j in range(0, len(S_dense_aux)):
                                for k in range (0, len(FE_pos)):
                                    points = (INNO_pos,N,S_orig,FE_pos)
                                    point = np.array([INNO_pos[i],N[n],S_dense_aux[j],FE_pos[k]])
                                    pol[i,n,j,k] = (interpn(points,PHIp_pol,point,method = "nearest"))
                                 
                    PHIp_dense[:,educ,:,:,:]      = pol
                        
                elif options.Fertility == 'Exo':
                    pol = np.zeros([len(INNO_pos),len(S_dense_aux),len(FE_pos)])
                    for i in range (0, len(INNO_pos)):
                        for j in range(0, len(S_dense_aux)):
                            for k in range (0, len(FE_pos)):
                                pol[i,j,k] = interp_3d(INNO_pos,S_orig,FE_pos,S_pol,INNO_pos[i],S_dense_aux[j],FE_pos[k])
                    S_dense[j_pos][:,educ,0,:,:]  = pol
                    
                    pol = np.zeros([len(INNO_pos),len(S_dense_aux),len(FE_pos)])
                    for i in range (0, len(INNO_pos)):
                        for j in range(0, len(S_dense_aux)):
                            for k in range (0, len(FE_pos)):
                                pol[i,j,k] = interp_3d(INNO_pos,S_orig,FE_pos,C_pol,INNO_pos[i],S_dense_aux[j],FE_pos[k])
                    C_dense[j_pos][:,educ,0,:,:] = pol
                    
                    pol = np.zeros([len(INNO_pos),len(S_dense_aux),len(FE_pos)])
                    for i in range (0, len(INNO_pos)):
                        for j in range(0, len(S_dense_aux)):
                            for k in range (0, len(FE_pos)):
                                pol[i,j,k] = interp_3d(INNO_pos,S_orig,FE_pos,Ck_pol,INNO_pos[i],S_dense_aux[j],FE_pos[k])
                    Ck_dense[j_pos][:,educ,0,:,:]   = pol
                    
                    Tp_pol                      = np.squeeze(Tp[j_pos][:,educ,:,:])
                    
                    pol = np.zeros([len(INNO_pos),len(S_dense_aux),len(FE_pos)])
                    for i in range (0, len(INNO_pos)):
                        for j in range(0, len(S_dense_aux)):
                            for k in range (0, len(FE_pos)):
                                points = (INNO_pos,S_orig,FE_pos)
                                point = np.array([INNO_pos[i],S_dense_aux[j],FE_pos[k]])
                                pol[i,j,k] = (interpn(points,Tp_pol,point,method = "nearest"))
                                
                    Tp_dense[j_pos][:,educ,0,:,:]   = pol                   

                    pol                         = np.zeros([len(INNO_pos),len(S_dense_aux),len(FE_pos)])
                    for i in range (0, len(INNO_pos)):
                        for j in range(0, len(S_dense_aux)):
                            for k in range (0, len(FE_pos)):
                                points = (INNO_pos,S_orig,FE_pos)
                                point = np.array([INNO_pos[i],S_dense_aux[j],FE_pos[k]])
                                pol[i,j,k] = (interpn(points,PHIp_pol,point,method = "nearest"))
                                
                    PHIp_dense[:,educ,0,:,:]        = pol

        
        ## 6. Interpolate policy function of savings for work old
        for j_pos       in range (par.Jc_pos+par.Je1_pos-2, par.Jr_pos-1):
            S_dense[j_pos]          = np.zeros([len(INNO_pos),len(EDUC),len(par.grids_dense[0][j_pos]),len(FE_pos)])
            S_dense[j_pos][:]       = np.nan
            C_dense[j_pos]          = np.zeros([len(INNO_pos),len(EDUC),len(par.grids_dense[0][j_pos]),len(FE_pos)])
            C_dense[j_pos][:]       = np.nan
            for educ in range  (0,len(EDUC)):
                S_dense_aux = par.grids_dense[educ][j_pos]
                S_orig      = par.grids[educ][j_pos]
                Sp_max      = np.max(par.grids_dense[educ][j_pos+1])
            
                S_pol       = np.squeeze(S[j_pos][:,educ,:,:])
                C_pol       = np.squeeze(C[j_pos][:,educ,:,:])
                
                if np.max(ravel(S_pol[:])) > Sp_max:
                    print(f'Attention: age = {par.age[j_pos]}, educ = {educ}, extrapolation in ergodic distribution, min extrap = {log(np.max(S_pol[:])/Sp_max)}')
                    S_pol           = np.where(S_pol>Sp_max, Sp_max, S_pol)
                
                    Sp_min           = np.min(par.grids_dense[educ][j_pos+1])
                if np.min(np.ravel(S_pol[:])) > Sp_min:
                    print(f'Attention: age = {par.age[j_pos]}, educ = {educ}, extrapolation in ergodic distribution, min extrap = {log(np.min(S_pol[:])/Sp_min)}')
                    S_pol           = np.where(S_pol<Sp_min, Sp_min, S_pol)
                                  
                    pol = np.zeros([len(INNO_pos),len(S_dense_aux),len(FE_pos)])
                    for i in range (0, len(INNO_pos)):
                        for j in range(0, len(S_dense_aux)):
                            for k in range (0, len(FE_pos)):
                                pol[i,j,k] = interp_3d(INNO_pos,S_orig,FE_pos,S_pol,INNO_pos[i],S_dense_aux[j],FE_pos[k])
                                
                    S_dense[j_pos][:,educ,:,:] = pol        

                    
                    pol = np.zeros([len(INNO_pos),len(S_dense_aux),len(FE_pos)])
                    for i in range (0, len(INNO_pos)):
                        for j in range(0, len(S_dense_aux)):
                            for k in range (0, len(FE_pos)):
                                pol[i,j,k] = interp_3d(INNO_pos,S_orig,FE_pos,C_pol,INNO_pos[i],S_dense_aux[j],FE_pos[k])
                    C_dense[j_pos][:,educ,:,:] = pol 


        
        ## 8. Retirement
        for j_pos       in range (par.Jr_pos-1,par.Jd_pos):
            S_dense[j_pos]          = np.zeros([len(EDUC),len(par.grids_dense[0][j_pos]),len(FE_pos)])
            S_dense[j_pos][:]       = np.nan
            C_dense[j_pos]          = np.zeros([len(EDUC),len(par.grids_dense[0][j_pos]),len(FE_pos)])
            C_dense[j_pos][:]       = np.nan
            
            for educ in range (0,len(EDUC)):
                S_dense_aux = par.grids_dense[educ][j_pos]
                S_orig      = par.grids[educ][j_pos]
                j_pos_1     = min(j_pos,par.Jd_pos-1)
                Sp_max      = np.max(par.grids_dense[educ][j_pos_1])
                S_pol       = np.squeeze(S[j_pos][educ,:,:])
                C_pol       = np.squeeze(C[j_pos][educ,:,:])
                
                if np.max(ravel(S_pol[:])) > Sp_max:
                    print(f'Attention: age = {par.age[j_pos]}, educ = {educ}, extrapolation in ergodic distribution, min extrap = {log(np.max(S_pol[:])/Sp_max)}')
                    S_pol           = np.where(S_pol>Sp_max, Sp_max, S_pol)
                
                    Sp_min           = np.min(par.grids_dense[educ][j_pos+1])
                if np.min(np.ravel(S_pol[:])) > Sp_min:
                    print(f'Attention: age = {par.age[j_pos]}, educ = {educ}, extrapolation in ergodic distribution, min extrap = {log(np.min(S_pol[:])/Sp_min)}')
                    S_pol           = np.where(S_pol<Sp_min, Sp_min, S_pol)
                             
                pol = np.zeros([len(S_dense_aux),len(FE_pos)])
                
                for j in range(0, len(S_dense_aux)):
                    for k in range (0, len(FE_pos)):
                        pol[j,k] = interp_2d(S_orig,FE_pos,S_pol,S_dense_aux[j],FE_pos[k])
                S_dense[j_pos][educ,:,:]  = pol
                
                pol = np.zeros([len(S_dense_aux),len(FE_pos)])
            
                for j in range(0, len(S_dense_aux)):
                    for k in range (0, len(FE_pos)):
                        pol[j,k] = interp_2d(S_orig,FE_pos,C_pol,S_dense_aux[j],FE_pos[k])
                C_dense[j_pos][educ,:,:]  = pol
            
            #del S
        
        ## Save interpolated policy functions
        
        # 9. Rename grid vectors
        par2                        = copy.deepcopy(par)
        del par2.grids
        del par2.grids_dense
        par2.grids                  = copy.deepcopy(par.grids_dense)
        
        # Policy functions
        Se                          = copy.deepcopy(Se_dense)
        S                           = copy.deepcopy(S_dense)
        Ce                          = copy.deepcopy(Ce_dense)
        C                           = copy.deepcopy(C_dense)
        Ck                          = copy.deepcopy(Ck_dense)
        Np                          = copy.deepcopy(Np_dense)
        PHIp                        = copy.deepcopy(PHIp_dense)
        # Grs                       = Grs_dense
        pol_dense.Se                = copy.deepcopy(Se)
        pol_dense.S                 = copy.deepcopy(S)
        pol_dense.Ce                = copy.deepcopy(Ce)
        pol_dense.C                 = copy.deepcopy(C)
        pol_dense.Ck                = copy.deepcopy(Ck)
        pol_dense.Np                = copy.deepcopy(Np)
        pol_dense.PHIp              = copy.deepcopy(PHIp)
        pol_dense.Tp                = copy.deepcopy(Tp_dense)
        #pol_dense.Grs               = Grs
        pol_dense.tau0              = copy.deepcopy(tau0)
    
    return pol_dense,par,par2