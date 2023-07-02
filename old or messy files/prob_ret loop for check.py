# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 13:35:23 2022

@author: Giorgia
"""


#####################
V_   = {}
for i in range (0, par.Jd_pos):
    V_[i] = [[]]

C_   = {}
for i in range (0, par.Jd_pos):
    C_[i] = [[]]
Ck   = {}
for i in range (0, par.Jd_pos):
    Ck[i] = [[]]
    
S   = {}
for i in range (0, par.Jd_pos):
    S[i] = [[]]
    
Tp   = {}
for i in range (0, par.Jd_pos):
    Tp[i] = [[]]

## Dead: V_Jd = 0
V_[par.Jd_pos-1] = np.zeros([len(par.educ),int(par.Ls[-1]),par.N_fe])
C_[par.Jd_pos-1] = np.zeros([len(par.educ),int(par.Ls[-1]),par.N_fe])
S[par.Jd_pos-1] = np.zeros([len(par.educ),int(par.Ls[-1]),par.N_fe])

for j_pos in range (par.Jd_pos-2, par.Jd_pos-3,-1):#in range (par.Jd_pos-2,par.Jr_pos-2,-1):   
    Cp = C_[j_pos + 1]
    Vp = V_[j_pos + 1]
    
    if j_pos == par.Jd_pos-2: #consume all, next period is dead
        S = par.grids[0][j_pos-1]## matlab uses {} because par.grids is a cell array which I do not manage to convert in python
        Sp = np.zeros([len(EDUC),len(S),len(FE_pos)])
        C = np.zeros([len(EDUC),len(S),len(FE_pos)])
        V = np.zeros([len(EDUC),len(S),len(FE_pos)])
        
        for educ in range(0,len(EDUC)):
            S = par.grids[educ][j_pos-1]
            C[educ,:,:,] = ((1+r_sav)*np.repeat([S],len(FE_pos),axis =0) + np.reshape(np.repeat(ret_rep(par,FE_pos,educ,options),len(S)),[len(FE_pos),len(S)])).T
    else:# solve Euler equation
        S = par.grids[0][j_pos-1]## matlab uses {} because par.grids is a cell array which I do not manage to convert in python
        Sp = np.zeros([len(EDUC),len(S),len(FE_pos)])
        C = np.zeros([len(EDUC),len(S),len(FE_pos)])
        V = np.zeros([len(EDUC),len(S),len(FE_pos)])

        boundgrid = np.zeros([len(FE_pos),len(EDUC)])

        for educ in range (0,len(EDUC)):
            S       = par.grids[educ][j_pos-1]
            Spgrid  = par.grids[educ][j_pos]
        
            for ife in range (0,len(FE_pos)):
                ##############################################
                #splVp = RegularGridInterpolator((Spgrid,),Vp[educ][ife][:]) 
                #splVp = spi.interp1d(Spgrid,Vp[educ,:,ife],axis=0, fill_value="extrapolate")
                #splVp = scipy.interpolate.interp1d(Spgrid,Vp[educ][ife][:])
                ##############################################
                Cpp   = Cp[educ,:,ife] #actually Cp, not C

                ucp     = beta*(1+r_sav)* Cpp**(-gammac)

                dispinc = ret_rep(par,FE_pos,educ,options)
                dispinc = dispinc[ife]

                C[educ,:,ife],Sp[educ,:,ife],boundgrid[ife,educ] = EGM(par,ucp,dispinc,Spgrid,S)

                V[educ,:,ife] = C[educ,:,ife]**(1-gammac)/(1-gammac) + beta*splVp(Sp[educ,:,ife])
                
    errors = sum(boundgrid[:])/np.size(C[:])