# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 12:47:17 2022

@author: Giorgia
"""
import numpy as np 

def prob_educ(par,options,Vp,Cp):
# Solve education choice

    # Age 16 (age of HS) 
    # Ve        = {}
    # for i in range(0,3):
    #     Ve[i] = [[],[],[]]
 ##########or   
    Ve =  {}
    for i in range (0, 3):
        a = [0]
        Ve[i]= [a]*3
        
    Ce =  {}
    for i in range (0, 3):
        a = [0]
        Ce[i]= [a]*3
        
    Se =  {}
    for i in range (0, 3):
        a = [0]
        Se[i]= [a]*3
    
            
    j_pos     = par.Je1_pos-1
    S         = par.grids[0][j_pos]
    FE_pos     = par.inc.fe_pos
    PSY       = par.psy_val_hs
    
    V0 = np.zeros([len(PSY),len(S),len(FE_pos)])
           
    tau = np.zeros([len(PSY),len(S),len(FE_pos)]) 
    
    ## Case 1: HS Dropout
    educ = 0
    Vhsd,Veaux,Ceaux,Seaux = prob_educ_hs_drop(par,options,Vp,Cp) ##ok<*ASGLU>
    
    for jp in range (0,3):
        Ve[educ][jp] = Veaux[0][jp]
        Ce[educ][jp] = Ceaux[0][jp]
        Se[educ][jp] = Seaux[0][jp]

    ## Case 2: HS Graduate
    educ = 1
    Veaux,Ceaux,Seaux = prob_educ_hs_grad(par,options,Vp,Cp)
    for jp in range (0,3):
        Ve[educ][jp] = Veaux[0][jp]
        Ce[educ][jp] = Ceaux[0][jp]
        Se[educ][jp] = Seaux[0][jp]
    
    ## Case 3: College Graduate
    educ = 2
    
    Veaux,Ceaux,Seaux = prob_educ_co_grad(par,options,Vp,Cp)
    for jp in range (0,3):
        Ve[educ][jp] = Veaux[0][jp]
        Ce[educ][jp] = Ceaux[0][jp]
        Se[educ][jp] = Seaux[0][jp]
    
    ## Optimal education:
    for ipsy in range (0,len(PSY)):
        for i_s in range (0,len(S)):
            for ife in range (0,len(FE_pos)):
                oo = Ce[1][0][99][i_s,ife]
                v_hsd = Vhsd[i_s][ife]
                v_hsg = Ve[1][0][ipsy][i_s,ife]
                v_cg  = Ve[2][0][ipsy][i_s,ife]
                
                v_aux = [v_hsd, v_hsg, v_cg]
                V0[ipsy,i_s,ife] = np.max(v_aux)
                tau[ipsy,i_s,ife] = np.argmax(v_aux)

    
    
    return Ve,Ce,Se,V0,tau