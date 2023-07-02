# -*- coding: utf-8 -*-
"""
Created on Sat Sep 17 16:53:16 2022

@author: Giorgia
"""
import types, copy
def prob_educ_co_grad(par,options,Vp,Cp):
    # Solve household problem when it is College Graduate, for age 12 and 16
    # College graduate:
    educ      = 2
    PSY       = par.psy_val_col
    r_sav     = par.r_sav
    beta      = par.beta
    gammac    = par.gammac
    Lambda    = par.Lambda
    
    V       = {}  
        
    for i in range (0, 1):
        a = [0]
        V[i]= [a]*3
        V[0][0] = [a]*len(PSY)
    
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
        
        

    FE_pos     = par.inc.fe_pos
    AGE_PROF   = par.inc.age_prof
    INNO_pos   = par.inc.inno_pos
    INNO_posT = np.reshape(INNO_pos,[1,5])

    FE         = par.inc.fe[educ]
    
    Vpaux     = np.squeeze(Vp[:,educ,:,:])
    Cpaux     = np.squeeze(Cp[:,educ,:,:])
    
    
    # Expectation about innovation
    j_pos      = par.Je2_pos+1
    #exp_prob   = np.reshape(np.squeeze(par.inc.z_prob[educ][:,j_pos,0]),[1,len(np.squeeze(par.inc.z_prob[educ][:,j_pos,0]))]) #############here we determine id xp_prob is empy or not which will determine Vaux in GEGM_college
    exp_prob   = np.squeeze(par.inc.z_prob[educ][:,j_pos,0])
    ## Age 14-20
    for j_pos in range (par.Je2_pos,par.Je1_pos-2,-1):
        INNOp_prob = par.inc.z_prob[educ][:,j_pos+1,:].T
        S          = par.grids[educ][j_pos] # This grids do not change by education (at this age)
        Spgrid     = par.grids[educ][j_pos+1]
        V2         = np.zeros([len(S),len(FE_pos)])
        C2         = np.zeros([len(S),len(FE_pos)])
        Sp2        = np.zeros([len(S),len(FE_pos)])
        boundgrid  = np.zeros([len(FE_pos)])
        
    
        if j_pos   == par.Je2_pos:
            r_debt_today        = 0
            r_debt_tomorrow     = copy.deepcopy(par.r_debt)
            col_fact_1          = (par.col_fact * (Spgrid[0]<0) + 1 * (Spgrid[0]>=0))
            par_temp            = copy.deepcopy(par)
            par_temp.r_debt     = 0
            FE                  = par.inc.fe[1] #Skills before getting educated
        elif j_pos   == par.Je2_pos-1:
            r_debt_today        = copy.deepcopy(par.r_debt)
            r_debt_tomorrow     = 0
            col_fact_1          = 1
            par_temp            = copy.deepcopy(par)
            par_temp.col_fact   = 1
            FE                  = par.inc.fe[1] #Skills before getting educated
        else:
            r_debt_today        = copy.deepcopy(par.r_debt)
            r_debt_tomorrow     = copy.deepcopy(par.r_debt)
            FE                  = par.inc.fe[0] #Skills before getting educated
        
        # r = np.zeros(len(Spgrid))
        # for i in range(0,len(Spgrid)):
        #     if Spgrid [i]>= 0:
        #         r[i] = r_sav
        #     else:
        #         r[i] = r_debt_tomorrow
                
        r         = (r_sav * (Spgrid>=0) + r_debt_tomorrow * (Spgrid<0)).T #( should be trnasposed but it is not) currently working but keep in  mind!!!
        rs        = (r_sav * (S>=0) + r_debt_today * (S<0))
        
        if j_pos > par.Je2_pos-1:
            inno3      = np.repeat(INNO_posT,len(S),axis=0)
            INNO2,S2  = np.meshgrid(INNO_pos,Spgrid) 

        
        for ife in range (0,len(FE_pos)):
            if j_pos > par.Je2_pos-1:##########################################################I used >par.Je2_pos-1 instead of > as in Matlab because all j_pos are shifted down to 1
                splVp = GlobalSpline2D(Spgrid, INNO_pos,np.squeeze(Vpaux[:,:,ife]) )
                Cpp           = np.squeeze(Cpaux[:,:,ife]).T  ###################################################################################### 
                if min(np.ravel(Cpp.any())<0 and abs(min(np.ravel(Cpp.any()))))>1e-6:
                    Cpp  = max(Cpp,0)
                        
                ucp = beta*np.reshape(((1+r)),[len(r),1])*np.dot(Cpp**(-gammac), np.reshape(exp_prob,[size(exp_prob),1]))
                
            else:
                splVp   = spi.interp1d(Spgrid,np.squeeze(Vpaux[:,ife]), fill_value="extrapolate")
                Cpp           = np.squeeze(Cpaux[:,ife]) ###################################################################################### 
                
                if min(np.ravel(Cpp.any())<0 and abs(min(np.ravel(Cpp.any()))))>1e-6:
                    Cpp  = max(Cpp,0)
                        
                ucp = beta*np.reshape(((1+r)),[len(r),])*Cpp**(-gammac)*exp_prob
                ucp = np.reshape(ucp,[size(ucp),1])                
            
            h   = exp(AGE_PROF[educ][0][j_pos]) * exp(FE[ife])
            
            # if ndim(Cpaux) == 3:
            #     Cpp           = np.squeeze(Cpaux[:,:,ife]).T
               
            # else:
            #     Cpp           = Cpaux[:,ife].T
           
            ############################################
            if j_pos >= par.Je2_pos-1:
                
                dispinc = par.w_college*h*(1-Lambda) - par.pe2 + par.init_trans[educ][j_pos]
            else:
                
                dispinc = -par.pe1 + par.init_trans[educ][j_pos]
            
            feasible = (dispinc + (1+rs)*S - Spgrid[0]*(1/col_fact_1)>0)
            
            not_feasible = (dispinc + (1+rs)*S - Spgrid[0]*(1/col_fact_1)<=0)
            
            
            if j_pos >= par.Je2_pos-1:###################also here I made adjustements which I need to this of as it is > instead of >= despite being par.Je2_pos-1 instead of par.Je2_pos
                C2[feasible,ife],Sp2[feasible,ife],boundgrid[ife] = GEGM_college(par_temp,ucp.T,dispinc,Spgrid,S[feasible],splVp,exp_prob)
                C2[not_feasible,ife]  = 0
                Sp2[not_feasible,ife] = Spgrid[0]
            else:
                C2[feasible,ife],Sp2[feasible,ife],boundgrid[ife] = GEGM(par,ucp.T,dispinc,Spgrid,S[feasible],splVp,exp_prob)
                C2[not_feasible,ife]  = 0
                Sp2[not_feasible,ife] = Spgrid[0]

            
            if j_pos > par.Je2_pos-1:
                Sp3 = np.repeat(np.reshape(Sp2[:,ife],[len(Sp2[:,ife]),1]),len(INNO_pos),axis = 1)
                vp = np.dot(splVp(Sp3[:,0],inno3[0,:]).T,exp_prob)
            else:
                vp          = splVp(Sp2[:,ife])

            
            V2[feasible,ife]=C2[feasible,ife]**(1-gammac)/(1-gammac) + beta*vp[feasible]
            V2[not_feasible,ife] = -(10**(5/gammac))**(1-gammac)/(1-gammac)

        
    
        if j_pos >= par.Je1_pos: # case j_pos = 10,9 I switch to >= and leave par.Je1_pos becasue we want to get '0'
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
        exp_prob = 1

    return V,C,Sp