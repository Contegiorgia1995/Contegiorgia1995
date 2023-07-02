# -*- coding: utf-8 -*-
"""
Created on Sat Sep 17 17:22:28 2022

@author: Giorgia
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 19:28:30 2022

@author: Giorgia
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 21:02:14 2022

@author: Giorgia
"""
import numpy as np

import code
#import matplotlib.pyplot as plt


def GEGM_college(par,ucp, dispinc, Spgrid, S, splVp, innop_prob):# inputs: ucp (1,80) (ucp.T compared to the ucp in prob_work_withchild), so c is (80,1)
# EGM
    #options = optimset('Display','off') #options = optimset(Name,Value) returns options with specified parameters set using one or more name-value pair arguments.
    #https://www.mathworks.com/help/matlab/ref/optimset.html
    r_sav  = par.r_sav
    r_debt = par.r_debt
    gammac = par.gammac
    beta   = par.beta
    INNO_pos   = range(0,size(innop_prob)) #Because only in last period need to take expectations
    
    ## Find no-concave region
    difmin = np.append(ucp[:,0:-1]- ucp[:,1::],0)
    indmin = (difmin<0)
    
    difmax = np.append(0, ucp[:,1::] - ucp[:,0:-1])
    indmax = (difmax>0);
    if sum(indmin)>0:
        #ucp = np.reshape(ucp,[len(ucp),1])
        vmin = np.min(ucp[:,indmin])
        vmax = np.max(ucp[:,indmax])
        #     vmin = min(ucp.*indmin+1e20.*(1-indmin));
        #     [vmax,~] = max(ucp.*indmax);
        
        imin = np.argmax(range(0,size(ucp))*(ucp>vmax))
        
        if vmin > np.min(ucp):
            imax = np.min(range(0,size(ucp))*(ucp<vmin)+1e20*(ucp>=vmin))
        elif vmin == np.min(ucp):
            imax = size(ucp)-1## hopefully -1 is ok, it was for GEGM
        
            #code.interact(local=locals())#https://docs.python.org/3/library/pdb.html

        
    else:
        imin = 0
        imax = 1

    
    # non concave region: [imin+1, imax-1]
    # figure
    # plot(Spgrid,ucp,Spgrid,vmin.*ones(size(Spgrid)),Spgrid,vmax.*ones(size(Spgrid)))
    
    
    ## Solve for a endogenous
    c = ( ucp.T)**(-1/gammac)
##    
    col_fact = (par.col_fact * (Spgrid<0) + 1 * (Spgrid>=0))
    col_fact = np.reshape(col_fact,[size(col_fact),1])
    lambdan = par.lambdan
    gamman  = par.gamman
    
    
    a_endo = (c +  np.reshape(Spgrid,[1,len(Spgrid)]).T*(1/col_fact) - dispinc )
    a_endo = a_endo/(1+r_sav)*(a_endo>=0) + a_endo/(1+r_debt)*(a_endo<0)
    
    # 
    # figure(2)
    ##############################################################
    # plot(a_endo,Spgrid)
    # fig, ax = plt.subplots()
    # ax.scatter(x=a_endo, y=Spgrid, marker='o', c='r', edgecolor='b')
    
    ## Find global solution:
    a_opt  = copy.deepcopy(a_endo)
    sp_opt = copy.deepcopy(Spgrid)
    c_opt  = copy.deepcopy(c)
    
    
    # Find global solution in non concave region
    if imin == 0:
        if imax != 1:
            for i in range(int(imin),int(imax)): #for i in range(int(imin+1),int(imax)):#################################IPERCAREFUL!!!!!!! i = 0 was excluded from the loop but chnaging loop starting from imin insteado of imin+1 it may mess tigs up!
                if i == imin: #if i == imin + 1:#################################IPERCAREFUL!!!!!!! i = 0 was excluded from the loop but it may mess tigs up!
                    if size(innop_prob)>1:
                        inno3,Sp3  = np.meshgrid(INNO_pos,Spgrid)
                       #Vpaux          = splVp(Sp3,inno3)*innop_prob.T
                        Vpaux          = np.dot(splVp(Sp3[:,0],inno3[0,:]).T,innop_prob)
                    else:
                        Vpaux          = splVp(Spgrid).T######################problem
                    
                    curv         = 2
                    Spgrid_dense_neg = -np.linspace(0,(-Spgrid[0])**(1/curv),300)**curv
                    Spgrid_dense = np.append(Spgrid_dense_neg[::-1], np.linspace(1e-6**(1/curv),np.max(Spgrid[:])**(1/curv),1000)**curv)
        ####            
                    col_fact_dense = (par.col_fact * (Spgrid_dense<0) + 1 * (Spgrid_dense>=0))
                    col_fact_dense = np.reshape(col_fact_dense,[size(col_fact_dense),1])
                    
                    
                    if size(innop_prob)>1:
                        inno3,Sp3  = np.meshgrid(INNO_pos,Spgrid_dense)
                        Vpaux_dense = np.dot(splVp(Sp3[:,0],inno3[0,:]).T,innop_prob)
                    else:
                        Vpaux_dense = splVp(Spgrid_dense).T
                a_candidate = a_endo[i]
               
                c = (1+r_sav)*a_candidate*(a_candidate>=0) + (1+r_debt)*a_candidate*(a_candidate<0)  +dispinc - np.reshape(Spgrid,[size(Spgrid),1])*(1./col_fact)
                
                objfunc = np.zeros(size(c))
                for j in range(0,size(c)):
                    if c[j]>=0:
                        objfunc[j] = c[j]**(1-gammac)/(1-gammac)+ beta*Vpaux[j]
                    else: 
                        objfunc[j] = (-(10**(5/gammac))**(1-gammac)/(1-gammac)) + beta*Vpaux[j]
                #objfunc = (c**(1-gammac)/(1-gammac))*(c>=0) + (-(10**(5/gammac))**(1-gammac)/(1-gammac)*(c<0)) + beta*Vpaux
               
                amax = np.argmax(objfunc)
               
                if amax != i: # Discard solution
                 
            #       fprintf('z:%3.2f,a''_candidate:%3.2f,a''_global:%3.2f, \n',a_candidate,Spgrid(i),Spgrid(amax))
            #       a_opt(i)  = nan;
            #       sp_opt(i) = nan;
            #       c_opt(i)  = nan;
                  
                    c = (1+r_sav)*a_candidate*(a_candidate>=0) + (1+r_debt)*a_candidate*(a_candidate<0)  +dispinc - np.reshape(Spgrid_dense,[size(Spgrid_dense),1])*(1./col_fact_dense)
                    
                    
                    objfunc = np.zeros(size(c))
                    for j in range(0,size(c)):
                        if c[j]>=0:
                            objfunc[j] = c[j]**(1-gammac)/(1-gammac)+ beta*Vpaux_dense[j]
                        else: 
                            objfunc[j] = (-(10**(5/gammac))**(1-gammac)/(1-gammac)) + beta*Vpaux_dense[j]
        
                    amax = np.argmax(objfunc)
                  
                    a_opt[i]  = a_candidate
                    sp_opt[i] = Spgrid_dense[amax]
                    c_opt[i]  = c[amax]
    else:
        for i in range(int(imin+1),int(imax)): # it may mess tigs up!
            if i == imin+1: 
                if size(innop_prob)>1:
                    inno3,Sp3  = np.meshgrid(INNO_pos,Spgrid)
                   #Vpaux          = splVp(Sp3,inno3)*innop_prob.T
                    Vpaux          = np.dot(splVp(Sp3[:,0],inno3[0,:]).T,innop_prob)
                else:
                    Vpaux          = splVp(Spgrid).T######################problem
                
                curv         = 2
                Spgrid_dense_neg = -np.linspace(0,(-Spgrid[0])**(1/curv),300)**curv
                Spgrid_dense = np.append(Spgrid_dense_neg[::-1], np.linspace(1e-6**(1/curv),np.max(Spgrid[:])**(1/curv),1000)**curv)
    ####            
                col_fact_dense = (par.col_fact * (Spgrid_dense<0) + 1 * (Spgrid_dense>=0))
                col_fact_dense = np.reshape(col_fact_dense,[size(col_fact_dense),1])
                
                
                if size(innop_prob)>1:
                    inno3,Sp3  = np.meshgrid(INNO_pos,Spgrid_dense)
                    Vpaux_dense = np.dot(splVp(Sp3[:,0],inno3[0,:]).T,innop_prob)
                else:
                    Vpaux_dense = splVp(Spgrid_dense).T
            a_candidate = a_endo[i]
           
            c = (1+r_sav)*a_candidate*(a_candidate>=0) + (1+r_debt)*a_candidate*(a_candidate<0)  +dispinc - np.reshape(Spgrid,[size(Spgrid),1])*(1./col_fact)
            
            objfunc = np.zeros(size(c))
            for j in range(0,size(c)):
                if c[j]>=0:
                    objfunc[j] = c[j]**(1-gammac)/(1-gammac)+ beta*Vpaux[j]
                else: 
                    objfunc[j] = (-(10**(5/gammac))**(1-gammac)/(1-gammac)) + beta*Vpaux[j]
            #objfunc = (c**(1-gammac)/(1-gammac))*(c>=0) + (-(10**(5/gammac))**(1-gammac)/(1-gammac)*(c<0)) + beta*Vpaux
           
            amax = np.argmax(objfunc)
           
            if amax != i: # Discard solution
             
        #       fprintf('z:%3.2f,a''_candidate:%3.2f,a''_global:%3.2f, \n',a_candidate,Spgrid(i),Spgrid(amax))
        #       a_opt(i)  = nan;
        #       sp_opt(i) = nan;
        #       c_opt(i)  = nan;
              
                c = (1+r_sav)*a_candidate*(a_candidate>=0) + (1+r_debt)*a_candidate*(a_candidate<0)  +dispinc - np.reshape(Spgrid_dense,[size(Spgrid_dense),1])*(1./col_fact_dense)
                
                
                objfunc = np.zeros(size(c))
                for j in range(0,size(c)):
                    if c[j]>=0:
                        objfunc[j] = c[j]**(1-gammac)/(1-gammac)+ beta*Vpaux_dense[j]
                    else: 
                        objfunc[j] = (-(10**(5/gammac))**(1-gammac)/(1-gammac)) + beta*Vpaux_dense[j]
    
                amax = np.argmax(objfunc)
              
                a_opt[i]  = a_candidate
                sp_opt[i] = Spgrid_dense[amax]
                c_opt[i]  = c[amax]
            
            
    
    # figure(3)
    # plot(a_endo,Spgrid,a_opt,sp_opt)
    
    ## Borrowing constraint
    # Case 1: Spgrid(1) is global solution
    col_fact_1 = (par.col_fact * (Spgrid[0]<0) + 1 * (Spgrid[0]>=0))
    if isnan(sp_opt[1]) == 0:
    #     fprintf('Borrowing constraint - global \n')
        a_bc    = a_opt[0] # Any current level of savings below this should be constrained 
        
        S_bc = []
        
        #a_endo = c + Spgrid.T- dispinc
        for i in range(0,len(S)):
            if np.any(S[i]<a_bc):
                S_bc.append(S[i])      
        S_bc = np.array(S_bc)# Borrowing constrained
         
        c_bc    = (1+r_sav) * S_bc * (S_bc >= 0) + (1+r_debt) * S_bc * (S_bc < 0) + dispinc - Spgrid[0]*(1/col_fact_1)
        
    else: # Approximate a_bc
        pos = np.argmin(sp_opt)
        a_opt1    = a_opt[pos]
        
        if a_opt1<S[0]:
    #         fprintf('Borrowing constraint - local, but never binding \n')
            a_bc    = a_opt[pos]
            S_bc = []
            
            #a_endo = c + Spgrid.T- dispinc
            for i in range(0,len(S)):
                if np.any(S[i]<a_bc):
                    S_bc.append(S[i])      
            S_bc = np.array(S_bc)# Borrowing constrained
            
            c_bc    = (1+r_sav) * S_bc * (S_bc >= 0) + (1+r_debt) * S_bc * (S_bc < 0) + dispinc - Spgrid[0]*(1/col_fact_1)
            
        else:
    #         fprintf('Borrowing constraint - local & potentially binding \n')
            S_bc = S[S<a_opt1]
            
            if size(innop_prob)>1:
               inno3,Sp3  = np.meshgrid(INNO_pos,Spgrid)
               #Vpaux          = splVp(Sp3,inno3)*innop_prob.T
               Vpaux          = np.dot(splVp(Sp3[:,0],inno3[0,:]).T,innop_prob)
            else:
               Vpaux          = splVp(Spgrid).T
          
            
            c = (np.repeat((1+r_sav) * S_bc * (S_bc >= 0) + (1+r_debt) * S_bc * (S_bc < 0),size(Spgrid.T)) +dispinc - np.repeat(Spgrid.T*(1/col_fact_1),size(S_bc)))
            c_k  = f_n * c
            
            objfunc = np.zeros(len(c))
            for j in range(0,len(c)):
                if c[j]>=0:
                    objfunc[j] = c[j]**(1-gammac)/(1-gammac) + beta*np.repeat(Vpaux,size(S_bc))[j]
                else: 
                    objfunc[j] = (-(10**(5/gammac))**(1-gammac)/(1-gammac)) + beta*np.repeat(Vpaux,size(S_bc))[j]
            
            
            amax = np.argmax(objfunc, axis =0)
            Sp_vfi   = Spgrid[amax]
            col_fact_vfi = (par.col_fact * (Sp_vfi<0) + 1 * (Sp_vfi>=0))
            c_bc    = (1+r_sav) * S_bc * (S_bc >= 0) + (1+r_debt) * S_bc * (S_bc < 0) + dispinc - Sp_vfi*(1./col_fact_vfi)
            

 ##########################################################################   
 ##################################################################  ##so fa exactly like matlab 
    # Discard local solutions and sort grid
    a_opt = a_opt[(isnan(a_opt)==0)]
    c_opt = c_opt[(isnan(c_opt)==0)]
    
    ord = np.argsort(a_opt)
    a_opt = np.sort(a_opt)
    c_opt = c_opt[ord]
    
    
    # Interpolate in originial grid S
    # (Since the borrowing limit is defined on strict inequality, all a_endo
    # should be included in interpolation)
    
    c_final  = approx_2d(np.append(S_bc, a_opt), np.append(c_bc, c_opt), S)
    
    sp_final = (1+r_sav)*S.T*(S.T>=0) + (1+r_debt)*S.T*(S.T<0) + dispinc - c_final 
    col_fact = (par.col_fact * (sp_final<0) + 1 * (sp_final>=0))
    sp_final = sp_final*col_fact # Back to original grid
    sp_final = sp_final * ((sp_final>=Spgrid[0]) +(sp_final<Spgrid[0]-1e-6)) + Spgrid[0] *(sp_final<Spgrid[0])*(sp_final>=Spgrid[0]-1e-6)
    
    ###if len(innop_prob)>1:
    if size(innop_prob)>1:
        inno3,Sp3       = np.meshgrid(INNO_pos,sp_final)
        Vpaux0          = np.dot(splVp(Sp3[:,0],inno3[0,:]).T,innop_prob)
    else:
        Vpaux0          = splVp(sp_final)
###################################################CAREFUL THAT objfunc is defined differently for c-final >00 or c_final <0!!!!!!!!!!!!!!!!!
    objfunc0 = np.zeros(len(c_final))

    for j in range(0,len(c_final)):
        if c_final[j]>=0:
            objfunc0[j] = c_final[j]**(1-gammac)/(1-gammac)+ beta*Vpaux0[j]
            
        else: 
            objfunc0[j] = (-(10**(5/gammac))**(1-gammac)/(1-gammac)) + beta*Vpaux0[j]
    
   

    
    boundgrid = sum(sp_final>1.025*Spgrid[-1])
    ##############################################c_final still as good as in matlab
    # Check if there are slow decreases in C, i.e., not jumps
    difc = np.append(0, c_final[1::] - c_final[0:-1])
    inddecc = (difc<0)
    inddecc = (inddecc*((np.append(1, inddecc[0:-1])==1) + (np.append(inddecc[1::], 1) == 1))>0)
    #########################################################################so far ok then messed up
    inddecc[0:-2] = np.max([inddecc[0:-2],inddecc[2::]],axis =0)
    inddecc[2::] = np.max([inddecc[2::],  inddecc[0:-2]],axis = 0)
    inddecc[objfunc0<0] = 1
    inddecc[np.append(0, objfunc0[0:-1])<0] = 1 #inddedd 78 but objfunc is 80
        
    grid_prob = np.arange(0,len(c_final))
    do_interp = 1
    for i in (grid_prob[inddecc]):
        if do_interp == 1:
            curv         = 2
            Spgrid_dense_neg = -np.linspace(0,(-Spgrid[0])**(1/curv),300)**curv
            Spgrid_dense = np.append(Spgrid_dense_neg[::-1], np.linspace(1e-6**(1/curv),max(Spgrid[:])**(1/curv),1000)**curv)
            col_fact_dense = (par.col_fact * (Spgrid_dense<0) + 1 * (Spgrid_dense>=0))
            
            
            if size(innop_prob)>1:
                inno3, Sp3  = np.meshgrid(INNO_pos,Spgrid_dense)
                Vpaux_dense = np.dot(splVp(Spgrid_dense,INNO_pos).T,innop_prob)
            else:
                Vpaux_dense          = splVp(Spgrid_dense)
            
            do_interp    = 0
        
        c = (1+r_sav)*S[i]*(S[i]>=0) + (1+r_debt)*S[i]*(S[i]<0)  +dispinc - Spgrid_dense.T*(1./col_fact_dense)
        
        objfunc = np.zeros(len(c))

        for j in range(0,len(c)):
            if c[j]>=0:
                objfunc[j] = c[j]**(1-gammac)/(1-gammac)+ beta*Vpaux_dense[j]
                
            else: 
                objfunc[j] = (-(10**(5/gammac))**(1-gammac)/(1-gammac)) + beta*Vpaux_dense[j]

       
        #objfunc0 = c**(1-gammac)/(1-gammac) * (c>=0) + -(10**(5/gammac))**(1-gammac)/(1-gammac)*(c<0) + beta*Vpaux_dense
    
        amax = np.argmax(objfunc)
        sp_final[i] = Spgrid_dense[amax]
        c_final[i]  = c[amax] #######HERE cfinal[5:6] is wrong!!!!!!!!!!!!! becasue objfunc is  different and so in argmax (objfunc)
        
    # check for bad extrapolation
    c_max   = (1+r_sav)*S.T*(S.T>=0) + (1+r_debt)*S.T*(S.T<0)+dispinc-Spgrid[0]*(1/col_fact_1)
   
    c_final = ( c_final<= c_max) * c_final + (c_final > c_max) * c_max
    
    sp_final = (1+r_sav)*S.T*(S.T>=0) + (1+r_debt)*S.T*(S.T<0) + dispinc - c_final 
    col_fact = (par.col_fact * (sp_final<0) + 1 * (sp_final>=0))
    sp_final = sp_final*col_fact # Back to original grid
    
    if ((np.argmin(sp_final)-Spgrid[0])/np.abs(Spgrid[0]) > 0.01) and (np.argmin(sp_final)-Spgrid[0]<0):
        print('EGM error: savings are below debt limit - wrong extrapolation, min(s''): {min(sp_final)} \n',100*(np.argmin(sp_final)-Spgrid[0])/np.abs(Spgrid[0]))


    return c_final, sp_final, boundgrid


#np.dot(np.reshape(splVp(Sp3[:,0]),[1,len(Sp3[:,0])]).T,np.reshape(innop_prob,[1,len(innop_prob)]))