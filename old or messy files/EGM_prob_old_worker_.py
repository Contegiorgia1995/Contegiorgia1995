# -*- coding: utf-8 -*-
"""
Created on Thu Aug 18 15:57:23 2022

@author: Giorgia
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 11:58:53 2022

@author: Giorgia
"""


from types import SimpleNamespace
import time
import itertools as it
import numpy as np
from scipy import optimize




def EGM_old_worker(par, ucp, dispinc, Spgrid, S):
    r_sav = par.r_sav
    r_debt = par.r_debt
    gammac = par.gammac
    
##check if ucp is decreasing
#Lag = (ucp(1:end-1)<=1e+10)*(ucp(2:end)-ucp(1:end-1))
#Lagmax = max(Lag)
#if Lagmax > 0:
    #fprint('ucp is not decreasing\n')
    
#solve for an endogenous
    c = (ucp)**(-1/gammac)# vector-- ucp, c, Spgrid and dipinc must have same dimension
    ################################################################
    # a_endo = np.zeros([len(c),len(Spgrid)])
    # for i in range(0,len(c)):
    #     for j in range(0,len(Spgrid)):
    #         a_endo[i,j] = c[i]+Spgrid[j] -dispinc
    # a_endo = a_endo.T
    a_endo = (c + np.reshape(Spgrid,[80,1]) - dispinc )############# careful here, it needs to be a one dimensional vector. c is (80,1) and Spgrid (1,80) which is why we use the transposed
    a_endo = a_endo/(1+r_sav)*(a_endo>=0) + a_endo/(1+r_debt)*(a_endo<0)
    #########################################################    
    #check borrowing limits
    #here we get scalars
    a_bc = a_endo[0] 
    S_bc = []
    
    #a_endo = c + Spgrid.T- dispinc
    for i in range(0,len(S)):
        if np.any(S[i]<a_bc):
            S_bc.append(S[i])      
    S_bc = np.array(S_bc)

    c_bc = (1 + r_sav)*S_bc * (S_bc >=0) + (1 +r_debt)*(S_bc)*(S_bc < 0) + dispinc - Spgrid[0] ######what is S_bc was not (2x1)?
    
    #Interpolate in original grid S (Since the borrowing limit is defined on strict inequality, all a_endo should be included in interpolation)
   
    c_final = approx_2d(np.append(S_bc,a_endo), np.append(c_bc,c), S)#### problem with dimensions
    #c_final = approx_2d(np.concatenate(np.reshape(S_bc,(len(S_bc),1)), a_endo), np.concatenate(np.reshape(c_bc,(len(c_bc),1)),c), S)#### problem with dimensions
    
    #c_final = approx_2d(np.concatenate( [S_bc,a_endo]), np.append(np.array([np.append(c_bc,c)]),np.array([np.append(c_bc,c)])), S)## works as in matlab for X, Y vector with same lenght and S vector even with different lenght (but while in Matlab it returns a value fro c_final, in python it gives an error even if thye deliver the same output y)
    
    sp_final = (1 +r_sav)*S*(S>=0) + (1 + r_debt)*S.T*(S.T<0) + dispinc - c_final #WORKS as in MATLAB!!! but same result as if I did not transpose S
    
    #sp_final = (1 +r_sav)*np.reshape(S,(len(S),1))*(np.reshape(S,(len(S),1))>=0) + (1 + r_debt)*np.reshape(S,(len(S),1))*(np.reshape(Spgrid,(len(Spgrid),1))<0) + dispinc #- c_final

    sp_final = sp_final * ((sp_final >= Spgrid[0]) + (sp_final<Spgrid[0]-1e-6)) + Spgrid[0] * (sp_final < Spgrid[0])*(sp_final >= Spgrid[0])*(sp_final>=Spgrid[0]-1e-6)### same result as matlab
    
    #if (sp_final < Spgrid[0]) > 0:
        #print('EGM error: savings are below debt limit - worng extrapolation, min(s''): %3.2f \n', min(sp_final))### must be fixed
    
    boundgrid = sum(sp_final >1.025*Spgrid[len(Spgrid)-1]) #should be ok, same result as Matlab

    
    #check for bad extrapolation
    c_max = (1 + r_sav)*S.T*(S.T >= 0) + (1 + r_debt)*S.T*(S.T<0) + dispinc-Spgrid[0] #working as in Matlab
    c_final = (c_final <= c_max)*c_final + (c_final > c_max) * c_max
    
    # check for bad extrapolation    
    sp_final = (1 +r_sav)*S.T*(S.T>=0) + (1 + r_debt)*S.T*(S.T<0) + dispinc - c_final
    
    return c_final, sp_final, boundgrid
        

                                                                                                          
                                                                                                 

    
    