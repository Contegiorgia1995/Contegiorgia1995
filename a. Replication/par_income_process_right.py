# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 17:28:23 2022

@author: Giorgia
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 18:23:35 2022

@author: Giorgia
"""
from types import SimpleNamespace
import time
import itertools as it
import numpy as np
from scipy import optimize
from scipy import stats
from __future__ import print_function
from numpy import *
from statistics import NormalDist



options.AdultRisk       = 'Y' #from start


def par_income_process(par,options):
## Income Process
# Based on Abbott, Gallipoli, Meghir, Violante (2013). Updated.
# agent i, education e = 1,2,3, age a:
# Earnings_i,e,a = h_e,a gamma_e,i eta_e,a,i
# h_e,a:     age profile of efficency units of labor by education group
# gamma_e,i: fixed effect
# eta_e,a,i: idiosyncratic shock, AR(1)
# clear; clc; close all; 



############################################################################### MAY BE UNNECESSARY
    # par.time_period     = 2
    # par.age             = np.linspace(0,80,41)
    # par.NJ              = len(par.age);
    # par.Je1             = 16                                    # HS age
    # par.Je1_pos         = list(par.age).index(par.Je1)+1        #np.where(par.age == par.Je1)
    # par.Je2             = 18                                    #College age
    # par.Je2_pos         = list(par.age).index(par.Je2)+1 
    # par.Jc              = 28                                    #Fertility
    # par.Jc_pos          = list(par.age).index(par.Jc)+1 
    # par.Jt              = 40                                    #OAS
    # par.Jt_pos          = list(par.age).index(par.Jt)+1
    # par.Jr              = 66                                    #Retirement
    # par.Jr_pos          =  list(par.age).index(par.Jr)+1
    # par.Jd              = 80                                    #Die, source: WB - life expectancy = 79
    # par.Jd_pos          =  list(par.age).index(par.Jd)+1
    
    par.educ = [1,2,3]
    EDUC         = par.educ
    par.time_period     = 2
    par.age             = np.linspace(0,80,41)
    age          = par.age[0:par.Jr_pos]
    par.N_markov = 5
    N_markov     = par.N_markov
    N_fe         = par.N_fe = 5
    inc.inno_pos = list(range(0,N_markov))
    inc.fe_pos   = range(0,N_fe)
    
    age_start       = [par.Je1_pos-1, par.Je2_pos-1, par.Je2_pos+1]
########################################################################################
    ## 1. age profile of efficency units of labor by education group: h_e,a
    # Source: PSID
    def age_prof(e,a):
        if e==0:
            return 1*(0.0507705*a-0.0005522*a**2)
        elif e==1:
            return 1*(0.0668012*a-0.0007312*a**2)
        elif e==2:
            return 1*(0.1221627*a-0.0013147*a**2)



    inc.age_prof = {}
    for i in range(0,len(EDUC)):
        inc.age_prof[i] = [[]]
        

    for ie in range(0,len(EDUC)):
        inc.age_prof[ie][0]      = np.zeros(par.Jr_pos)
        for ia in range (par.Je1_pos-1,par.Jr_pos):
            inc.age_prof[ie][0][ia] = age_prof(ie,age[ia])
        
        inc.age_prof[ie][0][range(par.Je1_pos-1,par.Jr_pos)]     = inc.age_prof[ie][0][range(par.Je1_pos-1,par.Jr_pos)] - inc.age_prof[ie][0][par.Je1_pos-1]


                 
       ## 2. Idiosyncratic shock: eta_e,a,i 
       # Source: NLSY79
    inc.inno_rho   = [0.875153, 0.961399, 0.966375] # autocorrelation of innovation for education groups
    inc.z0_var      = [0.22832, 0.10134, 0.0660616] # Variance of iid shocks to innovation for education groups
    inc.inno_var    = [0.0624742, 0.02341, 0.0291241]   # Variance of iid shocks to innovation for education groups
    if options.AdultRisk == 'N':
       inc.inno_var = 1e-6*inc.inno_var
       inc.z0_var   = 1e-6*inc.z0_var

    inc.z_var       = np.zeros([len(EDUC),len(age)])
    inc.z_var = reshape(inc.z_var,[len(EDUC),len(age)])
    
    for ie in range (0,len(EDUC)):
        rho      = inc.inno_rho[ie]
        z0_var   = inc.z0_var[ie]
        inno_var = inc.inno_var[ie]
        inc.z_var[ie,age_start[ie]] = z0_var
        for ia in range (age_start[ie]+1,par.Jr_pos):
            inc.z_var[ie,ia] = rho**2 * inc.z_var[ie, ia-1] + inno_var


    # aux_inc.inno_var     = zeros(length(inc.educ),length(age));
    inc.z_val        = {}
    for i in range(0,len(EDUC)):
        inc.z_val[i] = [[]]
    
    inc.z_prob    = {}
    for i in range(0,len(EDUC)):
         
         inc.z_prob[i] = [[]]
    # for i in range(0,1):
    #     a = [0]
    #     inc.z_prob[i] = [[a]*len(EDUC)]
    
    def normcdf(x,mu,sigma):
        result = np.zeros(len(x))
        for i in range(0,len(x)):
            result[i] = NormalDist(mu,sigma).cdf(x[i])
        return result
    
    
    
    # Z0 PROBABILITIES (AND VALUES) 
    for ie in range (0,len(EDUC)): 
        inc.z_val[ie]        = np.zeros([par.Jr_pos,N_markov])
        inc.z_prob[ie]     = np.repeat(np.zeros([par.Jr_pos,N_markov]),N_markov)
        inc.z_prob[ie] = np.reshape(inc.z_prob[ie],[N_markov,par.Jr_pos,N_markov])
                                 
        # Initial draw of z
        ia                      = age_start[ie]
        nu                      = sqrt(N_markov-1) * inc.z_var[ie,ia]**0.5
        etagrid                 = np.linspace(-nu, nu, N_markov).T
        inc.z_val[ie][ia,:]     = etagrid
   
        grid_med                 = 0.5*(etagrid[1:]+etagrid[0:-1])
        aux_cdf                  = normcdf(grid_med,0,inc.z_var[ie,ia]**0.5)
        aux_pdf                  = np.append(aux_cdf[0], [aux_cdf[1:]-aux_cdf[0:-1]])
        aux_pdf                  = np.append(aux_pdf, 1-sum(aux_pdf))
        inc.z_prob[ie][:,ia,:] = np.reshape(np.repeat(aux_pdf,N_markov),[5,5])
        
        for ia in range (age_start[ie]+1,par.Jr_pos):
            ############################################################
            P, etagrid  =  discretize_nonstationary_ar_rouwen(inc.inno_rho[ie],inc.z_var[ie,ia-1]**0.5,inc.z_var[ie,ia]**0.5, N_markov)           
            print(P)
            #inc_z_prob[ie][:,ia,:],inc_z_val[ie][ia,:] = discretize_nonstationary_ar_rouwen(inc_inno_rho[ie],inc_z_var[ie,ia-1]**0.5,inc_z_var[ie,ia]**0.5, N_markov)
            inc.z_prob[ie][:,ia,:] = P.T
            inc.z_val[ie][ia,:] = etagrid
            
    return inc
##need to create a SimpleNameSpace first

