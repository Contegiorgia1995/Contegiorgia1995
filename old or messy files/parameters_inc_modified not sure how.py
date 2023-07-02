# -*- coding: utf-8 -*-
"""
Created on Sun Aug  7 13:52:12 2022

@author: Giorgia
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jul 26 12:07:27 2022

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
import pyimsl
#import par_income_process_inc
#from par_income_process import par_income_process
#stats.lognorm(0.5, scale=np.exp(2)).ppf(0.005)


model = SimpleNamespace() #model
model.par = SimpleNamespace() #param
model.sol = SimpleNamespace() 
model.options = SimpleNamespace()
model.inc = SimpleNamespace()
model.par.inc = SimpleNamespace()
options = model.options
par = model.par
par_est = model.par
par = model.par
inc = model.inc
par.inc = model.par.inc
## Options

options.exerciseSS_On   = 'N'         # Default = N; if = Y: solve SS exercise (different wage and savings options)
options.Fertility       = 'Endo'       # 'Endo': endogenous fertility; 'Exo': constant fertility
options.Fertility_Exo_N = 2.15       # 'Exo' case: number of children
options.ParentTrans     = 'Endo'      # Endo: Benchmark case with parents transfers; 'Exo': exogenously given
options.ParentTransVal  = 0.05         # Value for exogenous trasnfers
options.ExPolRetirement = 'N'         # N: benchmark case; Y: Exercise Retirement Policy
options.GridH0orig      = 'N'        # Force initial grid of H0 -- To use when doing sensitivity analysis of parameters
options.patient         = 'N'         # 'N': Baseline Discount Factor; 'Y': Higher beta (lower discounting)
options.noborrowingwedge= 'N'         # 'N': Baseline Wedge between borrowing and lending; 'Y': Smaller wedge
options.fert_transfer   = 'N'
options.init_transfer   = 'N'
options.psy_avg         = 'N'
options.AdultRisk       = 'Y' #from start
options.timer_on    	= 'N';          # display timer

#inc.age_prof = cell(length(EDUC),1); #from file par_income_process
#inc.fe_pos   = np.arange(1,N_fe,N_fe)





x = np.matrix([[0, 0.247760398343029],
              [1, 0.756588865150413],
              [2, 0.253974483993104],
              [3, 25.7139793046902],
              [4, 1.73231100706877],
              [5, 3.98254310465696],
              [6, 0.557506847031282],
              [7, 0.236348281617006],
              [8, 0.163010113547491],
              [9, 0.109185709139714],
              [10, 0.668236388083336],
              [11, 0.664607727032972],
              [12, 0.185443539313036],
              [13, 0.470569749176502]])

x_val = np.array([ 0.247760398343029,
             0.756588865150413,
             0.253974483993104,
             25.7139793046902,
             1.73231100706877,
             3.98254310465696,
             0.557506847031282,
             0.236348281617006,
             0.163010113547491,
             0.109185709139714,
             0.668236388083336,
             0.664607727032972,
             0.185443539313036,
             0.470569749176502])


par_est.gamman      = x[0,1]
par_est.ln_level    = x[1,1]
par_est.mu          = 0.5
par_est.sigmah0     = x[2,1]
par_est.psy_max     = x[3,1]
par_est.psy_cor     = x[4,1]
par_est.psy_hs      = x[5,1]
par_est.R_hs_curv   = x[6,1]
par_est.R_col_curv  = x[7,1]
par_est.R_hs_level  = x[8,1]
par_est.R_col_level = x[9,1]
par_est.Nanny_oppC  = 1
par_est.Child_Ccurv = 0.6445
# par_est.Child_Ccurv = 1.000
par_est.Child_C_Tot = x[10,1] #Multiply by (1-tax) in child cost function to obtain parameter from paper
par_est.child_cost_inc_curv = x[11,1]
par_est.beta_inc    = x[12,1]
#par.meanh0 = 2
#par.sigmah0 = 0.5



def parameters(par_est,options):
    par = model.par
    par_est = model.par
## 1. Mean household income (unconditional of educ groups) in 2000 dollars
# * Source PSID, 2000 % 
    par.mean_inc        = 70179.45
    mean_inc_age42      = 75630                        # mean income at age 40-43 (to compute OAS transfers)

    options.exerciseSS_On  = 'N'         # Default = N; if = Y: solve SS exercise (different wage and savings options)
# Option to change w for SS exercises
    if options.exerciseSS_On == 'N':
        par.w       = 1
    else:
        par.w       = options.exerciseSS.w_ss

###############
    ## 2. Age structuretime_period     = 2;
    par.time_period     = 2
    par.age             = np.linspace(0,80,41)
    par.NJ              = len(par.age);
    par.Je1             = 16                                    # HS age
    par.Je1_pos         = list(par.age).index(par.Je1)+1        #np.where(par.age == par.Je1)
    par.Je2             = 18                                    #College age
    par.Je2_pos         = list(par.age).index(par.Je2)+1 
    par.Jc              = 28                                    #Fertility
    par.Jc_pos          = list(par.age).index(par.Jc)+1 
    par.Jt              = 40                                    #OAS
    par.Jt_pos          = list(par.age).index(par.Jt)+1
    par.Jr              = 66                                    #Retirement
    par.Jr_pos          =  list(par.age).index(par.Jr)+1
    par.Jd              = 80                                    #Die, source: WB - life expectancy = 79
    par.Jd_pos          =  list(par.age).index(par.Jd)+1

    

    ## 3. Income process
    par.educ     = np.array([1, 2, 3])
    par.dist_educ = np.array([0.094, 0.601, (1-.601-.094)])
    par.N_markov = 5
    par.N_fe     = 5
    par.sigmah0  = par_est.sigmah0 ##par_est.sigmah0     = x(3), x from function "calibrate"
    par.inc      = par_income_process(par,options) ##par_income_process fucntion defined in file par_income process
    if 1 == 0:
    # To calculate the persistence of the income process: XXX Need to Fix XXX
        solve_corr_model()### must define

    ##Education Returns
    par.Reducs      = np.zeros([2,2]) #Education stage, [constant, curvature]: h_(t+1) = h_t^curv * constant
    par.Reducs[0,:] = [par_est.R_hs_level, par_est.R_hs_curv] ## from calibratewith input x: par_est.R_hs_curv   = x(7); par_est.R_col_curv  = x(8); par_est.R_hs_level  = x(9); par_est.R_col_level = x(10); 
    par.Reducs[1,:] = [par_est.R_col_level, par_est.R_col_curv]
    
## 4. Initial H0 
########################################################################
# Distribution of parents: By Education 
###################################################################
    par.beta_inc        = par_est.beta_inc ## = x(13) from calibrate
    par.sigmah0         = par_est.sigmah0 ## x(3)
    par.meanh0          = par_est.meanh0 = 1 ##output of a loop in calibrate file and in par_grid_h0_old

    pp                  = np.linspace(0.025,0.975,par.N_fe)   
    
    def lognorminv(x,mu,sigma):
        lognorm = stats.lognorm(sigma, scale=math.exp(mu))
        lognorm.ppf(x)
        return lognorm.ppf(x)
    
    if options.GridH0orig == 'N':
        sigma_grid          = 1.5*par_est.sigmah0 #Larger grid of H0 so that it includes more points if beta_inc is large.
    else:
        sigma_grid          = 1.5*options.SigmaH0orig # Larger grid of H0 so that it includes more points if beta_inc is large.

    par.H0                  = lognorminv(pp,log(par.meanh0),sigma_grid) # eturns the inverse of the lognormal cdf with the distribution parameters mu (mean of logarithmic values) and sigma (standard deviation of logarithmic values), evaluated at the probability values in p.


#     def lognormcdf(x,mu,sigma):
#         for i in range(0,len(x_val)):
#             #for j in range(0,len(mu)):
#             lognorm = stats.lognorm(sigma, scale=math.exp(mu))
#             res = []
#             res = res.append(lognorm.cdf(x_val[i]))
#             i += 1
#             res_final = []
#             res_final = res_final.append(lognorm.cdf(x_val[i]))
#             return res_final
        
        
#     mat = np.zeros([len(mu),len(x)])  
#     def lognormcdf(x,mu,sigma):
#         for i in range(0,len(mu)):
#             mat[i,:] = mu[i]
#             for j in range (0, len(x)):
#                 mat[:,j] = x[j]
#                 return mat
            
# ####need to translate into function

    def logncdf(x,mu,sigma):  
        #mat = np.zeros([len(mu),len(x)])
        mat = np.empty([len(mu),len(x)])
        for i in range(0,len(mu)): 
            for j in range(0,len(x)):
                mat[i,j] = stats.lognorm(sigma, scale=math.exp(mu[i])).cdf(x[j])
            return mat

    # def lognormcdf(x_val,mu,sigma):
    #     #x_val = np.reshape(x_val,[1,len(x_val)])
    #     mat = np.zeros([len(mu),len(x_val)])
    #     for j in range(0,len(mu)):
    #             for i in range(0,len(x_val)):
    #                 lognorm = stats.lognorm(sigma, scale=math.exp(mu[j]))
    #                 mat[j,i] = lognorm.cdf(x_val[i])
    #                 return mat


######must be fixed
 #in the same code in  matlab they transpose the array of x so it is difficult to compare

###################################################################################
    #par.PG =  @(x) np.array([logncdf(np.reshape(np.append(np.array(par.H0[0:-1]),inf),[1,len(np.append(np.array(par.H0[0:-1]),inf))]),log(par.meanh0) + par_est.beta_inc * log(x_val),par_est.sigmah0) - lognormcdf(np.append(0, np.array(par.H0[0:-1])),log(par.meanh0) + par_est.beta_inc * log(x),par_est.sigmah0)]).T ##need to understand x()
    def PG(x):
        return np.array([logncdf(np.append(np.array(par.H0[0:-1]),inf),log(par.meanh0) + par_est.beta_inc * log(x_val),par_est.sigmah0) - logncdf(np.append(0, np.array(par.H0[0:-1])),log(par.meanh0) + par_est.beta_inc * log(x_val),par_est.sigmah0)]).T 
    par.PG = PG(x_val)##need to understand x()
########################################################################################                          
    # Transformed by education##
    ########### WARNING!!  par.inc.fe is a cell array...how to translate it?
    par.H0               = par.H0.T ### is it actually the right shape?
    H0                   = par.H0
    
    EDUC = par.educ
    par.inc.fe           = {}
    for i in range (0,len(EDUC)):
        par.inc.fe [i] = []
    par.inc.fe [0]      = log(H0)                                                #HS Drop
    par.inc.fe [1]      = log(   H0   + par.Reducs[0,0] * H0**par.Reducs[0,1])   #HS Grad
    par.inc.fe [2]      = log(   H0   + par.Reducs[1,0] * H0**par.Reducs[1,1])   #Col Grad

    ## Model-Data exchange rate: use mean income at age 40-43 for HS graduates, from average initial h0
    mean_inc_data = par.time_period*mean_inc_age42
    mean_inc_model = par.time_period # Solve to get h0 at age 42 = 2, 1 per year
    par.p_model_data = mean_inc_model / mean_inc_data
    
    ## 6. Preferences:
    # 6.1 Estimated from outside the model
    beta_annual = 0.975
    #options.patient = 'Y'
    if options.patient == 'Y':
        beta_annual = 0.99
    
    par.beta = beta_annual**par.time_period
    par.gammac = 0.5
    
    # 6.2 Estimated internally, altruism: b(n) = par.ln_level * n^(par.gamman)
    par.gamman = par_est.gamman # grid_gamman = [0.01, 0.8]
    par.lambdan = par_est.ln_level # grid_in = [0.01, 0.8]
    
    # 6.3 Education costs
    # 6.3.1 College: externally
    # Source: Delta Cost Project, average cost 2000 = 6588 per year.
    # In the model: 4 years of college and household ( two agents)
    par.pe2 = par.time_period * 2* 6588 * par.p_model_data
    # 6.3.2 High School: internally
    par.pe1 = par.pe2*0.09
    
    
    #options.exerciseSS_On = 'Y'
    options.exerciseSS_change_pe = 'Y'
    if options.exerciseSS_On == 'Y':
        if options.exerciseSS_change_pe == 'Y':
            par.pe2 = par.pe2*options.exerciseSS_p_mult
            par.pe1 = par.pe1*options.exerciseSS_p_mult
            
    pe2_total = (4/par.time_period)*par.pe2
    pe1_total = (4/par.time_period)*par.pe1
    
    # 6.4 Psychic cost: relative cost of education
    # par_est.psy         = 5;
    # par_est.psy         = 20;
    gridpsy = 100 # multiple of 5
    # par.psy_val_hs      = (linspace(0,par_est.psy_max,gridpsy).*par.pe1*par_est.psy_hs)';
    
    if options.exerciseSS_On == 'Y':
        if options.exercizeSS_change_psy == 'Y':
            par_est.psy_hs = par_est.psy_hs * options.exerciseSS_w_ss**options.psy_adj
            par_est.psy_max = par_est.psy_max * options.exerciseSS_w_ss**options.psy_adj
   
    par.psy_val_hs      = (np.linspace(0,par_est.psy_hs,gridpsy)).T
    par.psy_val_col     = (np.linspace(0,par_est.psy_max,gridpsy)).T
    # par.psy_val         = (linspace(0,par_est.sigmapsy,gridpsy).*par.pe2)';
    # ini                 = 100/(gridpsy) *.5;
    # psy_cdf             = linspace(ini,100-ini,gridpsy)';
    # psy_cdf             = psy_cdf/100;
    # par.psy_val         = norminv(psy_cdf,0,par_est.sigmapsy).*par.pe2;

    par.psy_prob        = np.zeros([3,gridpsy])
    Ngrid               = gridpsy+1
    # HS drop: more prob in higher values
    curv_psy            = par_est.psy_cor
    grid_cdf            = np.linspace(0**(1/curv_psy),1**(1/curv_psy),Ngrid)**curv_psy
    par.psy_prob[0,:]   = grid_cdf[1:]-grid_cdf[0:-1]
    
     # HS grad
    curv_psy            = 1
    grid_cdf            = np.linspace(0**(1/curv_psy),1**(1/curv_psy),Ngrid)**curv_psy
    par.psy_prob[1,:]   = grid_cdf[1:]-grid_cdf[0:-1]

    # CO grad: more prob in lower values
    curv_psy            = par_est.psy_cor
    grid_cdf            = np.linspace(0**(1/curv_psy),1**(1/curv_psy),Ngrid)**curv_psy
    grid_cdf_inv        = grid_cdf[1:]-grid_cdf[0:-1]
    par.psy_prob[2,:]   = grid_cdf_inv[::-1]

    prob_educ_data      = np.array([9.4, 60.1, 30.5])/100
    prob_psy            = np.dot(prob_educ_data,par.psy_prob)
    avg_psy             = np.array([np.dot(prob_psy,par.psy_val_hs), np.dot(prob_psy,par.psy_val_col)])
    
    if options.psy_avg == 'Y':
        gridpsy         = 2
        par.psy_val_hs  = np.array([avg_psy[0]*0.995,avg_psy[0]])
        par.psy_val_col = np.array([avg_psy[1]*0.995,avg_psy[1]])
            
        par.psy_prob    = np.ones([3,gridpsy])/2
            



       
    #6.5 Child cost parameters
    par.Nanny_oppC = par_est.Nanny_oppC
    par.Child_Ccurv = par_est.Child_Ccurv
    par.Child_C_Tot = par_est.Child_C_Tot
    par.child_cost_inc_curv = par_est.child_cost_inc_curv
    
    ##7. Partial Equilibrium prices and taxes
    ## 7. Partial Equilibrium prices and taxes
    ## 7.1 Wages: internally, to match mean income and education shares
    ## 7.2 Interest rate
    ## prime        = 1/beta_annual - 1;
    par.time_period     = 2 ## assumed
    prime        = 3/100; # Smets and Wouters AER 2007
    par.r        = (1+prime)**par.time_period-1
    par.r_sav     = (1+prime)**par.time_period-1            # Roys, Seshadri - 2014 XXX CONSIDER ADDING TAXES
    #par_est.int_iota = 0.1; # XXX Add to estimation? XXX
    par.int_iota = 0.1 ###### my modeification from par_est which is part of a more complicated function
    
    #options.noborrowingwedge = 'Y'
    if options.noborrowingwedge == 'Y':
        par_est.int_iota = 0.01
        
    par.r_debt   =  ((1+prime+par.int_iota)**par.time_period - 1)
    par.r_col    = ((1+prime+0.00925)**par.time_period - 1); #See /DK/data/Data_US/Education/Loans/College Loans.xls
    par.col_fact = par.r_col/(1-(1+par.r_col)**(-10))*(1-(1+par.r_debt)**(-10))/par.r_debt # Assume debt is repaid in 20 years, i.e. 10 periods.
    
    
    #Borrowing constraint: by education
    borr = par.w * (np.array([10000, 24000, 34000])*par.p_model_data) # self-reported limits on unsecured credit by family type from the SCF. Based on Abbot et al (2016).
    borr_sch = par.w * np.array([0, 0, 2*par.col_fact*23000])* par.p_model_data # Less borrowing if not finished college
    
    # 7.3 Government social security taxes
    if options.ExPolRetirement == 'N':
        par.labda = 0.124
    elif options.ExPolRetirement == 'Y':           # Social Security Ppayroll Tax, Krueger, Ludwig - 2015
        par.labda = options.PolRetTaxes
      
    par.pn = 0.54 * par.time_period* par.mean_inc*par.p_model_data #price of childcare, source: Folbre 2008, page 129. Wage of child care in 2000 = $7.43, mean wage = 13.74, 7.43/13.74=.54
    par.w_college = 0.56 * par.w # wage at college, eg inc = w_college * w * h; source: IPUMS
    
    if options.exerciseSS_On == 'Y':
        par.pn = par.pn * options.exerciseSS_p_mult ## actually commented in Matlab
        
    gdp_pc = 44308
    par.gov_debt = 0.2*gdp_pc*par.p_model_data
     
    # Source: AGMV 2013
    par.prod_alpha = 0.350
    par.prod_delta = 1-(1-0.0650)**par.time_period
    par.prod_s = np.array([0.160, 0.390, 0.450])
    par.prod_rho = 0.680
     
     ## 9. Grids for savings and discrete choices
     # 8.1 Number of children
     
    par.fam_size = 2
    if options.Fertility == 'Endo': #Endogenous fertility
        par.N = np.arange(0,4)
    elif options.Fertility == 'Exco': #Exogenous fertility
        par.N = options.Fertility_Exo_N/par.fam_size ####### from start_calibrate: options.Fertility_Exo_N = 2.15;         #'Exo' case: number of children
     
     # N = 0: no child
     # N = 1: fam_size children
     # N = 2: 2*fam_size children
     # etc.
     
     # 8.4 Savings
    par.Ls_0 = 20 #Length of savings grid for j = 12
    par.Ls_1 = 80 ##Length of savings grid for j = 16,20,24
     #% par.Ls              = [1 1 1 par.Ls_0 repmat(par.Ls_1,1,3) repmat(par.Ls_2,1,11) repmat(par.Ls_3,1,3)]; 
    par.Ls = np.append(np.append(np.ones(par.Je1_pos-1), par.Ls_0), np.repeat(par.Ls_1,par.Jd_pos - par.Je1_pos))

    curv = 3
        
    par.grids =  {}
    for i in range (0, len(EDUC)):
        a = [0]
        par.grids[i]= [a]*par.Jd_pos
    
        
    
    if options.ParentTrans == 'Endo':
        PHI_max = (pe1_total + pe2_total)*3
        par.PHI = np.append(np.array([0,pe1_total/3, 2/3*pe1_total]), (np.linspace((pe1_total*1.1)**(1/curv),PHI_max**(1/curv),int(par.Ls[par.Je1_pos-1])-3))**curv)        
        
    elif options.ParentTrans == 'Exo':
        # Computation of factor: 
        # 1. solve with factor = 1. 
        # 2. Solve for factor such that (res_MOM(4,100)*factor-res_MOM(1,100))/res_MOM(1,100) = 1
        if options.parentTransExo == 'MeanTrans':
            options.ParentTransVal = 2*30566*par.p_model_data
        elif options.parentTransExo == 'MeanHSCost':
            options.ParentTransVal = pe1_total*1.5
        elif options.parentTransExco == 'MeanColCost':
            options.ParentTransVal = (pe1_total + pe2_total)*1.5
        elif options.parentTransExo == 'search_val':
            options.ParentTransVal = 2*options.ParentTransVal*par.p_model_data
        par.PHI = options.ParentTransVal * np.array([0.999, 1.001])
               
    for ie in range (0,len(par.educ)):
        borr_limits = np.append(np.append(np.append(np.zeros(par.Je2_pos),borr_sch[ie]),np.append(borr_sch[ie],borr[ie]*np.ones(par.Jr_pos-1-(par.Je2_pos + 2)))), np.zeros(par.Jd_pos - (par.Jr_pos-1)))                      
        # j = Je1
        par.grids[ie][par.Je1_pos-1] = par.PHI ##### par.phi is a sequence
        #j = Je2
        maxE = par.PHI[-1]*2;

        for j in range(par.Je1_pos, par.Je2_pos+1):### fits only for j = 8
            #par.grids[j][ie]    = np.linspace(int(-borr_limits[j]**(1/curv)),int(maxE**(1/curv)),int(par.Ls[j]))**curv
            par.grids[ie][j] = np.linspace(-borr_limits[j]**(1/curv),maxE**(1/curv),int(par.Ls[j]))**curv

    # j = Je2+1 : Jc -1 CHECK! IT'S IMPORTANT THAT GRIDS STARTS CHANGING BY EDUCATION ONLY AFTER JE2+1
    
        for j in range(par.Je2_pos+1,par.Jc_pos-1):
            maxWY           = par.PHI[-1]*par.N[-1]*0.25

            max_inc         = par.w*exp(inc.age_prof[ie][0][j-1] + np.max(par.inc.fe[ie][:]) + np.max(inc.z_val[ie][j-1,:]))*2 # par.inc.age_prof to be changed in par_income_process to get rid of [0] as it is defined wrong and we should not need [0]
            maxS            = max(par.grids[ie][j-1][-1]*1.3,max(max_inc,maxWY))
            par.grids[ie][j] = np.linspace(-borr_limits[j]**(1/curv),maxS**(1/curv),int(par.Ls[j]))**curv
    
    # j = Jc: Jr
        for j in range (par.Jc_pos-1,par.Jr_pos):
            max_inc         = par.w*exp(inc.age_prof[ie][0][j-1] + np.max(par.inc.fe[ie][:]) + max(inc.z_val[ie][j-1,:]))
            maxS            = max(par.grids[ie][j-1][-1]*1.3,max_inc*3)
            par.grids[ie][j] = np.linspace(-borr_limits[j]**(1/curv),maxS**(1/curv),int(par.Ls[j]))**curv

    # j = Jr+1 : Jd
        maxR                = par.grids[ie][j][-1]*2
        for j in range (par.Jr_pos,par.Jd_pos):
            par.grids[ie][j] = np.linspace(-borr_limits[j]**(1/curv),maxR**(1/curv),int(par.Ls[j]))**curv


## Fertility Policies

    if options.fert_transfer == 'N':
        par.fert_trans = np.zeros(len(par.N))
    elif options.fert_transfer == 'only2':
        par.fert_trans = np.zeros(len(par.N))
        ind            = logical(par.N == 1);
        par.fert_trans[ind] = 2*options.fert_trans_size*par.p_model_data


    if options.init_transfer == 'N':
        par.init_trans = np.zeros([len(par.educ),par.Jd_pos])
        par.init_trans = np.reshape(par.init_trans,[len(par.educ),par.Jd_pos])
    elif options.init_transfer == 'Y':
        par.init_trans = np.zeros([len(par.educ),par.Jd_pos])
        par.init_trans = np.reshape(par.init_trans,[len(par.educ),par.Jd_pos])
###########################################################################################     missing options.init_trans_size
        par.init_trans[:,par.Je1_pos]   = 2*options.init_trans_size*par.p_model_data
    elif options.init_transfer == 'OnlyCollege':
        par.init_trans = np.zeros([len(par.educ),par.Jd_pos])
        par.init_trans = np.reshape(par.init_trans,[len(par.educ),par.Jd_pos])
        par.init_trans[2,(par.Je2_pos-1):(par.Je2_pos+1)]    = 2*options.init_trans_size*par.p_model_data
############################################################################################
#still need options.init_trans_size and to fix lognormcdf as results are wrong
    return par

parameters(par_est, options)

