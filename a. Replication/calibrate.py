# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 13:17:46 2022

@author: Giorgia
"""
import time
import itertools as it
import numpy as np
from scipy import optimize
from scipy import stats
from __future__ import print_function
from numpy import *
from statistics import NormalDist
import scipy.interpolate as spi
from scipy.io import loadmat
import os
from types import SimpleNamespace
import numpy as np
import random
import types, copy


model = SimpleNamespace() #model
model.par = SimpleNamespace() #param
model.sol = SimpleNamespace() 
model.options = SimpleNamespace()
model.inc = SimpleNamespace()
model.par.inc = SimpleNamespace()
model.par_temp = SimpleNamespace()
model.guess = SimpleNamespace()
model.options.guess = SimpleNamespace()
model.pol = SimpleNamespace()
model.flags = SimpleNamespace()

options = model.options
par = model.par
par_temp = model.par_temp
par_est = model.par
par = model.par
inc = model.inc
par.inc = model.par.inc
guess = model.guess
options.guess = model.options.guess 
pol = model.pol
flags = model.flags



def tic():
    #Homemade version of matlab tic and toc functions
    import time
    global startTime_for_tictoc
    startTime_for_tictoc = time.time()

def toc():
    import time
    if 'startTime_for_tictoc' in globals():
        print("Elapsed time is " + str(time.time() - startTime_for_tictoc) + " seconds.")
    else:
        print("Toc: start time not set")
tic()
toc()

# Load parameters
# Estimation:
x = np.matrix([[0, 0.24776039834302887588],
              [1, 0.75658886515041290366],  
              [2,  0.25397448399310357248],
              [3, 25.71397930469017012456],
              [4, 1.7323110070687675055],
              [5, 3.9825431046569619298],
              [6, 0.55750684703128172703],
              [7, 0.23634828161700616178],
              [8, 0.16301011354749064819],
              [9, 0.10918570913971439862],
              [10, 0.66823638808333596373],
              [11, 0.6646077270329715514],
              [12, 0.18544353931303594885],
              [13, 0.47056974917650218337]])
                   

x_val = np.array([ 0.24776039834302887588,
             0.75658886515041290366,
             0.25397448399310357248,
             25.71397930469017012456,
             1.7323110070687675055,
             3.9825431046569619298,
             0.55750684703128172703,
             0.23634828161700616178,
             0.16301011354749064819,
             0.10918570913971439862,
             0.66823638808333596373,
             0.6646077270329715514,
             0.18544353931303594885,
             0.47056974917650218337])

def calibrate(x, options):
    
##Solve model and compute moments  for estimated parameters x
    random.randint(0,666)                              #Set seed for shocks
    ticcalib = tic#########################fix
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
    #if options.solve_h0 = 'N'
    par_est.meanh0      = x[13,1]
    
    
    options.solve_ho = 'N' ###will need to change this
    options.soboln = 1
    options.n_fold = 1
    
    if options.solve_h0 == 'Y':
        h0_min              = 0.25
        h0_max              = 0.75
        dist_h0             = 100
        pace_h0             = 0.25
        iter_h0             = 1
        tol_h0              = 0.01
        max_iter_h0         = 50
        stop                = 0
            
            
        # while dist_h0 > tol_h0 && iter_h0 < max_iter_h0 && stop == 0
        #     par_est.meanh0  = 0.5*h0_min + 0.5*h0_max;
        #     % Parameters outside estimation:
        #     [par]                                   = parameters(par_est,options);    % Load parameters
        #     par.final_iter                          = 'N';
        #     %% Solve policy functions:
        #     switch options.calib_steps
        #         case {'Y'}
        #             fprintf('solve model, time: %3.2f sec\n',toc(ticcalib));
        #     end
            
        #     [par,pol,guess,flags_solve]             = solve_model(par,options); %#ok<*ASGLU,*NASGU>
        #     switch flags_solve.converge
        #         case 'Y'
        #             switch options.save_mat
        #                 case 'Y'
        #                     filename = strcat('guess_',num2str(options.soboln,'%i'));
        #                     save(sprintf('%s.mat',filename),'guess');
        #             end
                    
        #             switch options.save_guess
        #                 case {'Y'}
        #                     filename = strcat('guess');
        #                     save(sprintf('%s.mat',filename),'guess');
        #                     parentpath    = cd(cd('..'));
        #                     filename      = strcat(parentpath,'/RESULTS/','guess_',num2str(options.n_fold,'%i'));
        #                     save(sprintf('%s.mat',filename),'guess');
        #                     if options.soboln >=2
        #                         options.save_guess = 'N';
        #                     end
        #             end
                    
        #             clear guess
                
                
                
        #             % 1. Interpolate policy functions into dense grid
        #             [pol_dense,par,par2]                    = solve_ergodic_distribution1(par,pol,options);
        #             clear pol
        #             % 2. Find ergodic distribution across cohorts
        #             [mu,mu_ige0,muc,pop_final,Q_ergo]       = solve_ergodic_distribution2(par2,pol_dense,options);
    
        #             % 3. Find distribution for a given cohort
        #             %     [mu,mu_cs,dist_age]                     = solve_ergodic_distribution3(par2,pol_dense,mu,pop_final,Q_ergo,options);
    
        #             me_inc_42      = mean_income_age42( par2,mu );
        #             dif_h0         = me_inc_42 - 2;
        #             dist_h0        = abs(dif_h0);
        #             if dist_h0 > tol_h0
        #                 if dif_h0 > 0
        #                     h0_max = pace_h0*h0_max+(1-pace_h0)*par_est.meanh0;
        #                 else
        #                     h0_min = pace_h0*h0_min+(1-pace_h0)*par_est.meanh0;
        #                 end
        #                 fprintf('iter = %i: \t h0 = %3.3f dist = %3.3f \t new h0 = [%3.4f,%3.4f] \n ',iter_h0,par_est.meanh0,dist_h0,h0_min,h0_max);
        #             else
        #                 fprintf('CONVERGED H0 \n iter = %i: \t h0 = %3.3f dist = %3.3f \t \n ',iter_h0,par_est.meanh0,dist_h0);
        #             end
        #             case 'N'
        #                 h0_max  = h0_max + 0.05*rand;
        #                 h0_min  = h0_min - 0.05*rand;
        #                 fprintf('iter = %i: \t Solve Failed. h0 = %3.3f dist = %3.3f \t new h0 = [%3.4f,%3.4f] \n ',iter_h0,par_est.meanh0,dist_h0,h0_min,h0_max);
        #         end
        #         iter_h0 = iter_h0 + 1;
        #     end
        #     x(14)                =  par_est.meanh0;
        #     par2.final_iter      = 'Y';
            
    else:
        par_est.meanh0      = x[13,1]
        par                                  = parameters(par_est,options)# Load parameters
        par.final_iter                          = 'Y'
        ## Solve policy functions:
        if options.calib_steps == 'Y':
            print(f'solve model, time: %3.2f sec\n',toc(ticcalib))

        
        par,pol,guess,flags_solve            = solve_model(par,options) #ok<*ASGLU,*NASGU>
        if flags_solve.converge == 'Y':
            if options.save_mat = 'Y'
            filename = str('guess_'+ str(options.soboln))
            scipy.io.savemat(str(filename + '.mat'), {str(filename):guess_1}) ##########whilc file are we saveing???

                
                if options.save_guess =='Y':
                    filename = str('guess')
                    scipy.io.savemat(str(filename + '.mat'), {"guess":guess})#save variable guess as guess.mat
                    parentpath    = str('/Users/Giorgia/Dropbox/ReplicationProject/DairuchPaper/Replication/Files_replicated_Python')
                    filename      = str(parentpath + '/RESULTS/' +'guess_'+ str(options.n_fold))
                    scipy.io.savemat(str(filename + '.mat'), {str(filename):guess})
                    if options.soboln >=2
                        options.save_guess = 'N'

                
                del guess

        
        # 1. Interpolate policy functions into dense grid
        pol_dense,par,par2                    = solve_ergodic_distribution1(par,pol,options)
        del pol
