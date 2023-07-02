# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 13:17:46 2022

@author: Giorgia
"""
import time
import itertools as it
import numpy as np
import scipy
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

options.exerciseSS_On   = 'N'         # Default = N; if = Y: solve SS exercise (different wage and savings options)
options.Fertility       = 'Endo'       # 'Endo': endogenous fertility; 'Exo': constant fertility
options.Fertility_Exo_N = 2.15       # 'Exo' case: number of children
options.ParentTrans     = 'Endo'      # Endo: Benchmark case with parents transfers; 'Exo': exogenously given
options.ParentTransVal  = 0.05         # Value for exogenous trasnfers
options.ExPolRetirement = 'N'         # N: benchmark case; Y: Exercise Retirement  icy
options.GridH0orig      = 'N'        # Force initial grid of H0 -- To use when doing sensitivity analysis of parameters
options.patient         = 'N'         # 'N': Baseline Discount Factor; 'Y': Higher beta (lower discounting)
options.noborrowingwedge= 'N'         # 'N': Baseline Wedge between borrowing and lending; 'Y': Smaller wedge
options.fert_transfer   = 'N'
options.init_transfer   = 'N'
options.psy_avg         = 'N'
options.AdultRisk       = 'Y' #from start
options.timer_on    	= 'N';          # display timer
options.start_Calibrate = 'N'
options.AdultRisk       = 'Y'
options.Ck              = 'Yes'
options.maxiter         = 25
options.tolV            = 1e-3
options.solve_h0        = 'N'
options.maxiterD        = 100;  # max iter for demographics
options.AdultRisk       = 'Y'
options.calib_steps     = 'Y'
options.start_Calibrate = 'N'
options.timer_on        = 'Y'
options.save_mat        = 'N'
options.save_guess      ='N'
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
            print(f'solve model, time: {toc()} sec\n')

        
        pol,par,guess,flags_solve            = solve_model(par,options) #ok<*ASGLU,*NASGU>
        if flags_solve.converge == 'Y':
            if options.save_mat == 'Y':
                parentpath1    = str('/Users/Giorgia/Dropbox/ReplicationProject/DairuchPaper/Replication/Files_replicated_Python/')
                filename1    = str(parentpath1 + '/RESULTS/' +'guess_'+ str(options.soboln))
                scipy.io.savemat(str(filename1 + '.mat'), {str(filename1):guess_}) ##########whilc file are we saveing???

                
                if options.save_guess =='Y':
                    filename = str('guess')
                    scipy.io.savemat(str(filename + '.mat'), {"guess":guess})#save variable guess as guess.mat
                    parentpath1    = str('/Users/Giorgia/Dropbox/ReplicationProject/DairuchPaper/Replication/Files_replicated_Python')
                    filename1    = str(parentpath1 + '/RESULTS/' +'guess_'+ str(options.n_fold))
                    scipy.io.savemat(str(filename1 + '.mat'), {str(filename1):guess})
                    if options.soboln >=2:
                        options.save_guess = 'N'

                
                #del guess

        
        # 1. Interpolate policy functions into dense grid
        pol_dense,par,par2                    = solve_ergodic_distribution1(par,pol,options)
        #del pol
        
    if flags_solve.converge == 'N':
    #options.n_fold = 1;
    #options.soboln = 1; #in start_calibrate
        print(f'ERROR HERE: solve_model did not converge \n')
        error_type  = 2
        parentpath    = str('/Users/Giorgia/Dropbox/ReplicationProject/DairuchPaper/DairuchFiles')
        parentpath1    = str('/Users/Giorgia/Dropbox/ReplicationProject/DairuchPaper/Replication/Files_replicated_Python')
        filename    = str(parentpath + '/RESULTS/' +'errors_'+ str(options.n_fold))
        filename1    = str(parentpath1 + '/RESULTS/' +'errors_'+ str(options.n_fold))
        Y           = loadmat(str(filename + '.mat'))
        aux         = np.append(np.append(np.reshape(x_val,[1,len(x_val)]),options.n_fold), np.append(options.soboln, error_type))
        error_x     = np.append(Y["error_x"], np.reshape(aux,[1,len(aux)]))
        scipy.io.savemat(str(filename1 + '.mat'), {str(filename1):error_x})

    elif flags_solve.converge == 'Y':
        mu,mu_ige0,muc,pop_final,Q_ergo       = solve_ergodic_distribution2(par2,pol_dense,options)
    
    if options.save_mat == 'Y':
        model.pol_sel = SimpleNamespace()
        pol_sel = model.pol_sel
        pol_sel.Np      = copy.deepcopy(pol_dense.Np)
        pol_sel.mu_fert = copy.deepcopy(mu[par2.Jc_pos-1][0])#########??????
        filename1 = str(parentpath1 + '/RESULTS/' + 'PolicyFunctionsFertility_'+ str(options.soboln))
        PolicyFunctionsFertility = {"par2": par2, "pol_sel":pol_sel,"options" : options}
        scipy.io.savemat(str(filename1 + '.mat'), {str(filename1): PolicyFunctionsFertility})

    flag.pop_final       = pop_final
    
    if flag.pop_final == 0: # stop
        print(f'Population = 0, go to next set of parameters\n')
        m_data,~,~ = moments_data######################################to be coded
        Nm          = size(m_data,1)
        m_model     = np.zeros(Nm)
        m_model[:] = np.nan
    else :
        #3. Find distribution for a given cohort
        mu,mu_cs,dist_age                     = solve_ergodic_distribution3(par2,pol_dense,mu,pop_final,Q_ergo,options)
        # 4. Intergenerational distributions
        options.Nsample                 = 300000
#             switch options.ComputeOtherMus
#                 case 'Y'
#                     options.Nsample                 = 300000;
             
        IGE                           = gen_ige_sample(par2,pol_dense,muc,options)################################àto code
        
        mu_ige,mu_chetty,Inc_par_grid,mu_educ_ige,mu_par_ic = trans_ige(par2,mu_ige0,pol_dense.tau0,pol_dense.S,pol_dense.Se,options)
        
        if options.save_mat == 'Y':
            filename = str('dist_LE')
            dist_LE = {'par2': par2,'mu': mu}
            scipy.io.savemat(str(filename + '.mat'),'par2','mu')
                
        if options.ComputeMoments == 'Y':
            ## Compute moments:
            options.calib_steps == 'Y'
                print(f'compute moments, time: toc() sec\n')
                #             save temp.mat
                m_model,m_title,m_mod_id      = moments_model(par2,pol_dense,mu,mu_cs,mu_ige,mu_chetty,Inc_par_grid,mu_educ_ige,dist_age,IGE,options)###########################àto code

    
   # Pack-up output: Results:
    filename    = str(pwd,'/results/'+'results_'+ str(options.n_fold))
    Y           = loadmat(str(filename + '.mat'))
    aux         = np.append(np.append(np.append(np.reshape(x_val,[1,len(x_val)]),options.n_fold), np.append(options.soboln, error_type)), np.append(par.p_model_data flags_solve.iter))
    res_x       = np.append(Y["error_x"], np.reshape(aux,[1,len(aux)]))
    res_MOM     = np.append(Y["res_MOM"], np.reshape(m_model,[1,len(m_model)]))
    dict_ = {res_x : 'res_x', res_MOM :'res_MOM'}
    scipy.io.savemat(str(filename + '.mat'), {str(filename): dict_})
    
    print(f'Calibration time: toc() sec\n')

        
