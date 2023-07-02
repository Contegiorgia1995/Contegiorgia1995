
%% Commands to run in personal computer:
clear; close all; clc; tic;
dbstop if error
addpath([pwd,'/run_generic']);
options.n_fold = 111121;
ex_version = 1;
warning('off','MATLAB:nearlySingularMatrix')

%% Generic Options
options.print_moms      = 'Y';          % display moments
options.timer_on    	= 'Y';          % display timer
options.model_iter_on   = 'Y';          % display model solve iterations
options.calib_steps     = 'Y';          % display calibration steps
options.save_mat    	= 'N';          % save policy functions
options.save_guess	    = 'N';          % save guess
options.save_dist 	    = 'N';          % save distribution
options.start_Calibrate = 'N';          % Default  = N: use standard initial guess; Y: load specific initial guess
options.maxiter         = 50;           % max iter for policy functions
options.tolV            = 1e-3;         % tolerance for policy functions
options.maxiterD        = 100;           % max iter for demographics
options.tolD            = 1e-6;         % tolerance for demographics
options.Nsim            = 1e5;          % Number of simulations
options.exerciseSS.On   = 'N';          % Default = N; if = Y: solve SS exercise (different wage and savings options)
options.ComputeMoments  = 'Y';          % Y: compute moments
options.ExPolRetirement = 'N';          % N: benchmark case; Y: Exercise Retirement Policy
options.PolRetMult      = 1;            % multiplier of retirment benefit in first bracket. Benchmark = 1;
options.ComputeOtherMus = 'Y';          % Compute time-consuming moments
options.ComputeLE_age   = 'Y';          % Compute time-consuming Lifetime Earnings moments
options.transfers       = 'On';         % 'On': with OAS; 'Reduced'; 'Off'
options.ChildCost       = 'OppC';       % 'OC': opportunity cost or 'Constant'
options.Fertility       = 'Endo';       % 'Endo': endogenous fertility; 'Exo': constant fertility
options.Ck              = 'Yes';        % 'Yes': children consume at home; 'No'
options.ParentTrans     = 'Endo';       % Endo: Benchmark case with parents transfers; 'Exo': exogenously given
options.AdultRisk       = 'Y';          % Y: benchmark case; N: no adult risk
options.GridH0orig      = 'N';          % Force initial grid of H0 -- To use when doing sensitivity analysis of parameters
options.psy_avg         = 'N';          % To have only one level of school taste 
options.fert_transfer   = 'N';          % 'N': Baseline; 'only2': For policy analysis of giving money to families with two children
options.init_transfer   = 'N';          % 'N': Baseline; 'Y': For policy analysis of giving money to all agents when they are 16; 'OnlyCollege': For policy analysis of giving money to all agents when they are 16 if they go to college
options.patient         = 'N';          % 'N': Baseline Discount Factor; 'Y': Higher beta (lower discounting)
options.noborrowingwedge= 'N';          % 'N': Baseline Wedge between borrowing and lending; 'Y': Smaller wedge
options.solve_h0        = 'N';          % 'Y': Solve such that normalization of income holds; 'N': Do not solve (e.g., when it has already been solved or for policy exercises)

filename                = [pwd,'\results\xmin.mat'];
load(filename)

options.soboln  = 1;
options.continue  = 'N';
VLE_IC_byage(x,options);
