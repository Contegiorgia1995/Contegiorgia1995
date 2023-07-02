%% Start
clear; close all; clc; tic;
dbstop if error

profile off
% profile on
profile clear
% profile -memory on

warning('off','MATLAB:dispatcher:nameConflict')
warning('off','MATLAB:nearlySingularMatrix')
addpath([pwd,'/run_generic']);
options.n_fold = 1;
options.soboln = 1;
datestr(now)
pwd
warning('off','MATLAB:nearlySingularMatrix')

%% Options
options.print_moms      = 'Y';          % display moments
options.timer_on    	= 'Y';          % display timer
options.model_iter_on   = 'Y';          % display model solve iterations
options.calib_steps     = 'Y';          % display calibration steps
options.RecycleSobol    = 'N';          % If RecycleSobol = Y, start with sobolmin from previous run
options.save_mat    	= 'N';          % save policy functions + LE_dist
options.save_guess	    = 'N';          % save guess
options.save_dist 	    = 'N';          % save distribution
options.start_Calibrate = 'N';          % Default  = N: use standard initial guess; Y: load specific initial guess
options.maxiter         = 25;           % max iter for policy functions
options.tolV            = 1e-3;         % tolerance for policy functions
options.maxiterD        = 100;          % max iter for demographics
options.tolD            = 1e-6;         % tolerance for demographics
options.Nsim            = 1e5;          % Number of simulations
options.exerciseSS.On   = 'N';          % Default = N; if = Y: solve SS exercise (different wage and savings options)
options.ComputeMoments  = 'Y';          % Y: compute moments
options.transfers       = 'On';         % 'On': with OAS; 'Reduced'; 'Off'
options.ChildCost       = 'OppC'; 		% 'OC': opportunity cost or 'Constant'
options.Fertility       = 'Endo';       % 'Endo': endogenous fertility; 'Exo': constant fertility
options.Fertility_Exo_N = 2.15;         % 'Exo' case: number of children
options.Ck              = 'Yes';        % 'Yes': children consume at home; 'No'
options.ParentTrans     = 'Endo';       % Endo: Benchmark case with parents transfers; 'Exo': exogenously given
options.ParentTransVal  = .05;          % Value for exogenous trasnfers
options.AdultRisk       = 'Y';          % Y: benchmark case; N: no adult risk
options.ExPolRetirement = 'N';          % N: benchmark case; Y: Exercise Retirement Policy
options.PolRetMult      = 1;            % multiplier of retirment benefit in first bracket. Benchmark = 1;
options.ComputeOtherMus = 'N';          % Compute time-consuming moments
options.ComputeLE_age   = 'N';          % Compute time-consuming Lifetime Earnings moments
options.GridH0orig      = 'N';          % Force initial grid of H0 -- To use when doing sensitivity analysis of parameters
options.psy_avg         = 'N';          % To have only one level of school taste 
options.fert_transfer   = 'N';          % 'N': Baseline; 'only2': For policy analysis of giving money to families with two children
options.init_transfer   = 'N';          % 'N': Baseline; 'Y': For policy analysis of giving money to all agents when they are 16; 'OnlyCollege': For policy analysis of giving money to all agents when they are 16 if they go to college
options.patient         = 'N';          % 'N': Baseline Discount Factor; 'Y': Higher beta (lower discounting)
options.noborrowingwedge= 'N';          % 'N': Baseline Wedge between borrowing and lending; 'Y': Smaller wedge
options.solve_h0        = 'N';          % 'Y': Solve such that normalization of income holds; 'N': Do not solve (e.g., when it has already been solved or for policy exercises)

%1. Initiate storage
res_x                   = [];
res_MOM                 = [];
error_x                 = [];
% 2. Create 'result'.mat files
filename                = [pwd,'/results/','results_',num2str(options.n_fold,'%i')];
save(filename,'res_x','res_MOM');
clear res*

filename                = [pwd,'/results/','errors_',num2str(options.n_fold,'%i')];
save(filename,'error_x');
clear error_x

%% Solve calibration
filename                = [pwd,'\results\xmin.mat'];
load(filename)
calibrate(x,options);
