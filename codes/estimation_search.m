%% Estimation of model using fminsearch
datestr(now)
pwd
warning('off','MATLAB:nearlySingularMatrix')

%% Options
options.print_moms      = 'N';          % display moments
options.timer_on    	= 'N';          % display timer
options.model_iter_on   = 'N';          % display model solve iterations
options.calib_steps     = 'N';          % display calibration steps
options.save_mat    	= 'N';          % save policy functions + LE_dist
options.save_guess	    = 'Y';          % save guess
options.save_dist 	    = 'Y';          % save distribution
options.start_Calibrate = 'N';          % Default  = N: use standard initial guess; Y: load specific initial guess
options.maxiter         = 50;           % max iter for policy functions
options.tolV            = 1e-3;         % tolerance for policy functions
options.maxiterD        = 100;           % max iter for demographics
options.tolD            = 1e-6;         % tolerance for demographics
options.Nsim            = 1e5;          % Number of simulations
options.exerciseSS.On   = 'N';          % Default = N; if = Y: solve SS exercise (different wage and savings options)
options.ComputeMoments  = 'Y';          % Y: compute moments
options.transfers       = 'On';         % 'On': with OAS; 'Reduced'; 'Off'
options.ChildCost       = 'OppC'; 		% 'OC': opportunity cost or 'Constant'
options.Fertility       = 'Endo';       % 'Endo': endogenous fertility; 'Exo': constant fertility
options.Ck              = 'Yes';        % 'Yes': children consume at home; 'No'
options.ParentTrans     = 'Endo';       % Endo: Benchmark case with parents transfers; 'Exo': exogenously given
options.ParentTransVal  = .05;          % Value for exogenous trasnfers
options.AdultRisk       = 'Y';          % Y: benchmark case; N: no adult risk
options.ExPolRetirement = 'N';          % N: benchmark case; Y: Exercise Retirement Policy
options.PolRetMult      = 1;            % multiplier of retirment benefit in first bracket. Benchmark = 1;
options.ComputeOtherMus = 'N';
options.ComputeLE_age   = 'N';
options.equilibrium     = 'partial';    % 'partial' or 'general'
options.MeanRevEx       = 'Off';
options.GridH0orig      = 'N';
options.psy_avg         = 'N';
options.fert_transfer   = 'N';
options.init_transfer   = 'N';
options.patient         = 'N';
options.noborrowingwedge = 'N';
options.solve_h0        = 'Y';



options.data_year   = 'data_2000';      % data_2000 data_1960 data_fert_1960



%% Load x0:
options.xnames          = {'gamman', 'ln_level', 'sigmah0', 'psy_col' , 'cor_psy', 'psy_hs', 'R_hs_curv', 'R_col_curv'  , 'R_hs_level'  ,'R_col_level' , 'Child_C_Tot', 'child_cost_inc_curv'   , 'beta_inc'};
xlb                     = [ 0.010  , 0.100     , 0.010    , 0.010     , 1.000    , 0.010   , 0.010     , 0.010          , 0.010         , 0.01         , 0.010        , 0.010     , 0.05        ];
xub                     = [ 1.000  , 1.000     , 0.800    , 80.000    , 5.000    , 15.000   , 2.00      , 2.00           , 1.000         , 1.00         , 1.000        , 1.000     , 0.99       ];


parentpath  = fileparts(pwd);
load([parentpath,'/xmin_mat.mat']);
x0          = xmin_mat(options.n_fold,:);

%% Create Storage file
%1. Initiate storage
res_x                   = [];
res_MOM                 = [];
error_x                 = [];
% 2. Create 'result'.mat files
parentpath              = cd(cd('..'));
filename                = strcat(parentpath,'/RESULTS/','results_',num2str(options.n_fold,'%i'));
save(filename,'res_x','res_MOM');
clear res*

parentpath              = cd(cd('..'));
filename_diary          = strcat(parentpath,'/RESULTS/','diary_',num2str(options.n_fold,'%i'),'.txt');
diary(filename_diary);

parentpath              = cd(cd('..'));
filename                = strcat(parentpath,'/RESULTS/','errors_',num2str(options.n_fold,'%i'));
save(filename,'error_x');
clear error_x

%% Call fminsearch
fmin_opt    = optimset('Display','iter');
fval        = @(x) eval_model(x,xlb,xub,options);

[x,fval,exitflag,output] = fminsearch(fval,x0,fmin_opt);


%% Save Results
parentpath  = cd(cd('..'));
filename    = strcat(parentpath,'/RESULTS/','final_results_',num2str(options.n_fold,'%i'));
save(filename,'x','output');
