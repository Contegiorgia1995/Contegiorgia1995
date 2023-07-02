clear; close all; clc; tic;
dbstop if error
addpath([pwd,'/run_generic']);
options.n_fold = 111121;
ex_version = 1;

%% Generic Options
options.print_moms      = 'Y';          % display moments
options.timer_on    	= 'N';          % display timer
options.model_iter_on   = 'N';          % display model solve iterations
options.calib_steps     = 'Y';          % display calibration steps
options.save_mat    	= 'Y';          % save policy functions
options.save_guess	    = 'Y';          % save guess
options.save_dist 	    = 'Y';          % save distribution
options.start_Calibrate = 'N';          % Default  = N: use standard initial guess; Y: load specific initial guess
options.maxiter         = 50;           % max iter for policy functions
options.tolV            = 1e-3;         % tolerance for policy functions
options.maxiterD        = 50;           % max iter for demographics
options.tolD            = 1e-6;         % tolerance for demographics
options.Nsim            = 1e5;          % Number of simulations
% options.exerciseSS.On   = 'Y';          % Default = N; if = Y: solve SS exercise (different wage and savings options)
options.ComputeMoments  = 'Y';          % Y: compute moments
options.transfers       = 'On';         % 'On': with OAS; 'Reduced'; 'Off'
options.ChildCost       = 'OppC';       % 'OC': opportunity cost or 'Constant'
options.Fertility       = 'Endo';       % 'Endo': endogenous fertility; 'Exo': constant fertility
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

%% LOAD ESTIMATED PARAMETERS
filename                = [pwd,'\results\xmin.mat'];
load(filename)

options.exerciseSS.change_pe    = 'N';
options.exerciseSS.change_psy   = 'Y';

%% Find Psychic Cost 
if 1 == 1
    options.exerciseSS.On           = 'Y';          % Default = N; if = Y: solve SS exercise (different wage and savings options)
    options.save_guess              = 'Y';          % so that it is not replaced by 'N' in calibrate, line 51

    options.exerciseSS.file         = 1;
    options.soboln                  = 1;
    options.exerciseSS.w_ss         = 0.80;
    options.exerciseSS.p_mult       = options.exerciseSS.w_ss;
    HSdropout1960_change            = 0.1080564 - 0.094;
    HSdropout2000_model             = 0.1;

    xmin                            = 0.200;
    xmax                            = 0.800; 
    dist                            = 1;
    tol                             = 0.001;
    
    while dist > tol
        x0  = (xmin + xmax)/2;
        options.psy_adj = x0;
        fprintf('Solving X = %3.3f \n',x0);

        %% Initiate transitory storage
        res_x                   = [];
        res_MOM                 = [];
        error_x                 = [];
        filename                = strcat(pwd,'/results/','results_',num2str(options.n_fold,'%i'));
        filename_diary          = strcat(pwd,'/results/','diary_',num2str(options.n_fold,'%i'),'.txt');
        save(filename,'res_x','res_MOM');
        diary(filename_diary);
        clear res*

        filename                = strcat(pwd,'/results/','errors_',num2str(options.n_fold,'%i'));
        save(filename,'error_x');
        clear error_x
        calibrate(x,options);

        filename                = strcat(pwd,'/results/','results_',num2str(options.n_fold,'%i'));
        X                       = load(filename);
        HSdropout_model_change  = X.res_MOM(1,3) - HSdropout2000_model;


        fprintf('X = %3.3f, \t data = %3.3f, \t model = %3.3f \n',x0,HSdropout1960_change,HSdropout_model_change);
        dif                     = HSdropout_model_change - HSdropout1960_change;
        dist                    = abs(dif);

        if dif > tol 
            % Too many dropouts -> Need to reduce psychic cost more -> X needs to be larger
            xmin  = x0;
        elseif dif < - tol
            % Too few dropouts -> Need to increase psychic cost more -> X needs to be smaller
            xmax  = x0;
        end
        if dist > tol
            fprintf('dist = %3.3f, \t new xlim = [%3.3f,%3.3f] \n \n',dist,xmin,xmax);
        else
            fprintf('CONVERGED: dist = %3.3f, \t x = %3.3f \n \n',dist,x0);
        end
    end
    
else
    options.psy_adj                 = 0.477;    %243: .511, 169: .425, 109: .399
end

%% Initiate transitory storage
res_x                   = [];
res_MOM                 = [];
error_x                 = [];

filename                = strcat(pwd,'/results/','results_',num2str(options.n_fold,'%i'));
filename_diary          = strcat(pwd,'/results/','diary_',num2str(options.n_fold,'%i'),'.txt');
save(filename,'res_x','res_MOM');
diary(filename_diary);
clear res*

filename                = strcat(pwd,'/results/','errors_',num2str(options.n_fold,'%i'));
save(filename,'error_x');
clear error_x


%% SOLVE FOR DIFFERENT WAGES
options.exerciseSS.On           = 'Y';          % Default = N; if = Y: solve SS exercise (different wage and savings options)
% options.exerciseSS.change_pe    = 'N';
% options.exerciseSS.change_psy   = 'Y';

W   = 0.2:0.02:1.8;
GDP = zeros(size(W));

for iw = 1:length(W)
    fprintf('wage = %2.2f \n',W(iw));
	try
		%1. Solve policy functions and ergodic distribution for w
		options.save_guess	        = 'Y';          % so that it is not replaced by 'N' in calibrate, line 51
		options.soboln              = iw;
		options.exerciseSS.file     = iw;
		options.exerciseSS.w_ss     = W(iw);
		options.exerciseSS.p_mult   = W(iw);
		calibrate(x,options);
	catch
	end
end

%% PUT WAGES AND ELASTICITIES TOGETHER
filename      = strcat(pwd,'/results/','results_',num2str(options.n_fold,'%i'));
load(filename,'res_x','res_MOM');

filename_destination    = strcat(pwd,'/results/','results_compare_ss_',num2str(ex_version),'.mat');
save(filename_destination,'res_x','res_MOM','W','options');