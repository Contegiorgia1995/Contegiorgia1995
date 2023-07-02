%% Commands to run in personal computer:
clear; close all; clc; tic;
dbstop if error
addpath([pwd,'/run_generic']);
options.n_fold = 111131;
ex_version = 1;
warning('off','MATLAB:nearlySingularMatrix')


%% Solve models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.data_year   = 'data_2000';      % data_2000 data_1960 data_fert_1960

if 1 == 1
    %% Initiate storage
    res_x                   = [];
    res_MOM                 = [];
    error_x                 = [];
    % 2. Create 'result'.mat files
    filename_destination    = strcat(pwd,'/results/','results_model_options_',options.data_year,'_',num2str(ex_version),'.mat');
    save(filename_destination,'res_x','res_MOM');
    clear res*
    
    parentpath              = cd(cd('..'));
    filename                = strcat(pwd,'/results/','errors_',num2str(options.n_fold,'%i'));
    save(filename,'error_x');
    clear error_x
end

%% Generic Options
options.print_moms      = 'Y';          % display moments
options.timer_on    	= 'N';          % display timer
options.model_iter_on   = 'Y';          % display model solve iterations
options.calib_steps     = 'Y';          % display calibration steps
options.save_mat    	= 'N';          % save policy functions
options.save_guess	    = 'N';          % save guess
options.save_dist 	    = 'N';          % save distribution
options.start_Calibrate = 'N';          % Default  = N: use standard initial guess; Y: load specific initial guess
options.maxiter         = 200;           % max iter for policy functions
options.tolV            = 1e-3;         % tolerance for policy functions
options.maxiterD        = 100;           % max iter for demographics
options.tolD            = 1e-6;         % tolerance for demographics
options.Nsim            = 1e5;          % Number of simulations
options.exerciseSS.On   = 'N';          % Default = N; if = Y: solve SS exercise (different wage and savings options)
options.ComputeMoments  = 'Y';          % Y: compute moments
% options.transfers       = 'On';       % 'On': with OAS; 'Reduced'; 'Off'
% options.ChildCost       = 'OppC';     % 'OC': opportunity cost or 'Constant'
% options.Fertility       = 'Endo';       % 'Endo': endogenous fertility; 'Exo': constant fertility
% options.Ck              = 'Yes';        % 'Yes': children consume at home; 'No'
% options.ParentTrans     = 'Endo';       % Endo: Benchmark case with parents transfers; 'Exo': exogenously given
% options.ParentTransVal  = .05;          % Value for exogenous trasnfers
% options.AdultRisk       = 'Y';          % Y: benchmark case; N: no adult risk
options.ExPolRetirement = 'N';          % N: benchmark case; Y: Exercise Retirement Policy
options.PolRetMult      = 1;            % multiplier of retirment benefit in first bracket. Benchmark = 1;
options.ComputeOtherMus = 'Y';          % Compute time-consuming moments
options.ComputeLE_age   = 'N';          % Compute time-consuming Lifetime Earnings moments
options.fert_transfer   = 'N';          % 'N': Baseline; 'only2': For policy analysis of giving money to families with two children
options.init_transfer   = 'N';          % 'N': Baseline; 'Y': For policy analysis of giving money to all agents when they are 16; 'OnlyCollege': For policy analysis of giving money to all agents when they are 16 if they go to college
options.patient         = 'N';          % 'N': Baseline Discount Factor; 'Y': Higher beta (lower discounting)
options.noborrowingwedge= 'N';          % 'N': Baseline Wedge between borrowing and lending; 'Y': Smaller wedge
options.solve_h0        = 'N';          % 'Y': Solve such that normalization of income holds; 'N': Do not solve (e.g., when it has already been solved or for policy exercises)



for num_model = [1 4:5 7:14] %[ 1 5 12]: Evaluate Endogenous Forces (table 6); [ 1 7:10 13:14]: Evaluate Exogenous Forces (table 7); [ 1 15:16]: Evaluate Policies (table 8)
    filename                = [pwd,'\results\xmin.mat'];
    load(filename)
    
    if num_model == 1 % Baseline model
        exname                  = 'Base';
        options.transfers       = 'On';
        options.ChildCost       = 'OppC';       % 'OC': opportunity cost or 'Constant'
        options.Fertility       = 'Endo';       % 'Endo': endogenous fertility; 'Exo': constant fertility
        options.Ck              = 'Yes';        % 'Yes': children consume at home; 'No'
        options.ParentTrans     = 'Endo';       % Endo: Benchmark case with parents transfers; 'Exo': exogenously given
        options.AdultRisk       = 'Y';          % Y: benchmark case; N: no adult risk
        options.GridH0orig      = 'N';
        options.psy_avg         = 'N';
        options.fert_transfer   = 'N';
        options.init_transfer   = 'N';
    
    elseif num_model == 4 % Exogenous Fertility: 2
        exname                  = 'ExoFert2';
        options.transfers       = 'On';
        options.ChildCost       = 'OppC';       % 'OC': opportunity cost or 'Constant'
        options.Fertility       = 'Exo';        % 'Endo': endogenous fertility; 'Exo': constant fertility
        options.Fertility_Exo_N = 2.00;
        options.Ck              = 'Yes';        % 'Yes': children consume at home; 'No'
        options.ParentTrans     = 'Endo';       % Endo: Benchmark case with parents transfers; 'Exo': exogenously given
        options.AdultRisk       = 'Y';          % Y: benchmark case; N: no adult risk
        options.GridH0orig      = 'N';
        options.psy_avg         = 'N';
        options.fert_transfer   = 'N';
        options.init_transfer   = 'N';
        
    elseif num_model == 5 % Exogenous Fertility: 2.15
        exname                  = 'ExoFert215';
        options.transfers       = 'On';
        options.ChildCost       = 'OppC';       % 'OC': opportunity cost or 'Constant'
        options.Fertility       = 'Exo';        % 'Endo': endogenous fertility; 'Exo': constant fertility
        %%% SEE BELOW FOR N
        %         options.Fertility_Exo_N = 2.15;
        options.Ck              = 'Yes';        % 'Yes': children consume at home; 'No'
        options.ParentTrans     = 'Endo';       % Endo: Benchmark case with parents transfers; 'Exo': exogenously given
        options.AdultRisk       = 'Y';          % Y: benchmark case; N: no adult risk
        options.GridH0orig      = 'N';
        options.psy_avg         = 'N';
        options.fert_transfer   = 'N';
        options.init_transfer   = 'N';
        
        
    elseif num_model == 7  % Constant H0
        exname                  = 'ConstantH0';
        options.GridH0orig      = 'Y';
        options.SigmaH0orig     = x(3);
        x(3)                    = 0.0001; % H0 variance
        x(13)                   = 0; % H0 persistence
        options.transfers       = 'On';
        options.ChildCost       = 'OppC';       % 'OC': opportunity cost or 'Constant'
        options.Fertility       = 'Endo';       % 'Endo': endogenous fertility; 'Exo': constant fertility
        options.Ck              = 'Yes';        % 'Yes': children consume at home; 'No'
        options.ParentTrans     = 'Endo';       % Endo: Benchmark case with parents transfers; 'Exo': exogenously given
        options.AdultRisk       = 'Y';          % Y: benchmark case; N: no adult risk
        options.psy_avg         = 'N';
        options.fert_transfer   = 'N';
        options.init_transfer   = 'N';
        
    elseif num_model == 8  % H0 iid
        exname                  = 'H0iid';
        %         x(3)                    = 0.0001; % H0 variance
        options.GridH0orig      = 'N';
        x(13)                   = 0; % H0 persistence
        options.transfers       = 'On';
        options.ChildCost       = 'OppC';       % 'OC': opportunity cost or 'Constant'
        options.Fertility       = 'Endo';       % 'Endo': endogenous fertility; 'Exo': constant fertility
        options.Ck              = 'Yes';        % 'Yes': children consume at home; 'No'
        options.ParentTrans     = 'Endo';       % Endo: Benchmark case with parents transfers; 'Exo': exogenously given
        options.AdultRisk       = 'Y';          % Y: benchmark case; N: no adult risk
        options.GridH0orig      = 'N';
        options.psy_avg         = 'N';
        options.fert_transfer   = 'N';
        options.init_transfer   = 'N';
        
        
    elseif num_model == 9  % H0 without uncertainty
        exname                  = 'CertainH0';
        options.GridH0orig      = 'Y';
        options.SigmaH0orig     = x(3);
        x(3)                    = 0.0001; % H0 constant
        options.transfers       = 'On';
        options.ChildCost       = 'OppC';       % 'OC': opportunity cost or 'Constant'
        options.Fertility       = 'Endo';       % 'Endo': endogenous fertility; 'Exo': constant fertility
        options.Ck              = 'Yes';        % 'Yes': children consume at home; 'No'
        options.ParentTrans     = 'Endo';       % Endo: Benchmark case with parents transfers; 'Exo': exogenously given
        options.AdultRisk       = 'Y';          % Y: benchmark case; N: no adult risk
        options.psy_avg         = 'N';
        options.fert_transfer   = 'N';
        options.init_transfer   = 'N';
        
    elseif num_model == 10 % No adult risk
        exname                  = 'NoAdultRisk';
        options.GridH0orig      = 'N';
        options.transfers       = 'On';
        options.ChildCost       = 'OppC';       % 'OC': opportunity cost or 'Constant'
        options.Fertility       = 'Endo';       % 'Endo': endogenous fertility; 'Exo': constant fertility
        options.Ck              = 'Yes';        % 'Yes': children consume at home; 'No'
        options.ParentTrans     = 'Endo';       % Endo: Benchmark case with parents transfers; 'Exo': exogenously given
        options.AdultRisk       = 'N';          % Y: benchmark case; N: no adult risk
        options.GridH0orig      = 'N';
        options.psy_avg         = 'N';
        options.fert_transfer   = 'N';
        options.init_transfer   = 'N';
        
    elseif num_model == 11 % Constant Transfers = Mean Transfer
        exname                  = 'ConstantTransfers_MeanTrans';
        options.transfers       = 'On';
        options.ChildCost       = 'OppC';       % 'OC': opportunity cost or 'Constant'
        options.Fertility       = 'Endo';       % 'Endo': endogenous fertility; 'Exo': constant fertility
        options.Ck              = 'Yes';        % 'Yes': children consume at home; 'No'
        options.ParentTrans     = 'Exo';        % Endo: Benchmark case with parents transfers; 'Exo': exogenously given
        options.AdultRisk       = 'Y';          % Y: benchmark case; N: no adult risk
        options.GridH0orig      = 'N';
        options.ParentTransExo  = 'search_val';
        options.psy_avg         = 'N';
        options.fert_transfer   = 'N';
        options.init_transfer   = 'N';   
        
    elseif num_model == 12 % Constant Transfers = Mean Hs Cost
        exname                  = 'ConstantTransfers_Mean';
        options.transfers       = 'On';
        options.ChildCost       = 'OppC';       % 'OC': opportunity cost or 'Constant'
        options.Fertility       = 'Endo';       % 'Endo': endogenous fertility; 'Exo': constant fertility
        options.Ck              = 'Yes';        % 'Yes': children consume at home; 'No'
        options.ParentTrans     = 'Exo';        % Endo: Benchmark case with parents transfers; 'Exo': exogenously given
        options.ParentTransExo  = 'search_val';
        options.ParentTransVal  = 1;            % Value for exogenous trasnfers
        options.AdultRisk       = 'Y';          % Y: benchmark case; N: no adult risk
        options.GridH0orig      = 'N';
        options.psy_avg         = 'N';
        options.fert_transfer   = 'N';
        options.init_transfer   = 'N';
        
        %             Load Target Transfer
        filename_destination    = strcat(pwd,'/results/','results_model_options_',options.data_year,'_',num2str(ex_version),'.mat');
        Y           = load(filename_destination);
        options.ParentTransVal  = Y.res_MOM(1,98)*Y.res_MOM(1,1);
        
        
    elseif num_model == 13 % Psychic Cost: iid
        exname                  = 'Psy_iid';
        x(5)                    = 1;            % H0 persistence
        options.transfers       = 'On';
        options.ChildCost       = 'OppC';       % 'OC': opportunity cost or 'Constant'
        options.Fertility       = 'Endo';       % 'Endo': endogenous fertility; 'Exo': constant fertility
        options.Ck              = 'Yes';        % 'Yes': children consume at home; 'No'
        options.ParentTrans     = 'Endo';       % Endo: Benchmark case with parents transfers; 'Exo': exogenously given
        options.AdultRisk       = 'Y';          % Y: benchmark case; N: no adult risk
        options.GridH0orig      = 'N';
        options.psy_avg         = 'N';
        options.fert_transfer   = 'N';
        options.init_transfer   = 'N';
        
    elseif num_model == 14 % Psychic Cost: Average
        exname                  = 'Psy_avg';
        options.transfers       = 'On';
        options.ChildCost       = 'OppC';       % 'OC': opportunity cost or 'Constant'
        options.Fertility       = 'Endo';       % 'Endo': endogenous fertility; 'Exo': constant fertility
        options.Ck              = 'Yes';        % 'Yes': children consume at home; 'No'
        options.ParentTrans     = 'Endo';       % Endo: Benchmark case with parents transfers; 'Exo': exogenously given
        options.AdultRisk       = 'Y';          % Y: benchmark case; N: no adult risk
        options.GridH0orig      = 'N';
        options.psy_avg         = 'Y';
        options.fert_transfer   = 'N';
        options.init_transfer   = 'N';
        
    elseif num_model == 15 % Fertility Transfer Policy
        exname                  = 'FertPolicy';
        options.transfers       = 'On';
        options.ChildCost       = 'OppC';       % 'OC': opportunity cost or 'Constant'
        options.Fertility       = 'Endo';       % 'Endo': endogenous fertility; 'Exo': constant fertility
        options.Ck              = 'Yes';        % 'Yes': children consume at home; 'No'
        options.ParentTrans     = 'Endo';       % Endo: Benchmark case with parents transfers; 'Exo': exogenously given
        options.AdultRisk       = 'Y';          % Y: benchmark case; N: no adult risk
        options.GridH0orig      = 'N';
        options.psy_avg         = 'N';
        options.fert_transfer   = 'only2';
        options.fert_trans_size = 20000;
        options.init_transfer   = 'N';
    elseif num_model == 16 % Initial Transfer Policy
        exname                  = 'TransferPolicy';
        options.transfers       = 'On';
        options.ChildCost       = 'OppC';       % 'OC': opportunity cost or 'Constant'
        options.Fertility       = 'Endo';       % 'Endo': endogenous fertility; 'Exo': constant fertility
        options.Ck              = 'Yes';        % 'Yes': children consume at home; 'No'
        options.ParentTrans     = 'Endo';       % Endo: Benchmark case with parents transfers; 'Exo': exogenously given
        options.AdultRisk       = 'Y';          % Y: benchmark case; N: no adult risk
        options.GridH0orig      = 'N';
        options.psy_avg         = 'N';
        options.fert_transfer   = 'N';
        options.init_transfer   = 'Y';
        options.init_trans_size = 20000;
    end
    
    %% Initiate transitory storage
    res_x                   = [];
    res_MOM                 = [];
    filename                = strcat(pwd,'/results/','results_',num2str(options.n_fold,'%i'));
    filename_diary          = strcat(pwd,'/results/','diary_',num2str(options.n_fold,'%i'),'.txt');
    save(filename,'res_x','res_MOM');
    diary(filename_diary);
    clear res*
    
    fprintf('Start %s, time = %3.0f sec. \n',char(exname),toc)
    
    filename = strcat('PolicyFunctions_',exname,'.mat');
    if exist(filename,'file')==2
        delete(filename)
    end
    filename = strcat('dist_ss',exname,'.mat');
    if exist(filename,'file')==2
        delete(filename)
    end
    
    filename = strcat('guess_',num2str(num_model),'.mat');
    if exist(filename,'file')==2
        copyfile(filename,'guess.mat','f')
    end
    
    %% Config exofert
    switch exname
        case 'ExoFert215'
            filename_destination    = strcat(pwd,'/results/','results_model_options_',options.data_year,'_',num2str(ex_version),'.mat');
            Y                       = load(filename_destination);
            mean_fert_baseline      = Y.res_MOM(1,6);
            options.Fertility_Exo_N = mean_fert_baseline;
    end
    
    %% Config ParentTrans
    if num_model == 11
        options.model_iter_on   = 'N';          % display model solve iterations
        options.calib_steps     = 'N';          % display calibration steps
        options.print_moms      = 'N';
        options.ComputeOtherMus = 'N';
        options.save_guess	    = 'Y';          % save guess
        
        %             Load Target Transfer
        filename_destination    = strcat(pwd,'/results/','results_model_options_',options.data_year,'_',num2str(ex_version),'.mat');
        Y           = load(filename_destination);
        mean_transfer_baseline = Y.res_MOM(1,98);
        
        %   Start search
        phi_min     = 20000;
        phi_max     = 50000;
        dist0       = 100;
        
        
        % Evaluate
        iter0       = 0;
        fprintf('Looking for mean transfer point \n');
        while (dist0 > 0.01) && (iter0 < 20)
            phi0        = 0.5*phi_min + 0.5*phi_max;
            options.ParentTransVal  = phi0;
            
            options.soboln  = num_model;
            calibrate(x,options);
            
            parentpath              = cd(cd('..'));
            filename_source         = strcat(pwd,'/results/','results_',num2str(options.n_fold,'%i'));
            X                       = load(filename_source);
            mean_transfer_0         = X.res_MOM(1,98);
            
            dif0                    = mean_transfer_0 - mean_transfer_baseline;
            dist0                   = abs(dif0);
            
            if dif0 > 0
                phi_max = phi0;
            else
                phi_min = phi0;
            end
            
            fprintf('Iter = %i: \t Transfer = %3.2f, \t \n Transfer/Income: \t Baseline: %3.3f \t Case: %3.3f \t Difference = %3.4f \n',iter0,phi0,mean_transfer_baseline,mean_transfer_0,dif0);
            iter0  = iter0 + 1;
            if (dist0 > 0.01) && (iter0 < 20)
                fprintf('New Iter \t Min: %3.2f \t Max: %3.2f \n \n',phi_min,phi_max)
            end
            
            % Restart transitory storage
            res_x                   = [];
            res_MOM                 = [];
            parentpath              = cd(cd('..'));
            filename                = strcat(pwd,'/results/','results_',num2str(options.n_fold,'%i'));
            filename_diary          = strcat(pwd,'/results/','diary_',num2str(options.n_fold,'%i'),'.txt');
            save(filename,'res_x','res_MOM');
            diary(filename_diary);
            clear res*
            
        end
        
        options.model_iter_on   = 'Y';          % display model solve iterations
        options.calib_steps     = 'Y';          % display calibration steps
        options.print_moms      = 'Y';
        options.ComputeOtherMus = 'Y';
        options.save_guess	    = 'N';          % save guess
        filename = strcat('guess');
        if  exist(sprintf('%s.mat',filename),'file') == 2
            delete(sprintf('%s.mat',filename));
        end
        %             end
    end
    
    %% SOLVE
    
    diary off; diary on
    options.soboln  = num_model;
    calibrate(x,options);
    
    
    %% Pack-up output results
    parentpath              = cd(cd('..'));
    filename_source         = strcat(pwd,'/results/','results_',num2str(options.n_fold,'%i'));
    filename_destination    = strcat(pwd,'/results/','results_model_options_',options.data_year,'_',num2str(ex_version),'.mat');
    Y           = load(filename_destination);
    X           = load(filename_source);
    res_x       = [ Y.res_x;          reshape(X.res_x,1,numel(X.res_x))];
    res_MOM     = [ Y.res_MOM;        reshape(X.res_MOM,1,numel(X.res_MOM))];
    save(filename_destination,'res_x','res_MOM');
end