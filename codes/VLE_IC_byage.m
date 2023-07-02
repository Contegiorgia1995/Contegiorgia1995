function VLE_IC_byage(x,options)

switch options.continue
    case 'N'
        %% Solve model and compute moments  for estimated parameters x
        rng(666);                                             % Set seed for shocks
        ticcalib = tic;
        % Load parameters
        % Estimation:
        par_est.gamman      = x(1);
        par_est.ln_level    = x(2);
        par_est.mu          = 0.5;
        par_est.sigmah0     = x(3);
        par_est.psy_max     = x(4);
        par_est.psy_cor     = x(5);
        par_est.psy_hs      = x(6);
        par_est.R_hs_curv   = x(7);
        par_est.R_col_curv  = x(8);
        par_est.R_hs_level  = x(9);
        par_est.R_col_level = x(10);
        par_est.Nanny_oppC  = 1;
        par_est.Child_Ccurv = 0.6445;
        % par_est.Child_Ccurv = 1.000;
        par_est.Child_C_Tot = x(11);
        par_est.child_cost_inc_curv = x(12);
        % par_est.int_iota    = x(12);
        par_est.beta_inc    = x(13);
        par_est.meanh0      = x(14);
        
        
        % Parameters outside estimation:
        [par]               = parameters(par_est,options);    % Load parameters
        par.final_iter      = 'Y';
        
        %% Solve policy functions:
        switch options.calib_steps
            case {'Y'}
                fprintf('solve model, time: %3.2f sec\n',toc(ticcalib));
        end
        [par,pol,guess,flags_solve] = solve_model(par,options); %#ok<*ASGLU,*NASGU>
        clear guess
        
        switch flags_solve.converge
            case 'N'
                fprintf('ERROR HERE: solve_model did not converge \n')
                
            case 'Y'
                %% Ergodic distribution
                switch options.calib_steps
                    case {'Y'}
                        fprintf('\n Compute ergodic distribution, time: %3.2f sec\n \n',toc(ticcalib));
                end
                
                % 1. Interpolate policy functions into dense grid
                [pol_dense,par,par2]            = solve_ergodic_distribution1(par,pol,options);
                clear pol
                % 2. Find ergodic distribution across cohorts
                [mu,mu_ige0,muc,pop_final,Q_ergo]      = solve_ergodic_distribution2(par2,pol_dense,options);
                flag.pop_final       = pop_final;
                
                if flag.pop_final == 0 % stop
                    fprintf('Population = 0, go to next set of parameters\n')
                    [m_data,~,~] = moments_data;
                    Nm          = size(m_data,1);
                    m_model     = NaN(Nm,1);
                else % Compute next distributions
                    % 3. Find distribution for a given cohort
                    [mu,mu_cs,dist_age]             = solve_ergodic_distribution3(par2,pol_dense,mu,pop_final,Q_ergo,options);
                    
                    options.Nsample                 = 300000;
                    [IGE]                           = gen_ige_sample(par2,pol_dense,muc,options);
                    
                    % Empty array
                    count_m   = 1;
                    m_model   = NaN(1000,1);
                    m_mod_id  = NaN(1000,1);
                    m_title   = cell(1000,1);
                    VLE_IC_Age= NaN(par2.Jr_pos-1,1);
                    
                    % Age 16
                    j_pos  = par2.Je1_pos;
                    [~,aux,~,~ ] = mom_LE( par2,count_m,m_model,m_title,m_mod_id,IGE,options );
                    VLE_IC_Age(j_pos,1) = aux(5);
                    fprintf('Age %i, VLE-IC = %3.6f \n',par2.age(j_pos),aux(5));
                    clear mu_LE_age Grid_Avg_LE_age
                end
                
                j_pos_start = par2.Je2_pos;
                
                % Save to be able to continue from here
                parentpath  = cd(cd('..'));
                filename    = strcat(parentpath,'/RESULTS/','VLEIC_continue_',num2str(options.n_fold,'%i'),'.mat');
                save(filename);
        end
        
    case 'Y'
        % Load to be able to continue from here
        parentpath  = cd(cd('..'));
        filename    = strcat(parentpath,'/RESULTS/','VLEIC_continue_',num2str(options.n_fold,'%i'),'.mat');
        load(filename);
end

% keyboard


% Age 18+
for j_pos = j_pos_start:par2.Jr_pos-1
    [ ~,aux,~,~ ] = mom_LE_age( par2,count_m,m_model,m_title,m_mod_id,IGE,j_pos,options );
    VLE_IC_Age(j_pos,1) = aux(2);
    fprintf('Age %i, VLE-IC = %3.6f \n',par2.age(j_pos),aux(2));
    clear mu_LE_age Grid_Avg_LE_age
    
    j_pos_start = j_pos + 1;
    
    % Save to be able to continue from here
%     parentpath  = cd(cd('..'));
%     filename    = strcat(parentpath,'/RESULTS/','VLEIC_continue_',num2str(options.n_fold,'%i'),'.mat');
%     save(filename);
%     keyboard
end

%% Pack-up output: Results:
par         = par2;
parentpath  = cd(cd('..'));
filename    = strcat(parentpath,'/RESULTS/','LE_dist_age_',num2str(options.n_fold,'%i'));
save(filename,'x','par','VLE_IC_Age');

fprintf('Calibration time: %3.2f sec\n',toc(ticcalib));
end
