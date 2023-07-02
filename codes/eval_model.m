function out = eval_model(x,xlb,xub,options)
%% Check x within bounds
x       = min(max(x,xlb),xub);
fprintf('gamman = %3.4f, ln_level = %3.4f, mu = 1.00, sigmah0 = %3.4f, psy max = %3.4f, psy cor = %1.3f, psy hs = %1.3f, R hs curv = %3.4f, R col curv = %3.4f, R hs level = %3.4f, R col level = %3.4f,  Child_C_Tot= %1.3f,  child_cost_inc_curv= %1.3f,  BetaInc= %1.3f  \n',x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12),x(13));

%% Solve model and compute moments  for estimated parameters x
rng(666);                                             % Set seed for shocks
ticcalib = tic;
% Load parameters
% Estimation:
par_est.gamman      = x(1);
par_est.ln_level    = x(2);
par_est.mu          = 0.50;
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
par_est.beta_inc    = x(13);

switch options.solve_h0
    case 'Y'
        h0_min              = 0.05;
        h0_max              = 0.7;
        dist_h0             = 100;
        pace_h0             = 0.25;
        iter_h0             = 1;
        tol_h0              = 0.01;
        max_iter_h0         = 50;
        stop                = 0;
        
        
        while dist_h0 > tol_h0 && iter_h0 < max_iter_h0 && stop == 0
            par_est.meanh0  = 0.5*h0_min + 0.5*h0_max;
            % Parameters outside estimation:
            [par]                                   = parameters(par_est,options);    % Load parameters
            par.final_iter                          = 'N';
            %% Solve policy functions:
            switch options.calib_steps
                case {'Y'}
                    fprintf('solve model, time: %3.2f sec\n',toc(ticcalib));
            end
            
            [par,pol,guess,flags_solve]             = solve_model(par,options); %#ok<*ASGLU,*NASGU>
            switch flags_solve.converge
                case 'Y'
                    clear guess
            
            
            
                    % 1. Interpolate policy functions into dense grid
                    [pol_dense,par,par2]                    = solve_ergodic_distribution1(par,pol,options);
                    clear pol
                    % 2. Find ergodic distribution across cohorts
                    [mu,mu_ige0,muc,pop_final,Q_ergo]       = solve_ergodic_distribution2(par2,pol_dense,options);

                    % 3. Find distribution for a given cohort
                    %     [mu,mu_cs,dist_age]                     = solve_ergodic_distribution3(par2,pol_dense,mu,pop_final,Q_ergo,options);

                    me_inc_42      = mean_income_age42( par2,mu );
                    dif_h0         = me_inc_42 - 2;
                    dist_h0        = abs(dif_h0);
                    if dist_h0 > tol_h0
                        if dif_h0 > 0
                            h0_max = pace_h0*h0_max+(1-pace_h0)*par_est.meanh0;
                        else
                            h0_min = pace_h0*h0_min+(1-pace_h0)*par_est.meanh0;
                        end
                        fprintf('iter = %i: \t h0 = %3.3f dist = %3.3f \t new h0 = [%3.4f,%3.4f] \n ',iter_h0,par_est.meanh0,dist_h0,h0_min,h0_max);
                    else
                        fprintf('CONVERGED H0 \n iter = %i: \t h0 = %3.3f dist = %3.3f \t \n ',iter_h0,par_est.meanh0,dist_h0);
                    end
                case 'N'
                    h0_max  = h0_max + 0.05*rand;
                    h0_min  = h0_min - 0.05*rand;
                    fprintf('iter = %i: \t Solve Failed. h0 = %3.3f dist = %3.3f \t new h0 = [%3.4f,%3.4f] \n ',iter_h0,par_est.meanh0,dist_h0,h0_min,h0_max);
            end
            iter_h0 = iter_h0 + 1;
        end
        x(14)                =  par_est.meanh0;
        par2.final_iter      = 'Y';
        
    case 'N'
        par_est.meanh0      = x(14);
        [par]                                   = parameters(par_est,options);    % Load parameters
        par.final_iter                          = 'Y';
        %% Solve policy functions:
        switch options.calib_steps
            case {'Y'}
                fprintf('solve model, time: %3.2f sec\n',toc(ticcalib));
        end
        
        [par,pol,guess,flags_solve]             = solve_model(par,options); %#ok<*ASGLU,*NASGU>
        switch flags_solve.converge
            case 'Y'
                switch options.save_mat
                    case 'Y'
                        filename = strcat('guess_',num2str(options.soboln,'%i'));
                        save(sprintf('%s.mat',filename),'guess');
                end
                
                switch options.save_guess
                    case {'Y'}
                        filename = strcat('guess');
                        save(sprintf('%s.mat',filename),'guess');
                        parentpath    = cd(cd('..'));
                        filename      = strcat(parentpath,'/RESULTS/','guess_',num2str(options.n_fold,'%i'));
                        save(sprintf('%s.mat',filename),'guess');
                        if options.soboln >=2
                            options.save_guess = 'N';
                        end
                end
                
                clear guess
        end
        
        % 1. Interpolate policy functions into dense grid
        [pol_dense,par,par2]                    = solve_ergodic_distribution1(par,pol,options);
        clear pol
end

switch flags_solve.converge
    case 'N'
        fprintf('ERROR HERE: solve_model did not converge \n')
        out         = 1000;
        error_type  = 2;
        
        parentpath  = cd(cd('..'));
        filename    = strcat(parentpath,'/RESULTS/','errors_',num2str(options.n_fold,'%i'));
        Y           = load(sprintf('%s.mat',filename));
        aux         = [reshape(x,1,numel(x)) out options.n_fold error_type];
        error_x     = [ Y.error_x;        reshape(aux,1,numel(aux))];
        save(filename,'error_x');
        
    case 'Y'
        %% Ergodic distribution
        switch options.calib_steps
                case {'Y'}
                fprintf('\n Compute ergodic distribution, time: %3.2f sec\n \n',toc(ticcalib));
        end
        [mu,mu_ige0,muc,pop_final,Q_ergo]       = solve_ergodic_distribution2(par2,pol_dense,options);

       switch options.save_mat
            case {'Y'}
                pol_sel.Np      = pol_dense.Np;
                pol_sel.mu_fert = mu{par2.Jc_pos,1};
                filename = strcat('PolicyFunctionsFertility_',num2str(options.soboln,'%i'));
                save(sprintf('%s.mat',filename),'par2','pol_sel','options');
        end
        flag.pop_final       = pop_final;
        if flag.pop_final == 0 % stop
            fprintf('Population = 0, go to next set of parameters\n')
            [m_data,~,~] = moments_data;
            Nm          = size(m_data,1);
            m_model     = NaN(Nm,1);
        else 
            % 3. Find distribution for a given cohort
            [mu,mu_cs,dist_age]                     = solve_ergodic_distribution3(par2,pol_dense,mu,pop_final,Q_ergo,options);
            % 4. Intergenerational distributions
            options.Nsample                 = 300000;
            [IGE]                           = gen_ige_sample(par2,pol_dense,muc,options);
            
            [mu_ige,mu_chetty,Inc_par_grid,mu_educ_ige,mu_par_ic] = trans_ige(par2,mu_ige0,pol_dense.tau0,pol_dense.S,pol_dense.Se,options);
            
            switch options.save_mat
                case 'Y'
                    filename = strcat('dist_LE');
                    save(sprintf('%s.mat',filename),'par2','mu');
            end

            switch options.ComputeMoments
                case {'Y'}
                    %% Compute moments:
                    switch options.calib_steps
                        case {'Y'}
                            fprintf('compute moments, time: %3.2f sec\n',toc(ticcalib));
                    end
                     [m_model,m_title,m_mod_id]      = moments_model(par2,pol_dense,mu,mu_cs,mu_ige,mu_chetty,...
                        Inc_par_grid,mu_educ_ige,dist_age,IGE,options);            
            end
        end

        %% Pack-up output: use same as read_res to evaluate estimation
        %1. Moments from data
        switch options.data_year
            case 'data_2000'
                [m_data,mtitles,inc]= moments_data;
            case 'data_1960'
                [m_data,mtitles,inc]= moments_data_1960;
            case 'data_fert_1960'
                [m_data,mtitles,inc]= moments_data_fert_1960;
        end
        
%         [m_data,mtitles,inc]= moments_data;
%         m_data(7)           = -0.126;
        Nm                  = length(m_data); 
        inc                 = zeros(size(inc));

        moments             = [3 5 98 14 886 840 863 978 1001 121 114:118 ];
        inc(moments)        = 1;
		inc([114:118])      = 25;
        m_data              = m_data';
        Np                  = 13;

        %2. Moments from model
        results_model       = m_model;
        results_data        = m_data;

        % 3. Make weighting matrix with weights only on moments in inc
        d       = (1:Nm)';
        d       = d(inc>0);
        m_est   = results_data';
        W       = diag(inc.*(1./(max(0.1,abs(m_est)).^2)));         % Weight by 1/moment

        % 4. Compute objective function values for each of the model runs
        inc_not         = (inc==0);
        mhat            = results_model;
        mhat(inc_not)   = 0;
        m_est(inc_not)  = 0;
        out             = (m_est-mhat)'*W*(m_est-mhat);
        
        if isnan(out) == 1
            out     = 10000;
        end

        %% Pack-up output: Results:
        parentpath  = cd(cd('..'));
        filename    = strcat(parentpath,'/RESULTS/','results_',num2str(options.n_fold,'%i'));
        Y           = load(sprintf('%s.mat',filename));
        aux         = [reshape(x,1,numel(x)) out options.n_fold par2.p_model_data flags_solve.iter];
        res_x       = [ Y.res_x;          reshape(aux,1,numel(aux))];
        res_MOM     = [ Y.res_MOM;        reshape(m_model,1,numel(m_model))];
        save(filename,'res_x','res_MOM');
end

fprintf('Objective Function = %3.5f, time: %3.2f sec\n',out,toc(ticcalib));


end