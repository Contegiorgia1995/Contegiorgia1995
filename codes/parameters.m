function [ par] = parameters(par_est,options)
% keyboard
% Load parameters of the model. 

% Calibration to 2000 U.S.

% This version June 2015

% Parameters calibrated internally: 
% Altruism : par_est.gamman, par_est.ln_level
% Psy cost : par_est.hscost, par_est.psy_corr

% Other important functions with estimation: 
% par_income_process.m
% ret_rep.m
% ChildCost.m

%% 1. Mean household income (unconditional of educ groups) in 2000 dollars
% * Source PSID, 2000 % 
par.mean_inc        = 70179.45;
mean_inc_age42      = 75630;                        % mean income at age 40-43 (to compute OAS transfers)
% mean_inc_40_educ2   = 60198.72;
% Option to change w for SS exercises
switch options.exerciseSS.On
    case 'N'
        par.w       = 1;
    case {'Y'}
        par.w       = options.exerciseSS.w_ss;
end


%% 2. Age structure
par.time_period     = 2;
par.age             = 0:par.time_period:80;
par.NJ              = length(par.age);
par.Je1             = 16;                                       % HS age
par.Je1_pos         = find(par.age == par.Je1,1,'first');
par.Je2             = 18;                                       % College age
par.Je2_pos         = find(par.age == par.Je2,1,'first');
par.Jc              = 28;                                       % Fertility
par.Jc_pos          = find(par.age == par.Jc,1,'first');
% par.Jt              = 40;                                       % OAS
% par.Jt_pos          = find(par.age == par.Jt,1,'first');
par.Jr              = 66;                                       % Retirement
par.Jr_pos          = find(par.age == par.Jr,1,'first');
par.Jd              = 80;                                       % Die, source: WB - life expectancy = 79
par.Jd_pos          = find(par.age == par.Jd,1,'first');



%% 3. Income process
par.educ     = [1 2 3];
par.dist_educ = [0.094 0.601 (1-.601-.094)];
par.N_markov = 5;
par.N_fe     = 5;
par.sigmah0  = par_est.sigmah0;
par.inc      = par_income_process(par,options);
if 1 == 0
    % To calculate the persistence of the income process: XXX Need to Fix XXX
    solve_corr_model
end

%% Education Returns
par.Reducs      = zeros(2,2); %Education stage, [constant, curvature]: h_(t+1) = h_t^curv * constant
par.Reducs(1,:) = [par_est.R_hs_level par_est.R_hs_curv];
par.Reducs(2,:) = [par_est.R_col_level par_est.R_col_curv];


%% 4. Initial H0 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distribution of parents: By Education 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.beta_inc        = par_est.beta_inc;
par.sigmah0         = par_est.sigmah0;
par.meanh0          = par_est.meanh0;

pp                  = linspace(0.025,0.975,par.N_fe);
switch options.GridH0orig
    case 'N'
        sigma_grid          = 1.5*par_est.sigmah0; % Larger grid of H0 so that it includes more points if beta_inc is large.
    case 'Y'
        sigma_grid          = 1.5*options.SigmaH0orig; % Larger grid of H0 so that it includes more points if beta_inc is large.
end
par.H0              = logninv(pp,log(par.meanh0),sigma_grid);

par.PG              = @(x) [logncdf([par.H0(1:end-1) inf],log(par.meanh0) + par_est.beta_inc * log(x),par_est.sigmah0) ...
                            - logncdf([0 par.H0(1:end-1)],log(par.meanh0) + par_est.beta_inc * log(x),par_est.sigmah0)]';


% Transformed by education
par.H0               = par.H0';
H0                   = par.H0;
par.inc.fe           = cell(length(par.educ),1);
par.inc.fe{1,1}      = log(H0);                                                 %HS Drop
par.inc.fe{2,1}      = log(   H0   + par.Reducs(1,1) * H0 .^par.Reducs(1,2))  ;  %HS Grad
par.inc.fe{3,1}      = log(   H0   + par.Reducs(2,1) * H0 .^par.Reducs(2,2)) ;  %Col Grad

% par.inc.fe{3,1}      = log(par.inc.fe{3,1}  + par.Reducs(2,1) * par.inc.fe{3,1} .^par.Reducs(2,2)) ;    %Col Grad
% Incomes Cutoffs at age 40 ---> For Fertility Groups 
% ie               = 2;
% h0_pos           = ceil(par.N_fe/2);
% h0               = exp(par.inc.fe{ie,1}(h0_pos));
% % Increase from avg. labor market experience until age 40
% j_pos            = find(par.age == 42,1,'first');
% h0               = h0*exp(par.inc.age_prof{ie,1}(j_pos));
% 
% inc_age_42      = zeros(3,1);
% inc_age_42(2)   = par.w * h0;
% inc_age_42(1)   = inc_age_42(2) * 0.6426; 
% inc_age_42(3)   = inc_age_42(2) * 1.7872; 
% par.avg_inc     = inc_age_42;

% prob_educ_data  = [9.4 60.1 30.5]/100;
% inc_age_42(2)   = 1/(prob_educ_data(1)*0.6426 + prob_educ_data(2) + prob_educ_data(3)*1.7872); %Because avg inc at 42 is normalized to 1
% inc_age_42(1)   = inc_age_42(2) * 0.6426; 
% inc_age_42(3)   = inc_age_42(2) * 1.7872; 
% par.avg_inc     = inc_age_42;
% par.avg_inc_42  = prob_educ_data*par.avg_inc;




%% Model-Data exchange rate: use mean income at age 40-43 for HS graduates, from average initial h0
% Incomes Cutoffs at age 40 ---> For Fertility Groups 
% ie               = 2;
% h0_pos           = ceil(par.N_fe/2);
% h0               = exp(par.inc.fe{ie,1}(h0_pos));
% 
% % Increase from avg. labor market experience until age 40
% j_pos            = find(par.age == 40,1,'first');
% h0               = h0*exp(par.inc.age_prof{ie,1}(j_pos));
% 
mean_inc_data    = par.time_period*mean_inc_age42; % 
mean_inc_model   = par.time_period; % Solve to get h0 at age 42 = 2, 1 per year
par.p_model_data =  mean_inc_model / mean_inc_data;


%% 6. Preferences:
% 6.1 Estimated from outside the model
beta_annual         = .975;
switch options.patient
    case 'Y'
        beta_annual         = .99;   
end
par.beta            = beta_annual^par.time_period;                      % Roys, Seshadri - 2014
par.gammac          = 0.5;                                      % Risk aversion <1 to get positive utility, Jones, Schoonbroodt - 2010

% 6.2 Estimated internally, altruism: b(n) = par.ln_level * n^(par.gamman)
par.gamman          = par_est.gamman;                           % gird_gamman  = [0.01 0.8];
par.lambdan         = par_est.ln_level;                         % grid_ln      = [0.01 0.8];


% 6.3 Education costs
% 6.3.1 College: externally
% Source: Delta Cost Project, average cost 2000 = 6588 per year.
% In the model: 4 years of college and household ( two agents)
par.pe2             = par.time_period * 2 * 6588 * par.p_model_data;
% 6.3.2 High School: internally
par.pe1             = par.pe2*0.09;

switch options.exerciseSS.On
    case {'Y'}
        switch options.exerciseSS.change_pe
            case 'Y'
                par.pe2             = par.pe2 * options.exerciseSS.p_mult;
                par.pe1             = par.pe1 * options.exerciseSS.p_mult;
        end
end

pe2_total           = (4/par.time_period)*par.pe2;
pe1_total           = (4/par.time_period)*par.pe1;



% 6.4 Psychic cost: relative cost of education
% par_est.psy         = 5;
% par_est.psy         = 20;
gridpsy             = 100; % multiple of 5
% par.psy_val_hs      = (linspace(0,par_est.psy_max,gridpsy).*par.pe1*par_est.psy_hs)';

switch options.exerciseSS.On
    case 'Y'
        switch options.exerciseSS.change_psy
            case 'Y'
                par_est.psy_hs  = par_est.psy_hs * options.exerciseSS.w_ss^options.psy_adj;
                par_est.psy_max = par_est.psy_max * options.exerciseSS.w_ss^options.psy_adj;
        end     
end


par.psy_val_hs      = (linspace(0,par_est.psy_hs,gridpsy))';
par.psy_val_col     = (linspace(0,par_est.psy_max,gridpsy))';
% par.psy_val         = (linspace(0,par_est.sigmapsy,gridpsy).*par.pe2)';
% ini                 = 100/(gridpsy) *.5;
% psy_cdf             = linspace(ini,100-ini,gridpsy)';
% psy_cdf             = psy_cdf/100;
% par.psy_val         = norminv(psy_cdf,0,par_est.sigmapsy).*par.pe2;

par.psy_prob        = zeros(3,gridpsy);
Ngrid               = gridpsy+1;
% HS drop: more prob in higher values
curv_psy            = par_est.psy_cor;
grid_cdf            = linspace(0^(1/curv_psy),1^(1/curv_psy),Ngrid).^curv_psy;
par.psy_prob(1,:)   = grid_cdf(2:end)-grid_cdf(1:end-1);
% HS grad
curv_psy            = 1;
grid_cdf            = linspace(0^(1/curv_psy),1^(1/curv_psy),Ngrid).^curv_psy;
par.psy_prob(2,:)   = grid_cdf(2:end)-grid_cdf(1:end-1);
% CO grad: more prob in lower values
curv_psy            = par_est.psy_cor;
grid_cdf            = linspace(0^(1/curv_psy),1^(1/curv_psy),Ngrid).^curv_psy;
grid_cdf_inv        = grid_cdf(2:end)-grid_cdf(1:end-1);
par.psy_prob(3,:)   = grid_cdf_inv(gridpsy:-1:1);

prob_educ_data      = [9.4 60.1 30.5]/100;
prob_psy            = prob_educ_data * par.psy_prob;
avg_psy             = [prob_psy*par.psy_val_hs; prob_psy*par.psy_val_col];
switch options.psy_avg
    case 'Y'
        gridpsy         = 2;
        par.psy_val_hs  = [avg_psy(1)*0.995;avg_psy(1)];
        par.psy_val_col = [avg_psy(2)*0.995;avg_psy(2)];
        
        par.psy_prob    = repmat(0.5,3,gridpsy);
end



% 6.5 Child cost parameters
par.Nanny_oppC  = par_est.Nanny_oppC;
par.Child_Ccurv = par_est.Child_Ccurv;
par.Child_C_Tot = par_est.Child_C_Tot;
par.child_cost_inc_curv = par_est.child_cost_inc_curv;

%% 7. Partial Equilibrium prices and taxes
% 7.1 Wages: internally, to match mean income and education shares
% 7.2 Interest rate
% prime        = 1/beta_annual - 1;
prime        = 3/100; % Smets and Wouters AER 2007
par.r        = (1+prime)^par.time_period-1;
par.r_sav    = (1+prime)^par.time_period-1;            % Roys, Seshadri - 2014 XXX CONSIDER ADDING TAXES
par_est.int_iota = 0.1; % XXX Add to estimation? XXX
switch options.noborrowingwedge
    case 'Y'
        par_est.int_iota = 0.01;
end
par.r_debt   =  ((1+prime+par_est.int_iota)^par.time_period - 1);
par.r_col    = ((1+prime+0.00925)^par.time_period - 1); %See /DK/data/Data_US/Education/Loans/College Loans.xls
par.col_fact = par.r_col/(1-(1+par.r_col)^(-10))*(1-(1+par.r_debt)^(-10))/par.r_debt; % Assume debt is repaid in 20 years, i.e. 10 periods.
% switch options.noborrowingwedge
%     case 'Y'
%         par.r_col    = ((1+prime+0.0)^par.time_period - 1); %See /DK/data/Data_US/Education/Loans/College Loans.xls
%         par.col_fact = 1;
% end

% Borrowing contraint: By education
borr     = par.w * [10000 24000 34000].* par.p_model_data; %self-reported limits on unsecured credit by family type from the SCF. Based on Abbot et al (2016).
borr_sch = par.w * [0 0 2*par.col_fact*23000].* par.p_model_data; %Less borrowing if not finished college

% 7.3 Government social security taxes
switch options.ExPolRetirement
    case 'N'
        par.lambda  = .124;                                     % Social Security Ppayroll Tax, Krueger, Ludwig - 2015
    case 'Y'
        par.lambda  = options.PolRetTaxes;
end     
par.pn              = 0.54 * par.time_period * par.mean_inc * par.p_model_data;   % price of childcare, source: Folbre 2008, page 129. Wage of child care in 2000 = $7.43, mean wage = 13.74, 7.43/13.74=.54
par.w_college       = .56 * par.w;                              % wage at college, eg inc = w_college * w * h; source: IPUMS
switch options.exerciseSS.On
    case {'Y'}
%         par.pn				= par.pn  * options.exerciseSS.p_mult ;  
end

gdp_pc              = 44308;
par.gov_debt        = 0.2 * gdp_pc * par.p_model_data;
                   
% Source: AGMV 2013
par.prod_alpha      = 0.350;
par.prod_delta      = 1-(1-0.0650)^par.time_period;
par.prod_s          = [0.160 0.390 0.450];
par.prod_rho        = 0.680;

%% 9. Grids for savings and discrete choices
% 8.1 Number of children
par.fam_size        = 2;                                        % each child of the hh represents par.fam_size children
switch options.Fertility
    case 'Endo' % Endogenous fertility
        par.N       = [0:3];
    case 'Exo'  % Exogenous fertility
        par.N       = options.Fertility_Exo_N/par.fam_size;
end
% N = 0: no child
% N = 1: fam_size children
% N = 2: 2*fam_size children
% etc.


% 8.4 Savings
par.Ls_0            = 20;                                       % Length of savings grid for j = 12
par.Ls_1            = 80;                                      % Length of savings grid for j = 16,20,24
% par.Ls              = [1 1 1 par.Ls_0 repmat(par.Ls_1,1,3) repmat(par.Ls_2,1,11) repmat(par.Ls_3,1,3)]; 
par.Ls              = [ones(1,par.Je1_pos - 1) par.Ls_0 repmat(par.Ls_1,1,par.Jd_pos - par.Je1_pos)]; 

curv                = 3;
switch options.ParentTrans
    case 'Endo'
        PHI_max             = (pe1_total +  pe2_total)*3;
        par.PHI             = [0 pe1_total/3 2/3*pe1_total linspace((pe1_total*1.1)^(1/curv),PHI_max^(1/curv),par.Ls(par.Je1_pos)-3).^curv];
    case 'Exo'
        % Computation of factor: 
        % 1. solve with factor = 1. 
        % 2. Solve for factor such that (res_MOM(4,100)*factor-res_MOM(1,100))/res_MOM(1,100) = 1
        switch options.ParentTransExo
            case 'MeanTrans'
                options.ParentTransVal = 2*30566*par.p_model_data;
            
            case 'MeanHSCost'
                options.ParentTransVal = pe1_total*1.5;
                
            case 'MeanColCost'
                options.ParentTransVal = (pe1_total + pe2_total)*1.5;
            case 'search_val'
                options.ParentTransVal = 2*options.ParentTransVal*par.p_model_data;
        end
        par.PHI         = options.ParentTransVal * [0.999 1.001];
end

par.grids              = cell(length(par.educ),par.Jd_pos);

for ie = 1:length(par.educ)
    borr_limits    = [zeros(par.Je2_pos,1); borr_sch(ie); borr_sch(ie); borr(ie)*ones(par.Jr_pos-1-(par.Je2_pos+2),1); zeros(par.Jd_pos - (par.Jr_pos-1),1)];
    % j = Je1
    par.grids{ie,par.Je1_pos} = par.PHI;

    % j = Je2
    maxE = par.PHI(end)*2;
    for j       = par.Je1_pos + 1: par.Je2_pos + 1
        par.grids{ie,j}     = linspace(-borr_limits(j)^(1/curv),maxE^(1/curv),par.Ls(j)).^curv;
    end

    % j = Je2+1 : Jc -1 CHECK! IT'S IMPORTANT THAT GRIDS STARTS CHANGING BY EDUCATION ONLY AFTER JE2+1
    
    for j=par.Je2_pos+2:par.Jc_pos-1
        maxWY           = par.PHI(end)*par.N(end)*0.25;
        max_inc         = par.w*exp(par.inc.age_prof{ie,1}(j-1) + max(par.inc.fe{ie,1}(:)) + max(par.inc.z_val{ie,1}(j-1,:)))*2;
        maxS            = max(par.grids{ie,j-1}(end)*1.3,max(max_inc,maxWY));
        par.grids{ie,j} = linspace(-borr_limits(j)^(1/curv),maxS^(1/curv),par.Ls(j)).^curv;
    end

    % j = Jc: Jr
    for j = par.Jc_pos:1:par.Jr_pos
        max_inc         = par.w*exp(par.inc.age_prof{ie,1}(j-1) + max(par.inc.fe{ie,1}(:)) + max(par.inc.z_val{ie,1}(j-1,:)));
        maxS            = max(par.grids{ie,j-1}(end)*1.3,max_inc*3);
        par.grids{ie,j} = linspace(-borr_limits(j)^(1/curv),maxS^(1/curv),par.Ls(j)).^curv;
    end

    % j = Jr+1 : Jd
    maxR                = par.grids{ie,j}(end)*2;
    for j=par.Jr_pos+1:par.Jd_pos
        par.grids{ie,j}    = linspace(-borr_limits(j)^(1/curv),maxR^(1/curv),par.Ls(j)).^curv;
    end
end

%% Fertility Policies
switch options.fert_transfer
    case 'N'
        par.fert_trans = zeros(length(par.N),1);
    case 'only2'
        par.fert_trans = zeros(length(par.N),1);
        ind            = logical(par.N == 1);
        par.fert_trans(ind) = 2*options.fert_trans_size*par.p_model_data;
end

switch options.init_transfer
    case 'N'
        par.init_trans                  = zeros(length(par.educ),par.Jd_pos);
    case 'Y'
        par.init_trans                  = zeros(length(par.educ),par.Jd_pos);
        par.init_trans(:,par.Je1_pos)   = 2*options.init_trans_size*par.p_model_data;
    case 'OnlyCollege'
        par.init_trans                                  = zeros(length(par.educ),par.Jd_pos);
        par.init_trans(3,par.Je2_pos:par.Je2_pos+1)     = 2*options.init_trans_size*par.p_model_data;
end

end