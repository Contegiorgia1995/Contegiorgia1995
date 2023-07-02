%% Compute a simple idea of initial wages
clear; close all; clc;
% Estimation: (does not matter, only to get parameters)
options.Fertility       = 'Endo';       % 'Endo': endogenous fertility; 'Exo': constant fertility
options.ParentTrans     = 'Endo';       % Endo: Benchmark case with parents transfers; 'Exo': exogenously given
options.equilibrium     = 'partial';    % 'partial' or 'general'
options.transfers       = 'on';


par_est.gamman      = 0.200;
par_est.ln_level    = 0.200;
par_est.hscost      = 0.05;
par_est.psyc_mean   = 1;
par_est.psyc_var    = 15;
par_est.w           = [1 1 1];

[par]               = parameters_new(par_est,options);    % Load parameters


EDUC     = par.educ;
Agg_H    = zeros(3,1);
sum_mu_L = zeros(3,1);
start_L  = [par.Je1_pos par.Je2_pos par.Je2_pos + 1];

mu_educ  = [11.9 59.0 29.1]/100;
N_age    = 1/length(par.age);
prob_tot = 0;
for educ = 1:length(EDUC)
    for j_pos = start_L(educ):par.Jr_pos-1
        AGE_PROF  = par.inc.age_prof(educ,j_pos);
        h_prob     = mu_educ(educ)*N_age;
        h_val      = exp(AGE_PROF);
        Agg_H(educ)    = Agg_H(educ) + h_prob*h_val;
        sum_mu_L(educ) = sum_mu_L(educ) + h_prob;
        prob_tot   = prob_tot + h_prob;
    end
end

Y_tot = 1000;
K_tot = 3*Y_tot/par.time_period;

prod_alpha      = 0.350;
prod_delta      = 1-(1-0.0650)^par.time_period;
prod_s          = [0.160 0.390 0.450];
prod_rho        = 0.680;



H_tot               = (prod_s * Agg_H.^prod_rho)^1/prod_rho;
TFP                 = Y_tot/(K_tot^prod_alpha * H_tot^(1-prod_alpha));


w                   = TFP * prod_s' .* (1-prod_alpha) .* (K_tot/H_tot)^prod_alpha .* (H_tot./Agg_H).^(1-prod_rho);
r                   = TFP * prod_alpha*(H_tot/K_tot)^(1-prod_alpha)-prod_delta;
wage_premium_hs_drop =  w(1)/w(2);
wage_premium_college = w(3)/w(2);
