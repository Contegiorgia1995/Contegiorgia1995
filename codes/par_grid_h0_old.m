clear; close all; clc;

%% Solve Approximation for log-normal with mean = 1
if 1 == 0
% 1. Get standard normal with gauss hermite quadrature
n           = 5;
[x,w]       = GaussHermite(n);
prob        = w/sum(w);

% 2. Get normal
sigma       = 5;
exp_x_mean  = exp(sigma*2^.5*x') * prob;
mu          = log(1/exp_x_mean);
y           = mu + 2^.5 * sigma * x;

% 3. Get y
expy       = exp(y);

e_expy     = expy' * prob
end

%% Solve for h0
% Parents
n_hp                = 10;
grid_hp             = linspace(0.1,10,n_hp);
prob_hp             = repmat(1/n_hp,1,n_hp);
par.hp_mean         = log(grid_hp) * prob_hp';
par.h0_rho          = .5;

par.hp_rho_exp_mean = exp(par.h0_rho * (log(grid_hp)-par.hp_mean)) * prob_hp';

% Children
nsigma              = 11;
grid_sigma          = linspace(0,5,nsigma);
par.h0_n            = 5;

for is = 1:nsigma
    par.h0_sigma        = grid_sigma(is);
    fprintf('sigma = %3.2f \t',par.h0_sigma)
    
    
    
%     [zgrid, P]          = discretize_ar_rouwen(0, par.h0_sigma, par.h0_n);
    [zgrid,P]           = tauchen(0,0,1,par.h0_n);
    par.h0_prob         = P(1,:)';
    
    % Normal
    h0_mu               = log(1/par.hp_rho_exp_mean) - par.h0_sigma^2/2;
    par.h0_shock_grid	= h0_mu + par.h0_sigma * zgrid';
    
    % Log-normal
%     h0_mu               = - par.h0_sigma^2/2;
%     par.h0_shock_grid	= exp(h0_mu + par.h0_sigma * zgrid');

    mean_shock          = exp(par.h0_shock_grid) * par.h0_prob * par.hp_rho_exp_mean;
    fprintf('E[exp(x)] * E[exp(rho*y)] = %3.2f \t',mean_shock)

    
%     h0_mu               = - par.h0_sigma^2/2;
%     v                   = exp(par.h0_sigma^2) - 1;
%     h0_mu               = log(1/(1+v)^.5);
%     par.h0_shock_grid	= exp(h0_mu + par.h0_sigma*zgrid');
    
    
    
    %% Simulate h0 draws
    grid_h0             = nan(n_hp,par.h0_n);
    grid_h0log          = nan(n_hp,par.h0_n);
    grid_h0_prob        = nan(n_hp,par.h0_n);
    
    for ip = 1:n_hp
        grid_h0log(ip,:)	= par.h0_shock_grid + par.h0_rho * (par.hp_mean - log(grid_hp(ip)));
        grid_h0(ip,:)       = exp(grid_h0log(ip,:));
        grid_h0_prob(ip,:) = par.h0_prob .* prob_hp(ip);
    end
    
    
    % Compute mean h0
    mean_h0log          = grid_h0log(:)' * grid_h0_prob(:);
    fprintf('mean log(h0) = %3.2f \t',mean_h0log)
    mean_h0             = grid_h0(:)' * grid_h0_prob(:);
    fprintf('mean h0 = %3.2f \n',mean_h0)
    
    
end