clear; close all; clc;

%% Solve for h0
% Parents
par.h0_rho          = 0.5;

n_hp                = 5;
grid_hp             = linspace(0.1,10,n_hp);
prob_hp             = repmat(1/n_hp,1,n_hp);

% Moments from parents
par.hp_log_mean     = log(grid_hp) * prob_hp';
par.mean_exp_rho_hp = exp(par.h0_rho * (log(grid_hp) - par.hp_log_mean)) * prob_hp';

% Children
% nsigma              = 5;
% grid_sigma          = linspace(1,5,nsigma);
par.h0_n            = 5;
nsigma              = 1;
grid_sigma          = 1;

for is = 1:nsigma
    par.h0_sigma        = grid_sigma(is);
    fprintf('sigma = %3.2f \t',par.h0_sigma)
    
    % 1. Get standard normal
    [x,w]           = GaussHermite(par.h0_n);
    prob            = w/sum(w);
    par.h0_prob     = prob;
    
    % 2. Get Normal (mu,sigma)
    exp_x_mean      = exp(par.h0_sigma*2^.5*x') * prob;
    
    % Mu such that E(h0) = 1
    mu                  = log(1/(exp_x_mean*par.mean_exp_rho_hp));
    par.h0_shock_grid   = mu + 2^.5 * par.h0_sigma * x;
    
    mean_exp_shock      = exp(par.h0_shock_grid') * par.h0_prob;
    fprintf('mean exp shock = %3.2f \t',mean_exp_shock)
    
    
    
    %% Simulate h0 draws
    grid_h0             = nan(n_hp,par.h0_n);
    grid_h0log          = nan(n_hp,par.h0_n);
    grid_h0_prob        = nan(n_hp,par.h0_n);
    grid_h0_cond_prob	= nan(n_hp,par.h0_n);
    mean_h0_cond        = nan(n_hp,1);
    
    for ip = 1:n_hp
        grid_h0log(ip,:)	= par.h0_shock_grid + par.h0_rho * (log(grid_hp(ip))- par.hp_log_mean);
        grid_h0(ip,:)       = exp(grid_h0log(ip,:));
        grid_h0_prob(ip,:)  = par.h0_prob .* prob_hp(ip);
        grid_h0_cond_prob(ip,:)  = par.h0_prob ;
        mean_h0_cond(ip)    = grid_h0(ip,:) * grid_h0_cond_prob(ip,:)';
    end
    
    
    % Compute mean h0
    mean_h0             = grid_h0(:)' * grid_h0_prob(:);
    fprintf('mean h0 = %3.2f \t',mean_h0)
    
    if is == ceil(nsigma/6)
        figure(1)
        subplot(1,2,1)
        hold on
        leg  = {};
        for ip  = 1:n_hp
            h = plot(grid_h0(ip,:),cumsum(par.h0_prob));
            h.Marker = 'x';
            leg = [leg;{['hp = ',num2str(grid_hp(ip),'%3.1f')]}];
        end
        xlabel('h0')
        ylabel('CDF')
        legend(leg)
    end
    
    
    %% Move to same grid of h0
    h0_min                  = min(grid_h0(:));
    h0_max                  = max(grid_h0(:));
    % Choose points to minimze the distance between desired points and grid
    % points
    curv                    = 1;
    grid_h0_all             = linspace(h0_min^curv,h0_max^curv,par.h0_n).^(1/curv);
    desired_points          = grid_h0(:);
    desired_probs           = grid_h0_cond_prob(:);
    par.n_hp                = n_hp;
    par.grid_hp             = grid_hp;
    options.print           = 'N';
    obj                     = @(gridh) par_func_h0_make_grid(gridh,grid_h0,grid_h0_cond_prob,mean_h0_cond,prob_hp,par,options);
    
    opt                     = optimset('fminsearch');
    opt.Display             = 'notify';
    [x1,eq,exitflag]        = fminsearch(obj,grid_h0_all,opt);
    
    if exitflag <= 0
        fprintf('Error fminsearch grid h0')
    end
    par.grid_h0_all         = sort(x1);
    
    options.print           = 'Y';
    [obj,paraux]            = par_func_h0_make_grid(par.grid_h0_all,grid_h0,grid_h0_cond_prob,mean_h0_cond,prob_hp,par,options);
    par.grid_prob_h0        = paraux.grid_prob_h0;
    
    
    
    if is == ceil(nsigma/6)
        figure(1)
        subplot(1,2,2)
        hold on
        leg  = {};
        for ip  = 1:n_hp
            h = plot(par.grid_h0_all,cumsum(par.grid_prob_h0(ip,:)));
            h.Marker = 'x';
            %            leg = [leg;{['hp = ',num2str(grid_hp(ip),'%3.1f')]}];
        end
        xlabel('h0')
        ylabel('CDF')
        %        legend(leg)
    end
    
    
end



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
    
    e_expy     = expy' * prob;
end


