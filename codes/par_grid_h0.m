clear; close all; clc;

%% Get grid h0
% Parents
par.h0_rho          = 0.5;
n_hp                = 10;
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
    
    % Compute mean shock (for reference)
    mean_exp_shock      = exp(par.h0_shock_grid') * par.h0_prob;
    fprintf('mean exp shock = %3.2f \t',mean_exp_shock)
    
    
    
    %% Simulate h0 draws correlated with parents
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
    
    
    %% Move to same grid of h0
    h0_min                  = min(grid_h0(:));
    h0_max                  = max(grid_h0(:));
    % Choose points to minimze the distance between desired points and grid
    % points
    curv                    = 1;
    gridh0_0                = linspace(h0_min^curv,h0_max^curv,par.h0_n).^(1/curv);
    
    % Search grid h0 to minimze distance - not always the better option...
    if 1 == 0
        desired_points          = grid_h0(:);
        desired_probs           = grid_h0_cond_prob(:);
        
        obj                     = @(gridh) par_func_h0_make_grid(gridh,desired_points,desired_probs);
        
        opt                     = optimset('fminsearch');
        opt.Display             = 'notify';
        [x1,eq,exitflag]        = fminsearch(obj,gridh0_0,opt);
        
        if exitflag <= 0
            fprintf('Error fminsearch grid h0')
        end
        par.grid_h0_all         = sort(x1);
    else
        par.grid_h0_all         = gridh0_0;
    end
    
    % Compute probabilities in new grid
    grid_prob_h0            = zeros(n_hp,par.h0_n);
    total_prob_h0           = zeros(n_hp,par.h0_n);
    mean_h0_cond_approx     = nan(n_hp,1);
    
    for ip = 1:n_hp
        for ii = 1:par.h0_n
            % find left and right points
            if grid_h0(ip,ii) >= par.grid_h0_all(end)
                pos_r           = size(par.grid_h0_all,2);
                pos_l           = pos_r - 1;
                % compute distance to left and right
                weight_l        = 0;
                weight_r        = 1;
            elseif grid_h0(ip,ii) <= par.grid_h0_all(1)
                pos_r = 1;
                pos_l = 2;
                % compute distance to left and right
                weight_l        = 1;
                weight_r        = 0;
            else
                pos_r           = find(par.grid_h0_all>grid_h0(ip,ii),1,'First');
                pos_l           = pos_r - 1;
                % compute distance to left and right
                %                 keyboard
                dist_l          = grid_h0(ip,ii) - par.grid_h0_all(pos_l);
                dist_r          = par.grid_h0_all(pos_r) - grid_h0(ip,ii);
                weight_l        = dist_r/(dist_l + dist_r);
                weight_r        = dist_l/(dist_l + dist_r);
            end
            % Assign probabilities
            grid_prob_h0(ip,pos_l) = grid_prob_h0(ip,pos_l) + weight_l * grid_h0_cond_prob(ip,ii);
            grid_prob_h0(ip,pos_r) = grid_prob_h0(ip,pos_r) + weight_r * grid_h0_cond_prob(ip,ii);
        end
        mean_h0_cond_approx(ip)     = par.grid_h0_all * grid_prob_h0(ip,:)';
        total_prob_h0(ip,:)             =  grid_prob_h0(ip,:) .* prob_hp(ip);
    end
    
    
    aux_grid_h0         = repmat(par.grid_h0_all,n_hp,1);
    
    mean_h0_approx1     = aux_grid_h0(:)' * total_prob_h0(:);
    
    % Re-center grid to have same mean
    mean_distance       = mean_h0 - mean_h0_approx1;
    
    if abs(mean_distance) > 0
        fprintf('mean adjustment %3.2f\t',mean_distance)
    end
    
    %     keyboard
    if min(mean_distance + par.grid_h0_all) < 0
        fprintf('Need to truncate grid at 0 by %3.2f \t',min(mean_distance + par.grid_h0_all))
    end
    
    par.grid_h0_all     = max(0.001,mean_distance + par.grid_h0_all);
    aux_grid_h0         = repmat(par.grid_h0_all,n_hp,1);
    mean_h0_approx2     = aux_grid_h0(:)' * total_prob_h0(:);
    
    
    fprintf('mean h0 approx = %3.2f \n \n',mean_h0_approx2)
    for ip = 1:n_hp
        fprintf('hp = %3.2f \t distance = %3.2f \t mean h0 = %3.2f mean h0 approx = %3.2f\n',...
            grid_hp(ip),mean_h0_cond(ip)-mean_h0_cond_approx(ip),...
            mean_h0_cond(ip),mean_h0_cond_approx(ip));
        
    end
    
    %     fprintf('Compute mean for each ip, choose 5 points to minimize distance? \n\n')
    
    if is == ceil(nsigma/6)
        figure(1)
        nsp = ceil(n_hp^.5);
        for ip  = 1:n_hp
            subplot(nsp,nsp,ip)
            hold on
            %             h = plot(grid_h0(ip,:),cumsum(par.h0_prob));
            h = plot(grid_h0(ip,:),par.h0_prob);
            h.Marker = 'x';
            %             h = plot(par.grid_h0_all,cumsum(grid_prob_h0(ip,:)));
            h = plot(par.grid_h0_all,grid_prob_h0(ip,:));
            h.Marker = 'o';
            title(['hp = ',num2str(grid_hp(ip),'%3.1f')])
            if ip == 1
                legend('true','approx')
            end
            xlabel('h0')
            ylabel('pdf')
        end
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


