function [eq,par] = par_func_h0_make_grid(grid_h0_all,grid_h0,grid_h0_cond_prob,mean_h0_cond,prob_hp,par,options)
% keyboard
par.grid_h0_all  = sort(grid_h0_all);

%%

grid_prob_h0            = zeros(par.n_hp,par.h0_n);
total_prob_h0           = zeros(par.n_hp,par.h0_n);
mean_h0_cond_approx     = nan(par.n_hp,1);

for ip = 1:par.n_hp
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
            dist_l          = grid_h0(ip,ii) - par.grid_h0_all(pos_l);
            dist_r          = par.grid_h0_all(pos_r) + grid_h0(ip,ii);
            weight_l        = dist_r/(dist_l + dist_r);
            weight_r        = dist_l/(dist_l + dist_r);
        end
        % Assign probabilities
        grid_prob_h0(ip,pos_l) = grid_prob_h0(ip,pos_l) + weight_l * grid_h0_cond_prob(ip,ii);
        grid_prob_h0(ip,pos_r) = grid_prob_h0(ip,pos_r) + weight_r * grid_h0_cond_prob(ip,ii);
    end
    mean_h0_cond_approx(ip)         = par.grid_h0_all * grid_prob_h0(ip,:)';
    total_prob_h0(ip,:)             =  grid_prob_h0(ip,:) .* prob_hp(ip);
end


%% EQ aux

n_data  = length(mean_h0_cond_approx);
eqaux   = nan(1,n_data);
for ii = 1:n_data
    eqaux(ii)	= (mean_h0_cond_approx(ii)-mean_h0_cond(ii))^2;
end

% keyboard

eq = eqaux * prob_hp';

%%
par.grid_prob_h0    = grid_prob_h0;

switch options.print
    case 'Y'
        aux_grid_h0         = repmat(par.grid_h0_all,par.n_hp,1);
        
        mean_h0             = aux_grid_h0(:)' * total_prob_h0(:);
        fprintf('mean h0 approx = %3.2f \n',mean_h0)
        for ip = 1:par.n_hp
            fprintf('hp = %3.2f \t distance = %3.2f \t mean h0 = %3.2f mean h0 approx = %3.2f\n',...
                par.grid_hp(ip),mean_h0_cond(ip)-mean_h0_cond_approx(ip),...
                mean_h0_cond(ip),mean_h0_cond_approx(ip));
            
        end
end