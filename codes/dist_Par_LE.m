function [mu_Par_LE] = dist_Par_LE(par,mu_par_ic,Inc_par_grid,mu_LE,Grid_Avg_LE)
% save aux_debug_dist_Par_LE.mat

% mu_Par_LE(Par_educ,Par_inc,LE)
% mu_par_ic(educ_par,inc_par,S0_c,H0_c,psy_c) (grid: Inc_par_grid)
% mu_LE(S0,H0,psy,educ,avg_le) (grid: Grid_Avg_LE)
dist_par = sum(sum(sum(mu_par_ic,5),4),3);

EDUC            = par.educ;
j_pos           = par.Je1_pos;
PSY             = par.psy_val_hs;
S0              = par.grids{1,j_pos};
H0              = par.gridh{1,j_pos};

%% mu_Par_LE
mu_Par_LE  = zeros(length(EDUC),length(Inc_par_grid),length(EDUC)*length(Grid_Avg_LE(1,:)));

% Parent
for p_e = 1:length(EDUC)
    prob_e      = dist_par(p_e,:);
    ind_e       = (prob_e>0);
    aux_e       = 1:length(Inc_par_grid);
    for p_inc = aux_e(ind_e)
        prob_par    = dist_par(p_e,p_inc);
        
        prob_s      = squeeze(sum(sum(mu_par_ic(p_e,p_inc,:,:,:),5),4));
        ind_s       = (prob_s>0);
        aux_s       = 1:length(S0);
        % Initial conditions
        for is = aux_s(ind_s)
            prob_h      = squeeze(sum(mu_par_ic(p_e,p_inc,is,:,:),5));
            ind_h       = (prob_h>0);
            aux_h       = 1:length(H0);
            for ih = aux_h(ind_h)
                prob_p      = squeeze(mu_par_ic(p_e,p_inc,is,ih,:));
                ind_p       = (prob_p>0);
                aux_p       = 1:length(PSY);
                for ip = aux_p(ind_p)
                    prob_ic_par = mu_par_ic(p_e,p_inc,is,ih,ip)/prob_par;
                    aux         = mu_par_ic(:,:,is,ih,ip);
                    prob_ic     = sum(aux(:));
                    % Avg LE
                    prob_le_ic = squeeze(mu_LE(is,ih,ip,:))/prob_ic;
                    mu_Par_LE(p_e,p_inc,:) = squeeze(mu_Par_LE(p_e,p_inc,:)) + prob_par .* prob_ic_par .* prob_le_ic;
                end
            end
        end
    end
end

end

% mu_par_LE
%                         mu_all(p_e,p_inc,is,ih,ip,:) = prob_par .* prob_ic_par .* prob_le_ic;
%                     for i_avg_LE = 1:length(Grid_Avg_LE)
%                         prob_le_ic = mu_LE(is,ih,ip,i_avg_LE)/prob_ic;                        
%                         mu_all(p_e,p_inc,is,ih,ip,i_avg_LE) = prob_par * prob_ic_par * prob_le_ic;
%                     end
% mu_all          = zeros(length(EDUC),length(Inc_par_grid),length(S0),length(H0),length(PSY),length(Grid_Avg_LE));

% mu_Par_LE  = zeros(length(EDUC),length(Inc_par_grid),length(Grid_Avg_LE));
% mu_Par_LE  = squeeze(sum(sum(sum(mu_Par_LE,5),4),3));

% OLD
% % % % For parents type:
% % % for p_e = 1:length(EDUC)
% % %     for p_inc = 1:length(Inc_par_grid)
% % %         prob_ic_par                 = mu_par_ic(p_e,p_inc,:,:,:);
% % %         prob_ic_par                 = prob_ic_par(:);
% % %         prob_ic_par_sum             = sum(prob_ic_par);
% % % %         prob_ic_par                 = prob_ic_par./prob_ic_par_sum;
% % %         if prob_ic_par_sum > 0
% % %             % For Lifetime earnings value
% % %             for i_le = 1:length(Grid_Avg_LE)
% % %                 % Sum across initial conditions
% % %                 prob_le_ic                  = mu_LE(:,:,:,i_le);
% % %                 prob_le_ic                  = prob_le_ic(:);
% % %                 prob_le_ic_sum              = sum(prob_le_ic);
% % %                 if prob_le_ic_sum > 0
% % %                     mu_Par_LE(p_e,p_inc,i_le)   = sum(prob_le_ic.*prob_ic_par);
% % %                 end
% % %             end
% % %         end
% % %     end
% % %     %
% % % end
% % % 
% % % end
% % % 
% % % 
% % % %             for i_s0 = 1:length(S0)
% % % %                 for i_h0 = 1:length(H0)
% % % %                     for i_p0 = 1:length(PSY)
% % % %                         prob_le_ic                  = mu_LE(i_s0,i_h0,i_p0,i_le);
% % % %                         prob_ic_par                 = mu_par_ic(p_e,p_inc,i_s0,i_h0,i_p0);
% % % %                         mu_Par_LE(p_e,p_inc,i_le)   = mu_Par_LE(p_e,p_inc,i_le) + prob_le_ic*prob_ic_par;
% % % %                     end
% % % %                 end
% % % %             end
