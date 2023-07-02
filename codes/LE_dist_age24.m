function [mu_LE,Grid_Avg_LE] = LE_dist_age24(par,mu24,options)
% Lifetime earnings
% Enter with mu24: (S0,H0,educ)
switch options.timer_on
    case {'Y'}
        fprintf('\n Lifetime earnings \n \n');
end

% Create Grid for Lifetime Earnings
w               = par.w;
EDUC            = par.educ;
max_inc         = 0;
for j_pos = par.Je2_pos+2:par.Jr_pos-1;
    for educ = 1:length(EDUC)
        max_inc = max(max_inc,max(par.gridh{educ,j_pos}));
    end
end
max_inc        = w * max_inc;
curv           = 3;
Grid_Avg_LE    = linspace(0,max_inc^(1/curv),200).^curv;


%% 1.  Keep initial savings, h0 and EDUC constant and move hp and compute Lifetime Earnings.
% First year
j_pos           = par.Je2_pos+2;
S0              = par.grids{1,j_pos};
H0              = par.gridh{1,j_pos};
EDUC            = par.educ;
Ndeltas         = length(par.deltas{1,1}(j_pos,:));

% H policy function
Hpol            = zeros(length(H0),Ndeltas,length(EDUC));
Hpol_pr         = zeros(length(H0),Ndeltas,length(EDUC));
for educ = 1:length(EDUC)
    Hpol(:,:,educ)     = par.gridh{educ,j_pos}' * par.deltas{educ,1}(j_pos,:);
    Hpol_pr(:,:,educ)  = repmat(par.prob{educ,1}(j_pos,:),length(H0),1);
end
Hp              = par.gridh{1,j_pos}'; %In first period, Hp is H today.

mu_LE           = zeros(length(S0),length(H0),length(EDUC),length(Hp),length(Grid_Avg_LE));
for is = 1:length(S0)
    for ih = 1:length(H0)
        for ie = 1:length(EDUC)
            mu_LE(is,ih,ie,ih,1) = mu24(is,ih,ie);
        end
    end
end
clear mu24

%% 3. iterate on this until retirement
for j_pos  = par.Je2_pos+2:par.Jr_pos-1;
    period        = j_pos - (par.Je2_pos+1);
    Ndeltas       = length(par.deltas{1,1}(j_pos,:));
    H             = par.gridh{1,j_pos};
    Hp            = par.gridh{1,j_pos+1};
    
    Hpol            = zeros(length(H),Ndeltas,length(EDUC));
    Hpol_pr         = zeros(length(H),Ndeltas,length(EDUC));
    for educ = 1:length(EDUC)
        Hpol(:,:,educ)     = par.gridh{educ,j_pos}' * par.deltas{educ,1}(j_pos,:);
        Hpol_pr(:,:,educ)  = repmat(par.prob{educ,1}(j_pos,:),length(H),1);
    end
    
    mu0int          = zeros(length(S0),length(H0),length(EDUC),length(Hp),length(Grid_Avg_LE));
    
    prob_e   = squeeze(sum(sum(sum(sum(mu_LE,5),4),2),1));
    ind_e    = (prob_e>0);
    aux_e    = 1:length(EDUC);
    for educ = aux_e(ind_e)
        H_t              = par.gridh{educ,j_pos};
        % Transition of Hp:
        Hp              = par.gridh{educ,j_pos+1}';
        QH              = sparse(length(H_t),length(Hp));
        for delta = 1:Ndeltas
            % vectorize policy
            Hpol_vec                = Hpol(:,delta,educ);
            Hpol_vec                = Hpol_vec(:);
            Hp                      = par.gridh{educ,j_pos+1}';
            % create transition matrix of savings: linear interpolant for policy
            fspace_h                = fundef({'spli',Hp,0,1});
            Qh_aux                  = funbas(fspace_h,Hpol_vec);
            
            Hpol_pr_vec     = squeeze(Hpol_pr(:,delta,educ));
            Hpol_pr_vec     = Hpol_pr_vec(:);
            QH              = QH + Qh_aux .* repmat(Hpol_pr_vec,1,size(Qh_aux,2));
        end
        QH                  = kron(ones(length(Grid_Avg_LE),1),QH);
        clear Hpol_vec Hpol_pr_vec
        
        % Transition of avg_LE
        inc_t               = w * H_t';
        s_aux               = gridmake(inc_t,Grid_Avg_LE');
        avg_inc_p_vec       = (1/period) * s_aux(:,1) + ((period-1)/period) * s_aux(:,2);
        
        fspaceergeduc       = fundef({'spli',Grid_Avg_LE,0,1});
        QLE                 = funbas(fspaceergeduc,avg_inc_p_vec);
        
        % Total Q
        Q                   = dprod(QLE,QH);
        clear QLE QH
        
        prob_h   = squeeze(sum(sum(sum(mu_LE(:,:,educ,:,:),5),4),1));
        ind_h    = (prob_h>0);
        aux_h    = 1:length(H0);
        for h0 = aux_h(ind_h)
            prob_s   = squeeze(sum(sum(mu_LE(:,h0,educ,:,:),5),4));
            ind_s    = (prob_s>0);
            aux_s    = 1:length(S0);
            for s0 = aux_s(ind_s)
                mu0_vec             = squeeze(mu_LE(s0,h0,educ,:,:));
                mup                 = Q'* mu0_vec(:);
                mu0int(s0,h0,educ,:,:) = reshape(full(mup),1,1,1,length(Hp),length(Grid_Avg_LE));
            end
        end
    end
    mu_LE               = mu0int;
end
clear mu0int
% Remove Hp
mu_LE = squeeze(sum(mu_LE,4));

switch options.timer_on
    case {'Y'}
        fprintf('LE Age 24, time: %3.1f sec \n',toc)
end
end