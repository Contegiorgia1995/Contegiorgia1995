function [mu_LE,Grid_Avg_LE] = LE_dist(par,muc,tau0,options)
% Lifetime earnings
% Enter with muc: (Sc,Hc)
switch options.timer_on
    case {'Y'}
        fprintf('\n Lifetime earnings \n \n');
end

% Create Grid for Lifetime Earnings
w               = par.w;
EDUC            = par.educ;
max_inc         = 0;
for j_pos = par.Je2_pos:par.Jr_pos-1;
    for educ = 1:length(EDUC)
        max_inc = max(max_inc,max(par.gridh{educ,j_pos}));
    end
end
max_inc        = w * max_inc;
curv           = 3;
Grid_Avg_LE    = linspace(0,max_inc^(1/curv),200).^curv;

j_pos           = par.Je2_pos;
PSY             = par.psy_val_hs;
S0              = par.grids{1,j_pos};
H0              = par.gridh{1,j_pos};
EDUC            = par.educ;

%% 2.  Keep initial savings, h0 and EDUC constant and move hp and compute Lifetime Earnings.
% First year
j_pos           = par.Je2_pos;
Ndeltas         = length(par.deltas{1,1}(j_pos,:));

% H policy function
Hpol            = zeros(length(H0),Ndeltas,length(EDUC));
Hpol_pr         = zeros(length(H0),Ndeltas,length(EDUC));
Hpol(:,:,1)     = par.gridh{1,j_pos}' * par.deltas{1,1}(j_pos,:);
Hpol_pr(:,:,1)  = repmat(par.prob{1,1}(j_pos,:),length(H0),1);
Haux            = repmat(par.gridh{2,j_pos}',1,Ndeltas);
Hpol(:,:,2)     = Haux + par.Reducs(1,1) * Haux.^par.Reducs(1,2);
Hpol_pr(:,:,2)  = repmat([1 0 0],length(H0),1);
Haux            = repmat(par.gridh{3,j_pos}',1,Ndeltas);
Hpol(:,:,3)     = Haux + par.Reducs(1,1) * Haux.^par.Reducs(1,2);
Hpol_pr(:,:,3)  = repmat([1 0 0],length(H0),1);

Hp              = par.gridh{1,j_pos+1}';
  

mu0int          = zeros(length(S0),length(H0),length(PSY),length(EDUC),length(Hp),length(Grid_Avg_LE));

prob_e   = squeeze(sum(sum(sum(mu_LE,3),2),1));
ind_e    = (prob_e>0);
aux_e    = 1:length(EDUC);
for educ = aux_e(ind_e)
    prob_h   = squeeze(sum(sum(mu_LE(:,:,:,educ),3),1));
    ind_h    = (prob_h>0);
    aux_h    = 1:length(H0);
    for h0 = aux_h(ind_h)
        H0              = par.gridh{educ,j_pos};
        % Transition of Hp:
        Hp              = par.gridh{educ,j_pos+1}';
        QH              = sparse(1,length(Hp));
        for delta = 1:Ndeltas
            % vectorize policy
            Hpol_vec                = Hpol(h0,delta,educ);
            Hpol_vec                = Hpol_vec(:);
            Hp                      = par.gridh{educ,j_pos+1}';
            % create transition matrix of savings: linear interpolant for policy
            fspace_h                = fundef({'spli',Hp,0,1});
            Qh_aux                  = funbas(fspace_h,Hpol_vec);
            
            Hpol_pr_vec     = squeeze(Hpol_pr(h0,delta,educ));
            Hpol_pr_vec     = Hpol_pr_vec(:);
            QH              = QH + Qh_aux .* repmat(Hpol_pr_vec,1,size(Qh_aux,2));
        end
        clear Hpol_vec Hpol_pr_vec
        % Transition of avg_LE
        if educ == 1
            w_aux = w;
        else
            w_aux = 0;
        end
        inc_t               = w_aux * H0(h0);
        avg_inc_p_vec       = inc_t(:);
        fspaceergeduc       = fundef({'spli',Grid_Avg_LE,0,1});
        QLE                 = funbas(fspaceergeduc,avg_inc_p_vec);
        
        % Total Q
        Q                   = dprod(QLE,QH);
        clear QLE QH
        
        prob_s   = squeeze(sum(mu_LE(:,h0,:,educ),3));
        ind_s    = (prob_s>0);
        aux_s    = 1:length(S0);
        for s0 = aux_s(ind_s)
            prob_p   = squeeze(mu_LE(s0,h0,:,educ));
            ind_p    = (prob_p>0);
            aux_p    = 1:length(PSY);
            for psy = aux_p(ind_p)
                mu0_vec             = mu_LE(s0,h0,psy,educ);
                mup                 = Q'* mu0_vec;
                mu0int(s0,h0,psy,educ,:,:) = reshape(full(mup),1,1,1,1,length(Hp),length(Grid_Avg_LE));
            end
        end
    end
end
mu_LE               = mu0int;
clear Q mu0int mup


%% 3. iterate on this until retirement

for j_pos  = par.Je1_pos+1:par.Jr_pos-1;
    period        = j_pos - par.Je1_pos;
    Ndeltas       = length(par.deltas{1,1}(j_pos,:));
    H             = par.gridh{1,j_pos};
    Hp            = par.gridh{1,j_pos+1};
    
    if j_pos == par.Je2_pos
        % H policy function
        Hpol            = zeros(length(H),Ndeltas,length(EDUC));
        Hpol_pr         = zeros(length(H),Ndeltas,length(EDUC));
        
        Hpol(:,:,1)     = par.gridh{1,j_pos}' * par.deltas{1,1}(j_pos,:);
        Hpol_pr(:,:,1)  = repmat(par.prob{1,1}(j_pos,:),length(H),1);
        
        Hpol(:,:,2)     = par.gridh{2,j_pos}' * par.deltas{2,1}(j_pos,:);
        Hpol_pr(:,:,2)  = repmat(par.prob{2,1}(j_pos,:),length(H),1);
        
        Haux            = repmat(par.gridh{3,j_pos}',1,Ndeltas);
        Hpol(:,:,3)     = Haux + par.Reducs(2,1)* Haux.^par.Reducs(2,2);
        Hpol_pr(:,:,3)  = repmat([1 0 0],length(H),1);
        
    else
        Hpol            = zeros(length(H),Ndeltas,length(EDUC));
        Hpol_pr         = zeros(length(H),Ndeltas,length(EDUC));
        for educ = 1:length(EDUC)
            Hpol(:,:,educ)     = par.gridh{educ,j_pos}' * par.deltas{educ,1}(j_pos,:);
            Hpol_pr(:,:,educ)  = repmat(par.prob{educ,1}(j_pos,:),length(H),1);
        end
    end
    
    mu0int          = zeros(length(S0),length(H0),length(PSY),length(EDUC),length(Hp),length(Grid_Avg_LE));
    
    prob_e   = squeeze(sum(sum(sum(sum(sum(mu_LE,6),5),3),2),1));
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
        if j_pos == par.Je2_pos && educ == 3
            w_aux = par.w_college;
        else
            w_aux = w;
        end
        inc_t               = w_aux * H_t';
        s_aux               = gridmake(inc_t,Grid_Avg_LE');
        avg_inc_p_vec       = (1/period) * s_aux(:,1) + ((period-1)/period) * s_aux(:,2);
        
        fspaceergeduc       = fundef({'spli',Grid_Avg_LE,0,1});
        QLE                 = funbas(fspaceergeduc,avg_inc_p_vec);
        
        % Total Q
        Q                   = dprod(QLE,QH);
        clear QLE QH
        
        prob_h   = squeeze(sum(sum(sum(sum(mu_LE(:,:,:,educ,:,:),6),5),3),1));
        ind_h    = (prob_h>0);
        aux_h    = 1:length(H0);
        for h0 = aux_h(ind_h)
            prob_s   = squeeze(sum(sum(sum(mu_LE(:,h0,:,educ,:,:),6),5),3));
            ind_s    = (prob_s>0);
            aux_s    = 1:length(S0);
            for s0 = aux_s(ind_s)
                prob_p   = squeeze(sum(sum(mu_LE(s0,h0,:,educ,:,:),6),5));
                ind_p    = (prob_p>0);
                aux_p    = 1:length(PSY);
                for psy = aux_p(ind_p)
                    mu0_vec             = squeeze(mu_LE(s0,h0,psy,educ,:,:));
                    mup                 = Q'* mu0_vec(:);
                    mu0int(s0,h0,psy,educ,:,:) = reshape(full(mup),1,1,1,1,length(Hp),length(Grid_Avg_LE));
                end
            end
        end
    end
    mu_LE               = mu0int;
end
clear mu0int
% Remove Hp and Educ
mu_LE = squeeze(sum(sum(mu_LE,5),4));

switch options.timer_on
    case {'Y'}
        fprintf('LE, time: %3.1f sec \n',toc)
end
end