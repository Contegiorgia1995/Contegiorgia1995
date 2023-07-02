function [mu_IC_Inc,grid_IC_Inc] = dist_IniCond_IncJc(par,muc,tau0,options)
% Type: Initial conditions (H0,S0,PSY)
% Outcome: Income at age Jc

% Enter with muc: (Sc,Hc)
switch options.timer_on
    case {'Y'}
        fprintf('\n Initial conditions and income at Jc \n \n');
end

% Create Grid for Earnings
w               = par.w;
EDUC            = par.educ;
FE_pos          = par.inc.fe_pos;
INNO_pos        = par.inc.inno_pos;

max_inc         = 0;
j_pos           = par.Jc_pos;
for educ = 1:length(EDUC)
    max_inc = max(max_inc,max(par.gridh{educ,j_pos}));
end
max_inc        = w * max_inc;
curv           = 3;
grid_IC_Inc    = linspace(0,max_inc^(1/curv),200).^curv;


%% 1. Choice of education
j_pos           = par.Je1_pos;
PSY             = par.psy_val_hs;
S0              = par.grids{1,j_pos};
H0              = par.gridh{1,j_pos};
EDUC            = par.educ;
mu0int          = zeros(length(S0),length(H0),length(PSY),length(EDUC));
mu0int(:,:,:,1) = muc;
mu0int_vec      = mu0int(:);
% extend policy for education tau
tau0p           = ones(length(S0),length(H0),length(PSY),length(EDUC));
tau0p(:,:,:,1)  = tau0;
tau0int_vec     = tau0p(:);
clear tau0p

% create transition matrix of educ: linear interpolant for policy of education
fspaceergeduc   = fundef({'spli',EDUC,0,1});
Qeduc           = funbas(fspaceergeduc,tau0int_vec);
clear tau0int_vec
% create transition matrix of fixed states: Qfixed
Qfixed          = kron(ones(length(EDUC),1),speye(length(S0)*length(H0)*length(PSY)));
% create aggregate transition: Q
Q               = dprod(Qeduc,Qfixed);
clear Qfixed Qeduc

% new distribution:
mu_int_vec      = Q' * mu0int_vec;
clear Q
mu_IC_Inc       = reshape(mu_int_vec,length(S0),length(H0),length(PSY),length(EDUC));
clear mu_int_vec
%% 2.  Keep initial savings, h0 and EDUC constant and move hp and compute Earnings.
% First year
j_pos           = par.Je1_pos;
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
  

mu0int          = zeros(length(S0),length(H0),length(PSY),length(EDUC),length(Hp),length(grid_IC_Inc));

prob_e   = squeeze(sum(sum(sum(mu_IC_Inc,3),2),1));
ind_e    = (prob_e>0);
aux_e    = 1:length(EDUC);
for educ = aux_e(ind_e)
    prob_h   = squeeze(sum(sum(mu_IC_Inc(:,:,:,educ),3),1));
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
        
        fspaceergeduc       = fundef({'spli',grid_IC_Inc,0,1});
        QLE                 = funbas(fspaceergeduc,avg_inc_p_vec);
        
        % Total Q
        Q                   = dprod(QLE,QH);
        clear QLE QH
        
        prob_s   = squeeze(sum(mu_IC_Inc(:,h0,:,educ),3));
        ind_s    = (prob_s>0);
        aux_s    = 1:length(S0);
        for s0 = aux_s(ind_s)
            prob_p   = squeeze(mu_IC_Inc(s0,h0,:,educ));
            ind_p    = (prob_p>0);
            aux_p    = 1:length(PSY);
            for psy = aux_p(ind_p)
                mu0_vec             = mu_IC_Inc(s0,h0,psy,educ);
                mup                 = Q'* mu0_vec;
                mu0int(s0,h0,psy,educ,:,:) = reshape(full(mup),1,1,1,1,length(Hp),length(grid_IC_Inc));
            end
        end
    end
end
mu_IC_Inc               = mu0int;
clear Q mu0int mup


%% 3. iterate on this until Jc

for j_pos  = par.Je1_pos+1:par.Jc_pos
    period        = j_pos - par.Je1_pos;
    Ndeltas       = length(par.deltas{1,1}(j_pos,:));
    H             = par.gridh{1,j_pos};
    Hp            = par.gridh{1,j_pos+1};
    
    if j_pos == par.Je1_pos + 1
        % H policy function
        Hpol            = zeros(length(H),Ndeltas,length(EDUC));
        Hpol_pr         = zeros(length(H),Ndeltas,length(EDUC));
        
        Hpol(:,:,1)     = par.gridh{1,j_pos}' * par.deltas{1,1}(j_pos,:);
        Hpol_pr(:,:,1)  = repmat(par.prob{1,1}(j_pos,:),length(H),1);
        
        Haux            = repmat(par.gridh{2,j_pos}',1,Ndeltas);
        Hpol(:,:,2)     = Haux + par.Reducs(1,1)* Haux.^par.Reducs(1,2);
        Hpol_pr(:,:,2)  = repmat([1 0 0],length(H),1);
        
        Haux            = repmat(par.gridh{3,j_pos}',1,Ndeltas);
        Hpol(:,:,3)     = Haux + par.Reducs(1,1)* Haux.^par.Reducs(1,2);
        Hpol_pr(:,:,3)  = repmat([1 0 0],length(H),1);
    elseif j_pos <= par.Je2_pos +1
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
    
    mu0int          = zeros(length(S0),length(H0),length(PSY),length(EDUC),length(Hp),length(grid_IC_Inc));
    
    prob_e   = squeeze(sum(sum(sum(sum(sum(mu_IC_Inc,6),5),3),2),1));
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
        QH                  = kron(ones(length(grid_IC_Inc),1),QH);
        clear Hpol_vec Hpol_pr_vec
        
        % Transition of avg_LE:  
        % NOTE: All weight in par.Jc_pos
        if j_pos <= par.Je1_pos+1 && educ >= 2
            w_aux = 0;
        elseif j_pos <= par.Je2_pos + 1 && educ == 3
            w_aux = par.w_college;
        else
            w_aux = w;
        end
        inc_t               = w_aux * H_t';
        s_aux               = gridmake(inc_t,grid_IC_Inc');
        avg_inc_p_vec       = (1/period) * s_aux(:,1) + ((period-1)/period) * s_aux(:,2);
        
        if j_pos == par.Jc_pos
            avg_inc_p_vec       =  s_aux(:,1);
        end
        
        fspaceergeduc       = fundef({'spli',grid_IC_Inc,0,1});
        QLE                 = funbas(fspaceergeduc,avg_inc_p_vec);
        
        % Total Q
        Q                   = dprod(QLE,QH);
        clear QLE QH
        
        prob_h   = squeeze(sum(sum(sum(sum(mu_IC_Inc(:,:,:,educ,:,:),6),5),3),1));
        ind_h    = (prob_h>0);
        aux_h    = 1:length(H0);
        for h0 = aux_h(ind_h)
            prob_s   = squeeze(sum(sum(sum(mu_IC_Inc(:,h0,:,educ,:,:),6),5),3));
            ind_s    = (prob_s>0);
            aux_s    = 1:length(S0);
            for s0 = aux_s(ind_s)
                prob_p   = squeeze(sum(sum(mu_IC_Inc(s0,h0,:,educ,:,:),6),5));
                ind_p    = (prob_p>0);
                aux_p    = 1:length(PSY);
                for psy = aux_p(ind_p)
                    mu0_vec             = squeeze(mu_IC_Inc(s0,h0,psy,educ,:,:));
                    mup                 = Q'* mu0_vec(:);
                    mu0int(s0,h0,psy,educ,:,:) = reshape(full(mup),1,1,1,1,length(Hp),length(grid_IC_Inc));
                end
            end
        end
    end
    mu_IC_Inc               = mu0int;
end
clear mu0int
% Remove Hp and Educ
mu_IC_Inc = squeeze(sum(sum(mu_IC_Inc,5),4));

switch options.timer_on
    case {'Y'}
        fprintf('mu_IC_Inc, time: %3.1f sec \n',toc)
end
end