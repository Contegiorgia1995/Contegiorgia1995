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
FE_pos    = par.inc.fe_pos;
INNO_pos  = par.inc.inno_pos;
Grid_Avg_LE     = zeros(3,200);
curv           = 3;
for educ = 1:length(EDUC)
    max_inc         = 0;
    for j_pos = par.Je1_pos:par.Jr_pos-1;
        INNO       = repmat(reshape(par.inc.z_val{educ,1}(j_pos,:),1,length(INNO_pos)),length(FE_pos),1);
        FE         = repmat(reshape(par.inc.fe{educ,1},length(FE_pos),1),1,length(INNO_pos));
        AGE_PROF   = par.inc.age_prof{educ,1}(j_pos);
        disc       = (1/(1+par.r_sav))^(j_pos-par.Je1_pos);
        max_inc = max(max_inc,max(disc*exp(AGE_PROF) .* exp(FE(:)) .* exp(INNO(:))));
    end
    max_inc        = w * max_inc;
    Grid_Avg_LE(educ,:) = linspace(0,max_inc^(1/curv),200).^curv;
end


%% 1. Choice of education
j_pos           = par.Je1_pos;
PSY             = par.psy_val_hs;
S0              = par.grids{1,j_pos};
FE_pos    = par.inc.fe_pos;
INNO_pos  = par.inc.inno_pos;
EDUC            = par.educ;
mu0int          = zeros(length(S0),length(FE_pos),length(PSY),length(EDUC));
mu0int(:,:,:,1) = muc;
mu0int_vec      = mu0int(:);
% extend policy for education tau
tau0p           = ones(length(S0),length(FE_pos),length(PSY),length(EDUC));
tau0p(:,:,:,1)  = tau0;
tau0int_vec     = tau0p(:);
clear tau0p

% create transition matrix of educ: linear interpolant for policy of education
fspaceergeduc   = fundef({'spli',EDUC,0,1});
Qeduc           = funbas(fspaceergeduc,tau0int_vec);
clear tau0int_vec
% create transition matrix of fixed states: Qfixed
Qfixed          = kron(ones(length(EDUC),1),speye(length(S0)*length(FE_pos)*length(PSY)));
% create aggregate transition: Q
Q               = dprod(Qeduc,Qfixed);
clear Qfixed Qeduc

% new distribution:
mu_int_vec      = Q' * mu0int_vec;
clear Q
mu_LE           = reshape(mu_int_vec,length(S0),length(FE_pos),length(PSY),length(EDUC));
clear mu_int_vec

% 2.2 draw of innovation z0
% load distribution of z0 conditional on education
z_prob     = par.inc.z_prob;
% z0_prob    = [squeeze(z_prob{1,1}(par.Je2_pos,1,:))'; squeeze(z_prob{2,1}(par.Je2_pos+1,1,:))'; squeeze(z_prob{3,1}(par.Je2_pos+1,1,:))'];      % (educ,prob)
z0_prob    = [squeeze(z_prob{1,1}(par.Je1_pos,1,:))'; squeeze(z_prob{2,1}(par.Je2_pos,1,:))'; squeeze(z_prob{3,1}(par.Je2_pos+2,1,:))'];      % (educ,prob)
INNO_pos   = par.inc.inno_pos;

% create transition matrix of Z0 conditional on education
Qz                = kron(z0_prob,ones(length(S0)*length(FE_pos)*length(PSY),1));

% create transition matrix of fixed states: Qfixed
Qfixed              = speye(length(S0)*length(FE_pos)*length(PSY)*length(EDUC)); 

% create aggregate transition: Q
Q   = dprod(Qz,Qfixed); 

mup_vec = Q' * mu_LE(:);
mu_LE    = reshape(mup_vec,length(S0),length(FE_pos),length(PSY),length(EDUC),length(INNO_pos));


%% 2.  Keep Initial savings, FE and EDUC constant and move INNO to compute Lifetime Earnings.
% First year
j_pos           = par.Je1_pos;
mu0int          = zeros(length(S0),length(FE_pos),length(PSY),length(EDUC),length(INNO_pos),length(Grid_Avg_LE(1,:)));

prob_e   = squeeze(sum(sum(sum(sum(mu_LE,5),3),2),1));
ind_e    = (prob_e>0);
aux_e    = 1:length(EDUC);
for educ = aux_e(ind_e)
    prob_fe   = squeeze(sum(sum(sum(mu_LE(:,:,:,educ,:),5),3),1));
    ind_fe    = (prob_fe>0);
    aux_fe    = 1:length(FE_pos);
    for ife = aux_fe(ind_fe)
        % Transition of avg_LE
        if educ == 1
            w_aux = w;
            QZ    = squeeze(par.inc.z_prob{1,1}(j_pos,:,:));
        else
            w_aux = 0;
            QZ    = speye(length(INNO_pos));
        end
        AGE_PROF            = par.inc.age_prof{educ,1}(j_pos);
        FE                  = par.inc.fe{educ,1}(ife);
        INNO                = par.inc.z_val{educ,1}(j_pos,:);

        inc_t               = w_aux * exp(AGE_PROF) .* exp(FE) .* exp(INNO(:)) ;
        avg_inc_p_vec       = inc_t(:);
        fspaceergeduc       = fundef({'spli',Grid_Avg_LE(educ,:),0,1});
        QLE                 = funbas(fspaceergeduc,avg_inc_p_vec);
        
        % Total Q
        Q                   = dprod(QLE,QZ);
        clear QLE QZ
        
        prob_s   = squeeze(sum(sum(mu_LE(:,ife,:,educ,:),5),3));
        ind_s    = (prob_s>0);
        aux_s    = 1:length(S0);
        for s0 = aux_s(ind_s)
            prob_p   = squeeze(sum(mu_LE(s0,ife,:,educ,:),5));
            ind_p    = (prob_p>0);
            aux_p    = 1:length(PSY);
            for psy = aux_p(ind_p)
                mu0_vec             = mu_LE(s0,ife,psy,educ,:);
                mup                 = Q'* mu0_vec(:);
                mu0int(s0,ife,psy,educ,:,:) = reshape(full(mup),1,1,1,1,length(INNO_pos),length(Grid_Avg_LE(educ,:)));
            end
        end
    end
end
mu_LE               = mu0int;
clear Q mu0int mup


%% 3. iterate on this until retirement

for j_pos  = par.Je1_pos+1:par.Jr_pos-1
    period        = j_pos - (par.Je1_pos-1); % First period should be 2 here.
    
    mu0int          = zeros(length(S0),length(FE_pos),length(PSY),length(EDUC),length(INNO_pos),length(Grid_Avg_LE(1,:)));
    
    prob_e   = squeeze(sum(sum(sum(sum(sum(mu_LE,6),5),3),2),1));
    ind_e    = (prob_e>0);
    aux_e    = 1:length(EDUC);
    for educ = aux_e(ind_e)
        
        prob_fe   = squeeze(sum(sum(sum(sum(mu_LE(:,:,:,educ,:,:),6),5),3),1));
        ind_fe    = (prob_fe>0);
        aux_fe    = 1:length(FE_pos);
        for ife = aux_fe(ind_fe)
        
            % Transition of avg_LE
            if j_pos <= par.Je1_pos+1 && educ <= 2
                w_aux = 0;
                QZ    = squeeze(par.inc.z_prob{educ,1}(j_pos,:,:));
            elseif j_pos <= par.Je2_pos + 1 && educ == 3
                w_aux = par.w_college;
                QZ    = speye(length(INNO_pos));
            else
                w_aux = w;
                QZ    = squeeze(par.inc.z_prob{educ,1}(j_pos,:,:));
            end
            QZ                  = kron(ones(length(Grid_Avg_LE(educ,:)),1),QZ);

            AGE_PROF            = par.inc.age_prof{educ,1}(j_pos);
            FE                  = par.inc.fe{educ,1}(ife);
            INNO                = par.inc.z_val{educ,1}(j_pos,:);

            disc                = (1/(1+par.r_sav))^(j_pos-par.Je1_pos);
            inc_t               = disc*w_aux * exp(AGE_PROF) .* exp(FE) .* exp(INNO(:)) ;
            s_aux               = gridmake(inc_t,Grid_Avg_LE(educ,:)');
            avg_inc_p_vec       = (1/period) * s_aux(:,1) + ((period-1)/period) * s_aux(:,2);

            fspaceergeduc       = fundef({'spli',Grid_Avg_LE(educ,:),0,1});
            QLE                 = funbas(fspaceergeduc,avg_inc_p_vec);

            % Total Q
            Q                   = dprod(QLE,QZ);
            clear QLE QH
        
        
            prob_s   = squeeze(sum(sum(sum(mu_LE(:,ife,:,educ,:,:),6),5),3));
            ind_s    = (prob_s>0);
            aux_s    = 1:length(S0);
            for s0 = aux_s(ind_s)
                prob_p   = squeeze(sum(sum(mu_LE(s0,ife,:,educ,:,:),6),5));
                ind_p    = (prob_p>0);
                aux_p    = 1:length(PSY);
                for psy = aux_p(ind_p)
                    mu0_vec             = mu_LE(s0,ife,psy,educ,:);
                    mup                 = Q'* mu0_vec(:);
                    mu0int(s0,ife,psy,educ,:,:) = reshape(full(mup),1,1,1,1,length(INNO_pos),length(Grid_Avg_LE(educ,:)));
                end
            end
        end
    end
    mu_LE               = mu0int;
end
clear mu0int
% Remove InnoP,Educ
mu_LE = squeeze(sum(mu_LE,5));

switch options.timer_on
    case {'Y'}
        fprintf('LE, time: %3.1f sec \n',toc)
end
end