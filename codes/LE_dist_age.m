function [mu_LE,Grid_Avg_LE] = LE_dist_age(par,mu0,j_pos_0,options)
% Lifetime earnings
% Enter with mu24: (S0,H0,educ)
switch options.timer_on
    case {'Y'}
        fprintf('\n Lifetime earnings \n \n');
end

% Create Grid for Lifetime Earnings
w               = par.w;
EDUC            = par.educ;
FE_pos          = par.inc.fe_pos;
INNO_pos        = par.inc.inno_pos;
curv            = 3;
if j_pos_0  <= par.Jr_pos-1
    Grid_Avg_LE     = zeros(3,200);
else
    Grid_Avg_LE     = zeros(3,length(INNO_pos)*length(FE_pos));
end
for educ = 1:length(EDUC)
    if j_pos_0  <= par.Jr_pos-1
        max_inc         = 0;
        for j_pos = j_pos_0:par.Jr_pos-1
            INNO       = repmat(reshape(par.inc.z_val{educ,1}(j_pos,:),1,length(INNO_pos)),length(FE_pos),1);
            FE         = repmat(reshape(par.inc.fe{educ,1},length(FE_pos),1),1,length(INNO_pos));
            AGE_PROF   = par.inc.age_prof{educ,1}(j_pos);
            disc       = (1/(1+par.r_sav))^(j_pos-par.Je1_pos);
            
            max_inc = max(max_inc,max(disc*exp(AGE_PROF) .* exp(FE(:)) .* exp(INNO(:))));
        end
        max_inc             = w * max_inc;
        Grid_Avg_LE(educ,:) = linspace(0,max_inc^(1/curv),200).^curv;
    else
        Grid_Avg_LE(educ,:) = exp(AGE_PROF) .* exp(FE(:)) .* exp(INNO(:)); % CHECK XXX
    end
end



%% 1.  Keep initial savings, h0 and EDUC constant and move hp and compute Lifetime Earnings.
% First year
j_pos           = j_pos_0;
S0              = par.grids{1,j_pos};
EDUC            = par.educ;

mu_LE           = zeros(length(S0),length(FE_pos),length(EDUC),length(INNO_pos),length(INNO_pos),length(Grid_Avg_LE(1,:)));
for is = 1:length(S0)
    for ife = 1:length(FE_pos)
        for ie = 1:length(EDUC)
            for inno = 1:length(INNO_pos)
                if j_pos_0 >= par.Je1_pos+1
                    mu_LE(is,ife,ie,inno,inno,1) = squeeze(sum(mu0(is,ife,ie,inno,:),5));
                else
                    mu_LE(is,ife,ie,inno,inno,1) = squeeze(sum(mu0(is,ife,:,ie,inno),3));
                end
            end
        end
    end
end
clear mu24

%% 3. iterate on this until retirement
prob_e   = squeeze(sum(sum(sum(sum(sum(mu_LE,6),5),4),2),1));
ind_e    = (prob_e>0);
aux_e    = 1:length(EDUC);



for j_pos  = j_pos_0:par.Jr_pos-1
    period        = j_pos - (j_pos_0 - 1); % First period should be 1 here.
    
    mu0int          = zeros(length(S0),length(FE_pos),length(EDUC),length(INNO_pos),length(INNO_pos),length(Grid_Avg_LE(1,:)));
    
    for educ = aux_e(ind_e)
        
        prob_fe   = squeeze(sum(sum(sum(sum(mu_LE(:,:,educ,:,:,:),6),5),4),1));
        ind_fe    = (prob_fe>0);
        aux_fe    = 1:length(FE_pos);
        for ife = aux_fe(ind_fe)
            for inno0 = 1:length(INNO_pos)
            

                % Transition of avg_LE
                if j_pos <= par.Je2_pos+1 && educ == 3
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
                clear QLE QZ

                prob_s   = squeeze(sum(sum(mu_LE(:,ife,educ,inno0,:,:),6),5));
                ind_s    = (prob_s>0);
                aux_s    = 1:length(S0);
                for s0 = aux_s(ind_s)
                    mu0_vec             = mu_LE(s0,ife,educ,inno0,:,:);
                    mup                 = Q'* mu0_vec(:);
                    mu0int(s0,ife,educ,inno0,:,:) = reshape(full(mup),1,1,1,1,length(INNO_pos),length(Grid_Avg_LE(educ,:)));
                end
            end
        end
    end
    mu_LE               = mu0int;
end
clear mu0int
% Remove Hp
mu_LE = squeeze(sum(mu_LE,5));

switch options.timer_on
    case {'Y'}
        fprintf('LE Age 24, time: %3.1f sec \n',toc)
end
end