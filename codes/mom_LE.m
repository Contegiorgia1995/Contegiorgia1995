function [ count_m,m_model,m_title,m_mod_id ] = mom_LE( par,count_m,m_model,m_title,m_mod_id,IGE,options )
EDUC      = par.educ;
FE_pos    = par.inc.fe_pos;
INNO_pos  = par.inc.inno_pos;

switch options.ComputeOtherMus
    case 'Y'
        fprintf('Doing mom_LE \n');
        
        j_pos_min = find(par.age == 22,1,'first');
        j_pos_max = find(par.age == 64,1,'first');

        LE_nosav     = zeros(size(IGE.lab{j_pos_min,2},1),1);
        LE_wsav      = zeros(size(IGE.lab{j_pos_min,2},1),1);
        for j_p_pos = j_pos_min:j_pos_max
            LE_nosav     = LE_nosav + (1/(1+par.r_sav))^(j_p_pos-j_pos_min)*IGE.lab{j_p_pos,2};
            LE_wsav      = LE_wsav  + (1/(1+par.r_sav))^(j_p_pos-j_pos_min)*(IGE.lab{j_p_pos,2}+IGE.sav{j_p_pos,2});
        end
        ind         = ~isnan(LE_nosav);

        LE_nosav    = log(LE_nosav(ind));
        LE_wsav     = log(LE_wsav(ind));
        
        
        S           = IGE.States{par.Je1_pos-1,2}(ind,5);
        FE          = IGE.States{par.Je1_pos-1,2}(ind,2);
        PSY         = IGE.States{par.Je1_pos-1,2}(ind,3);
        IC_num      = findgroups(S,FE,PSY);
        EDUC        = IGE.States{par.Je1_pos-1,2}(ind,4);
        LE          = [S FE PSY  IC_num LE_nosav EDUC];
        
%         j_pos = par.Je1_pos;
%         S     = par.grids{1,j_pos};
%         PSY   = par.psy_val_hs;
%         LE    = nan(length(S)*length(FE_pos)*length(PSY)*length(EDUC)*length(Grid_Avg_LE(1,:)),7);
%         count = 1;
%         ic_num = 1;
%         for is = 1:length(S)
%             for ife = 1:length(FE_pos)
%                 for ipsy = 1:length(PSY)
%                     tot_prob = sum(sum(squeeze(mu_LE(is,ife,ipsy,:,:))));
%                     if tot_prob>0
%                         for ie = 1:length(EDUC)
%                             LE(count:count+length(Grid_Avg_LE(ie,:))-1,1) = is;
%                             LE(count:count+length(Grid_Avg_LE(ie,:))-1,2) = ife;
%                             LE(count:count+length(Grid_Avg_LE(ie,:))-1,3) = ipsy;
%                             LE(count:count+length(Grid_Avg_LE(ie,:))-1,4) = ic_num;
%                             LE(count:count+length(Grid_Avg_LE(ie,:))-1,5) = Grid_Avg_LE(ie,:)';
%                             LE(count:count+length(Grid_Avg_LE(ie,:))-1,6) = squeeze(mu_LE(is,ife,ipsy,ie,:));
%                             LE(count:count+length(Grid_Avg_LE(ie,:))-1,7) = ie;
%                             count = count+length(Grid_Avg_LE(ie,:));
%                         end
%                         ic_num = ic_num + 1;
%                     end
%                 end
%             end
%         end
% 
%         clear mu_LE Grid_Avg_LE
%         
%         ind_nan    = logical(1-isnan(LE(:,1)));
%         LE         = LE(ind_nan,:);
% 
%         ind_pos    = (LE(:,6)>0);
%         LE         = LE(ind_pos,:);

        [Savings,~] = grp2idx(LE(:,1));
        [Hs,~]      = grp2idx(LE(:,2));
        [Psys,~]    = grp2idx(LE(:,3));
        [ICs,~]     = grp2idx(LE(:,4));
        Educs       = LE(:,6);

%         mu0        = LE(:,6);
        mu0        = ones(length(Savings),1)*1/length(Savings);
        LE_Y       = LE(:,5);
        mu2        = mu0.^0.5;
        clear LE
        
        % Variance of Log Lifetime earnings
        mean_log_LE   = mu0'*log(LE_Y);
        logY          = log(LE_Y) - mean_log_LE;
        Var_log_LE    = mu0'*logY.^2;


        % num_phi    = size(D_phi,2);
        % num_hk     = size(D_fe,2);
        % Do regression wrt dummy
        mean_LE   = mu0'*LE_Y;
        Y         = LE_Y - mean_LE;
        
        % Variance of Lifetime Earnings
        Var_LE    = mu0'*Y.^2;

        mean_LE_educ    = zeros(length(EDUC),1);
        Var_LE_educ     = zeros(length(EDUC),1);
        CV_LE_educ      = zeros(length(EDUC),1);
        for ie   = 1:length(EDUC)
            ind  = logical(Educs == ie);
            
            mu0_ie              = mu0(ind)/sum(mu0(ind));
            LE_Y_ie             = LE_Y(ind);
            mean_LE_educ(ie,1)  = mu0_ie'*LE_Y_ie;
            Y_ie                = LE_Y_ie - mean_LE_educ(ie,1);
            
            Var_LE_educ(ie,1)   = mu0_ie'*Y_ie.^2;
            
            CV_LE_educ(ie,1)    = Var_LE_educ(ie,1).^0.5 /mean_LE_educ(ie,1);
        end
        
        Yw                  = mu0.^0.5 .*Y;

        % From Transfers
%         Xw                   = D_phiw(:,2:end);
        D_phi      = dummyvar(Savings);
        clear Savings
        D_phiw      = repmat(mu2,1,size(D_phi,2)).*D_phi;
        betas               = (D_phiw(:,2:end)'*D_phiw(:,2:end))\ (D_phiw(:,2:end)'*Yw);
        clear D_phiw
        ind_beta            = isnan(betas);
        ind_beta            = logical(1- ind_beta);
        betas_ok            = zeros(length(betas),1);
        betas_ok(ind_beta)  = betas(ind_beta);
        err_phi             = Y - D_phi(:,2:end)*betas_ok;
        Var_er_phi          = mu0'*err_phi.^2;
        clear D_phi

        % From H
%         Xw                  = D_hw(:,2:end);
        D_h                 = dummyvar(Hs);
        clear Hs
        D_hw        = repmat(mu2,1,size(D_h,2)).*D_h;
        betas               = (D_hw(:,2:end)'*D_hw(:,2:end))\ (D_hw(:,2:end)'*Yw);
        clear D_hw
        ind_beta            = isnan(betas);
        ind_beta            = logical(1- ind_beta);
        betas_ok            = zeros(length(betas),1);
        betas_ok(ind_beta)  = betas(ind_beta);
        err_h              = Y - D_h(:,2:end)*betas_ok;
        Var_er_h           = mu0'*err_h.^2;
        clear D_h

        % From PSY
%         Xw                  = D_psyw(:,2:end);
        D_psy       = dummyvar(Psys);
        clear Psys
        D_psyw      = repmat(mu2,1,size(D_psy,2)).*D_psy;
        betas               = (D_psyw(:,2:end)'*D_psyw(:,2:end))\ (D_psyw(:,2:end)'*Yw);
        clear D_psyw
        ind_beta            = isnan(betas);
        ind_beta            = logical(1- ind_beta);
        betas_ok            = zeros(length(betas),1);
        betas_ok(ind_beta)  = betas(ind_beta);
        err_psy             = Y - D_psy(:,2:end)*betas_ok;
        Var_er_psy          = mu0'*err_psy.^2;
        clear D_psy
        
        % All Initial Conditions (IC)
        D_ic        = dummyvar(ICs);
        clear ICs
        D_icw       = repmat(mu2,1,size(D_ic,2)).*D_ic;
%         Xw                  = D_icw(:,2:end);
        betas               = (D_icw(:,2:end)'*D_icw(:,2:end))\ (D_icw(:,2:end)'*Yw);
        clear D_icw
        ind_beta            = isnan(betas);
        ind_beta            = logical(1- ind_beta);
        betas_ok            = zeros(length(betas),1);
        betas_ok(ind_beta)  = betas(ind_beta);
        err_ic              = Y - D_ic(:,2:end)*betas_ok;
        Var_er_ic           = mu0'*err_ic.^2;
        clear D_ic
        
        clear ind_nan ind_pos mu0 LE_Y mu2 Yw betas ind_beta betas_ok
        
    case 'N'
        Var_LE     = NaN;
        mean_LE    = NaN;
        Var_er_ic    = NaN;
        Var_er_h    = NaN;
        Var_er_phi    = NaN;
        Var_er_psy    = NaN;
        CV_LE_educ    = NaN(3,1);
        Var_log_LE      = NaN;
        
end
        
% Rename variables:
m_model(count_m)    = Var_LE^.5/mean_LE;
m_title{count_m}    = 'CV Lifetime Earnings           ';
m_mod_id(count_m)   = 283;
    
count_m             = count_m+1;

for ie = 1:3
    m_model(count_m)    = CV_LE_educ(ie,1);
    m_title{count_m}    = 'CV Lifetime Earnings Educ          ';
    m_mod_id(count_m)   = 338 + ie;

    count_m             = count_m+1;
end

Var_LE_ic   = 100*(Var_LE-Var_er_ic)/Var_LE;
m_model(count_m)    = Var_LE_ic;
m_title{count_m}    = '     $\%$  expl. by initial conds    ';
m_mod_id(count_m)   = 284;
count_m             = count_m+1;

Var_LE_FE    = 100*(Var_LE-Var_er_h)/Var_LE;
m_model(count_m)    = Var_LE_FE;
m_title{count_m}    = '     $\%$  expl. by FE         ';
m_mod_id(count_m)   = 285;
count_m             = count_m+1;

Var_LE_trans   = 100*(Var_LE-Var_er_phi)/Var_LE;
m_model(count_m)    = Var_LE_trans;
m_title{count_m}    = '     $\%$  expl. by transfers    ';
m_mod_id(count_m)   = 286;
count_m             = count_m+1;

Var_LE_psy   =  100*(Var_LE-Var_er_psy)/Var_LE;
% 100*(Var_LE-Var_er_psy)/Var_LE;
m_model(count_m)    = Var_LE_psy;
m_title{count_m}    = '     $\%$  expl. by school taste    ';
m_mod_id(count_m)   = 287;
count_m             = count_m+1;

Var_LE_inter= Var_LE_ic - Var_LE_trans - Var_LE_FE -Var_LE_psy;
m_model(count_m)    = Var_LE_inter;
m_title{count_m}    = '     $\%$  expl. by interaction  ';
m_mod_id(count_m)   = 288;
count_m             = count_m+1;


m_model(count_m)    = Var_log_LE;
m_title{count_m}    = 'Variance of log lifetime earnings';
m_mod_id(count_m)   = 342;
count_m             = count_m+1;


end