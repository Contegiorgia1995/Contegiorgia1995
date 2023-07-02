function [ count_m,m_model,m_title,m_mod_id ] = mom_LE_age24( par,count_m,m_model,m_title,m_mod_id,mu_LE,Grid_Avg_LE,options )
EDUC      = par.educ;

switch options.ComputeOtherMus
    case 'Y'
        j_pos = par.Je2_pos+2;
        S     = par.grids{1,j_pos};
        H     = par.gridh{1,j_pos};
        LE    = nan(length(S)*length(H)*length(EDUC)*length(Grid_Avg_LE),6);
        count = 1;
        ic_num = 1;
        for is = 1:length(S)
            for h = 1:length(H)
                for ie = 1:length(EDUC)
                    tot_prob = sum(squeeze(mu_LE(is,h,ie,:)));
                    if tot_prob>0
                        LE(count:count+length(Grid_Avg_LE)-1,1) = is;
                        LE(count:count+length(Grid_Avg_LE)-1,2) = h;
                        LE(count:count+length(Grid_Avg_LE)-1,3) = ie;
                        LE(count:count+length(Grid_Avg_LE)-1,4) = ic_num;
                        LE(count:count+length(Grid_Avg_LE)-1,5) = Grid_Avg_LE';
                        LE(count:count+length(Grid_Avg_LE)-1,6) = squeeze(mu_LE(is,h,ie,:));
                        count = count+length(Grid_Avg_LE);
                        ic_num = ic_num + 1;
                    end
                end
            end
        end

        ind_nan    = logical(1-isnan(LE(:,1)));
        LE         = LE(ind_nan,:);

        ind_pos    = (LE(:,6)>0);
        LE         = LE(ind_pos,:);

        [LE(:,1),~] = grp2idx(LE(:,1));
        [LE(:,2),~] = grp2idx(LE(:,2));
        [LE(:,3),~] = grp2idx(LE(:,3));
        [LE(:,4),~] = grp2idx(LE(:,4));

        mu0        = LE(:,6);
        LE_Y       = LE(:,5);
        mu2        = mu0.^0.5;
        D_phi      = dummyvar(LE(:,1));
        D_h        = dummyvar(LE(:,2));
        D_educ       = dummyvar(LE(:,3));
        D_ic        = dummyvar(LE(:,4));
        D_phiw      = repmat(mu2,1,size(D_phi,2)).*D_phi;
        D_hw        = repmat(mu2,1,size(D_h,2)).*D_h;
        D_educw      = repmat(mu2,1,size(D_educ,2)).*D_educ;
        D_icw       = repmat(mu2,1,size(D_ic,2)).*D_ic;

        % num_phi    = size(D_phi,2);
        % num_hk     = size(D_fe,2);
        % Do regression wrt dummy
        mean_LE   = mu0'*LE_Y;
        Y         = LE_Y - mean_LE;

        % Variance of Lifetime Earnings
        Var_LE    = mu0'*Y.^2;

        Yw                  = mu0.^0.5 .*Y;

        % All Initial Conditions (IC)
%         Xw                  = D_icw(:,2:end);
        betas               = (D_icw(:,2:end)'*D_icw(:,2:end))\ (D_icw(:,2:end)'*Yw);
        ind_beta            = isnan(betas);
        ind_beta            = logical(1- ind_beta);
        betas_ok            = zeros(length(betas),1);
        betas_ok(ind_beta)  = betas(ind_beta);
        err_ic              = Y - D_ic(:,2:end)*betas_ok;
        Var_er_ic           = mu0'*err_ic.^2;

        % From Transfers
%         Xw                   = D_phiw(:,2:end);
        betas               = (D_phiw(:,2:end)'*D_phiw(:,2:end))\ (D_phiw(:,2:end)'*Yw);
        ind_beta            = isnan(betas);
        ind_beta            = logical(1- ind_beta);
        betas_ok            = zeros(length(betas),1);
        betas_ok(ind_beta)  = betas(ind_beta);
        err_phi             = Y - D_phi(:,2:end)*betas_ok;
        Var_er_phi          = mu0'*err_phi.^2;

        % From H
%         Xw                  = D_hw(:,2:end);
        betas               = (D_hw(:,2:end)'*D_hw(:,2:end))\ (D_hw(:,2:end)'*Yw);
        ind_beta            = isnan(betas);
        ind_beta            = logical(1- ind_beta);
        betas_ok            = zeros(length(betas),1);
        betas_ok(ind_beta)  = betas(ind_beta);
        err_h              = Y - D_h(:,2:end)*betas_ok;
        Var_er_h           = mu0'*err_h.^2;

        % From Educ
%         Xw                  = D_psyw(:,2:end);
        betas               = (D_educw(:,2:end)'*D_educw(:,2:end))\ (D_educw(:,2:end)'*Yw);
        ind_beta            = isnan(betas);
        ind_beta            = logical(1- ind_beta);
        betas_ok            = zeros(length(betas),1);
        betas_ok(ind_beta)  = betas(ind_beta);
        err_educ             = Y - D_educ(:,2:end)*betas_ok;
        Var_er_educ          = mu0'*err_educ.^2;
        
        clear LE ind_nan ind_pos mu0 LE_Y mu2 D_phi D_h D_educ D_phiw D_hw D_educw D_icw Yw betas ind_beta betas_ok
        
    case 'N'
        Var_LE     = NaN;
        mean_LE    = NaN;
        Var_er_ic    = NaN;
        Var_er_h    = NaN;
        Var_er_phi    = NaN;
        Var_er_educ    = NaN;
end
        
% Rename variables:
m_model(count_m)    = Var_LE^.5/mean_LE;
m_title{count_m}    = 'CV Lifetime Earnings           ';
m_mod_id(count_m)   = 531;
    
count_m             = count_m+1;

Var_LE_ic   = 100*(Var_LE-Var_er_ic)/Var_LE;
m_model(count_m)    = Var_LE_ic;
m_title{count_m}    = '     $\%$  expl. by initial conds    ';
m_mod_id(count_m)   = 532;
count_m             = count_m+1;

Var_LE_FE    = 100*(Var_LE-Var_er_h)/Var_LE;
m_model(count_m)    = Var_LE_FE;
m_title{count_m}    = '     $\%$  expl. by FE         ';
m_mod_id(count_m)   = 533;
count_m             = count_m+1;

Var_LE_trans   = 100*(Var_LE-Var_er_phi)/Var_LE;
m_model(count_m)    = Var_LE_trans;
m_title{count_m}    = '     $\%$  expl. by transfers    ';
m_mod_id(count_m)   = 534;
count_m             = count_m+1;

Var_LE_educ   =  100*(Var_LE-Var_er_educ)/Var_LE;
% 100*(Var_LE-Var_er_psy)/Var_LE;
m_model(count_m)    = Var_LE_educ;
m_title{count_m}    = '     $\%$  expl. by education    ';
m_mod_id(count_m)   = 535;
count_m             = count_m+1;

Var_LE_inter= Var_LE_ic - Var_LE_trans - Var_LE_FE -Var_LE_educ;
m_model(count_m)    = Var_LE_inter;
m_title{count_m}    = '     $\%$  expl. by interaction  ';
m_mod_id(count_m)   = 536;
count_m             = count_m+1;


end