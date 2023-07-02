function [ count_m,m_model,m_title,m_mod_id ] = mom_educ_var_decomp( par,count_m,m_model,m_title,m_mod_id,mu_ige,options,pol )
FE_pos    = par.inc.fe_pos;
EDUC      = par.educ;
INNO_pos  = par.inc.inno_pos;

% Educ Parents - Educ Children
switch options.ComputeOtherMus
    case 'Y'
        j_pos     = par.Je1_pos;
        S         = par.grids{1,j_pos};
        PSY       = par.psy_val_hs;
        years     = zeros(length(S),length(FE_pos),length(PSY));
        phi_c     = zeros(length(S),length(FE_pos),length(PSY));
        fe_c       = zeros(length(S),length(FE_pos),length(PSY));
        psy_c     = zeros(length(S),length(FE_pos),length(PSY));
        mu0       = zeros(length(S),length(FE_pos),length(PSY));
        
        YEARS     = [8 12 16];

        for phi  = 1:length(S)
            for ife = 1:length(FE_pos)
                for psy = 1:length(PSY)
                    pr                  = mu_ige{j_pos}(:,:,:,:,phi,ife,psy,:);
                    mu0(phi,ife,psy)      = sum(pr(:));
                    phi_c(phi,ife,psy)    = phi;
                    fe_c(phi,ife,psy)      = ife;
                    psy_c(phi,ife,psy)    = psy;
                    educ                  = pol.tau0(phi,ife,psy);
                    years(phi,ife,psy)    = YEARS(educ);
                end
            end
        end

        years               = years(:);
        fe_c                 = fe_c(:);
        phi_c               = phi_c(:);
        psy_c               = psy_c(:);
        mu0                 = mu0(:).^0.5;
        ind                 = (mu0>0);
        mu0                 = mu0(ind);
        years               = years(ind);
        fe_c                 = fe_c(ind);
        phi_c               = phi_c(ind);
        psy_c               = psy_c(ind);

        [fe_c, ~]            = grp2idx(fe_c);
        [phi_c, ~]          = grp2idx(phi_c);
        [psy_c, ~]          = grp2idx(psy_c);

        D_phi               = dummyvar(phi_c);
        D_h                 = dummyvar(fe_c);
        D_psy               = dummyvar(psy_c);

        yearsw              = mu0.*years;
        D_phiw              = repmat(mu0,1,size(D_phi,2)).*D_phi;
        D_hw                = repmat(mu0,1,size(D_h,2)).*D_h;
        D_psyw              = repmat(mu0,1,size(D_psy,2)).*D_psy;
        
        % 
        % ind2                = sum(D_phiw,1);
        % ind2                = (ind2>0);
        % D_phi               = D_phi(:,ind2);
        % D_phiw              = D_phiw(:,ind2);
        % ind3                = sum(D_hw,1);
        % ind3                = (ind3>0);
        % D_h                = D_h(:,ind3);
        % D_hw               = D_hw(:,ind3);
        
        % Variance of years
        Var_years           = var(years,mu0.^2);
        if Var_years < 1e-8
            Var_years = 0;
            Var_years_educ      = Var_years;
            m_model(count_m)    = Var_years_educ;
            m_title{count_m}    = 'Variance years of education           ';
            m_mod_id(count_m)   = 56;
            count_m             = count_m+1;

            Var_years_educ_h0 = 100;
            m_model(count_m)    = Var_years_educ_h0;
            m_title{count_m}    = '     $\%$  expl. by FE           ';
            m_mod_id(count_m)   = 57;
            count_m             = count_m+1;

            Var_years_educ_trans= 100;
            m_model(count_m)    = Var_years_educ_trans;
            m_title{count_m}    = '     $\%$  expl. by transfers     ';
            m_mod_id(count_m)   = 58;
            count_m             = count_m+1;

            Var_years_educ_psy  = 100;
            m_model(count_m)    = Var_years_educ_psy;
            m_title{count_m}    = '     $\%$  expl. by school taste  ';
            m_mod_id(count_m)   = 59;
            count_m             = count_m+1;

            Var_years_educ_int  = 0;
            m_model(count_m)    = Var_years_educ_int;
            m_title{count_m}    = '     $\%$  expl. by interaction       ';
            m_mod_id(count_m)   = 60;
            count_m             = count_m+1;

        else

            % Do regression wrt dummy for PHI
%             Y                   = yearsw;
%             X                   = D_phiw;
            betas               = (D_phiw'*D_phiw)\ (D_phiw'*yearsw);
            ind_beta            = isnan(betas);
            ind_beta            = logical(1- ind_beta);
            betas_ok            = zeros(length(betas),1);
            betas_ok(ind_beta)  = betas(ind_beta);
            err_phi             = years - D_phi*betas;
            Var_er_phi          = var(err_phi,mu0.^2);
            
            % Do regression wrt dummy for H0
%             Y                   = yearsw;
%             X                   = D_hw;
            betas               = (D_hw'*D_hw)\ (D_hw'*yearsw);
            ind_beta            = isnan(betas);
            ind_beta            = logical(1- ind_beta);
            betas_ok            = zeros(length(betas),1);
            betas_ok(ind_beta)  = betas(ind_beta);
            err_hk              = years - D_h*betas;
            Var_er_hk           = var(err_hk,mu0.^2);

            % Do regression wrt dummy for school taste
%             Y                   = yearsw;
%             X                   = D_psyw;
            betas               = (D_psyw'*D_psyw)\ (D_psyw'*yearsw);
            ind_beta            = isnan(betas);
            ind_beta            = logical(1- ind_beta);
            betas_ok            = zeros(length(betas),1);
            betas_ok(ind_beta)  = betas(ind_beta);
            err_psy             = years - D_psy*betas;
            Var_er_psy          = var(err_psy,mu0.^2);
        %     Y                   = years2;
        %     X                   = D_psy2;
        %     fe_psy              = (X'*X)\ (X'*Y);
        %     err_psy             = years - D_psy*fe_psy;
        %     Var_er_psy          = var(err_psy,mu0);

            % Do regression wrt dummy for PHI and H0
        %     num_phi             = size(D_phi2,2);
        %     num_hk              = size(D_hk2,2);
        %     Y                   = years2;
        %     X                   = D_phi2 D_hk2];
        %     fe                  = (X'*X)\ (X'*Y);
        %     fe_phi2             = fe(1:num_phi);
        %     fe_hk2              = fe(num_phi+1:num_phi+num_hk);
        %     err_hk              = years - D_fe*fe_hk2;
        %     Var_er_hk2          = var(err_hk,mu0);
        %     err_phi             = years - D_phi*fe_phi2;
        %     Var_er_phi2         = var(err_phi,mu0);

            % Rename variables:
            Var_years_educ      = Var_years;
            m_model(count_m)    = Var_years_educ;
            m_title{count_m}    = 'Variance years of education           ';
            m_mod_id(count_m)   = 56;
            count_m             = count_m+1;

            Var_years_educ_h0 = 100*(Var_years-Var_er_hk)/Var_years;
            m_model(count_m)    = Var_years_educ_h0;
            m_title{count_m}    = '     $\%$  expl. by FE           ';
            m_mod_id(count_m)   = 57;
            count_m             = count_m+1;

            Var_years_educ_trans    = 100*(Var_years-Var_er_phi)/Var_years;
            m_model(count_m)    = Var_years_educ_trans;
            m_title{count_m}    = '     $\%$  expl. by transfers     ';
            m_mod_id(count_m)   = 58;
            count_m             = count_m+1;

            Var_years_educ_psy  = 100*(Var_years-Var_er_psy)/Var_years;
        %     100*(Var_years-Var_er_psy)/Var_years;
            m_model(count_m)    = Var_years_educ_psy;
            m_title{count_m}    = '     $\%$  expl. by school taste   ';
            m_mod_id(count_m)   = 59;
            count_m             = count_m+1;

            Var_years_educ_int  = 100 - Var_years_educ_trans - Var_years_educ_h0 - Var_years_educ_psy;
        %     100 - Var_years_educ_trans - Var_years_educ_h0 - Var_years_educ_psy;
            m_model(count_m)    = Var_years_educ_int;
            m_title{count_m}    = '     $\%$  expl. by interaction       ';
            m_mod_id(count_m)   = 60;
            count_m             = count_m+1;
        end
        clear phi_c years h_c psy_c mu0 yearsw D_phiw D_hw D_psyw betas ind_beta betas_ok
        
    case 'N'
            Var_years_educ      = NaN;
            m_model(count_m)    = Var_years_educ;
            m_title{count_m}    = 'Variance years of education           ';
            m_mod_id(count_m)   = 56;
            count_m             = count_m+1;

            Var_years_educ_h0 = NaN;
            m_model(count_m)    = Var_years_educ_h0;
            m_title{count_m}    = '     $\%$  expl. by FE           ';
            m_mod_id(count_m)   = 57;
            count_m             = count_m+1;

            Var_years_educ_trans= NaN;
            m_model(count_m)    = Var_years_educ_trans;
            m_title{count_m}    = '     $\%$  expl. by transfers     ';
            m_mod_id(count_m)   = 58;
            count_m             = count_m+1;

            Var_years_educ_psy  = NaN;
            m_model(count_m)    = Var_years_educ_psy;
            m_title{count_m}    = '     $\%$  expl. by school taste  ';
            m_mod_id(count_m)   = 59;
            count_m             = count_m+1;

            Var_years_educ_int  = NaN;
            m_model(count_m)    = Var_years_educ_int;
            m_title{count_m}    = '     $\%$  expl. by interaction       ';
            m_mod_id(count_m)   = 60;
            count_m             = count_m+1;
end


end