function [ count_m,m_model,m_title,m_mod_id ] = mom_LE_age( par,count_m,m_model,m_title,m_mod_id,IGE,j_pos_0,options )
EDUC      = par.educ;
FE_pos    = par.inc.fe_pos;
INNO_pos  = par.inc.inno_pos;

switch options.ComputeLE_age
    case 'Y'
        fprintf('Doing mom_LE_age \n');
        
        j_pos_min = j_pos_0;
        j_pos_max = find(par.age == 64,1,'first');

        LE_nosav     = zeros(size(IGE.lab{j_pos_min,2},1),1);
        LE_wsav      = zeros(size(IGE.lab{j_pos_min,2},1),1);
        for j_p_pos = j_pos_min:j_pos_max
            LE_nosav     = LE_nosav + IGE.lab{j_p_pos,2};
            LE_wsav      = LE_wsav + IGE.lab{j_p_pos,2}+IGE.sav{j_p_pos,2};
        end
        ind         = ~isnan(LE_nosav);

        LE_nosav    = log(LE_nosav(ind));
        LE_wsav     = log(LE_wsav(ind));
        
        S           = IGE.States{j_pos_0,2}(ind,1);
        FE          = IGE.States{j_pos_0,2}(ind,2);
        EDUC        = IGE.States{j_pos_0,2}(ind,3);
        if j_pos_0 < find(par.age == 22,1,'first')
            ind_coll = logical(EDUC == 3);
            IGE.States{j_pos_0,2}(ind_coll,4) = 1;
        end
        
        if size(IGE.States{j_pos_0,2}(ind,:),2) == 3
            IC_num      = findgroups(S,FE,EDUC);
        elseif size(IGE.States{j_pos_0,2}(ind,:),2) == 4
            A4      = IGE.States{j_pos_0,2}(ind,4);
            IC_num  = findgroups(S,FE,EDUC,A4);
        elseif size(IGE.States{j_pos_0,2}(ind,:),2) == 5
            A4      = IGE.States{j_pos_0,2}(ind,4);
            A5      = IGE.States{j_pos_0,2}(ind,5); 
            IC_num  = findgroups(S,FE,EDUC,A4,A5);
        end
        LE          = [S FE EDUC  IC_num LE_nosav ones(length(S),1)/length(S)];
        
%         j_pos = j_pos_0;
%         S     = par.grids{1,j_pos};
%         LE    = nan(length(S)*length(FE_pos)*length(INNO_pos)*length(EDUC)*length(Grid_Avg_LE(1,:)),6);
%         count = 1;
%         ic_num = 1;
%         for ie = 1:length(EDUC)
%             S  = par.grids{ie,j_pos};
%             for is = 1:length(S)
%                 for ife = 1:length(FE_pos)
%                     for inno = 1:length(INNO_pos)
%                         tot_prob = sum(squeeze(mu_LE(is,ife,ie,inno,:)));
%                         if tot_prob>0
%                             LE(count:count+length(Grid_Avg_LE(ie,:))-1,1) = S(is);
%                             LE(count:count+length(Grid_Avg_LE(ie,:))-1,2) = ife;
%                             LE(count:count+length(Grid_Avg_LE(ie,:))-1,3) = ie;
%                             LE(count:count+length(Grid_Avg_LE(ie,:))-1,4) = ic_num;
%                             LE(count:count+length(Grid_Avg_LE(ie,:))-1,5) = Grid_Avg_LE(ie,:)';
%                             LE(count:count+length(Grid_Avg_LE(ie,:))-1,6) = squeeze(mu_LE(is,ife,ie,inno,:));
%                             count = count+length(Grid_Avg_LE(ie,:));
%                             ic_num = ic_num + 1;
%                         end
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
        
%         [Savings,~] = grp2idx(LE(:,1));
%         [Hs,~] = grp2idx(LE(:,2));
%         [Educs,~] = grp2idx(LE(:,3));
        [ICs,~] = grp2idx(LE(:,4));
        
        mu0        = LE(:,6);
        LE_Y       = LE(:,5);
        mu2        = mu0.^0.5;
        
        clear LE
        
        % num_phi    = size(D_phi,2);
        % num_hk     = size(D_fe,2);
        % Do regression wrt dummy
        mean_LE   = mu0'*LE_Y;
        Y         = LE_Y - mean_LE;
        
        % Variance of Lifetime Earnings
        Var_LE    = mu0'*Y.^2;
        
        Yw                  = mu0.^0.5 .*Y;
        
        % From Transfers
        %         Xw                   = D_phiw(:,2:end);
%         D_phi      = dummyvar(Savings);
%         clear Savings
%         D_phiw      = repmat(mu2,1,size(D_phi,2)).*D_phi;
%         betas               = (D_phiw(:,2:end)'*D_phiw(:,2:end))\ (D_phiw(:,2:end)'*Yw);
%         clear D_phiw
%         ind_beta            = isnan(betas);
%         ind_beta            = logical(1- ind_beta);
%         betas_ok            = zeros(length(betas),1);
%         betas_ok(ind_beta)  = betas(ind_beta);
%         err_phi             = Y - D_phi(:,2:end)*betas_ok;
%         Var_er_phi          = mu0'*err_phi.^2;
%         clear D_phi
        
        % From H
        %         Xw                  = D_hw(:,2:end);
%         D_h                 = dummyvar(Hs);
%         clear Hs
%         D_hw        = repmat(mu2,1,size(D_h,2)).*D_h;
%         betas               = (D_hw(:,2:end)'*D_hw(:,2:end))\ (D_hw(:,2:end)'*Yw);
%         clear D_hw
%         ind_beta            = isnan(betas);
%         ind_beta            = logical(1- ind_beta);
%         betas_ok            = zeros(length(betas),1);
%         betas_ok(ind_beta)  = betas(ind_beta);
%         err_h              = Y - D_h(:,2:end)*betas_ok;
%         Var_er_h           = mu0'*err_h.^2;
%         clear D_h
        
        % From PSY
        %         Xw                  = D_psyw(:,2:end);
%         D_educ       = dummyvar(Educs);
%         clear Psys
%         D_educw      = repmat(mu2,1,size(D_educ,2)).*D_educ;
%         betas               = (D_educw(:,2:end)'*D_educw(:,2:end))\ (D_educw(:,2:end)'*Yw);
%         clear D_psyw
%         ind_beta            = isnan(betas);
%         ind_beta            = logical(1- ind_beta);
%         betas_ok            = zeros(length(betas),1);
%         betas_ok(ind_beta)  = betas(ind_beta);
%         err_educ             = Y - D_educ(:,2:end)*betas_ok;
%         Var_er_educ          = mu0'*err_educ.^2;
%         clear D_psy
        
        % All Initial Conditions (IC)
        %         D_ic        = dummyvar(ICs);
        %         clear ICs
        %         D_icw       = repmat(mu2,1,size(D_ic,2)).*D_ic;
        
        maxIC   = max(ICs);
        D_ic    = sparse(size(ICs,1),maxIC);
        D_icw    = sparse(size(ICs,1),maxIC);
        for iic = 1:maxIC
            ind = logical(ICs == iic);
            D_icw(ind,iic) = mu2(ind);
            D_ic(ind,iic) = 1;
        end
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

% Var_LE_FE    = 100*(Var_LE-Var_er_h)/Var_LE;
% m_model(count_m)    = Var_LE_FE;
% m_title{count_m}    = '     $\%$  expl. by FE         ';
% m_mod_id(count_m)   = 533;
% count_m             = count_m+1;
% 
% Var_LE_trans   = 100*(Var_LE-Var_er_phi)/Var_LE;
% m_model(count_m)    = Var_LE_trans;
% m_title{count_m}    = '     $\%$  expl. by transfers    ';
% m_mod_id(count_m)   = 534;
% count_m             = count_m+1;
% 
% Var_LE_educ   =  100*(Var_LE-Var_er_educ)/Var_LE;
% % 100*(Var_LE-Var_er_psy)/Var_LE;
% m_model(count_m)    = Var_LE_educ;
% m_title{count_m}    = '     $\%$  expl. by education    ';
% m_mod_id(count_m)   = 535;
% count_m             = count_m+1;
% 
% Var_LE_inter= Var_LE_ic - Var_LE_trans - Var_LE_FE -Var_LE_educ;
% m_model(count_m)    = Var_LE_inter;
% m_title{count_m}    = '     $\%$  expl. by interaction  ';
% m_mod_id(count_m)   = 536;
% count_m             = count_m+1;


end