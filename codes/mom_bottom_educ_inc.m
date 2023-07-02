function [ count_m,m_model,m_title,m_mod_id ] = mom_bottom_educ_inc( par,count_m,m_model,m_title,m_mod_id,mu_chetty_large,Inc_par_grid,options)
FE_pos    = par.inc.fe_pos;
EDUC      = par.educ;
INNO_pos  = par.inc.inno_pos;
r_sav     = par.r_sav;
r_debt    = par.r_debt;

switch options.ComputeOtherMus
    case 'Y'
        fprintf('Doing mom_theil_ParInc \n');
        %% mu_chetty_large(educ_p,Inc_par,Spchild,FeChild,EDUCchild,InnoChild)
        tot_quant   = 100;
        j_c_pos     = par.Jc_pos;
        S_ch        = par.grids{1,j_c_pos};
        j_p_pos     = find(par.age == 40,1,'first');
        
        N_inc        =  max([length(Inc_par_grid{1,1}) length(Inc_par_grid{2,1}) length(Inc_par_grid{3,1})]);
        income_ige   = nan(length(EDUC)*N_inc*length(S_ch)*length(FE_pos)*length(EDUC)*length(INNO_pos),6);
        
        count    = 1;
        aux_p      = sum(sum(sum(sum(sum(mu_chetty_large,2),3),4),5),6);
        ind_educp  = (aux_p > 0);
        for educ_p = EDUC(ind_educp)
            for s_c = 1:length(S_ch)
                for educ_c = 1:length(EDUC)
                    S_ch        = par.grids{educ_c,j_c_pos};
                    AGE_PROF   = par.inc.age_prof{educ_c,1}(j_c_pos);
                    for ife_c = 1:length(FE_pos)
                        FE         = par.inc.fe{educ_c,1}(ife_c);
                        for inno_c = 1:length(INNO_pos)
                            INNO = par.inc.z_val{educ_c,1}(j_c_pos,inno_c);
                            pr    = squeeze(mu_chetty_large(educ_p,:,s_c,ife_c,educ_c,inno_c));
                            inc_p = Inc_par_grid{educ_p,1};
                            %						  r     = (r_sav .* (S_ch(s_c)>=0) + r_debt .* (S_ch(s_c)<0))';
                            %                         inc_c = r. * S_ch(s_c) + par.w * H_ch(h_c);
                            inc_c = par.w * exp(AGE_PROF) .* exp(FE) .* exp(INNO) ;
                            income_ige(count:count+N_inc-1,1) = inc_p;
                            income_ige(count:count+N_inc-1,2) = inc_c;
                            income_ige(count:count+N_inc-1,3) = educ_p;
                            income_ige(count:count+N_inc-1,4) = ife_c;
                            income_ige(count:count+N_inc-1,5) = educ_c;
                            income_ige(count:count+N_inc-1,6) = pr;
                            count = count + N_inc;
                        end
                    end
                end
            end
        end
        
        %% Rank Parents
        ind        = (income_ige(:,6)>0);
        income_ige = income_ige(ind,:);
        [~,pos]    = sort(income_ige(:,1));
        income_ige = income_ige(pos,:);
        income_ige = [income_ige cumsum(income_ige(:,6))];
        
        % Quantile thresholds:
        nquant              = tot_quant;
        for i = 1:nquant-1
            bot  = (1/nquant)*(i-1);
            top  = (1/nquant)*i;
            
            if income_ige(1,7)<top
                imax = find((income_ige(:,7)<top),1,'last');
                %         prob_tot = income_ige(imax,6);
                %         prob_new = top - prob_tot;
            else
                imax = 0;
                %         prob_new = top;
            end
            
            ind        = (income_ige(:,1) == income_ige(imax+1,1));
            last_ind   = find((income_ige(:,1) < income_ige(imax+1,1)),1,'last');
            if isempty(last_ind) == 1
                prob_tot = 0;
            else
                prob_tot   = income_ige(last_ind,7);
            end
            prob_new   = top - prob_tot;
            old_prob   = sum(income_ige(ind,6));
            new_prob     = old_prob - prob_new;
            new_prob_rat = new_prob/old_prob;
            num        = sum(ind);
            
            income_ige = [income_ige(:,1:6); income_ige(ind,1:6)];
            
            income_ige(end-num+1:end,6) = income_ige(ind,6).*(1-new_prob_rat);
            income_ige(end-num+1:end,1) = income_ige(ind,1).*(1 - 0.00001);
            income_ige(ind,6) = income_ige(ind,6).*new_prob_rat;
            
            [~,pos]    = sort(income_ige(:,1));
            income_ige = income_ige(pos,:);
            income_ige = [income_ige cumsum(income_ige(:,6))];
        end
        
        income_ige = [income_ige zeros(size(income_ige,1),1)];
        
        % Quantile thresholds:
        nquant              = tot_quant;
        for i = 1:nquant
            bot  = (1/nquant)*(i-1);
            top  = (1/nquant)*i;
            if i < nquant
                ind  =  logical((income_ige(:,7)>bot) .* (income_ige(:,7) <= top));
            else
                ind  =  logical((income_ige(:,7)>bot) .* (income_ige(:,7) <= top));
            end
            income_ige(ind,8) = i;
        end
        
        %% Rank Children
        [~,pos]   = sort(income_ige(:,2));
        income_ige = income_ige(pos,:);
        income_ige = [income_ige cumsum(income_ige(:,6))];
        
        % Quantile thresholds:
        nquant              = tot_quant;
        for i = 1:nquant-1
            bot  = (1/nquant)*(i-1);
            top  = (1/nquant)*i;
            
            if income_ige(1,9)<top
                imax = find((income_ige(:,9)<top),1,'last');
                %         prob_tot = income_ige(imax,6);
                %         prob_new = top - prob_tot;
            else
                imax = 0;
                %         prob_new = top;
            end
            
            ind        = (income_ige(:,2) == income_ige(imax+1,2));
            last_ind   = find((income_ige(:,2) < income_ige(imax+1,2)),1,'last');
            if isempty(last_ind) == 1
                prob_tot = 0;
            else
                prob_tot   = income_ige(last_ind,9);
            end
            prob_new   = top - prob_tot;
            old_prob   = sum(income_ige(ind,6));
            new_prob   = old_prob - prob_new;
            new_prob_rat = new_prob/old_prob;
            num        = sum(ind);
            
            income_ige = [income_ige(:,1:8); income_ige(ind,1:8)];
            
            income_ige(end-num+1:end,6) = income_ige(ind,6).*(1-new_prob_rat);
            income_ige(end-num+1:end,2) = income_ige(ind,2).*(1 - 0.00001);
            income_ige(ind,6)           = income_ige(ind,6).*new_prob_rat;
            
            [~,pos]    = sort(income_ige(:,2));
            income_ige = income_ige(pos,:);
            income_ige = [income_ige cumsum(income_ige(:,6))];
        end
        
        income_ige = [income_ige zeros(size(income_ige,1),1)];
        
        % Quantile thresholds:
        nquant              = tot_quant;
        for i = 1:nquant
            bot  = (1/nquant)*(i-1);
            top  = (1/nquant)*i;
            if i < nquant
                ind  =  logical((income_ige(:,9)>bot) .* (income_ige(:,9) <= top));
            else
                ind  =  logical((income_ige(:,9)>bot) .* (income_ige(:,9) <= top));
            end
            income_ige(ind,10) = i;
        end
        
        %% Parents "poor" Educ =1 , Q<=10
        ind_p               = logical((income_ige(:,3) == 1).*(income_ige(:,8) <= 10));
        pos                 = ind_p;
        
        prob_par_ind        =  sum(income_ige(pos,6));
        
        par_max_inc         = max(income_ige(pos,1))/par.p_model_data;
        m_model(count_m)    = par_max_inc;
        m_title{count_m}    = 'Parent educ = 1 inc <= Q10, Max parent income';
        m_mod_id(count_m)   = 112;
        count_m             = count_m+1;
        
        child_inc_avg       = (income_ige(pos,2)' * income_ige(pos,6) ) /prob_par_ind;
        m_model(count_m)    = child_inc_avg/par.p_model_data;
        m_title{count_m}    = 'Parent educ = 1 inc <= Q10, child avg inc';
        m_mod_id(count_m)   = 107;
        count_m             = count_m+1;
        
        child_incrank_avg   = (income_ige(pos,10)' * income_ige(pos,6) ) /prob_par_ind;
        m_model(count_m)    = child_incrank_avg;
        m_title{count_m}    = 'Parent educ = 1 inc <= Q10, child avg inc rank';
        m_mod_id(count_m)   = 108;
        count_m             = count_m+1;
        
        mean_educ           = zeros(3,1);
        for educ_c = 1:3
            ind_c               = logical((income_ige(:,5) == educ_c));
            pos                 = logical(ind_p .* ind_c);
            mean_educ(educ_c)   = sum(income_ige(pos,6))/prob_par_ind;
            
            
            m_model(count_m)    = mean_educ(educ_c);
            m_title{count_m}    = ['Parent educ = 1 inc <= Q10, child educ ',num2str(educ_c)];
            m_mod_id(count_m)   = 108+educ_c;
            count_m             = count_m+1;
        end
        
        
        
        
        
    case 'N'
        
        m_model(count_m)    = nan;
        m_title{count_m}    = 'xxx';
        m_mod_id(count_m)   = 107;
        count_m             = count_m+1;
        
        m_model(count_m)    = nan;
        m_title{count_m}    = 'xxx';
        m_mod_id(count_m)   = 108;
        count_m             = count_m+1;
        
        
        m_model(count_m)    = nan;
        m_title{count_m}    = 'xxx';
        m_mod_id(count_m)   = 109;
        count_m             = count_m+1;
        
        m_model(count_m)    = nan;
        m_title{count_m}    = 'xxx';
        m_mod_id(count_m)   = 110;
        count_m             = count_m+1;
        
        m_model(count_m)    = nan;
        m_title{count_m}    = 'xxx';
        m_mod_id(count_m)   = 111;
        count_m             = count_m+1;
        
        m_model(count_m)    = nan;
        m_title{count_m}    = 'xxx';
        m_mod_id(count_m)   = 112;
        count_m             = count_m+1;
        
end
end