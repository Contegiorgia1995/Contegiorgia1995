function [ count_m,m_model,m_title,m_mod_id ] = mom_theil_ParInc( par,count_m,m_model,m_title,m_mod_id,mu_chetty_large,Inc_par_grid,options)
EDUC      = par.educ;
FE_pos    = par.inc.fe_pos;
INNO_pos  = par.inc.inno_pos;

r_sav     = par.r_sav;
r_debt    = par.r_debt;

switch options.ComputeOtherMus
    case 'Y'
        fprintf('Doing mom_theil_ParInc \n');
        % mu_chetty_large(educ_p,Inc_par,Spchild,FeChild,EDUCchild,InnoChild)
        tot_quant   = 100;
        j_c_pos     = par.Jc_pos;
        S_ch        = par.grids{1,j_c_pos};
        FE_pos      = par.inc.fe_pos;
        j_p_pos     = find(par.age == 40,1,'first');
        
        N_inc        =  max([length(Inc_par_grid{1,1}) length(Inc_par_grid{2,1}) length(Inc_par_grid{3,1})]);
        income_ige   = nan(length(EDUC)*N_inc*length(S_ch)*length(FE_pos)*length(EDUC)*length(INNO_pos),5);
        
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
        						  r     = (r_sav .* (S_ch(s_c)>=0) + r_debt .* (S_ch(s_c)<0))';
        %                         inc_c = r. * S_ch(s_c) + par.w * H_ch(h_c);
                                inc_c = par.w * exp(AGE_PROF) .* exp(FE) .* exp(INNO) ;
                                income_ige(count:count+N_inc-1,1) = inc_p;
                                income_ige(count:count+N_inc-1,2) = inc_c;
                                income_ige(count:count+N_inc-1,3) = educ_p;
                                income_ige(count:count+N_inc-1,4) = ife_c;
                                income_ige(count:count+N_inc-1,5) = pr;
                                count = count + N_inc;
                        end
                    end
                end
            end
        end
        
        % Rank Parents
        ind        = (income_ige(:,5)>0);
        income_ige = income_ige(ind,:);
        [~,pos]    = sort(income_ige(:,1));
        income_ige = income_ige(pos,:);
        income_ige = [income_ige cumsum(income_ige(:,5))];
        
        % Quantile thresholds:
        nquant              = tot_quant;
        for i = 1:nquant-1
            bot  = (1/nquant)*(i-1);
            top  = (1/nquant)*i;
            
            if income_ige(1,6)<top
                imax = find((income_ige(:,6)<top),1,'last');
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
                prob_tot   = income_ige(last_ind,6);
            end
            prob_new   = top - prob_tot;
            old_prob   = sum(income_ige(ind,5));
            new_prob     = old_prob - prob_new;
            new_prob_rat = new_prob/old_prob;
            num        = sum(ind);
            
            income_ige = [income_ige(:,1:5); income_ige(ind,1:5)];
            
            income_ige(end-num+1:end,5) = income_ige(ind,5).*(1-new_prob_rat);
            income_ige(end-num+1:end,1) = income_ige(ind,1).*(1 - 0.00001);
            income_ige(ind,5) = income_ige(ind,5).*new_prob_rat;
            
            [~,pos]    = sort(income_ige(:,1));
            income_ige = income_ige(pos,:);
            income_ige = [income_ige cumsum(income_ige(:,5))];
        end
        
        income_ige = [income_ige zeros(size(income_ige,1),1)];
        
        % Quantile thresholds:
        nquant              = tot_quant;
        for i = 1:nquant
            bot  = (1/nquant)*(i-1);
            top  = (1/nquant)*i;
            if i < nquant
                ind  =  logical((income_ige(:,6)>bot) .* (income_ige(:,6) <= top));
            else
                ind  =  logical((income_ige(:,6)>bot) .* (income_ige(:,6) <= top));
            end
            income_ige(ind,7) = i;
        end
        
        % 12.1  Theil-L index 1:
        % Outcome: child's income
        % Types: Parent income decile & Parent educ
        mean_outcome        = income_ige(:,2)' * income_ige(:,5);
        theil_tot           = (log(mean_outcome) - log(income_ige(:,2)))' * income_ige(:,5);
        
        nquant = 10;
        norig  = tot_quant;
        theil  = nan(nquant*length(EDUC),2);
        group = 1;
        for educ_p = 1:length(EDUC)
            for ip = 1:nquant
                bot_p  = (1/nquant)*(ip-1)*norig;
                top_p  = (1/nquant)*ip*norig;
                ind_p  = logical((income_ige(:,7)>bot_p) .* (income_ige(:,7) <= top_p) .* (income_ige(:,3) == educ_p));
                if sum(ind_p)>0
                    theil(group,1) = sum(income_ige(ind_p,5));
                    theil(group,2) = income_ige(ind_p,2)' * income_ige(ind_p,5) ./ theil(group,1);
                    group = group + 1;
                end
            end
        end
        
        theil = theil(1:group-1,:);
        
        theil_abs = (log(mean_outcome) - log(theil(:,2)))' * theil(:,1);
        m_model(count_m)    = theil_abs;
        m_title{count_m}    = 'Theil-L index absolute (par inc&educ) ';
        m_mod_id(count_m)   = 17;
        count_m             = count_m+1;
        
        theil_rel = theil_abs/theil_tot;
        m_model(count_m)    = theil_rel ;
        m_title{count_m}    = 'Theil-L index relative (par inc&educ) ';
        m_mod_id(count_m)   = 18;
        count_m             = count_m+1;
        
    case 'N'
        theil_tot   = NaN;
        theil_abs   = NaN;
        
        m_model(count_m)    = theil_abs;
        m_title{count_m}    = 'Theil-L index absolute (par inc&educ) ';
        m_mod_id(count_m)   = 17;
        count_m             = count_m+1;
        
        theil_rel = theil_abs/theil_tot;
        m_model(count_m)    = theil_rel ;
        m_title{count_m}    = 'Theil-L index relative (par inc&educ) ';
        m_mod_id(count_m)   = 18;
        count_m             = count_m+1;     
end
end