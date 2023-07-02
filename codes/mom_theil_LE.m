function [ count_m,m_model,m_title,m_mod_id ] = mom_theil_LE( par,count_m,m_model,m_title,m_mod_id,mu_LE,Grid_Avg_LE,options )

switch options.ComputeOtherMus
    case 'Y'
        fprintf('Doing mom_theil_LE \n');
        
        j_c_pos     = par.Je1_pos;
        S0          = par.grids{1,j_c_pos};
        FE_pos      = par.inc.fe_pos;
        Outcome     = Grid_Avg_LE(:);
        PSY         = par.psy_val_hs;
        income_ige  = nan(length(S0)*length(FE_pos)*length(PSY)*length(Outcome),5);
        N_inc       = length(Outcome);
        count       = 1;
        for s_c = 1:length(S0)
            for ife = 1:length(FE_pos)
                for ipsy = 1:length(PSY)
                    pr   = squeeze(mu_LE(s_c,ife,ipsy,:));
                    income_ige(count:count+N_inc-1,1) = ipsy;
                    income_ige(count:count+N_inc-1,2) = Outcome;
                    income_ige(count:count+N_inc-1,3) = S0(s_c);
                    income_ige(count:count+N_inc-1,4) = ife;
                    income_ige(count:count+N_inc-1,5) = pr(:);
                    count = count + N_inc;
                end
            end
        end

        % Types: Initial conditions: transfers, FE and PSY
        types_num       = length(S0)*length(FE_pos)*length(PSY);
        theil           = zeros(types_num,2);
        type_pos = 1;
        for iphi = 1:length(S0)
            for ife = 1:length(FE_pos)
                for ipsy = 1:length(PSY)
                    ind = logical((income_ige(:,3) == S0(iphi)) .* (income_ige(:,4) == ife)  .* (income_ige(:,1) == ipsy));
                    pr = sum(income_ige(ind,5));
                    if pr>0
                        theil(type_pos,1) = pr;
                        theil(type_pos,2) = income_ige(ind,2)' * income_ige(ind,5) ./ theil(type_pos,1);
                        type_pos = type_pos + 1;
                    end
                end
            end
        end
        theil = theil(1:type_pos-1,:);
        mean_outcome        = income_ige(:,2)' * income_ige(:,5);
        theil_tot           = (log(mean_outcome) - log(max(0.00001,income_ige(:,2))))' * income_ige(:,5);
        theil_abs           = (log(mean_outcome) - log(theil(:,2)))' * theil(:,1);
    
        clear income_ige pr theil Outcome
        
    case 'N'
        theil_tot   = NaN;
        theil_abs   = NaN;
end
m_model(count_m)    = theil_abs;
m_title{count_m}    = 'Theil-L index absolute (initial cond) ';
m_mod_id(count_m)   = 96;
count_m             = count_m+1;

theil_rel = theil_abs/theil_tot;
m_model(count_m)    = theil_rel ;
m_title{count_m}    = 'Theil-L index relative (initial cond) ';
m_mod_id(count_m)   = 97;
count_m             = count_m+1;

end