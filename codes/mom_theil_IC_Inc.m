function [ count_m,m_model,m_title,m_mod_id ] = mom_theil_IC_Inc( par,count_m,m_model,m_title,m_mod_id,mu_IC_Inc,grid_IC_Inc,options )

switch options.ComputeOtherMus
    case 'Y'
        fprintf('Doing mom_theil_IC_Inc \n');
        
        j_c_pos     = par.Je1_pos;
        S0          = par.grids{1,j_c_pos};
        H0          = par.gridh{1,j_c_pos};
        PSY         = par.psy_val_hs;
        Outcome     = grid_IC_Inc;
        income_ige  = nan(length(S0)*length(H0)*length(PSY)*length(Outcome),5);
        N_inc       = length(Outcome);
        count       = 1;
        for s_c = 1:length(S0)
            for h_c = 1:length(H0)
                for ipsy = 1:length(PSY)
                    pr   = squeeze(mu_IC_Inc(s_c,h_c,ipsy,:));
                    income_ige(count:count+N_inc-1,1) = ipsy;
                    income_ige(count:count+N_inc-1,2) = Outcome;
                    income_ige(count:count+N_inc-1,3) = S0(s_c);
                    income_ige(count:count+N_inc-1,4) = H0(h_c);
                    income_ige(count:count+N_inc-1,5) = pr;
                    count = count + N_inc;
                end
            end
        end

        % Types: Initial conditions: transfers, H0 and PSY
        types_num       = length(S0)*length(H0)*length(PSY);
        theil           = zeros(types_num,2);
        type_pos = 1;
        for iphi = 1:length(S0)
            for h0 = 1:length(H0)
                for ipsy = 1:length(PSY)
                    ind = logical((income_ige(:,3) == S0(iphi)) .* (income_ige(:,4) == H0(h0))  .* (income_ige(:,1) == ipsy));
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
m_title{count_m}    = 'Theil-L index absolute (InCond, inc) ';
m_mod_id(count_m)   = 22;
count_m             = count_m+1;

theil_rel = theil_abs/theil_tot;
m_model(count_m)    = theil_rel ;
m_title{count_m}    = 'Theil-L index relative (InCond, inc) ';
m_mod_id(count_m)   = 23;
count_m             = count_m+1;

end