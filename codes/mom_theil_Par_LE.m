function [ count_m,m_model,m_title,m_mod_id ] = mom_theil_Par_LE( par,count_m,m_model,m_title,m_mod_id,mu_Par_LE,Inc_par_grid,Grid_Avg_LE,options )
% save aux_debug_mom_theil_Par_LE.mat

switch options.ComputeOtherMus
    case 'Y'
        fprintf('Doing mom_theil_Par_LE \n');
        
        % Type: Parents educ and Income
        EDUC        = par.educ;
        % Outcome: Average lifetime earnings
        Outcome     = Grid_Avg_LE(:);
        income_ige  = nan(length(EDUC)*length(Inc_par_grid)*length(Outcome),5);
        N_inc       = length(Outcome);
        count       = 1;
        for e_p = 1:length(EDUC)
            for inc_p = 1:length(Inc_par_grid)
                pr   = squeeze(mu_Par_LE(e_p,inc_p,:));
                income_ige(count:count+N_inc-1,1) = e_p;
                income_ige(count:count+N_inc-1,2) = Outcome;
                income_ige(count:count+N_inc-1,3) = Inc_par_grid(inc_p);
                income_ige(count:count+N_inc-1,5) = pr;
                count = count + N_inc;
            end
        end

        % Types: Parents educ and Income
        types_num       = length(EDUC)*length(Inc_par_grid);
        theil           = zeros(types_num,2);
        type_pos = 1;
        for e_p = 1:length(EDUC)
            for inc_p = 1:length(Inc_par_grid)
                ind = logical((income_ige(:,3) == Inc_par_grid(inc_p))  .* (income_ige(:,1) == e_p));
                pr = sum(income_ige(ind,5));
                if pr>0
                    theil(type_pos,1) = pr;
                    theil(type_pos,2) = income_ige(ind,2)' * income_ige(ind,5) ./ theil(type_pos,1);
                    type_pos = type_pos + 1;
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
m_title{count_m}    = 'Theil-L index absolute (parents,LE)    ';
m_mod_id(count_m)   = 24;
count_m             = count_m+1;

theil_rel = theil_abs/theil_tot;
m_model(count_m)    = theil_rel ;
m_title{count_m}    = 'Theil-L index relative (parents,LE)    ';
m_mod_id(count_m)   = 25;
count_m             = count_m+1;

end