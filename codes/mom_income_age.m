function [ count_m,m_model,m_title,m_mod_id ] = mom_income_age( par,count_m,m_model,m_title,m_mod_id,mu_cs )
FE_pos    = par.inc.fe_pos;
EDUC      = par.educ;
INNO_pos  = par.inc.inno_pos;

mean_income_age = zeros(par.Jr_pos-1,1);
var_income_age  = zeros(par.Jr_pos-1,1);
N_age = length(par.Je2_pos+1:par.Jr_pos-1);
mom_mean_age = 675;
mom_var_age  = 698;
mom_var_log_age = 882;
for j_pos = par.Je2_pos+1:par.Jr_pos-1
    S           = par.grids{1,j_pos};
    inc_val     = zeros(length(S)*length(FE_pos)*length(INNO_pos),length(EDUC));
    inc_prob    = zeros(length(S)*length(FE_pos)*length(INNO_pos),length(EDUC));

    for educ = 1:length(EDUC)
        S          = par.grids{educ,j_pos};
        INNO       = squeeze(repmat(reshape(par.inc.z_val{educ,1}(j_pos,:),1,1,length(INNO_pos)),length(S),length(FE_pos)));
        FE         = squeeze(repmat(reshape(par.inc.fe{educ,1},1,length(FE_pos)),length(S),1,length(INNO_pos)));
        AGE_PROF   = par.inc.age_prof{educ,1}(j_pos);
        
        
        if j_pos == par.Je2_pos+1 && educ == 3
            inc_val(:,educ)	= par.w_college	* exp(AGE_PROF) .* exp(FE(:)) .* exp(INNO(:)) ./ par.time_period;
        else
            inc_val(:,educ)   = par.w	* exp(AGE_PROF) .* exp(FE(:)) .* exp(INNO(:)) ./ par.time_period;
        end
        if j_pos == par.Je1_pos
            prob_aux          = squeeze(sum(mu_cs{j_pos}(:,:,:,educ,:),3));
        else
            prob_aux          = mu_cs{j_pos}(:,:,educ,:,:);
        end
        prob_aux            = squeeze(sum(prob_aux,5));
        inc_prob(:,educ)    = prob_aux(:);
    end
    inc_prob                = inc_prob./sum(inc_prob(:));
    
    mean_income_age(j_pos) = inc_prob(:)'*inc_val(:);
    m_model(count_m) = mean_income_age(j_pos)/par.p_model_data;
    m_title{count_m} = ['Mean Income at Age ',num2str(par.age(j_pos))];
    m_mod_id(count_m)= mom_mean_age;
    count_m  = count_m +1;
    mom_mean_age = mom_mean_age + 1;
    
    var_income_age(j_pos)  = inc_prob(:)'*(inc_val(:)-mean_income_age(j_pos)).^2;
    m_model(count_m) = (var_income_age(j_pos))^0.5 /mean_income_age(j_pos);
    m_title{count_m} = ['CV Income at Age ',num2str(par.age(j_pos))]; %#ok<*AGROW>
    m_mod_id(count_m)  = mom_var_age;
    count_m  = count_m +1;
    mom_var_age = mom_var_age + 1;
    
    mean_log_income_age    = inc_prob(:)'*log(inc_val(:));
    var_log_income_age     = inc_prob(:)'*(log(inc_val(:))-mean_log_income_age).^2;
    m_model(count_m) = var_log_income_age;
    m_title{count_m} = ['Var log Income at Age ',num2str(par.age(j_pos))]; %#ok<*AGROW>
    m_mod_id(count_m)  = mom_var_log_age;
    count_m  = count_m +1;
    mom_var_log_age = mom_var_log_age + 1;
end

end