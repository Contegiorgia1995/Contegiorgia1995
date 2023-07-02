function [ count_m,m_model,m_title,m_mod_id,mean_income_educ ] = mom_income_educ( par,count_m,m_model,m_title,m_mod_id,mu_cs,dist_age )
% keyboard
EDUC      = par.educ;
FE_pos    = par.inc.fe_pos;
INNO_pos  = par.inc.inno_pos;
AGE_PROF   = par.inc.age_prof;
mean_income_educ        = zeros(length(EDUC),1);
mean_income_educ_age    = zeros(length(EDUC),par.Jr_pos-1);
mean_log_income_educ_age= zeros(length(EDUC),par.Jr_pos-1);
var_income_educ_age     = zeros(length(EDUC),par.Jr_pos-1);
N_educ_age              = length(EDUC)*length(par.Je2_pos+1:par.Jr_pos-1);
mom_mean_educ_age       = 537;
mom_var_educ_age        = 606;
mom_mean_educ           = 182;
mom_var_educ            = 185;
mom_var_log_educ_age    = 813;
mom_var_log_educ        = 372;
mom_mean_log_educ_age   = 905;
mom_mean_log_educ       = 459;

for educ = 1:length(EDUC)
   income      = [];
   for j_pos = par.Je2_pos+1:par.Jr_pos-1
        S           = par.grids{educ,j_pos};
        INNO       = squeeze(repmat(reshape(par.inc.z_val{educ,1}(j_pos,:),1,1,length(INNO_pos)),length(S),length(FE_pos)));
        FE         = squeeze(repmat(reshape(par.inc.fe{educ,1},1,length(FE_pos)),length(S),1,length(INNO_pos)));
        AGE_PROF   = par.inc.age_prof{educ,1}(j_pos);
        
       if j_pos == par.Je2_pos+1 && educ == 3
            inc_val           = par.w_college	* exp(AGE_PROF) .* exp(FE) .* exp(INNO)./ par.time_period;
       else 
            inc_val           = par.w	* exp(AGE_PROF) .* exp(FE) .* exp(INNO)./ par.time_period;
       end
       if j_pos == par.Je1_pos
            prob_aux          = squeeze(sum(mu_cs{j_pos}(:,:,:,educ,:,:),3));
       else
           prob_aux          = mu_cs{j_pos}(:,:,educ,:,:);
       end
       inc_prob          = squeeze(sum(prob_aux,5));
       inc_prob          = inc_prob(:);
       
       income            = [income; inc_val(:) inc_prob];
       inc_prob          = inc_prob./sum(inc_prob);
       
        mean_income_educ_age(educ,j_pos) = inc_prob(:)'*inc_val(:);
        m_model(count_m) = mean_income_educ_age(educ,j_pos)/par.p_model_data;
        m_title{count_m} = ['Mean Income at Age ',num2str(par.age(j_pos)),', Education Level ',num2str(educ)]; %#ok<*AGROW>
        m_mod_id(count_m)= mom_mean_educ_age;
        count_m = count_m + 1;
        
        var_income_educ_age(educ,j_pos)  = inc_prob(:)'*(inc_val(:)-mean_income_educ_age(educ,j_pos)).^2;
        m_model(count_m) = (var_income_educ_age(educ,j_pos)).^(0.5)/mean_income_educ_age(educ,j_pos);
        m_title{count_m} = ['CV Income at Age ',num2str(par.age(j_pos)),', Education Level ',num2str(educ)]; %#ok<*AGROW>
        m_mod_id(count_m)= mom_var_educ_age;
        count_m = count_m + 1;
        
        mean_log_income_educ_age(educ,j_pos) = inc_prob(:)'*log(inc_val(:)./par.p_model_data);
        m_model(count_m) = mean_log_income_educ_age(educ,j_pos);
        m_title{count_m} = ['Mean of Income at Age ',num2str(par.age(j_pos)),', Education Level ',num2str(educ)]; %#ok<*AGROW>
        m_mod_id(count_m)= mom_mean_log_educ_age;
        count_m = count_m + 1;
        
        var_log_income_educ_age  = inc_prob(:)'*(log(inc_val(:)./par.p_model_data)-mean_log_income_educ_age(educ,j_pos)).^2;
        m_model(count_m) = var_log_income_educ_age;
        m_title{count_m} = ['Var of Income at Age ',num2str(par.age(j_pos)),', Education Level ',num2str(educ)]; %#ok<*AGROW>
        m_mod_id(count_m)= mom_var_log_educ_age;
        count_m = count_m + 1;
        
        mom_mean_educ_age = mom_mean_educ_age +1;
        mom_var_educ_age  = mom_var_educ_age +1;
        mom_mean_log_educ_age= mom_mean_log_educ_age + 1;
        mom_var_log_educ_age = mom_var_log_educ_age + 1;
   end
   income(:,2) = income(:,2)./sum(income(:,2)); 
   
   mean_income_educ(educ) = income(:,1)'*income(:,2);
   m_model(count_m) = mean_income_educ(educ)/par.p_model_data;
   m_title{count_m} = ['Mean Income with Education Level ',num2str(educ)]; %#ok<*AGROW>
   m_mod_id(count_m)= mom_mean_educ;
   count_m  = count_m +1;
   mom_mean_educ = mom_mean_educ + 1;
   
   var_income_educ(educ) = (income(:,1)'-mean_income_educ(educ)).^2*income(:,2);
   m_model(count_m) = (var_income_educ(educ))^0.5 /mean_income_educ(educ);
   m_title{count_m} = ['CV Income with Education Level ',num2str(educ)]; %#ok<*AGROW>
   m_mod_id(count_m)= mom_var_educ;
   count_m  = count_m +1;
   mom_var_educ = mom_var_educ + 1;
   
   mean_log_income_educ = log(income(:,1)./par.p_model_data)'*income(:,2);
   m_model(count_m) = mean_log_income_educ;
   m_title{count_m} = ['Mean Log Income with Education Level ',num2str(educ)]; %#ok<*AGROW>
   m_mod_id(count_m)= mom_mean_log_educ;
   count_m  = count_m +1;
   mom_mean_log_educ = mom_mean_log_educ + 1;
   
   var_log_income_educ = (log(income(:,1)./par.p_model_data)'-mean_log_income_educ).^2*income(:,2);
   m_model(count_m) = var_log_income_educ;
   m_title{count_m} = ['Var Log Income with Education Level ',num2str(educ)]; %#ok<*AGROW>
   m_mod_id(count_m)= mom_var_log_educ;
   count_m  = count_m +1;
   mom_var_log_educ = mom_var_log_educ + 1;
   
end
HS_drop_HS_grad     = mean_income_educ(1)/mean_income_educ(2);
m_model(count_m) = HS_drop_HS_grad     ;
m_title{count_m} = 'Mean inc HS Dropout / HS Grad          ';
m_mod_id(count_m)= 304;
count_m  = count_m +1;

CO_grad_HS_grad     = mean_income_educ(3)/mean_income_educ(2);
m_model(count_m) = CO_grad_HS_grad     ;
m_title{count_m} = 'Mean inc College Graduates / HS Grad   ';
m_mod_id(count_m)= 305;
count_m  = count_m +1;


% Moments 339 to 350
moms   = count_m:count_m+length(par.Je2_pos+1:par.Jr_pos-1)-1;
m_model(moms) = mean_income_educ_age(1,par.Je2_pos+1:par.Jr_pos-1)./mean_income_educ_age(2,par.Je2_pos+1:par.Jr_pos-1);
m_mod_id(moms)= 767:1:789;
count_m  = count_m+length(par.Je2_pos+1:par.Jr_pos-1)-1 +1;

% Moments 351 to 362
moms   = count_m:count_m+length(par.Je2_pos+1:par.Jr_pos-1)-1;
m_model(moms) = mean_income_educ_age(3,par.Je2_pos+1:par.Jr_pos-1)./mean_income_educ_age(2,par.Je2_pos+1:par.Jr_pos-1);
m_mod_id(moms)= 790:1:812;
count_m  = count_m+length(par.Je2_pos+1:par.Jr_pos-1)-1 +1;

% Moments 462 to 473
moms   = count_m:count_m+length(par.Je2_pos+1:par.Jr_pos-1)-1;
m_model(moms) = mean_log_income_educ_age(1,par.Je2_pos+1:par.Jr_pos-1) - mean_log_income_educ_age(2,par.Je2_pos+1:par.Jr_pos-1);
m_mod_id(moms)= 974:1:996;
count_m  = count_m+length(par.Je2_pos+1:par.Jr_pos-1)-1 +1;

% Moments 474 to 485
moms   = count_m:count_m+length(par.Je2_pos+1:par.Jr_pos-1)-1;
m_model(moms) = mean_log_income_educ_age(3,par.Je2_pos+1:par.Jr_pos-1) - mean_log_income_educ_age(2,par.Je2_pos+1:par.Jr_pos-1);
m_mod_id(moms)= 997:1:1019;
count_m  = count_m+length(par.Je2_pos+1:par.Jr_pos-1)-1 +1;

end

