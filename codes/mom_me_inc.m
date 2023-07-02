function [ count_m,m_model,m_title,m_mod_id, me_inc ] = mom_me_inc( par,count_m,m_model,m_title,m_mod_id,mean_income_educ,dist_educ )

ind    = logical(1-isnan(mean_income_educ));
me_inc = mean_income_educ(ind)'*dist_educ(ind);

m_model(count_m)    = me_inc/par.p_model_data;
m_title{count_m}    = 'Mean Income                        ';
m_mod_id(count_m)   = 1;
count_m             = count_m+1;

end