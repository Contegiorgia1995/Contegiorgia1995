function [ count_m,m_model,m_title,m_mod_id ] = mom_childcare( par,count_m,m_model,m_title,m_mod_id,mu,pol )
j_pos       = par.Jc_pos;
% Fertility mean
Nan_pol     = pol.Tp{j_pos}(:);
prob        = mu{j_pos}(:);
ind         = logical(1-isnan(Nan_pol));

Nan_pol     = Nan_pol(ind);
prob        = prob(ind);

ind2        = (Nan_pol == 1);
prob_nan    = sum(prob(ind2));
m_model(count_m)    = prob_nan;
m_title{count_m}    = 'Hire Childcare (Share)              ';
m_mod_id(count_m)   = 55;
count_m             = count_m + 1;


end