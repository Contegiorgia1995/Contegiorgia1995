function [ count_m,m_model,m_title,m_mod_id ] = mom_borr_col( par,count_m,m_model,m_title,m_mod_id,mu,Sp )
j_pos  = par.Je2_pos+1;
educ    = 3;
mom_col_borr = 367;

sav_val  = Sp(:);
sav_prob = squeeze(sum(mu{j_pos}(:,:,educ,:),4));
sav_prob = sav_prob(:)./(sum(sav_prob(:)));

ind         = logical((sav_val<0));
share_debt  = sum(sav_prob(ind));
m_model(count_m)     = share_debt;
m_title{count_m}     = 'Share of Coll Students with Loans';
m_mod_id(count_m)    = mom_col_borr;
count_m = count_m +1;

end