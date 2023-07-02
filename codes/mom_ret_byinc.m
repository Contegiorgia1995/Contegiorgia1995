function [ count_m,m_model,m_title,m_mod_id ] = mom_ret_byinc( par,count_m,m_model,m_title,m_mod_id,mu,options )
j_pos               = par.Jr_pos;
FE_pos    = par.inc.fe_pos;
EDUC      = par.educ;

NH                  = length(par.gridh{1,par.Je2_pos+2});
NS                  = length(par.grids{1,par.Je2_pos+2});
income              = nan(length(EDUC)*NH*NS,2);
count               = 1;

for educ = 1:length(EDUC)
    H       = par.gridh{educ,j_pos};
    S       = par.grids{educ,j_pos};
    for s = 1:length(S)
        aux                         = squeeze(mu{j_pos}(s,:,educ));
        pr                          = squeeze(sum(aux,1));
        income(count:count+NH-1,1)  = H;
        income(count:count+NH-1,2)  = pr;
        income(count:count+NH-1,3)  = S(s);
        income(count:count+NH-1,4)  = ret_rep(par,H,educ,options);
        income(count:count+NH-1,5)  = 0;
        count                       = count + NH;
    end
end

income(:,2) = income(:,2)./sum(income(:,2));
ind        = (income(:,2)>0);
income     = income(ind,:);
[~,pos]    = sort(income(:,1));
income     = income(pos,:);
income     = [income cumsum(income(:,2)) zeros(length(income),1)];

nquant     = 5;
for i = 1:nquant
    bot  = (1/nquant)*(i-1);
    top  = (1/nquant)*i;
    if i < nquant
        ind  =  logical((income(:,end-1)>bot) .* (income(:,end-1) <= top));
    else
        ind  =  logical((income(:,end-1)>bot));
    end
    income(ind,end) = i;
end

npv = 0;
for j = 0:par.Jd_pos-1-par.Jr_pos
    npv              = npv + (1/(1+par.r))^j;
end
mom_ss  = 245;
mom_oas = 254;
mom_sav = 263;
for inc_group = 1:nquant
    ind = (income(:,end) == inc_group);
    income_sample = income(ind,:);
    prob          = income_sample(:,2)/sum(income_sample(:,2));
    ret_val       = income_sample(:,4);
    ret_mean      = npv*ret_val'*prob;
    
    tau_val       = income_sample(:,5);
    tau_mean      = npv*tau_val'*prob;
    
    sav_val       = income_sample(:,3);
    sav_mean      = npv*sav_val'*prob;
    
    
    tot_assets = ret_mean+tau_mean+sav_mean;
    
    m_model(count_m)    = ret_mean/tot_assets;
    m_title{count_m}    = ['Retirement Benefits Share, Income group',num2str(inc_group)];
    m_mod_id(count_m)   = mom_ss;
    mom_ss              = mom_ss + 1;
    count_m             = count_m+1;
    
    m_model(count_m)    = tau_mean/tot_assets ;
    m_title{count_m}    = ['Old Age Support Share, Income group',num2str(inc_group)];
    m_mod_id(count_m)   = mom_oas;
    mom_oas             = mom_oas + 1;
    count_m             = count_m+1;
    
    m_model(count_m)    = sav_mean/tot_assets ;
    m_title{count_m}    = ['Savings Share, Income group',num2str(inc_group)];
    m_mod_id(count_m)   = mom_sav;
    mom_sav             = mom_sav + 1;
    count_m             = count_m+1;
end


end