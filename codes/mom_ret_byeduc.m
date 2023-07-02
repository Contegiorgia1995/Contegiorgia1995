function [ count_m,m_model,m_title,m_mod_id ] = mom_ret_byeduc( par,count_m,m_model,m_title,m_mod_id,mu,options )
FE_pos    = par.inc.fe_pos;
EDUC      = par.educ;
INNO_pos  = par.inc.inno_pos;

j_pos               = par.Jr_pos;
% S                   = par.grids{1,j_pos};
% H                   = par.gridh{1,j_pos};
assets_educ         = nan(length(EDUC),3);

npv = 0;
for j = 0:par.Jd_pos-1-par.Jr_pos
   npv              = npv + (1/(1+par.r))^j; 
end
mom_ss   = 40;
mom_oas  = 45;
mom_sav  = 50;
for educ = 1:length(EDUC)
    S           = par.grids{educ,j_pos};
    ret_val     = zeros(length(FE_pos),1);
    ret_prob    = zeros(length(FE_pos),1);
    for ife = 1:length(FE_pos)
        aux = mu{j_pos}(:,ife,educ,:);
        ret_prob(ife)   = sum(aux(:));
        ret_val(ife)    = ret_rep(par,ife,educ,options);
    end
    
    ret_mean  = npv*ret_val'*ret_prob;
    
    tau_mean = 0;
    
    aux      = squeeze(mu{j_pos}(:,:,educ,:));
    sav_prob = squeeze(sum(sum(aux,2),3));
    sav_mean = S*sav_prob;
    
    tot_assets = ret_mean+tau_mean+sav_mean;
    
    assets_educ(educ,1) = ret_mean/tot_assets;
    m_model(count_m)    = assets_educ(educ,1) ;
    m_title{count_m}    = ['Retirement Benefits Share, Education Level ',num2str(educ)]; %#ok<*AGROW>
    m_mod_id(count_m)   = mom_ss;
    mom_ss              = mom_ss + 1;
    count_m             = count_m+1;
    
    assets_educ(educ,2) = tau_mean/tot_assets;
    m_model(count_m)    = assets_educ(educ,2) ;
    m_title{count_m}    = ['Old Age Support Share, Education Level ',num2str(educ)]; %#ok<*AGROW>
    m_mod_id(count_m)   = mom_oas;
    mom_oas             = mom_oas + 1;
    count_m             = count_m+1;
    
    assets_educ(educ,3) = sav_mean/tot_assets;
    m_model(count_m)    = assets_educ(educ,3) ;
    m_mod_id(count_m)   = mom_sav;
    m_title{count_m}    = ['Savings Share, Education Level ',num2str(educ)]; %#ok<*AGROW>
    count_m             = count_m+1;
    mom_sav             = mom_sav + 1;
end


end