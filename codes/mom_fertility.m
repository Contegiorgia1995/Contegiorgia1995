function [ count_m,m_model,m_title,m_mod_id ] = mom_fertility( par,count_m,m_model,m_title,m_mod_id,mu,options )

j_pos               = par.Jc_pos+1;
% Fertility mean
N                   = par.N;
dist_fert           = zeros(size(N));
for iN = 1:length(N)
    aux             = mu{j_pos}(:,:,:,:,iN);
    dist_fert(iN)   = sum(aux(:));
end
m_model(count_m)    = dist_fert * (par.fam_size * par.N');
m_title{count_m}    = 'Mean fertility                        ';
m_mod_id(count_m)   = 6;
count_m             = count_m+1;

switch options.Fertility
    case 'Endo'
        mom_fert            = 36;
        for iN = 1:3
            m_model(count_m)    = dist_fert(iN);
            m_title{count_m} = ['Prob(Fertility =  ',num2str(par.fam_size * par.N(iN)),')']; %#ok<*AGROW>
            m_mod_id(count_m)   = mom_fert;
            count_m             = count_m+1;
            mom_fert            = mom_fert+1;
        end
        
        m_model(count_m)    = sum(dist_fert(4:end));
        m_title{count_m} = ['Prob(Fertility >=  ',num2str(par.fam_size * par.N(4)),')']; %#ok<*AGROW>
        m_mod_id(count_m)   = mom_fert;
        count_m             = count_m+1;
end


%% Fertility mean by education groups
j_pos               = par.Jc_pos+1;
% Fertility mean
N                   = par.N;
EDUC                = par.educ;
mom_fert_educ       = 019;

for educ = 1:length(EDUC)
    dist_fert           = zeros(size(N));
    mu_aux              = mu{j_pos}(:,:,educ,:,:);
    pr_educ             = sum(mu_aux(:));
    for iN = 1:length(N)
        aux             = mu{j_pos}(:,:,educ,:,iN);
        dist_fert(iN)   = sum(aux(:))/ pr_educ;
    end
    m_model(count_m)    = dist_fert * (par.fam_size * par.N');
    m_title{count_m}    = ['Mean Fertility educ ',num2str(educ),'               '];
    m_mod_id(count_m)   = mom_fert_educ;
    count_m             = count_m+1;
    mom_fert_educ = mom_fert_educ + 1;
end
end