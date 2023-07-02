function [ count_m,m_model,m_title,m_mod_id,dist_educ ] = mom_dist_educ( par,count_m,m_model,m_title,m_mod_id,mu )
EDUC      = par.educ;

j_pos = par.Je1_pos;
dist_educ = zeros(length(EDUC),1);
for educ = 1:length(EDUC)
    aux  = mu{j_pos}(:,:,:,educ,:);
    dist_educ(educ) = sum(aux(:));
end

m_model(count_m)    = dist_educ(1);
m_title{count_m}    = 'HS Dropouts                        ';
m_mod_id(count_m)   = 3;
count_m = count_m +1;

m_model(count_m)   = dist_educ(2);
m_title{count_m}   = 'HS Graduates                        ';
m_mod_id(count_m)  = 4;
count_m = count_m +1;

m_model(count_m)   = dist_educ(3);
m_title{count_m}   = 'College Grads                        ';
m_mod_id(count_m)  = 5;
count_m             = count_m+1;


aux  = squeeze(sum(sum(mu{8,1},1),3));
m_model(count_m)   = par.H0' * aux';
m_title{count_m}   = 'Mean H0                        ';
m_mod_id(count_m)  = 113;
count_m             = count_m+1;
end