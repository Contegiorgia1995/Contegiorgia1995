function [count_m,m_model,m_title,m_mod_id,Agg_K] = Aggregates(par,mu_cs,options,count_m,m_model,m_title,m_mod_id)
% keyboard
FE_pos    = par.inc.fe_pos;
EDUC      = par.educ;
INNO_pos  = par.inc.inno_pos;

% Aggregate Capital
Agg_K = 0;
sum_mu_K = 0;
for j_pos = par.Je1_pos:par.Jd_pos-1
    for educ = 1:length(EDUC)
        S     = par.grids{educ,j_pos};
        if j_pos <= par.Je1_pos
            mu_j          = squeeze(sum(mu_cs{j_pos}(:,:,:,educ,:),3));
        else
            mu_j          = squeeze(mu_cs{j_pos}(:,:,educ,:));
        end
        mu_j  = squeeze(sum(mu_j,3));
        mu_j  = squeeze(sum(mu_j,2));

        Agg_K = Agg_K + S*mu_j;
        sum_mu_K = sum_mu_K + sum(mu_j);
        
    end
end

% Aggregate Labor
Agg_L    = zeros(3,1);
sum_mu_L = zeros(3,1);
start_L  = [par.Je1_pos par.Je2_pos par.Je2_pos + 1];

for educ = 1:length(EDUC)
    for j_pos = start_L(educ):par.Jr_pos-1
        INNO       = squeeze(repmat(reshape(par.inc.z_val{educ,1}(j_pos,:),1,length(INNO_pos)),length(FE_pos),1));
        FE         = squeeze(repmat(reshape(par.inc.fe{educ,1},length(FE_pos),1),1,length(INNO_pos)));
        AGE_PROF   = par.inc.age_prof{educ,1}(j_pos);
        if j_pos <= par.Je1_pos
            mu_j  = squeeze(sum(mu_cs{j_pos}(:,:,:,educ,:),3));
        else
            mu_j  = squeeze(mu_cs{j_pos}(:,:,educ,:,:));
        end
        mu_j  = squeeze(sum(mu_j,4));
        mu_j  = squeeze(sum(mu_j,1));
        H     = exp(AGE_PROF) .* exp(FE(:)) .* exp(INNO(:)) ;
        Agg_L(educ)    = Agg_L(educ) + H'*mu_j(:);
        sum_mu_L(educ) = sum_mu_L(educ) + sum(mu_j(:));
    end
end

Agg_L_D  = Agg_L(1);
Agg_L_HS = Agg_L(2);
Agg_L_C  = Agg_L(3);

% Government expenses in Retirement Benefits
Agg_SS = 0;
sum_mu_SS = 0;
for j_pos = par.Jr_pos:par.Jd_pos-1
    for educ = 1:length(EDUC)
        val = ret_rep(par,FE_pos,educ,options);
        mu_j  = squeeze(mu_cs{j_pos}(:,:,educ,:));
        mu_j  = squeeze(sum(mu_j,3));
        mu_j  = squeeze(sum(mu_j,1));
        Agg_SS = Agg_SS + val*mu_j';
        sum_mu_SS = sum_mu_SS + sum(mu_j);
    end
end

m_model(count_m)    = Agg_K;
m_title{count_m}    = 'Aggregate Capital                   ';
m_mod_id(count_m)   = 289;
count_m             = count_m+1;

m_model(count_m)    = Agg_L_D;
m_title{count_m}    = 'Aggregate H HS drop                 ';
m_mod_id(count_m)   = 290;
count_m             = count_m+1;

m_model(count_m)    = Agg_L_HS;
m_title{count_m}    = 'Aggregate H HS grad                 ';
m_mod_id(count_m)   = 291;
count_m             = count_m+1;

m_model(count_m)    = Agg_L_C;
m_title{count_m}    = 'Aggregate H College grad             ';
m_mod_id(count_m)   = 292;
count_m             = count_m+1;

m_model(count_m)    = Agg_SS;
m_title{count_m}    = 'Gov exp SS                           ';
m_mod_id(count_m)   = 293;
count_m             = count_m+1;


Agg_H               = [Agg_L_D Agg_L_HS Agg_L_C];
H_tot               = (sum(par.prod_s .* Agg_H.^par.prod_rho))^1/par.prod_rho;
K_tot               = max(Agg_K - par.gov_debt,1e-6);

r_new               = par.prod_alpha*(H_tot/K_tot)^(1-par.prod_alpha)-par.prod_delta;
w_new               = par.prod_s .* (1-par.prod_alpha) .* (K_tot/H_tot)^par.prod_alpha .* (H_tot./Agg_H).^(1-par.prod_rho);

m_model(count_m)    = r_new/par.r;
m_title{count_m}    = 'Implied/guess  r                     ';
m_mod_id(count_m)   = 294;
count_m             = count_m+1;

m_model(count_m)    = w_new(1);
m_title{count_m}    = 'Implied w(1)                 ';
m_mod_id(count_m)   = 295;
count_m             = count_m+1;

m_model(count_m)    = w_new(2);
m_title{count_m}    = 'Implied w(2)                 ';
m_mod_id(count_m)   = 296;
count_m             = count_m+1;

m_model(count_m)    = w_new(3);
m_title{count_m}    = 'Implied w(3)                 ';
m_mod_id(count_m)   = 297;
count_m             = count_m+1;


end
