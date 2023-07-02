function [ count_m,m_model,m_title,m_mod_id ] = mom_educ_ret(par,count_m,m_model,m_title,m_mod_id,mu_cs)

% keyboard
mu_h0       = mu_cs{1};
mu_h0       = squeeze(sum(sum(mu_h0,3),1))./sum(mu_h0(:));
FE_pos    = par.inc.fe_pos;
EDUC        = par.educ;
h_dist      = zeros(length(FE_pos),length(EDUC),length(par.age));
inc_dist    = zeros(length(FE_pos),length(EDUC),length(par.age));

%% Educ = 1
educ = 1;
for j_pos = par.Je1_pos:par.Jr_pos
    AGE_PROF                 = par.inc.age_prof{educ,1}(j_pos);
    h_dist(:,educ,j_pos)     = exp(par.inc.fe{educ,1} + AGE_PROF);
end
inc_dist(:,educ,:) 	= par.w * h_dist(:,educ,:) .* repmat(reshape((1/(1+par.r)).^(par.age/par.time_period-par.Je1_pos),1,1,length(par.age)),length(FE_pos),1,1);

%% Educ = 2
educ = 2;
j_pos = par.Je1_pos;
AGE_PROF                 = par.inc.age_prof{educ,1}(j_pos);
h_dist(:,educ,j_pos)     = exp(par.inc.fe{1,1} + AGE_PROF);
for j_pos = par.Je2_pos:par.Jr_pos
    AGE_PROF                 = par.inc.age_prof{educ,1}(j_pos);
    h_dist(:,educ,j_pos)     = exp(par.inc.fe{educ,1} + AGE_PROF);
end
inc_dist(:,educ,:) 	= par.w * h_dist(:,educ,:) .* repmat(reshape((1/(1+par.r)).^(par.age/par.time_period-par.Je1_pos),1,1,length(par.age)),length(FE_pos),1,1);
inc_dist(:,educ,par.Je1_pos) 	= 0;

%% Educ = 3
educ = 3;
j_pos = par.Je1_pos;
AGE_PROF                 = par.inc.age_prof{educ,1}(j_pos);
h_dist(:,educ,j_pos)     = exp(par.inc.fe{1,1} + AGE_PROF);
j_pos = par.Je2_pos;
AGE_PROF                 = par.inc.age_prof{educ,1}(j_pos);
h_dist(:,educ,j_pos)     = exp(par.inc.fe{2,1} + AGE_PROF);
j_pos = par.Je2_pos+1;
AGE_PROF                 = par.inc.age_prof{educ,1}(j_pos);
h_dist(:,educ,j_pos)     = exp(par.inc.fe{2,1} + AGE_PROF);

for j_pos = par.Je2_pos+3:par.Jr_pos
    AGE_PROF                 = par.inc.age_prof{educ,1}(j_pos);
    h_dist(:,educ,j_pos)     = exp(par.inc.fe{educ,1} + AGE_PROF);
end
inc_dist(:,educ,:) 	= par.w * h_dist(:,educ,:) .* repmat(reshape((1/(1+par.r)).^(par.age/par.time_period-par.Je1_pos),1,1,length(par.age)),length(FE_pos),1,1);
inc_dist(:,educ,par.Je1_pos) 	= 0;
inc_dist(:,educ,par.Je2_pos:par.Je2_pos+1) 	= par.w_college * squeeze(h_dist(:,educ,par.Je2_pos:par.Je2_pos+1)) .* (1/(1+par.r)).^(par.age(par.Je2_pos:par.Je2_pos+1)./par.time_period-par.Je1_pos);


%% Return to Educ: Life time earnings
inc_dist_lte 		= sum(inc_dist,3);
% Costs:
inc_dist_lte(:,2:3)  = inc_dist_lte(:,2:3)-par.pe1;
inc_dist_lte(:,3)    = inc_dist_lte(:,3)-par.pe2*(1/(1+par.r));

ret_educ             = [inc_dist_lte(:,2)./inc_dist_lte(:,1) inc_dist_lte(:,3)./inc_dist_lte(:,2)];
avg_ret_educ         = ret_educ' *mu_h0';

m_model(count_m)    = avg_ret_educ(1);
m_title{count_m}    = 'Net return to HS: LE                          ';
m_mod_id(count_m)   = 363;
count_m             = count_m + 1;

m_model(count_m)    = avg_ret_educ(2);
m_title{count_m}    = 'Net return to College: LE                     ';
m_mod_id(count_m)   = 364;
count_m             = count_m + 1;

%% Return to Educ: age 40
inc_dist_Jt 		= squeeze(inc_dist(:,:,find(par.age == 40,1,'first')));

ret_educ_Jt             = [inc_dist_Jt(:,2)./inc_dist_Jt(:,1) inc_dist_Jt(:,3)./inc_dist_Jt(:,2)];
avg_ret_educ_Jt         = ret_educ_Jt' *mu_h0';

m_model(count_m)    = avg_ret_educ_Jt(1);
m_title{count_m}    = 'Return to HS: age 40                          ';
m_mod_id(count_m)   = 365;
count_m             = count_m + 1;

m_model(count_m)    = avg_ret_educ_Jt(2);
m_title{count_m}    = 'Return to College: age 40                      ';
m_mod_id(count_m)   = 366;
% count_m             = count_m + 1;

end