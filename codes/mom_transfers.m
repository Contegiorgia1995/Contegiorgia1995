function [ count_m,m_model,m_title,m_mod_id ] = mom_transfers( par,count_m,m_model,m_title,m_mod_id,mu,me_inc,IGE )

j_pos = par.Je1_pos-1;
phi_prob   = sum(sum(mu{j_pos},2),3);
phi_val    = par.PHI/2;
mean_transfer       = phi_val*phi_prob;
cv_transfer         = ((phi_val-mean_transfer).^2 *phi_prob)^.5 / mean_transfer;
m_model(count_m)    = mean_transfer/me_inc;
m_title{count_m}    = 'Mean Transfer to Children            ';
m_mod_id(count_m)   = 98;
count_m = count_m +1;
m_model(count_m)    = cv_transfer;
m_title{count_m}    = 'CV Transfer to Children            ';
m_mod_id(count_m)   = 99;
count_m = count_m +1;

ind  = logical((par.PHI>0));
mean_trans_cond         = phi_val(ind)*phi_prob(ind)/sum(phi_prob(ind));
cv_transfer_cond        = ((phi_val(ind)-mean_trans_cond).^2 *phi_prob(ind)/sum(phi_prob(ind)))^.5 / mean_trans_cond;
m_model(count_m) = mean_trans_cond/me_inc;
m_title{count_m}    = 'Mean Transfer to Children (Cond >0)';
m_mod_id(count_m)   = 100;
count_m = count_m +1;
m_model(count_m) = cv_transfer_cond;
m_title{count_m}    = 'CV Transfer to Children (Cond >0) ';
m_mod_id(count_m)   = 101;
count_m = count_m + 1;

% Transfers by Groups
EDUC        = par.educ;
j_pos_par   = find(par.age == 42,1,'first');
S           = IGE.States{par.Je1_pos-1,2}(:,1);
FE          = IGE.States{par.Je1_pos-1,2}(:,2);
EDUC_CH     = IGE.States{par.Je1_pos-1,2}(:,4);
EDUC_PAR    = IGE.States{j_pos_par,1}(:,3);
INC_PAR     = IGE.lab{j_pos_par,1}(:,1);
S_PAR       = IGE.States{j_pos_par,1}(:,1);

tot_quant   = 5;
quants      = 100*(1:tot_quant)/tot_quant;

PrctileInc  = prctile(INC_PAR,quants);
LeftInc     = [0 PrctileInc(1:end-1)];
RightInc    = PrctileInc;

PrctileSav  = prctile(S_PAR,quants);
LeftSav     = [-Inf PrctileSav(1:end-1)];
RightSav    = PrctileSav;

RankParInc  = nan(length(S_PAR),1);
RankParSav   = nan(length(S_PAR),1);
for in = 1:length(INC_PAR)
    RankParInc(in)     = find((INC_PAR(in)> LeftInc).*(INC_PAR(in)<=RightInc));
    RankParSav(in)     = find((S_PAR(in)> LeftSav).*(S_PAR(in)<=RightSav));
end

% By education of parents
mom_tr_par_educ  = 126;
for ie = 1:length(EDUC)
    ind         = logical(EDUC_PAR == ie);
    me_trans    = mean(S(ind)/2);
    
    m_model(count_m)    = me_trans/me_inc;
    m_title{count_m}    = 'Mean Transfer';
    m_mod_id(count_m)   = mom_tr_par_educ;
    count_m             = count_m +1;
    mom_tr_par_educ     = mom_tr_par_educ + 1;
end

% By education of children
mom_tr_par_educ  = 129;
for ie = 1:length(EDUC)
    ind         = logical(EDUC_CH == ie);
    me_trans    = mean(S(ind)/2)/me_inc;
    
    m_model(count_m)    = me_trans;
    m_title{count_m}    = 'Mean Transfer';
    m_mod_id(count_m)   = mom_tr_par_educ;
    count_m             = count_m +1;
    mom_tr_par_educ     = mom_tr_par_educ + 1;
end

% By Income of Parent
mom_tr_par_educ  = 132;
for ig = 1:tot_quant
    ind         = logical(RankParInc == ig);
    me_trans    = mean(S(ind)/2)/me_inc;
    
    m_model(count_m)    = me_trans;
    m_title{count_m}    = 'Mean Transfer';
    m_mod_id(count_m)   = mom_tr_par_educ;
    count_m             = count_m +1;
    mom_tr_par_educ     = mom_tr_par_educ + 1;
end

% By Savings of Parent
mom_tr_par_educ  = 137;
for ig = 1:tot_quant
    ind         = logical(RankParSav == ig);
    me_trans    = mean(S(ind)/2);
    
    m_model(count_m)    = me_trans/me_inc;
    m_title{count_m}    = 'Mean Transfer';
    m_mod_id(count_m)   = mom_tr_par_educ;
    count_m             = count_m +1;
    mom_tr_par_educ     = mom_tr_par_educ + 1;
end

end