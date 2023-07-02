function [ count_m,m_model,m_title,m_mod_id ] = mom_rank_rank( par,count_m,m_model,m_title,m_mod_id,mu_chetty,Inc_par_grid )
FE_pos    = par.inc.fe_pos;
EDUC      = par.educ;
INNO_pos  = par.inc.inno_pos;

r_sav     = par.r_sav;
r_debt    = par.r_debt;

% mu_chetty(Inc_par,Spchild,Hpchild,EDUCchild)
tot_quant   = 50;
j_c_pos     = par.Jc_pos;
S_ch        = par.grids{1,j_c_pos};
j_p_pos     = find(par.age == 40,1,'first');

N_inc        =  max([length(Inc_par_grid{1,1}) length(Inc_par_grid{2,1}) length(Inc_par_grid{3,1})]);
income_ige   = nan(length(EDUC)*N_inc*length(S_ch)*length(FE_pos)*length(EDUC)*length(INNO_pos),5);
count        = 1;

aux_p      = sum(sum(sum(sum(sum(mu_chetty,2),3),4),5),6);
ind_educp  = (aux_p > 0);
for educ_p = EDUC(ind_educp)
    for s_c = 1:length(S_ch)
        for educ_c = 1:length(EDUC)
            AGE_PROF   = par.inc.age_prof{educ_c,1}(j_c_pos);
            S_ch        = par.grids{educ_c,j_c_pos};
            for ife_c = 1:length(FE_pos)
                FE         = par.inc.fe{educ_c,1}(ife_c);
                for inno_c = 1:length(INNO_pos)
                    INNO = par.inc.z_val{educ_c,1}(j_c_pos,inno_c);
                    
                    pr   = squeeze(mu_chetty(educ_p,:,s_c,ife_c,educ_c,inno_c));
                    inc_p = Inc_par_grid{educ_p,1};
                    r     = (r_sav .* (S_ch(s_c)>=0) + r_debt .* (S_ch(s_c)<0))';
                    inc_c = r * S_ch(s_c) + par.w * exp(AGE_PROF) .* exp(FE) .* exp(INNO) ;
                    
                    income_ige(count:count+N_inc-1,1) = inc_p;
                    income_ige(count:count+N_inc-1,2) = inc_c;
                    income_ige(count:count+N_inc-1,3) = s_c;
                    income_ige(count:count+N_inc-1,4) = ife_c;
                    income_ige(count:count+N_inc-1,5) = pr;
                    count = count + N_inc;
                end
            end
        end
    end
end

% Rank Parents
ind        = (income_ige(:,5)>0);
income_ige = income_ige(ind,:);
[~,pos]    = sort(income_ige(:,1));
income_ige = income_ige(pos,:);
income_ige = [income_ige cumsum(income_ige(:,5))];

% Quantile thresholds:
nquant              = tot_quant;
for i = 1:nquant-1
    bot  = (1/nquant)*(i-1);
    top  = (1/nquant)*i;
    
    if income_ige(1,6)<top
        imax = find((income_ige(:,6)<top),1,'last');
        %         prob_tot = income_ige(imax,6);
        %         prob_new = top - prob_tot;
    else
        imax = 0;
        %         prob_new = top;
    end
    
    ind        = (income_ige(:,1) == income_ige(imax+1,1));
    last_ind   = find((income_ige(:,1) < income_ige(imax+1,1)),1,'last');
    if isempty(last_ind) == 1
        prob_tot = 0;
    else
        prob_tot   = income_ige(last_ind,6);
    end
    prob_new        = top - prob_tot;
    old_prob        = sum(income_ige(ind,5));
    new_prob        = old_prob - prob_new;
    new_prob_rat    = new_prob/old_prob;
    num             = sum(ind);
    
    income_ige = [income_ige(:,1:5); income_ige(ind,1:5)];
    
    income_ige(end-num+1:end,5) = income_ige(ind,5).*(1-new_prob_rat);
    income_ige(end-num+1:end,1) = income_ige(ind,1).*(1 - 0.00001);
    income_ige(ind,5) = income_ige(ind,5).*new_prob_rat;
    
    [~,pos]    = sort(income_ige(:,1));
    income_ige = income_ige(pos,:);
    income_ige = [income_ige cumsum(income_ige(:,5))];
end

income_ige = [income_ige zeros(size(income_ige,1),1)];

% Quantile thresholds:
nquant              = tot_quant;
for i = 1:nquant
    bot  = (1/nquant)*(i-1);
    top  = (1/nquant)*i;
    if i < nquant
        ind  =  logical((income_ige(:,6)>bot) .* (income_ige(:,6) <= top));
    else
        ind  =  logical((income_ige(:,6)>bot) .* (income_ige(:,6) <= top));
    end
    income_ige(ind,7) = i;
end

% Rank Children
[~,pos]   = sort(income_ige(:,2));
income_ige = income_ige(pos,:);
income_ige = [income_ige cumsum(income_ige(:,5))];

% Quantile thresholds:
nquant              = tot_quant;
for i = 1:nquant-1
    bot  = (1/nquant)*(i-1);
    top  = (1/nquant)*i;
    
    if income_ige(1,8)<top
        imax = find((income_ige(:,8)<top),1,'last');
        %         prob_tot = income_ige(imax,6);
        %         prob_new = top - prob_tot;
    else
        imax = 0;
        %         prob_new = top;
    end
    
    ind        = (income_ige(:,2) == income_ige(imax+1,2));
    last_ind   = find((income_ige(:,2) < income_ige(imax+1,2)),1,'last');
    if isempty(last_ind) == 1
        prob_tot = 0;
    else
        prob_tot   = income_ige(last_ind,8);
    end
    prob_new   = top - prob_tot;
    old_prob   = sum(income_ige(ind,5));
    new_prob   = old_prob - prob_new;
    new_prob_rat = new_prob/old_prob;
    num        = sum(ind);
    
    income_ige = [income_ige(:,1:7); income_ige(ind,1:7)];
    
    income_ige(end-num+1:end,5) = income_ige(ind,5).*(1-new_prob_rat);
    income_ige(end-num+1:end,2) = income_ige(ind,2).*(1 - 0.00001);
    income_ige(ind,5)           = income_ige(ind,5).*new_prob_rat;
    
    [~,pos]    = sort(income_ige(:,2));
    income_ige = income_ige(pos,:);
    income_ige = [income_ige cumsum(income_ige(:,5))];
end

income_ige = [income_ige zeros(size(income_ige,1),1)];

% Quantile thresholds:
nquant              = tot_quant;
for i = 1:nquant
    bot  = (1/nquant)*(i-1);
    top  = (1/nquant)*i;
    if i < nquant
        ind  =  logical((income_ige(:,8)>bot) .* (income_ige(:,8) <= top));
    else
        ind  =  logical((income_ige(:,8)>bot) .* (income_ige(:,8) <= top));
    end
    income_ige(ind,9) = i;
end

% Now do Weighted OLS regression
X = [income_ige(:,5).^0.5.*ones(size(income_ige,1),1) income_ige(:,5).^0.5.*income_ige(:,7)];
Y = income_ige(:,5).^0.5 .* income_ige(:,9);

beta_coeff = (X'*X)\(X'*Y);
rank_rank  = beta_coeff(2);

m_model(count_m)    = rank_rank;
m_title{count_m}    = 'Rank-Rank IGE (100)                 ';
m_mod_id(count_m)   = 16;
count_m             = count_m+1;

X = [income_ige(:,5).^0.5.*ones(size(income_ige,1),1) income_ige(:,5).^0.5.*log(income_ige(:,1))];
Y = income_ige(:,5).^0.5 .* log(income_ige(:,2));

beta_coeff = (X'*X)\(X'*Y);
log_log     = beta_coeff(2);

m_model(count_m)    = log_log;
m_title{count_m}    = 'Log-Log IGE                    ';
m_mod_id(count_m)   = 9;
count_m             = count_m+1;

% 1. IGE Transtion Matrix
nquant = 5;
norig  = tot_quant;
mom_trans = 61;
% tot_p = 0;
% tot_c = zeros(nquant,1);
for ip = 1:nquant
    bot_p  = round((1/nquant)*(ip-1)*norig+1);
    top_p  = round((1/nquant)*ip*norig);
    ind_p  = logical((income_ige(:,7)>=bot_p) .* (income_ige(:,7) <= top_p));
    for ic = 1:nquant
        bot_c  = round((1/nquant)*(ic-1)*norig+1);
        top_c  = round((1/nquant)*ic*norig);
        %         fprintf(' Parent [%3.2f %3.2f] \t child [%3.2f %3.2f]\n',bot_p,top_p,bot_c,top_c)
        ind_c  = logical((income_ige(:,9)>=bot_c) .* (income_ige(:,9) <= top_c));
        ind_pc = logical(ind_p .* ind_c);
        m_model(count_m) = sum(income_ige(ind_pc,5))/sum(income_ige(ind_p,5));
        m_title{count_m} = ['Transition Par Q',num2str(ip),'-','Child Q',num2str(ic),'               ']; %#ok<*AGROW>
        m_mod_id(count_m) = mom_trans;
        mom_trans        = mom_trans + 1;
        count_m          = count_m+1;
        %         tot_c(ip,1) = tot_c(ip,1) + sum(income_ige(ind_pc,5))/sum(income_ige(ind_p,5));
    end
    %     tot_p = tot_p + sum(income_ige(ind_p,5));
end

% 12.1  Theil-L index 1:
% Outcome: child's income
% Types: Parent income decile
mean_outcome        = income_ige(:,2)' * income_ige(:,5);
theil_tot           = (log(mean_outcome) - log(income_ige(:,2)))' * income_ige(:,5);
m_model(count_m)    = theil_tot;
m_title{count_m}    = 'Theil-L index (Total) ';
m_mod_id(count_m)   = 93;
count_m             = count_m+1;

nquant = 10;
norig  = tot_quant;
theil  = zeros(nquant,2);
for ip = 1:nquant
    bot_p  = (1/nquant)*(ip-1)*norig+1;
    top_p  = (1/nquant)*ip*norig;
    ind_p  = logical((income_ige(:,7)>=bot_p) .* (income_ige(:,7) <= top_p));
    
    theil(ip,1) = sum(income_ige(ind_p,5));
    theil(ip,2) = income_ige(ind_p,2)' * income_ige(ind_p,5) ./ theil(ip,1);
end

theil_abs = (log(mean_outcome) - log(theil(:,2)))' * theil(:,1);
m_model(count_m)    = theil_abs;
m_title{count_m}    = 'Theil-L index absolute (parents inc) ';
m_mod_id(count_m)   = 94;
count_m             = count_m+1;

theil_rel = theil_abs/theil_tot;
m_model(count_m)    = theil_rel ;
m_title{count_m}    = 'Theil-L index relative (parents inc) ';
m_mod_id(count_m)   = 95;
count_m             = count_m+1;

end