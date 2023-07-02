function [ count_m,m_model,m_title,m_mod_id ] = mom_rank_rank_ige( par,count_m,m_model,m_title,m_mod_id,IGE )
par_age_min = [find(par.age == 42,1,'first') find(par.age == 40,1,'first') find(par.age == 28,1,'first') find(par.age == 28,1,'first') find(par.age == 24,1,'first')];
par_age_max = [find(par.age == 42,1,'first') find(par.age == 44,1,'first') find(par.age == 42,1,'first') find(par.age == 42,1,'first') find(par.age == 62,1,'first')];

ch_age_min = [find(par.age == 28,1,'first') find(par.age == 30,1,'first') find(par.age == 28,1,'first') find(par.age == 28,1,'first') find(par.age == 24,1,'first')];
ch_age_max = [find(par.age == 28,1,'first') find(par.age == 30,1,'first') find(par.age == 28,1,'first') find(par.age == 42,1,'first') find(par.age == 62,1,'first')];

mom_rank_rank = 119;
for igr = 1:length(par_age_min)
    j_p_pos_min = par_age_min(igr);
    j_p_pos_max = par_age_max(igr);
    j_c_pos_min = ch_age_min(igr);
    j_c_pos_max = ch_age_max(igr);

    IncPar      = zeros(size(IGE.lab{j_p_pos_min,1},1),1);
    IncCh       = zeros(size(IGE.lab{j_p_pos_min,1},1),1);

    for j_p_pos = j_p_pos_min:j_p_pos_max
        IncPar      = IncPar + IGE.lab{j_p_pos,1}+IGE.sav{j_p_pos,1}; %
    end
    for j_c_pos = j_c_pos_min:j_c_pos_max
        IncCh      = IncCh + IGE.lab{j_c_pos,2}+IGE.sav{j_c_pos,2}; %
    end
    ind         = ~isnan(IncCh);
    IncPar      = IncPar(ind);
    IncCh       = IncCh(ind);
    
    tot_quant   = 50;
    quants      = 100*(1:tot_quant)/tot_quant;

    PrctilePar  = prctile(IncPar,quants);
    LeftPar     = [0 PrctilePar(1:end-1)];
    RightPar    = PrctilePar;

    PrctileCh   = prctile(IncCh,quants);
    LeftCh     = [0 PrctileCh(1:end-1)];
    RightCh    = PrctileCh;

    RankPar  = nan(length(IncPar),1);
    RankCh   = nan(length(IncPar),1);
    for in = 1:length(IncPar)
        RankPar(in)     = find((IncPar(in)> LeftPar).*(IncPar(in)<=RightPar));
        RankCh(in)      = find((IncCh(in)> LeftCh).*(IncCh(in)<=RightCh));
    end

    X = [ones(size(RankPar,1),1) RankPar];
    Y = RankCh;

    beta_coeff = (X'*X)\(X'*Y);
    rank_rank  = beta_coeff(2);

    m_model(count_m)    = rank_rank;
    m_title{count_m}    = 'Rank-Rank IGE (100)                 ';
    m_mod_id(count_m)   = mom_rank_rank;
    count_m             = count_m+1;
    mom_rank_rank       = mom_rank_rank + 1;
end

%%
j_pos_min = find(par.age == 16,1,'first');
j_pos_max = find(par.age == 78,1,'first');

LE_nosav     = zeros(size(IGE.lab{j_pos_min,1},1),1);
LE_wsav      = zeros(size(IGE.lab{j_pos_min,1},1),1);
for j_p_pos = j_pos_min:j_pos_max
%     LE_nosav     = LE_nosav + IGE.lab{j_p_pos,2};
%     LE_wsav      = LE_wsav + IGE.lab{j_p_pos,2}+IGE.sav{j_p_pos,2};
    LE_nosav     = LE_nosav + (1/(1+par.r_sav))^(j_p_pos-j_pos_min)*IGE.lab{j_p_pos,2};
    LE_wsav      = LE_wsav  + (1/(1+par.r_sav))^(j_p_pos-j_pos_min)*(IGE.lab{j_p_pos,2}+IGE.sav{j_p_pos,2});
end
ind         = ~isnan(LE_nosav);

LE_nosav    = log(LE_nosav(ind));
LE_wsav     = log(LE_wsav(ind));

m_model(count_m)    = var(LE_nosav) ;
m_title{count_m}    = 'Var of Lifetime Labor Income ';
m_mod_id(count_m)   = 124;
count_m             = count_m+1;

m_model(count_m)    = var(LE_wsav) ;
m_title{count_m}    = 'Var of Lifetime Labor Income and Savings Returns';
m_mod_id(count_m)   = 125;
count_m             = count_m+1;

end