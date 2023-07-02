function [ count_m,m_model,m_title,m_mod_id ] = mom_fert_elast_labinc( par,count_m,m_model,m_title,m_mod_id,mu,pol )
FE_pos    = par.inc.fe_pos;
EDUC      = par.educ;
INNO_pos  = par.inc.inno_pos;

% Income groups of parent
j_pos               = par.Jc_pos;
NH                  = length(FE_pos)*length(INNO_pos);
inc_groups          = nan(length(EDUC)*NH,3);

pos    = 1;
for educ = 1:length(EDUC)
    INNO       = squeeze(repmat(reshape(par.inc.z_val{educ,1}(j_pos,:),1,length(INNO_pos)),length(FE_pos),1));
    FE         = squeeze(repmat(reshape(par.inc.fe{educ,1},length(FE_pos),1),1,length(INNO_pos)));
    AGE_PROF   = par.inc.age_prof{educ,1}(j_pos);
    
    inc                         = par.w * exp(AGE_PROF) .* exp(FE(:)) .* exp(INNO(:)) ;
    mu_aux                      = mu{j_pos}(:,:,educ,:);
    pr                          = sum(mu_aux,1);
    fert                        = par.fam_size*par.N(pol.Np(:,:,educ,:));
    avg_fert                    = squeeze(sum(fert .* mu_aux,1) ./ pr);
%     ch_cost                     = ChildCost(par,inc,avg_fert(:),0,options);
    inc_groups(pos:pos+NH-1,1)  = inc;
    inc_groups(pos:pos+NH-1,2)  = avg_fert(:);    
    inc_groups(pos:pos+NH-1,3)  = pr(:);    
    pos                         = pos + NH;
end


ind        = (inc_groups(:,3)>0);
inc_groups = inc_groups(ind,:);
[~,pos]    = sort(inc_groups(:,1));
inc_groups = inc_groups(pos,:);
inc_groups = [inc_groups cumsum(inc_groups(:,3))];

% Quantile thresholds:
nquant              = 10;
for i = 1:nquant-1
    bot  = (1/nquant)*(i-1);
    top  = (1/nquant)*i;  
    
    if inc_groups(1,4)<top
        imax = find((inc_groups(:,4)<top),1,'last');
        prob_tot = inc_groups(imax,4);
        prob_new = top - prob_tot;
    else
        imax = 0;
        prob_new = top;
    end
    
    inc_groups = [inc_groups(:,1:3); inc_groups(imax+1,1:3)];
    inc_groups(imax+1,3) = inc_groups(imax+1,3) - prob_new;
    inc_groups(end,1) = inc_groups(end,1)*(1 - 0.00001);
    inc_groups(end,3) = prob_new;

    [~,pos]    = sort(inc_groups(:,1));
    inc_groups = inc_groups(pos,:);
    inc_groups = [inc_groups cumsum(inc_groups(:,3))]; 
end

inc_groups = [inc_groups zeros(size(inc_groups,1),1)];

inc_quant           = zeros(nquant,1);
fert_quant          = zeros(nquant,1);
mom_fert            = 26;
for i = 1:nquant
    bot  = (1/nquant)*(i-1);
    top  = (1/nquant)*i;
    if i < nquant
        ind  =  logical((inc_groups(:,4)>bot) .* (inc_groups(:,4) <= top));
    else
        ind  =  logical((inc_groups(:,4)>bot) .* (inc_groups(:,4) <= top));
    end
    inc_groups(ind,5) = i; 
    inc_quant(i)     =  inc_groups(ind,1)'*inc_groups(ind,3)/sum(inc_groups(ind,3));
    fert_quant(i)    =  inc_groups(ind,2)'*inc_groups(ind,3)/sum(inc_groups(ind,3));
    m_model(count_m) = fert_quant(i);
    m_title{count_m} = ['Fertility quintile ',num2str(i),'               '];
    m_mod_id(count_m)   = mom_fert;
    count_m             = count_m+1;
    mom_fert            = mom_fert +1;
end


% Keep only obs with fert_quant >0 
ind                 = (fert_quant>1e-5);
if sum(ind)>2
    Y               = log(fert_quant(ind));
    X               = [ones(size(fert_quant(ind))), log(inc_quant(ind))];
    b               = regress(Y,X);
    m_model(count_m)= b(2);
else
    m_model(count_m)= 0;
end
m_title{count_m}    = 'Fertility elasticity                  ';
m_mod_id(count_m)   = 7;
count_m             = count_m+1;

% Quantile thresholds:
inc_groups          = inc_groups(:,1:4);
nquant              = 5;
for i = 1:nquant-1
    bot  = (1/nquant)*(i-1);
    top  = (1/nquant)*i;  
    
    if inc_groups(1,4)<top
        imax = find((inc_groups(:,4)<top),1,'last');
        prob_tot = inc_groups(imax,4);
        prob_new = top - prob_tot;
    else
        imax = 0;
        prob_new = top;
    end
    
    inc_groups = [inc_groups(:,1:3); inc_groups(imax+1,1:3)];
    inc_groups(imax+1,3) = inc_groups(imax+1,3) - prob_new;
    inc_groups(end,1) = inc_groups(end,1)*(1 - 0.00001);
    inc_groups(end,3) = prob_new;

    [~,pos]    = sort(inc_groups(:,1));
    inc_groups = inc_groups(pos,:);
    inc_groups = [inc_groups cumsum(inc_groups(:,3))]; 
end

inc_groups = [inc_groups zeros(size(inc_groups,1),1)];

inc_quant           = zeros(nquant,1);
fert_quant          = zeros(nquant,1);
mom_fert            = 114;
for i = 1:nquant
    bot  = (1/nquant)*(i-1);
    top  = (1/nquant)*i;
    if i < nquant
        ind  =  logical((inc_groups(:,4)>bot) .* (inc_groups(:,4) <= top));
    else
        ind  =  logical((inc_groups(:,4)>bot) .* (inc_groups(:,4) <= top));
    end
    inc_groups(ind,5) = i; 
    inc_quant(i)     =  inc_groups(ind,1)'*inc_groups(ind,3)/sum(inc_groups(ind,3));
    fert_quant(i)    =  inc_groups(ind,2)'*inc_groups(ind,3)/sum(inc_groups(ind,3));
    m_model(count_m) = fert_quant(i);
    m_title{count_m} = ['Fertility quintile ',num2str(i),'               '];
    m_mod_id(count_m)   = mom_fert;
    count_m             = count_m+1;
    mom_fert            = mom_fert +1;
end

end