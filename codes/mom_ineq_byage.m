function [ count_m,m_model,m_title,m_mod_id ] = mom_ineq_byage( par,count_m,m_model,m_title,m_mod_id,mu_cs)
FE_pos    = par.inc.fe_pos;
EDUC      = par.educ;
INNO_pos  = par.inc.inno_pos;

mom_gini  = 721;
mom_tb    = 744;
for j_pos = par.Je2_pos+1:par.Jr_pos-1
    EDUC        = par.educ;
    NH          = length(FE_pos)*length(INNO_pos);
    income   = zeros(length(EDUC)*NH,2);
    count    = 1;
    for educ = 1:length(EDUC)
        INNO       = squeeze(repmat(reshape(par.inc.z_val{educ,1}(j_pos,:),1,length(INNO_pos)),length(FE_pos),1));
        FE         = squeeze(repmat(reshape(par.inc.fe{educ,1},length(FE_pos),1),1,length(INNO_pos)));
        AGE_PROF   = par.inc.age_prof{educ,1}(j_pos);
        
        H = exp(AGE_PROF) .* exp(FE(:)) .* exp(INNO(:));
        aux          = mu_cs{j_pos}(:,:,educ,:,:);
        if j_pos == par.Je2_pos+1 && educ == 3
            inc_val           = par.w_college	* H./ par.time_period;
        else
            inc_val           = par.w	* H./ par.time_period;
        end
        pr   = squeeze(sum(sum(aux,1),5));
        income(count:count+NH-1,1) = inc_val(:);
        income(count:count+NH-1,2) = pr(:);
        count = count + NH;
    end
    
    
    income(:,2) = income(:,2)./sum(income(:,2));
    [income(:,1),ord] = sort(income(:,1));
    income(:,2) = income(ord,2);
    income      = [income cumsum(income(:,2).*income(:,1))]; % Si = sum_j=1^i [prob(yj)*yj], S0 = 0;
    Gini        = 1-sum(income(:,2).*(income(:,3) + [0; income(1:end-1,3)]))/income(end,3);
    m_model(count_m)    = Gini ;
    m_title{count_m} = ['Gini (Ages ',num2str(par.age(j_pos)),'-',num2str(par.age(j_pos+1)-1),')']; %#ok<*AGROW>
    m_mod_id(count_m) = mom_gini;
    mom_gini            = mom_gini + 1;
    count_m          = count_m+1;
    
    % Top-Bottom
    income     = income(:,1:2);
    ind        = (income(:,2)>0);
    income     = income(ind,:);
    [~,pos]    = sort(income(:,1));
    income     = income(pos,:);
    income     = [income cumsum(income(:,2))];
    
    % Quantile thresholds:
    nquant              = 20;
    for i = 1:nquant-1
        bot  = (1/nquant)*(i-1);
        top  = (1/nquant)*i;
        
        if income(1,3)<top
            imax = find((income(:,3)<top),1,'last');
            prob_tot = income(imax,3);
            prob_new = top - prob_tot;
        else
            imax = 0;
            prob_new = top;
        end
        
        income = [income(:,1:2); income(imax+1,1:2)];
        income(imax+1,2) = income(imax+1,2) - prob_new;
        income(end,1)    = income(end,1).*(1 - 0.00001);
        income(end,2)    = prob_new;
        
        [~,pos]    = sort(income(:,1));
        income = income(pos,:);
        income = [income cumsum(income(:,2))];
    end
    
    income = [income zeros(size(income,1),1)];
    ind_b      = logical((income(:,3)<=0.2).*(income(:,3)>0.05));
    ind_t      = logical((income(:,3)<0.95).*(income(:,3)>=0.8));
    
    income_bot = income(ind_b,1:2);
    prob_bot   = sum(income_bot(:,2));
    me_inc_bot = income_bot(:,2)'*income_bot(:,1)/prob_bot;
    
    income_top = income(ind_t,1:2);
    prob_top   = sum(income_top(:,2));
    me_inc_top = income_top(:,2)'*income_top(:,1)/prob_top;
    
    top_bot    = me_inc_top/me_inc_bot;
    m_model(count_m)    = top_bot ;
    m_title{count_m} = ['Top-Bottom (Ages ',num2str(par.age(j_pos)),'-',num2str(par.age(j_pos+1)-1),')']; %#ok<*AGROW>
    m_mod_id(count_m) = mom_tb;
    mom_tb           = mom_tb + 1;
    count_m          = count_m+1;
end

end