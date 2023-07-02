function [ count_m,m_model,m_title,m_mod_id ] = mom_parents_initialh( par,count_m,m_model,m_title,m_mod_id,mu_ige,options,pol )
% mu_ige   = zeros(length(Sp),length(Hp),length(EDUC),length(S0c),length(H0c),length(PSY));
EDUC      = par.educ;
Jp        = par.Jc_pos + 3;
Hp        = par.gridh{1,Jp};
Jc        = par.Je1_pos;
Hc       = par.gridh{1,Jc};


data   = nan(length(EDUC)*length(Hp)*length(Hc),4); %Educ,Hp,H0c,prob
count    = 1;
for educ_p = 1:length(EDUC)
    for iHp = 1:length(Hp)
        hp  = par.gridh{educ_p,Jp}(iHp);
        for iHc = 1:length(Hc)
            hc  = par.gridh{1,Jc}(iHc);
            pr   = sum(sum(mu_ige(:,iHp,educ_p,iHc,:)));
            
            data(count,1) = educ_p;
            data(count,2) = par.w *hp;
            data(count,3) = hc;
            data(count,4) = pr;
            count = count + 1;
        end
    end
end

% Rank Children
ind        = (data(:,4)>0);
data = data(ind,:);
[~,pos]    = sort(data(:,3));
data = data(pos,:);
data = [data cumsum(data(:,4))];

% Quantile thresholds:
nquant              = 5;
for i = 1:nquant-1
    bot  = (1/nquant)*(i-1);
    top  = (1/nquant)*i;  
    
    if data(1,5)<top
        imax = find((data(:,5)<top),1,'last');
    else
        imax = 0;
    end
    
    ind        = (data(:,3) == data(imax+1,3));
    last_ind   = find((data(:,3) < data(imax+1,3)),1,'last');
    if isempty(last_ind) == 1
        prob_tot = 0;
    else
        prob_tot   = data(last_ind,5);
    end
    prob_new   = top - prob_tot;
    old_prob   = sum(data(ind,4));
    new_prob     = old_prob - prob_new;
    new_prob_rat = new_prob/old_prob;
    num        = sum(ind);
    
    data = [data(:,1:4); data(ind,1:4)];
    
    data(end-num+1:end,4) = data(ind,4).*(1-new_prob_rat);
    data(end-num+1:end,3) = data(ind,3).*(1 - 0.00001);
    data(ind,4) = data(ind,4).*new_prob_rat;
   
    [~,pos]    = sort(data(:,3));
    data = data(pos,:);
    data = [data cumsum(data(:,4))]; 
end

data = [data zeros(size(data,1),1)];

% Quantile thresholds:
for i = 1:nquant
    bot  = (1/nquant)*(i-1);
    top  = (1/nquant)*i;
    if i < nquant
        ind  =  logical((data(:,5)>bot) .* (data(:,5) <= top));
    else
        ind  =  logical((data(:,5)>bot) .* (data(:,5) <= top));
    end
    data(ind,6) = i; 
end

data = [data(:,1:6) zeros(size(data,1),3)]; % Add probability in education, cumsum(prob), income group

% Rank Parents Within Education
for educ_p = 1:length(EDUC)
   ind_educ =  (data(:,1) == educ_p);
   % Create probability within education
   data(ind_educ,7) = data(ind_educ,4)./sum(data(ind_educ,4));
   
    [~,pos]    = sort(data(:,2));
    data = data(pos,:);
    data = [data(:,1:7) cumsum(data(:,7)) data(:,9)];
    
    % Quantile thresholds:
    nquant              = 3;
    for i = 1:nquant-1
        bot  = (1/nquant)*(i-1);
        top  = (1/nquant)*i;  

        if data(1,8)<top
            imax = find((data(:,8)<top),1,'last');
        else
            imax = 0;
        end

        ind        = (data(:,2) == data(imax+1,2));
        last_ind   = find((data(:,2) < data(imax+1,2)),1,'last');
        if isempty(last_ind) == 1
            prob_tot = 0;
        else
            prob_tot   = data(last_ind,8);
        end
        prob_new   = top - prob_tot;
        old_prob   = sum(data(ind,7));
        new_prob     = old_prob - prob_new;
        new_prob_rat = new_prob/old_prob;
        num        = sum(ind);

        data = [data(:,1:9); data(ind,1:9)];

        data(end-num+1:end,7) = data(ind,7).*(1-new_prob_rat);
        data(end-num+1:end,2) = data(ind,2).*(1 - 0.00001);
        data(ind,7) = data(ind,7).*new_prob_rat;

        [~,pos]    = sort(data(:,2));
        data = data(pos,:);
        data = [data(:,1:7) cumsum(data(:,7)) data(:,9)]; 
    end

    % Quantile thresholds:
    for i = 1:nquant
        bot  = (1/nquant)*(i-1);
        top  = (1/nquant)*i;
        if i < nquant
            ind  =  logical((data(:,1)==educ_p).*(data(:,8)>bot) .* (data(:,8) <= top));
        else
            ind  =  logical((data(:,1)==educ_p).*(data(:,8)>bot) .* (data(:,8) <= top));
        end
        data(ind,9) = i; 
    end
    
    % Reset probabilities for next education group
    data(:,7:8) = 0;
end


% Organize info
mom_trans = 486;
trans_mat = nan(length(EDUC)*3*5,4); %Educ,Hp,H0c,prob
count = 1;
for educ_p = 1:length(EDUC)
    for igroup_p = 1:3
        ind_p    = logical((data(:,1)==educ_p).*(data(:,9)==igroup_p));
        for igroup_c = 1:5
            ind_c    = logical((data(:,1)==educ_p).*(data(:,9)==igroup_p).*(data(:,6)==igroup_c));
            trans_mat(count,1) = educ_p;
            trans_mat(count,2) = igroup_p;
            trans_mat(count,3) = igroup_c;
            trans_mat(count,4) = sum(data(ind_c,4))/sum(data(ind_p,4));
            
            m_model(count_m) = trans_mat(count,4);
            m_title{count_m} = ['Transition Par Educ ',num2str(educ_p),'-Inc Group ',num2str(igroup_p),'-','Child Q',num2str(igroup_c),'               ']; %#ok<*AGROW>
            m_mod_id(count_m) = mom_trans;
            mom_trans         = mom_trans + 1;
            
            count_m          = count_m+1;
            count = count + 1;
        end
    end
end


end