function [ count_m,m_model,m_title,m_mod_id ] = mom_borr( par,count_m,m_model,m_title,m_mod_id,mu_cs )
S          = par.grids{1,find(par.age == 40,1,'first')};
FE_pos    = par.inc.fe_pos;
EDUC      = par.educ;
INNO_pos  = par.inc.inno_pos;
savings   = nan(length(EDUC)*length(S)*length(par.Je2_pos+3:par.Jr_pos-3),2);

count    = 1;
% Years after educ, until fertility
NH       = length(S);
for j_pos = par.Je2_pos+3:par.Jc_pos
    for educ = 1:length(EDUC)
        S          = par.grids{educ,j_pos}';
        prob_aux   = sum(sum(mu_cs{j_pos}(:,:,educ,:),2),4);
        
        savings(count:count+NH-1,1) = S(:);
        savings(count:count+NH-1,2) = prob_aux(:);
        count = count + NH;
    end
end

% Years with children
NH       = length(S);
for j_pos = par.Jc_pos+1:find(par.age == 40,1,'first')-1;
    for educ = 1:length(EDUC)
        S          = par.grids{educ,j_pos};
        prob_aux   = sum(sum(mu_cs{j_pos}(:,:,educ,:),4),2);
        
        savings(count:count+NH-1,1) = S(:);
        savings(count:count+NH-1,2) = prob_aux(:);
        count = count + NH;
    end
end

% Years after children
NH       = length(S);
for j_pos = find(par.age == 40,1,'first'):par.Jr_pos-3;
    for educ = 1:length(EDUC)
        S          = par.grids{educ,j_pos};
        prob_aux   = sum(sum(mu_cs{j_pos}(:,:,educ,:),4),2);
        
        savings(count:count+NH-1,1) = S(:);
        savings(count:count+NH-1,2) = prob_aux(:);
        count = count + NH;
    end
end
ind  = ~isnan(savings(:,1));
savings = [savings(ind,1) savings(ind,2)];

savings(:,2) = savings(:,2)./sum(savings(:,2)); 

ind_debt = logical((savings(:,1)<0));
share_debt = sum(savings(ind_debt,2));

m_model(count_m)    = share_debt ;
m_title{count_m}    = 'Share of HH with negative assets ';
m_mod_id(count_m)   = 368;
count_m             = count_m+1;

% Savings Average and CV
sav_me      = savings(:,1)'*savings(:,2);
sav_var     = (savings(:,1)'-sav_me).^2*savings(:,2);
sav_cv      = sav_var^0.5/sav_me;

m_model(count_m)    = sav_me/par.p_model_data ;
m_title{count_m}    = 'Average assets of HH (age 24-60) ';
m_mod_id(count_m)   = 369;
count_m             = count_m+1;

m_model(count_m)    = sav_cv ;
m_title{count_m}    = 'CV assets of HH (age 24-60) ';
m_mod_id(count_m)   = 370;
count_m             = count_m+1;

end