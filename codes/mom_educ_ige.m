function [ count_m,m_model,m_title,m_mod_id ] = mom_educ_ige( par,count_m,m_model,m_title,m_mod_id,mu_educ_ige,options )
EDUC      = par.educ;


educ_years = [8 12 16];
educ_ige   = nan(length(EDUC)*length(EDUC),5);
educ_trans = nan(length(EDUC),length(EDUC));
count      = 1;
mom_trans  = 212;
for educ_p = 1:length(EDUC)
    pr_p   = mu_educ_ige(educ_p,:);
    pr_p   = sum(pr_p(:));
    for educ_c = 1:length(EDUC)
        educ_ige(count,1) = educ_p;
        educ_ige(count,4) = educ_years(educ_p);
        educ_ige(count,2) = educ_c;
        educ_ige(count,5) = educ_years(educ_c);
        pr   = mu_educ_ige(educ_p,educ_c);
        pr   = pr(:);
        educ_ige(count,3) = sum(pr);
        
        m_model(count_m) = educ_ige(count,3)/pr_p;
        m_title{count_m} = ['Transition Par Educ ',num2str(educ_p),'-','Child Educ',num2str(educ_c)]; %#ok<*AGROW>
        m_mod_id(count_m)  = mom_trans;
        if pr_p>0
            educ_trans(educ_p,educ_c) =  m_model(count_m);
        end
        mom_trans        = mom_trans + 1;
        count_m          = count_m+1;
        count            = count + 1;
        
    end
end


% Log Log regression

X = [educ_ige(:,3).^0.5.*ones(size(educ_ige,1),1) educ_ige(:,3).^0.5.*log(educ_ige(:,4))];
Y = educ_ige(:,3).^0.5 .* log(educ_ige(:,5));

beta_coeff = (X'*X)\(X'*Y);
log_log    = beta_coeff(2);
m_model(count_m)    = log_log ;
m_title{count_m}    = 'IGE Education Persistence: Log-Log ';
m_mod_id(count_m)   = 12;
count_m          = count_m+1;

clear educ_ige X Y

educ_missing = sum(isnan(educ_trans(:,1)));
if educ_missing == 1
    keep_educ       = ~isnan(educ_trans(:,1));
    educ_trans      = educ_trans(keep_educ,keep_educ);
    EDUC            = EDUC(keep_educ);
end


% Second biggest eigenvalue

if sum(isnan(educ_trans(:)))==0
    eig_educ = eig(educ_trans);
    eig_educ = eig_educ(end-1);
    m_model(count_m)    = eig_educ ;
    m_title{count_m}    = 'IGE Education Persistence: 2nd Eig ';
    m_mod_id(count_m)   = 13;
    count_m          = count_m+1;
    
    aux                 = (length(EDUC) - trace(educ_trans))/(length(EDUC)-1);
    m_model(count_m)    = aux ;
    m_title{count_m}    = 'IGE Education Persistence: Trace ';
    m_mod_id(count_m)   = 14;
    count_m          = count_m+1;
    
    aux                 = 1-abs(det(educ_trans))^(1/(length(EDUC)-1));
    m_model(count_m)    = aux ;
    m_title{count_m}    = 'IGE Education Persistence: det ';
    m_mod_id(count_m)   = 15;
    count_m          = count_m+1;
    
else
    m_model(count_m)    = NaN ;
    m_title{count_m}    = 'IGE Education Persistence: 2nd Eig ';
    m_mod_id(count_m)   = 13;
    count_m          = count_m+1;
    
    m_model(count_m)    = NaN ;
    m_title{count_m}    = 'IGE Education Persistence: Trace ';
    m_mod_id(count_m)   = 14;
    count_m          = count_m+1;
    
    m_model(count_m)    = NaN ;
    m_title{count_m}    = 'IGE Education Persistence: det ';
    m_mod_id(count_m)   = 15;
    count_m          = count_m+1;
end

clear educ_trans



end