function [ count_m,m_model,m_title,m_mod_id ] = mom_agg( par,count_m,m_model,m_title,m_mod_id,mu_cs,pol,Agg_K,options )
FE_pos    = par.inc.fe_pos;
EDUC      = par.educ;
INNO_pos  = par.inc.inno_pos;

Agg_C = 0;
Agg_S = 0;
Agg_H = 0;
Agg_R = 0;

r_sav     = par.r_sav;
r_debt    = par.r_debt;

N         = 1:length(par.N);

% 22.1 C and S
% Education
start_L  = [par.Je1_pos par.Je2_pos par.Je2_pos + 2];
for j_pos = par.Je1_pos:par.Je2_pos+1;
    for educ = 1:length(par.educ)
        S     = par.grids{educ,j_pos};
%         if j_pos >= start_L(educ) & j_pos
        if j_pos == par.Je1_pos && educ == 1
            prob        = squeeze(sum(mu_cs{j_pos}(:,:,:,educ,:),3));
        elseif j_pos == par.Je1_pos && educ > 1
            prob        = squeeze(sum(mu_cs{j_pos}(:,:,:,educ,:),5));
        elseif j_pos > par.Je1_pos && educ > 2
%             prob        = squeeze(mu_cs{j_pos}(:,:,:,educ));
            prob        = squeeze(sum(mu_cs{j_pos}(:,:,educ,:),4));
        else
            prob        = squeeze(mu_cs{j_pos}(:,:,educ,:));
        end
        
        cons        = pol.Ce{educ,j_pos-(par.Je1_pos-1)}(:);        
        Agg_C_aux   = cons' * prob(:);
        Agg_C       = Agg_C + Agg_C_aux;

		r           = (r_sav .* (S>=0) + r_debt .* (S<0));
        sav         = r .* S;
        prob_s      = sum(prob(:,:),2);
        Agg_S_aux   = sav * prob_s;
        Agg_S       = Agg_S + Agg_S_aux;

    end
end
% Next periods
for j_pos = par.Je2_pos + 2:par.Jd_pos - 1
    for educ = 1:length(EDUC)
        S           = par.grids{educ,j_pos};
        prob        = squeeze(mu_cs{j_pos}(:,:,educ,:));
        cons        = squeeze(pol.C{j_pos}(:,:,educ,:));
        Agg_C_aux   = cons(:)' * prob(:);
        Agg_C       = Agg_C + Agg_C_aux;

        r           = (r_sav .* (S>=0) + r_debt .* (S<0));
        sav         = r .* S;
        prob_s      = sum(prob(:,:),2);
        Agg_S_aux   = sav * prob_s;
        Agg_S       = Agg_S + Agg_S_aux;
    end
end
% Children's consumption
for j_pos = par.Jc_pos:par.Jc_pos+par.Je1_pos-2
    prob        = mu_cs{j_pos}(:);
    cons        = pol.Ck{j_pos}(:);
    Agg_C_aux   = cons' * prob;
    Agg_C       = Agg_C + Agg_C_aux;
end
% 22.2 Aggregate H
educ = 1;
for j_pos = par.Je1_pos:par.Je2_pos+1
    if j_pos == par.Je1_pos
        prob        = squeeze(sum(mu_cs{j_pos}(:,:,:,educ,:),3));
    else
        prob        = squeeze(mu_cs{j_pos}(:,:,educ,:));
    end
    prob        = sum(prob,1);
     
    INNO       = squeeze(repmat(reshape(par.inc.z_val{educ,1}(j_pos,:),1,length(INNO_pos)),length(FE_pos),1));
    FE         = squeeze(repmat(reshape(par.inc.fe{educ,1},length(FE_pos),1),1,length(INNO_pos)));
    AGE_PROF   = par.inc.age_prof{educ,1}(j_pos);

    inc         = (1-par.lambda) * par.w * exp(AGE_PROF) .* exp(FE(:)) .* exp(INNO(:));
    Agg_H_aux   = inc' * prob(:);
    Agg_H       = Agg_H + Agg_H_aux;
end
educ = 2;
for j_pos = par.Je2_pos:par.Je2_pos+1
    if j_pos == par.Je1_pos
        prob        = squeeze(sum(sum(mu_cs{j_pos}(:,:,:,educ,:),3),5));
        INNO       = 0;
    else
        prob        = squeeze(mu_cs{j_pos}(:,:,educ,:));
        INNO       = squeeze(repmat(reshape(par.inc.z_val{educ,1}(j_pos,:),1,length(INNO_pos)),length(FE_pos),1));
    end
    prob        = sum(prob,1);

    FE         = squeeze(repmat(reshape(par.inc.fe{educ,1},length(FE_pos),1),1,length(INNO_pos)));
    AGE_PROF   = par.inc.age_prof{educ,1}(j_pos);

     
     inc         = (1-par.lambda) * par.w * exp(AGE_PROF) .* exp(FE(:)) .* exp(INNO(:));
     Agg_H_aux   = inc' * prob(:);
     Agg_H       = Agg_H + Agg_H_aux;
end
educ = 3;
for j_pos = par.Je2_pos:par.Je2_pos+1
    prob        = squeeze(sum(sum(mu_cs{j_pos}(:,:,:,educ,:),5),3));
    prob        = sum(prob,1);

    FE         = squeeze(reshape(par.inc.fe{educ-1,1},length(FE_pos),1));
    AGE_PROF   = par.inc.age_prof{educ,1}(j_pos);

     
     inc         = (1-par.lambda) * par.w_college * exp(AGE_PROF) .* exp(FE(:));
     Agg_H_aux   = inc' * prob(:);
     Agg_H       = Agg_H + Agg_H_aux;
end
% Next periods until fertility
for j_pos = par.Je2_pos + 2:par.Jc_pos - 1
    for educ = 1:length(par.educ)
        INNO       = squeeze(repmat(reshape(par.inc.z_val{educ,1}(j_pos,:),1,length(INNO_pos)),length(FE_pos),1));
        FE         = squeeze(repmat(reshape(par.inc.fe{educ,1},length(FE_pos),1),1,length(INNO_pos)));
        AGE_PROF   = par.inc.age_prof{educ,1}(j_pos);

        inc         = (1-par.lambda) * par.w * exp(AGE_PROF) .* exp(FE(:)) .* exp(INNO(:));
        prob        = sum(sum(mu_cs{j_pos}(:,:,educ,:,:),1),5);
        
        Agg_H_aux   = inc' * prob(:);
        Agg_H       = Agg_H + Agg_H_aux;
    end
end
% Next periods until fertility
for j_pos = par.Jc_pos:par.Jc_pos
    for educ = 1:length(par.educ)
        AGE_PROF   = par.inc.age_prof{educ,1}(j_pos);     
        for ife = 1:length(FE_pos)     
            for inno = 1:length(INNO_pos)
                H           = exp(AGE_PROF) .* exp(par.inc.fe{educ,1}(ife)) .* exp(par.inc.z_val{educ,1}(j_pos,inno));

                Nchoice     = pol.Np(:,ife,educ,inno);
                for in      = 1:length(N)
                    ind     = logical(Nchoice == in);

                    inc         = (1-par.lambda) * par.w * H;
                    if in == 1
                        ChildCost1 = 0;
                    else
                        ChildCost1 = ChildCost(par,inc,par.N(in),0,options);
                    end

                    inc         = (1-par.lambda) * par.w * H - ChildCost1;

                    prob        = sum(mu_cs{j_pos}(ind,ife,educ,inno));
                    Agg_H_aux   = inc * prob;
                    Agg_H       = Agg_H + Agg_H_aux;
                end
            end
        end
    end
end
% Next periods with children
for j_pos = par.Jc_pos+1:par.Jc_pos+par.Je1_pos-2
    for educ = 1:length(par.educ)
        INNO       = squeeze(repmat(reshape(par.inc.z_val{educ,1}(j_pos,:),1,length(INNO_pos)),length(FE_pos),1));
        FE         = squeeze(repmat(reshape(par.inc.fe{educ,1},length(FE_pos),1),1,length(INNO_pos)));
        AGE_PROF   = par.inc.age_prof{educ,1}(j_pos);

        for in = 1:length(N)
            inc         = (1-par.lambda) * par.w * exp(AGE_PROF) .* exp(FE(:)) .* exp(INNO(:));
            if in == 1
                ChildCost1 = 0;
            else
                ChildCost1 = ChildCost(par,inc,par.N(in),0,options);
            end
            inc         = inc - ChildCost1;
            prob        = sum(mu_cs{j_pos}(:,:,educ,:,in),1);
           
            Agg_H_aux   = inc' * prob(:);
            Agg_H       = Agg_H + Agg_H_aux;
        end
    end
end
% Next periods after children leave the house
for j_pos = par.Jc_pos+par.Je1_pos-1:par.Jr_pos - 1
    for educ = 1:length(par.educ)
        INNO       = squeeze(repmat(reshape(par.inc.z_val{educ,1}(j_pos,:),1,length(INNO_pos)),length(FE_pos),1));
        FE         = squeeze(repmat(reshape(par.inc.fe{educ,1},length(FE_pos),1),1,length(INNO_pos)));
        AGE_PROF   = par.inc.age_prof{educ,1}(j_pos);
        
        inc         = (1-par.lambda) * par.w * exp(AGE_PROF) .* exp(FE(:)) .* exp(INNO(:));
        prob        = sum(mu_cs{j_pos}(:,:,educ,:),1);

        Agg_H_aux   = inc' * prob(:);
        Agg_H       = Agg_H + Agg_H_aux;
    end
end
Agg_Tax   = par.lambda * Agg_H/(1-par.lambda);
% Income from retirement
for j_pos = par.Jr_pos :par.Jd_pos - 1
    for educ = 1:length(par.educ)
        inc         = ret_rep(par,FE_pos,educ,options);
        prob        = sum(mu_cs{j_pos}(:,:,educ,:),1);

        Agg_R_aux   = inc * prob(:);
        Agg_R       = Agg_R + Agg_R_aux;
    end
end

m_model(count_m)    = Agg_R - Agg_Tax;
m_title{count_m}    = 'Ret. Benefits - Tax Revenue          ';
m_mod_id(count_m)   = 091;
count_m             = count_m + 1;

m_model(count_m)    = Agg_R/Agg_Tax;
m_title{count_m}    = 'Ret. Benefits/Tax Revenue          ';
m_mod_id(count_m)   = 092;
count_m             = count_m + 1;

m_model(count_m)    = Agg_C/par.p_model_data;
m_title{count_m}    = 'GDP pc Cons              ';
m_mod_id(count_m)   = 102;
count_m             = count_m+1;
m_model(count_m)    = (Agg_H + Agg_S + Agg_R)/par.p_model_data;
m_title{count_m}    = 'GDP pc Inc              ';
m_mod_id(count_m)   = 103;
count_m             = count_m+1;

m_model(count_m)    = Agg_K/(Agg_H + Agg_S + Agg_R);
m_title{count_m}    = 'K/Y              ';
m_mod_id(count_m)   = 104;
count_m             = count_m+1;


end