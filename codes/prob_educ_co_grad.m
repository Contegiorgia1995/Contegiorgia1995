function [V,C,Sp] = prob_educ_co_grad(par,options,Vp,Cp)
% Solve household problem when it is College Graduate, for age 12 and 16
% College graduate:
educ      = 3;
PSY       = par.psy_val_col;
r_sav     = par.r_sav;
beta      = par.beta;
gammac    = par.gammac;
lambda    = par.lambda;

V         = cell(1,3);
C         = cell(1,3);
Sp        = cell(1,3);

Vpaux     = squeeze(Vp(:,:,educ,:));
Cpaux     = squeeze(Cp(:,:,educ,:));

FE_pos     = par.inc.fe_pos;
AGE_PROF   = par.inc.age_prof;
INNO_pos   = par.inc.inno_pos;
FE        = par.inc.fe;


% Expectation about innovation
j_pos      = par.Je2_pos+2;
exp_prob   = squeeze(par.inc.z_prob{educ,1}(j_pos,1,:))';
                         

%% Age 14-20
for j_pos = par.Je2_pos+1:-1:par.Je1_pos
    S          = par.grids{educ,j_pos};
    
    V2         = zeros(length(S),length(FE_pos));
    C2         = zeros(length(S),length(FE_pos));
    Sp2        = zeros(length(S),length(FE_pos));
    boundgrid  = zeros(length(FE_pos));
    Spgrid     = par.grids{educ,j_pos+1};
    if j_pos   == par.Je2_pos+1
        r_debt_today        = 0;
        r_debt_tomorrow     = par.r_debt;
        col_fact_1          = (par.col_fact .* (Spgrid(1)<0) + 1 .* (Spgrid(1)>=0));
        par_temp            = par;
        par_temp.r_debt     = 0;
        FE                  = par.inc.fe{2,1}; %Skills before getting educated
    elseif j_pos   == par.Je2_pos
        r_debt_today        = par.r_debt;
        r_debt_tomorrow     = 0;
        col_fact_1          = 1;
        par_temp            = par;
        par_temp.col_fact   = 1;
        FE         = par.inc.fe{2,1}; %Skills before getting educated
    else
        r_debt_today        = par.r_debt;
        r_debt_tomorrow     = par.r_debt;
        FE                  = par.inc.fe{1,1}; %Skills before getting educated
    end
    r          = (r_sav .* (Spgrid>=0) + r_debt_tomorrow .* (Spgrid<0))';
    rs         = (r_sav .* (S>=0)      + r_debt_today .* (S<0));
    
    if j_pos > par.Je2_pos
        inno3 = repmat(INNO_pos,length(S),1);
        [S2,INNO2] = ndgrid(Spgrid,INNO_pos);
    end
    
    for ife = 1:length(FE_pos)
        if j_pos > par.Je2_pos
            splVp   = griddedInterpolant(S2,INNO2,squeeze(Vpaux(:,ife,:)));
        else
            splVp   = griddedInterpolant(Spgrid,squeeze(Vpaux(:,ife)));
        end
        
        h       = exp(AGE_PROF{educ,1}(j_pos)) * exp(FE(ife));
        
        Cpp           = squeeze(Cpaux(:,ife,:));
        if min(Cpp(:))<0 && abs(min(Cpp(:)))<1e-6
            Cpp  = max(Cpp,0);
        end
        ucp = beta*(1+r).* (Cpp.^(-gammac)* exp_prob');
        
        if j_pos >= par.Je2_pos
            dispinc = par.w_college*h*(1-lambda) - par.pe2 + par.init_trans(educ,j_pos);
        else
            dispinc = -par.pe1 + par.init_trans(educ,j_pos);
        end
        feasible = (dispinc + (1+rs).*S - Spgrid(1)*(1/col_fact_1)>0);
        not_feasible = logical(1-feasible);
        if j_pos >= par.Je2_pos
            [C2(feasible,ife),Sp2(feasible,ife),boundgrid(ife)] = ...
                GEGM_college(par_temp,ucp',dispinc,Spgrid,S(feasible),splVp,exp_prob);
            C2(not_feasible,ife)  = 0;
            Sp2(not_feasible,ife) = Spgrid(1);
        else
            [C2(feasible,ife),Sp2(feasible,ife),boundgrid(ife)] = ...
                GEGM(par,ucp',dispinc,Spgrid,S(feasible),splVp,exp_prob);
            C2(not_feasible,ife)  = 0;
            Sp2(not_feasible,ife) = Spgrid(1);
        end
        
        if j_pos > par.Je2_pos
            Sp3         = repmat(Sp2(:,ife),1,length(INNO_pos));
            vp          = splVp(Sp3,inno3)*exp_prob';
        else
            vp          = splVp(Sp2(:,ife));
        end
        
        V2(feasible,ife) = C2(feasible,ife).^(1-gammac)/(1-gammac) ...
            + beta*vp(feasible);
        
        V2(not_feasible,ife) = -(10^(5/gammac)).^(1-gammac)/(1-gammac);
        
    end
    
    if j_pos > par.Je1_pos
        V{1,j_pos - (par.Je1_pos-1)} = V2;
        C{1,j_pos - (par.Je1_pos-1)} = C2;
        Sp{1,j_pos - (par.Je1_pos-1)} = Sp2;
    else
        V{1,j_pos - (par.Je1_pos-1)} = repmat(V2,1,1,length(PSY))-repmat(reshape(PSY,1,1,length(PSY)),size(V2,1),size(V2,2));
        C{1,j_pos - (par.Je1_pos-1)} =  repmat(C2,1,1,length(PSY));
        Sp{1,j_pos - (par.Je1_pos-1)} =  repmat(Sp2,1,1,length(PSY));
    end
    
    Vpaux   = V2;
    Cpaux   = C2;
    exp_prob = 1;
end


end