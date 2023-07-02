function [V,C,Sp] = prob_workT (par,options,Vp,Cp,j_pos)
% Solve household problem last period of working

r_sav     = par.r_sav;
r_debt    = par.r_debt;
beta      = par.beta;
gammac    = par.gammac;
w         = par.w;
lambda    = par.lambda;

S         = par.grids{1,j_pos};
EDUC      = par.educ;
AGE_PROF  = par.inc.age_prof;
INNO_pos  = par.inc.inno_pos;
FE_pos    = par.inc.fe_pos;


V   = zeros(length(S),length(FE_pos),length(EDUC),length(INNO_pos));
C   = zeros(length(S),length(FE_pos),length(EDUC),length(INNO_pos));
Sp  = zeros(length(S),length(FE_pos),length(EDUC),length(INNO_pos));
boundgrid = zeros(length(FE_pos),length(EDUC),length(INNO_pos));

for educ = 1:length(EDUC)
    S         = par.grids{educ,j_pos};
    FE        = par.inc.fe{educ,1};
    INNO      = par.inc.z_val{educ,1}(j_pos,:);
    Spgrid    = par.grids{educ,j_pos+1};
    r         = (r_sav .* (Spgrid>=0) + r_debt .* (Spgrid<0))';
    rs        = (r_sav .* (S>=0) + r_debt .* (S<0));
    for ife=1:length(FE_pos)
        splVp     = griddedInterpolant(Spgrid,squeeze(Vp(:,ife,educ)));
        % People don't work next period so only disutility of consumption. 
        Cpp = Cp(:,ife,educ);
        ucp = beta*(1+r).* Cpp.^(-gammac);
        
        for inno = 1:length(INNO_pos)
            dispinc = (1-lambda) * exp(AGE_PROF{educ,1}(j_pos)) * exp(FE(ife)) * exp(INNO(inno));

            feasible = (dispinc + (1+rs).*S - Spgrid(1)>0);
            not_feasible = logical(1-feasible);

            [C(feasible,ife,educ,inno),Sp(feasible,ife,educ,inno),boundgrid(ife,educ,inno)] = ...
                EGM(par,ucp,dispinc,Spgrid,S(feasible));

            vp = splVp(Sp(feasible,ife,educ,inno));
            V(feasible,ife,educ,inno) = C(feasible,ife,educ,inno).^(1-gammac)/(1-gammac) ...
                + beta*vp;
            C(not_feasible,ife,educ,inno) = 0;
            V(not_feasible,ife,educ,inno) = -(10^(5/gammac)).^(1-gammac)/(1-gammac);
            Sp(not_feasible,ife,educ,inno)= Spgrid(1);
        end
        
    end
end

errors = sum(boundgrid(:))/length(C(:));
switch options.timer_on
    case {'Y'}
        if sum(boundgrid(:)) >= length(EDUC)*length(FE_pos)*length(INNO_pos)
            fprintf('j : %i, Share of errors (increase grid) =%3.3f \n',par.age(j_pos),errors);
        end
end

end