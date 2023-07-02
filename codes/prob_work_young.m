function [V,C,Sp] = prob_work_young(par,options,Vp,Cp,j_pos)
% Solve household problem when young
r_sav     = par.r_sav;
r_debt    = par.r_debt;
beta      = par.beta;
gammac    = par.gammac;
w         = par.w;
lambda    = par.lambda;

S         = par.grids{1,j_pos};
EDUC      = par.educ;

AGE_PROF   = par.inc.age_prof;
INNO_pos   = par.inc.inno_pos;
FE_pos     = par.inc.fe_pos;

V         = zeros(length(S),length(FE_pos),length(EDUC),length(INNO_pos));
C         = zeros(length(S),length(FE_pos),length(EDUC),length(INNO_pos));
Sp        = zeros(length(S),length(FE_pos),length(EDUC),length(INNO_pos));
boundgrid = zeros(length(FE_pos),length(EDUC),length(INNO_pos));

for educ = 1:length(EDUC)
    INNOp_prob = squeeze(par.inc.z_prob{educ,1}(j_pos+1,:,:)); % Note I need j_pos + 1 here!
    
    S          = par.grids{educ,j_pos};
    Spgrid     = par.grids{educ,j_pos+1};
    r          = (r_sav .* (Spgrid>=0) + r_debt .* (Spgrid<0))';
    rs         = (r_sav .* (S>=0) + r_debt .* (S<0));
    
    inno3 = repmat(INNO_pos,length(S),1);
    [S2,INNO2] = ndgrid(Spgrid,INNO_pos);
    FE        = par.inc.fe{educ,1};
    INNO      = par.inc.z_val{educ,1}(j_pos,:);
    for ife=1:length(FE_pos)
        splVp = griddedInterpolant(S2,INNO2,squeeze(Vp(:,ife,educ,:)));
        for inno = 1:length(INNO_pos)
            innop_prob      = INNOp_prob(inno,:);
            h   = exp(AGE_PROF{educ,1}(j_pos)) * exp(FE(ife)) * exp(INNO(inno));
            
            Cpp = squeeze(Cp(:,ife,educ,:));
            if min(Cpp(:))<0 && abs(min(Cpp(:)))<1e-6
                Cpp  = max(Cpp,0);
            end
            
            ucp = beta*(1+r).* (Cpp.^(-gammac) * innop_prob');
            dispinc = w*h*(1-lambda);
            
            feasible = (dispinc + (1+rs).*S - Spgrid(1)>0);
            not_feasible = logical(1-feasible);
            
            [C(feasible,ife,educ,inno),Sp(feasible,ife,educ,inno),boundgrid(ife,educ,inno)] = ...
                GEGM(par,ucp',dispinc,Spgrid,S(feasible),splVp,innop_prob);
            %             EGM(par,ucp,dispinc,Spgrid,S);
            
            Sp3         = repmat(Sp(:,ife,educ,inno),1,length(INNO_pos));
            vp          = splVp(Sp3,inno3)*innop_prob';
            
            V(feasible,ife,educ,inno) =C(feasible,ife,educ,inno).^(1-gammac)/(1-gammac) ...
                + beta*vp(feasible);
            
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
