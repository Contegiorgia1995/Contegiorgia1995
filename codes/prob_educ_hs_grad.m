function [V,C,Sp] = prob_educ_hs_grad(par,options,Vp,Cp)
% Solve household problem when it is HS Graduate, for age 12 and 16

% HS graduate:
educ      = 2;
PSY    = par.psy_val_hs;
r_sav     = par.r_sav;
r_debt    = par.r_debt;
beta      = par.beta;
gammac    = par.gammac;
w         = par.w;
lambda    = par.lambda;

FE_pos     = par.inc.fe_pos;
AGE_PROF   = par.inc.age_prof;
INNO_pos   = par.inc.inno_pos;


V         = cell(1,3);
C         = cell(1,3);
Sp        = cell(1,3);

Vpaux     = squeeze(Vp(:,:,educ,:));
Cpaux     = squeeze(Cp(:,:,educ,:));

for j_pos = par.Je2_pos+1:-1:par.Je2_pos
    FE         = par.inc.fe{educ,1};
    INNOp_prob = squeeze(par.inc.z_prob{educ,1}(j_pos+1,:,:)); % Note I need j_pos + 1 here!
    S          = par.grids{educ,j_pos}; % This grids do not change by education (at this age)
    inno3 = repmat(INNO_pos,length(S),1);
    
    V2         = zeros(length(S),length(FE_pos),length(INNO_pos));
    C2         = zeros(length(S),length(FE_pos),length(INNO_pos));
    Sp2        = zeros(length(S),length(FE_pos),length(INNO_pos));
    boundgrid  = zeros(length(FE_pos),length(INNO_pos));
    
    Spgrid     = par.grids{educ,j_pos+1};
    r          = (r_sav .* (Spgrid>=0) + r_debt .* (Spgrid<0))';
    rs         = (r_sav .* (S>=0)      + r_debt .* (S<0));
    
    [S2,INNO2] = ndgrid(Spgrid,INNO_pos);
    INNO      = par.inc.z_val{educ,1}(j_pos,:);
    
    for ife = 1:length(FE_pos)
        splVp = griddedInterpolant(S2,INNO2,squeeze(Vpaux(:,ife,:)));
        for inno = 1:length(INNO_pos)
            innop_prob      = INNOp_prob(inno,:);
            h   = exp(AGE_PROF{educ,1}(j_pos)) * exp(FE(ife))* exp(INNO(inno));

            Cpp           = squeeze(Cpaux(:,ife,:));
            if min(Cpp(:))<0 && abs(min(Cpp(:)))<1e-6
                Cpp  = max(Cpp,0);
            end

            ucp = beta*(1+r).* (Cpp.^(-gammac) * innop_prob');
            dispinc = w*h*(1-lambda) + par.init_trans(educ,j_pos);

            feasible = (dispinc + (1+rs).*S - Spgrid(1)>0);
            not_feasible = logical(1-feasible);

            [C2(feasible,ife,inno),Sp2(feasible,ife,inno),boundgrid(ife,inno)] = ...
                GEGM(par,ucp',dispinc,Spgrid,S(feasible),splVp,innop_prob);

            C2(not_feasible,ife,inno) = 0;
            Sp2(not_feasible,ife,inno) = -Spgrid(1);

            Sp3         = repmat(Sp2(:,ife,inno),1,length(INNO_pos));
            vp          = splVp(Sp3,inno3)*innop_prob';

            V2(feasible,ife,inno) =C2(feasible,ife,inno).^(1-gammac)/(1-gammac) ...
                + beta*vp(feasible);
            V2(not_feasible,ife,inno) = -(10^(5/gammac)).^(1-gammac)/(1-gammac);
        end
    end
    
    errors = sum(boundgrid(:))/length(C2(:));
    switch options.timer_on
        case {'Y'}
            if sum(boundgrid(:)) >= length(FE_pos)*length(INNO_pos)
                fprintf('j : %i, Share of errors (increase grid) =%3.3f \n',par.age(j_pos),errors);
            end
    end
    
    V{1,j_pos - (par.Je1_pos-1)} = V2;
    C{1,j_pos - (par.Je1_pos-1)} = C2;
    Sp{1,j_pos - (par.Je1_pos-1)} = Sp2;
    
    Vpaux = V2;
    Cpaux = C2;
end




%% Age 16 (age of HS)
% Expectation about innovation
exp_prob   = squeeze(par.inc.z_prob{educ,1}(j_pos,1,:))';

for j_pos = par.Je1_pos:-1:par.Je1_pos
    S         = par.grids{educ,j_pos};
    inno3 = repmat(INNO_pos,length(S),1);
    FE         = par.inc.fe{1,1}; %Skills before getting educated
    
    V2         = zeros(length(S),length(FE_pos));
    C2         = zeros(length(S),length(FE_pos));
    Sp2        = zeros(length(S),length(FE_pos));
    boundgrid  = zeros(length(FE_pos));
    
    Spgrid     = par.grids{educ,j_pos+1};
    r          = (r_sav .* (Spgrid>=0) + r_debt .* (Spgrid<0))';
    rs         = (r_sav .* (S>=0) + r_debt .* (S<0));
    
    
    for ife = 1:length(FE_pos)
        splVp = griddedInterpolant(S2,INNO2,squeeze(Vpaux(:,ife,:)));
        h     = exp(AGE_PROF{educ,1}(j_pos)) * exp(FE(ife))* exp(INNO(inno));
        
        Cpp           = squeeze(Cpaux(:,ife,:));
        if min(Cpp(:))<0 && abs(min(Cpp(:)))<1e-6
            Cpp  = max(Cpp,0);
        end
        
        ucp     = beta*(1+r).* (Cpp.^(-gammac)* exp_prob');
        dispinc = -par.pe1 + par.init_trans(educ,j_pos);
        
        feasible = (dispinc + (1+rs).*S - Spgrid(1)>0);
        not_feasible = logical(1-feasible);
        
        [C2(feasible,ife),Sp2(feasible,ife),boundgrid(ife)] = ...
            GEGM(par,ucp',dispinc,Spgrid,S(feasible),splVp,exp_prob);
        
        C2(not_feasible,ife)  = 0;
        Sp2(not_feasible,ife) = 0;
        
        Sp3         = repmat(Sp2(:,ife),1,length(INNO_pos));
        vp          = splVp(Sp3,inno3)*exp_prob';
        
        V2(feasible,ife) = C2(feasible,ife).^(1-gammac)/(1-gammac) ...
            + beta*vp(feasible);
        
        V2(not_feasible,ife) = -(10^(5/gammac)).^(1-gammac)/(1-gammac);
    end
    
    if j_pos > par.Je1_pos+1
        V{1,j_pos - (par.Je1_pos-1)} = V2;
        C{1,j_pos - (par.Je1_pos-1)} = C2;
        Sp{1,j_pos - (par.Je1_pos-1)} = Sp2;
    else
        V{1,j_pos - (par.Je1_pos-1)} = repmat(V2,1,1,length(PSY))-repmat(reshape(PSY,1,1,length(PSY)),size(V2,1),size(V2,2));
        C{1,j_pos - (par.Je1_pos-1)} =  repmat(C2,1,1,length(PSY));
        Sp{1,j_pos - (par.Je1_pos-1)} =  repmat(Sp2,1,1,length(PSY));
    end
    
    Vpaux = V2;
    Cpaux = C2;
end


end