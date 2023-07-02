function [V_2,C_2,Ck_2,Sp_2,PHIp_2,Tp_2] = prob_work_trans(par,options,Vp,Cp,Vc0,j_pos)
% Solve household problem when j = Jc+2: children consume at home and has to set a found with transfers to children
% keyboard
r_sav     = par.r_sav;
r_debt    = par.r_debt;
beta      = par.beta;
gammac    = par.gammac;
w         = par.w;
lambda    = par.lambda;
lambdan   = par.lambdan;
gamman    = par.gamman;

S         = par.grids{1,j_pos};
EDUC      = par.educ;

AGE_PROF  = par.inc.age_prof;
INNO_pos  = par.inc.inno_pos;
FE_pos    = par.inc.fe_pos;

N         = 1:length(par.N);
PHI_pos   = 1:length(par.PHI);

V   = NaN(length(S),length(FE_pos),length(EDUC),length(INNO_pos),length(N),length(PHI_pos));
C   = NaN(length(S),length(FE_pos),length(EDUC),length(INNO_pos),length(N),length(PHI_pos));
Sp  = NaN(length(S),length(FE_pos),length(EDUC),length(INNO_pos),length(N),length(PHI_pos));
Ck  = NaN(length(S),length(FE_pos),length(EDUC),length(INNO_pos),length(N),length(PHI_pos));
Tp  = NaN(length(S),length(FE_pos),length(EDUC),length(INNO_pos),length(N),length(PHI_pos));
boundgrid = zeros(length(FE_pos),length(EDUC),length(INNO_pos),length(N),length(PHI_pos));


%% Without children and phi = 1
switch options.Fertility
    case 'Endo' % Only need to solve case without children in benchmark case where we allow them to choose endogenous fertility
        in   = 1;
        iphi = 1;
        ig   = 1; %Group does not matter in value function tomorrow if no children
        
        for educ = 1:length(EDUC)
            INNOp_prob  = squeeze(par.inc.z_prob{educ,1}(j_pos+1,:,:)); % Note I need j_pos + 1 here!
            S          = par.grids{educ,j_pos};
            Spgrid     = par.grids{educ,j_pos+1};
            r         = (r_sav .* (Spgrid>=0) + r_debt .* (Spgrid<0))';
            rs        = (r_sav .* (S>=0) + r_debt .* (S<0));
            
            inno3       = repmat(INNO_pos,length(S),1);
            FE          = par.inc.fe{educ,1};
            INNO        = par.inc.z_val{educ,1}(j_pos,:);
            [S2,INNO2]  = ndgrid(Spgrid,INNO_pos);
            
            for ife=1:length(FE_pos)
                splVp = griddedInterpolant(S2,INNO2,squeeze(Vp(:,ife,educ,:)));
                for inno = 1:length(INNO_pos)
                    innop_prob      = INNOp_prob(inno,:);
                    
                    Cpp           = squeeze(Cp(:,ife,educ,:));
                    if min(Cpp(:))<0 && abs(min(Cpp(:)))<1e-6
                        Cpp  = max(Cpp,0);
                    end
                    ucp = beta*(1+r).* (Cpp.^(-gammac) * innop_prob');
                    
                    h   = exp(AGE_PROF{educ,1}(j_pos)) * exp(FE(ife))* exp(INNO(inno));
                    dispinc = w*h*(1-lambda);
                    
                    feasible = (dispinc + (1+rs).*S - Spgrid(1)>0);
                    not_feasible = logical(1-feasible);
                    
                    [C(feasible,ife,educ,inno,in,iphi),Sp(feasible,ife,educ,inno,in,iphi),boundgrid(ife,educ,inno,in,iphi)] = ...
                        EGM(par,ucp,dispinc,Spgrid,S(feasible));
                    Ck(feasible,ife,educ,inno,in,iphi) = zeros(size(C(feasible,ife,educ,inno,in,iphi)));
                    
                    Sp3         = repmat(Sp(:,ife,educ,inno),1,length(INNO_pos));
                    vp          = splVp(Sp3,inno3)*innop_prob';
                    
                    V(feasible,ife,educ,inno,in,iphi) = C(feasible,ife,educ,inno,in,iphi).^(1-gammac)/(1-gammac) ...
                        + beta*vp(feasible);
                    C(not_feasible,ife,educ,inno,in,iphi) = 0;
                    V(not_feasible,ife,educ,inno,in,iphi) = -(10^(5/gammac)).^(1-gammac)/(1-gammac);
                    Sp(not_feasible,ife,educ,inno,in,iphi)= Spgrid(1);
                    
                    Tp(:,ife,educ,inno,in,iphi) = zeros(size(C(:,ife,educ,inno,in,iphi)));
                end
            end
        end
        
        % Fill cell for all PHI;
        for iphi = 2:length(PHI_pos)
            for educ = 1:length(EDUC)
                for ife=1:length(FE_pos)
                    for inno = 1:length(INNO_pos)
                        C(:,ife,educ,inno,in,iphi)  = C(:,ife,educ,inno,in,1);
                        Ck(:,ife,educ,inno,in,iphi) = Ck(:,ife,educ,inno,in,1);
                        Sp(:,ife,educ,inno,in,iphi) = Sp(:,ife,educ,inno,in,1);
                        V(:,ife,educ,inno,in,iphi)  = V(:,ife,educ,inno,in,1);
                        Tp(:,ife,educ,inno,in,iphi) = Tp(:,ife,educ,inno,in,1);
                        boundgrid(ife,educ,inno,in,iphi) = boundgrid(ife,educ,inno,in,1);
                    end
                end
            end
        end
end

%% Case with Children
switch options.Fertility
    case 'Endo' % Only need to solve case without children in benchmark case where we allow them to choose endogenous fertility
        in_1    = 2;
    case 'Exo'
        in_1    = 1;
end

for educ = 1:length(EDUC)
    INNOp_prob  = squeeze(par.inc.z_prob{educ,1}(j_pos+1,:,:)); % Note I need j_pos + 1 here!
    S          = par.grids{educ,j_pos};
    Spgrid     = par.grids{educ,j_pos+1};
    r          = (r_sav .* (Spgrid>=0) + r_debt .* (Spgrid<0))';
    rs        = (r_sav .* (S>=0) + r_debt .* (S<0));
    
    inno3       = repmat(INNO_pos,length(S),1);
    FE          = par.inc.fe{educ,1};
    INNO        = par.inc.z_val{educ,1}(j_pos,:);
    [S2,INNO2]  = ndgrid(Spgrid,INNO_pos);
    
    for in = in_1:length(N)
        n             = par.N(in);
        for iphi = 1:length(PHI_pos)
            phi        = par.PHI(iphi);
            for ife=1:length(FE_pos)
                splVp = griddedInterpolant(S2,INNO2,squeeze(Vp(:,ife,educ,:)));
                for inno = 1:length(INNO_pos)
                    innop_prob      = INNOp_prob(inno,:);
                    h               = exp(AGE_PROF{educ,1}(j_pos)) * exp(FE(ife)) * exp(INNO(inno))*(1-lambda);
                    
                    Cpp = squeeze(Cp(:,ife,educ,:));
                    if min(Cpp(:))<0 && abs(min(Cpp(:)))<1e-6
                        Cpp  = max(Cpp,0);
                    end
                    
                    ucp     = beta*(1+r).* (Cpp.^(-gammac) * innop_prob');
                    labor_inc = w*h*(1-lambda);
                    
                    n_final = n* (par.fam_size/2);
                    ChildCost_0 = ChildCost(par,labor_inc,n_final,0,options);
                    ChildCost_1 = ChildCost(par,labor_inc,n_final,1,options);
                    ChildCost_opt = min(ChildCost_0,ChildCost_1);
                    Tp(:,ife,educ,inno,in,iphi)   = (ChildCost_1>ChildCost_0);
                    
                    dispinc = labor_inc- ChildCost_opt -  n_final  *phi;
                    feasible = (dispinc + (1+rs).*S - Spgrid(1)>0);
                    not_feasible = logical(1-feasible);
                    
                    % Altruism
                    hp              = exp(AGE_PROF{educ,1}(j_pos)) * exp(FE(ife)) * exp(INNO(inno));
                    Vc0_aux     = reshape(squeeze(Vc0(iphi,:,:)),par.N_fe*length(par.psy_val_hs),1);
                    Vc0_prob_h0 = squeeze(par.PG(hp));
                    Vc0_prob_psy= par.psy_prob(educ,:)';
                    Vc0_prob    = gridmake(Vc0_prob_h0,Vc0_prob_psy);
                    Vc0_prob    = Vc0_prob(:,1) .* Vc0_prob(:,2);
                    
                    Gn = beta * lambdan * ( n_final )^(gamman) * Vc0_aux' * Vc0_prob;
                    
                    switch options.Ck
                        case 'Yes'
                            [C(feasible,ife,educ,inno,in,iphi),Sp(feasible,ife,educ,inno,in,iphi),boundgrid(ife,educ,inno,in,iphi)] = ...
                                EGM_withchild(par,ucp,dispinc,Spgrid,S(feasible),n_final);
                            
                            f_n     = (lambdan/n_final^(1-gamman))^(1/gammac);
                            Ck(feasible,ife,educ,inno,in,iphi)     = f_n * C(feasible,ife,educ,inno,in,iphi);
                                                      
                            Sp3         = repmat(Sp(:,ife,educ,inno,in,iphi),1,length(INNO_pos));
                            vp          = splVp(Sp3,inno3)*innop_prob';
                            
                            
                            V(feasible,ife,educ,inno,in,iphi) = C(feasible,ife,educ,inno,in,iphi).^(1-gammac)/(1-gammac).*(C(feasible,ife,educ,inno,in,iphi)>=0) ...
                                + lambdan* n_final ^(gamman) * Ck(feasible,ife,educ,inno,in,iphi).^(1-gammac)/(1-gammac).*(C(feasible,ife,educ,inno,in,iphi)>=0) ...
                                + -(10^(5/gammac)).^(1-gammac)/(1-gammac).*(C(feasible,ife,educ,inno,in,iphi)<0) ...
                                + beta*vp(feasible) .*(C(feasible,ife,educ,inno,in,iphi)>=0)...
                                + Gn .*(C(feasible,ife,educ,inno,in,iphi)>=0);
                            % Fix bad extrapolation:
                            Sp(feasible,ife,educ,inno,in,iphi) = Sp(feasible,ife,educ,inno,in,iphi).*(C(feasible,ife,educ,inno,in,iphi) >=0) + Spgrid(1).*(C(feasible,ife,educ,inno,in,iphi)<0);
                            C(feasible,ife,educ,inno,in,iphi)  = C(feasible,ife,educ,inno,in,iphi) .*(C(feasible,ife,educ,inno,in,iphi) >= 0) + 0.*(C(feasible,ife,educ,inno,in,iphi) <0);
                            Ck(feasible,ife,educ,inno,in,iphi) = Ck(feasible,ife,educ,inno,in,iphi) .*(Ck(feasible,ife,educ,inno,in,iphi)>= 0) + 0.*(Ck(feasible,ife,educ,inno,in,iphi)<0);
                            
                            C(not_feasible,ife,educ,inno,in,iphi) = 0;
                            Ck(not_feasible,ife,educ,inno,in,iphi) = 0;
                            V(not_feasible,ife,educ,inno,in,iphi) = -(10^(5/gammac)).^(1-gammac)/(1-gammac);
                            Sp(not_feasible,ife,educ,inno,in,iphi)= Spgrid(1);
                            
                        case 'No'
                            [C(feasible,ife,educ,inno,in,iphi),Sp(feasible,ife,educ,inno,in,iphi),boundgrid(ife,educ,inno,in,iphi)] = ...
                                EGM(par,ucp,dispinc,Spgrid,S(feasible));
                            Ck(feasible,ife,educ,inno,in,iphi)=zeros(size(C(feasible,ife,educ,inno,in,iphi)));
                                                        
                            Sp3         = repmat(Sp(:,ife,educ,inno,in,iphi),1,length(INNO_pos));
                            vp          = splVp(Sp3,inno3)*innop_prob';
                            
                            V(feasible,ife,educ,inno,in,iphi) = C(feasible,ife,educ,inno,in,iphi).^(1-gammac)/(1-gammac).*(C(feasible,ife,educ,inno,in,iphi)>=0) ...
                                + -(10^(5/gammac)).^(1-gammac)/(1-gammac).*(C(feasible,ife,educ,inno,in,iphi)<0) ...
                                + beta*vp(feasible) .*(C(feasible,ife,educ,inno,in,iphi)>=0) ...
                                + Gn .*(C(feasible,ife,educ,inno,in,iphi)>=0);
                            
                            % Fix bad extrapolation:
                            Sp(feasible,ife,educ,inno,in,iphi) = Sp(feasible,ife,educ,inno,in,iphi).*(C(feasible,ife,educ,inno,in,iphi) >=0) + Spgrid(1).*(C(feasible,ife,educ,inno,in,iphi)<0);
                            C(feasible,ife,educ,inno,in,iphi)  = C(feasible,ife,educ,inno,in,iphi) .*(C(feasible,ife,educ,inno,in,iphi) >=0) + 0.*(C(feasible,ife,educ,inno,in,iphi)<0);
                            
                            C(not_feasible,ife,educ,inno,in,iphi) = 0;
                            Ck(not_feasible,ife,educ,inno,in,iphi) = 0;
                            V(not_feasible,ife,educ,inno,in,iphi) = -(10^(5/gammac)).^(1-gammac)/(1-gammac);
                            Sp(not_feasible,ife,educ,inno,in,iphi)= Spgrid(1);
                    end
                end
            end
        end
    end
end


%% Search Grid:

V_2    = NaN(length(S),length(FE_pos),length(EDUC),length(INNO_pos),length(N));
C_2    = NaN(length(S),length(FE_pos),length(EDUC),length(INNO_pos),length(N));
Ck_2   = NaN(length(S),length(FE_pos),length(EDUC),length(INNO_pos),length(N));
Sp_2   = NaN(length(S),length(FE_pos),length(EDUC),length(INNO_pos),length(N));
PHIp_2 = NaN(length(S),length(FE_pos),length(EDUC),length(INNO_pos),length(N));
Tp_2   = NaN(length(S),length(FE_pos),length(EDUC),length(INNO_pos),length(N));

posP = zeros(length(S),length(FE_pos),length(EDUC),length(INNO_pos));

for is = 1:length(S)
    for educ = 1:length(EDUC)
        for ife=1:length(FE_pos)
            for inno = 1:length(INNO_pos)
                for in = 1:length(N)
                    Vaux                       = squeeze(V(is,ife,educ,inno,in,:));
                    [~,posauxP]                = max(Vaux);  % Max wrt phi
                    posP(is,ife,educ,inno,in)   = posauxP;
                    V_2(is,ife,educ,inno,in)     = V(is,ife,educ,inno,in,posP(is,ife,educ,inno,in));
                    C_2(is,ife,educ,inno,in)     = C(is,ife,educ,inno,in,posP(is,ife,educ,inno,in));
                    Ck_2(is,ife,educ,inno,in)    = Ck(is,ife,educ,inno,in,posP(is,ife,educ,inno,in));
                    Sp_2(is,ife,educ,inno,in)    = Sp(is,ife,educ,inno,in,posP(is,ife,educ,inno,in));
                    Tp_2(is,ife,educ,inno,in)    = Tp(is,ife,educ,inno,in,posP(is,ife,educ,inno,in));
                    PHIp_2(is,ife,educ,inno,in)  = posP(is,ife,educ,inno,in);
                end
            end
        end
    end
end


errors = sum(boundgrid(:))/length(C(:));
switch options.timer_on
    case {'Y'}
        if sum(boundgrid(:)) >= length(EDUC)*length(FE_pos)*length(INNO_pos)*length(PHI_pos)*length(N)
            fprintf('j : %i, Share of errors (increase grid) =%3.3f \n',par.age(j_pos),errors);
        end
end
end