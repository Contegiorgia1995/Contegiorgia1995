function [V,C,Ck,Sp,Tp] = prob_work_with_child (par,options,Vp,Cp,j_pos)
% Solve household problem when j = Jc+1: children consume at home
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

AGE_PROF   = par.inc.age_prof;
INNO_pos   = par.inc.inno_pos;
FE_pos     = par.inc.fe_pos;

N         = 1:length(par.N);

V   = zeros(length(S),length(FE_pos),length(EDUC),length(INNO_pos),length(N));
C   = zeros(length(S),length(FE_pos),length(EDUC),length(INNO_pos),length(N));
Sp  = zeros(length(S),length(FE_pos),length(EDUC),length(INNO_pos),length(N));
Ck  = zeros(length(S),length(FE_pos),length(EDUC),length(INNO_pos),length(N));
Tp  = zeros(length(S),length(FE_pos),length(EDUC),length(INNO_pos),length(N));
boundgrid  = zeros(length(FE_pos),length(EDUC),length(INNO_pos),length(N));

%% Without children and phi = 1
switch options.Fertility
    case 'Endo' % Only need to solve case without children in benchmark case where we allow them to choose endogenous fertility
        in   = 1;
        
        for educ = 1:length(EDUC)
            INNOp_prob = squeeze(par.inc.z_prob{educ,1}(j_pos+1,:,:)); % Note I need j_pos + 1 here!
            S          = par.grids{educ,j_pos};
            Spgrid     = par.grids{educ,j_pos+1};
            r          = (r_sav .* (Spgrid>=0) + r_debt .* (Spgrid<0))';
            rs         = (r_sav .* (S>=0) + r_debt .* (S<0));
            
            inno3       = repmat(INNO_pos,length(S),1);
            [S2,INNO2] = ndgrid(Spgrid,INNO_pos);
            FE        = par.inc.fe{educ,1};
            INNO      = par.inc.z_val{educ,1}(j_pos,:);
            
            for ife=1:length(FE_pos)
                splVp = griddedInterpolant(S2,INNO2,squeeze(Vp(:,ife,educ,:,in)));
                for inno = 1:length(INNO_pos)
                    innop_prob      = INNOp_prob(inno,:);
                    
                    h   = exp(AGE_PROF{educ,1}(j_pos)) * exp(FE(ife)) * exp(INNO(inno));
                    
                    
                    Cpp           = squeeze(Cp(:,ife,educ,:,in));
                    if min(Cpp(:))<0 && abs(min(Cpp(:)))<1e-6
                        Cpp  = max(Cpp,0);
                    end
                    ucp = beta*(1+r).* (Cpp.^(-gammac) * innop_prob');
                    
                    dispinc = w*h*(1-lambda);
                    feasible = (dispinc + (1+rs).*S - Spgrid(1)>0);
                    not_feasible = logical(1-feasible);
                    
                    [C(feasible,ife,educ,inno,in),Sp(feasible,ife,educ,inno,in),boundgrid(ife,educ,inno,in)] = ...
                        EGM(par,ucp,dispinc,Spgrid,S(feasible));
                    Ck(feasible,ife,educ,inno,in) = zeros(size(C(feasible,ife,educ,inno,in)));
                    
                    
                    Sp3         = repmat(Sp(:,ife,educ,inno,in),1,length(INNO_pos));
                    vp          = splVp(Sp3,inno3)*innop_prob';
                    
                    V(feasible,ife,educ,inno,in) = C(feasible,ife,educ,inno,in).^(1-gammac)/(1-gammac) ...
                        + beta*vp(feasible);
                    
                    C(not_feasible,ife,educ,inno,in) = 0;
                    V(not_feasible,ife,educ,inno,in) = -(10^(5/gammac)).^(1-gammac)/(1-gammac);
                    Sp(not_feasible,ife,educ,inno,in)= Spgrid(1);
                    
                    Tp(:,ife,educ,inno,in) = zeros(size(C(:,ife,educ,inno,in)));
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
    INNOp_prob = squeeze(par.inc.z_prob{educ,1}(j_pos+1,:,:)); % Note I need j_pos + 1 here!
    S          = par.grids{educ,j_pos};
    Spgrid     = par.grids{educ,j_pos+1};
    r          = (r_sav .* (Spgrid>=0) + r_debt .* (Spgrid<0))';
    rs         = (r_sav .* (S>=0) + r_debt .* (S<0));
    
    [S2,INNO2] = ndgrid(Spgrid,INNO_pos);
    FE        = par.inc.fe{educ,1};
    INNO      = par.inc.z_val{educ,1}(j_pos,:);
    inno3     = repmat(INNO_pos,length(S),1);

    for in = in_1:length(N)
        n             = par.N(in);
        for ife=1:length(FE_pos)
            splVp = griddedInterpolant(S2,INNO2,squeeze(Vp(:,ife,educ,:,in)));
            for inno = 1:length(INNO_pos)
                innop_prob      = INNOp_prob(inno,:);
                
                h   = exp(AGE_PROF{educ,1}(j_pos)) * exp(FE(ife)) * exp(INNO(inno));
                
                Cpp = squeeze(Cp(:,ife,educ,:,in));
                if min(Cpp(:))<0 && abs(min(Cpp(:)))<1e-6
                    Cpp  = max(Cpp,0);
                end
                
                ucp       = beta*(1+r).* (Cpp.^(-gammac) * innop_prob');
                labor_inc = w*h*(1-lambda);
                
                n_final = n* (par.fam_size/2);
                
                ChildCost_0 = ChildCost(par,labor_inc,n_final,0,options);
                ChildCost_1 = ChildCost(par,labor_inc,n_final,1,options);
                ChildCost_opt = min(ChildCost_0,ChildCost_1);
                Tp(:,ife,educ,inno,in)   = (ChildCost_1>ChildCost_0);
                
                dispinc = labor_inc - ChildCost_opt;
                feasible = (dispinc + (1+rs).*S - Spgrid(1)>0);
                not_feasible = logical(1-feasible);
                
                
                switch options.Ck
                    case 'Yes'
                        [C(feasible,ife,educ,inno,in),Sp(feasible,ife,educ,inno,in),boundgrid(ife,educ,inno,in)] = ...
                            GEGM_withchild(par,ucp',dispinc,Spgrid,S(feasible),splVp,innop_prob,n_final);
                        %                         EGM_withchild(par,ucp,dispinc,Spgrid,S,n_final); %%%% To do: Generalized EGM with child
                        f_n     = (lambdan/n_final^(1-gamman))^(1/gammac);
                        Ck(feasible,ife,educ,inno,in)     = f_n * C(feasible,ife,educ,inno,in);
                        
                        
                        Sp3         = repmat(Sp(:,ife,educ,inno,in),1,length(INNO_pos));
                        vp          = splVp(Sp3,inno3)*innop_prob';
                        
                        V(feasible,ife,educ,inno,in) = C(feasible,ife,educ,inno,in).^(1-gammac)/(1-gammac).*(C(feasible,ife,educ,inno,in)>=0) ...
                            + lambdan* n_final ^(gamman) * Ck(feasible,ife,educ,inno,in).^(1-gammac)/(1-gammac).*(C(feasible,ife,educ,inno,in)>=0) ...
                            + -(10^(5/gammac)).^(1-gammac)/(1-gammac).*(C(feasible,ife,educ,inno,in)<0) ...
                            + beta*vp(feasible);
                        % Fix bad extrapolation:
                        Sp(feasible,ife,educ,inno,in) = Sp(feasible,ife,educ,inno,in).*(C(feasible,ife,educ,inno,in) >=0) + Spgrid(1).*(C(feasible,ife,educ,inno,in)<0);
                        C(feasible,ife,educ,inno,in)  = C(feasible,ife,educ,inno,in) .*(C(feasible,ife,educ,inno,in) >= 0) + 0.*(C(feasible,ife,educ,inno,in) <0);
                        Ck(feasible,ife,educ,inno,in) = Ck(feasible,ife,educ,inno,in) .*(Ck(feasible,ife,educ,inno,in)>= 0) + 0.*(Ck(feasible,ife,educ,inno,in)<0);
                        
                        C(not_feasible,ife,educ,inno,in) = 0;
                        Ck(not_feasible,ife,educ,inno,in) = 0;
                        V(not_feasible,ife,educ,inno,in) = -(10^(5/gammac)).^(1-gammac)/(1-gammac);
                        Sp(not_feasible,ife,educ,inno,in)= Spgrid(1);
                        
                    case 'No'
                        [C(feasible,ife,educ,inno,in),Sp(feasible,ife,educ,inno,in),boundgrid(ife,educ,inno,in)] = ...
                            GEGM(par,ucp',dispinc,Spgrid,S(feasible),splVp,innop_prob);
                        % EGM(par,ucp,dispinc,Spgrid,S);
                        Ck(feasible,ife,educ,inno,in) = zeros(size(C(feasible,ife,educ,inno,in)));
                        
                        Sp3         = repmat(Sp(:,ife,educ,inno,in),1,length(INNO_pos));
                        vp          = splVp(Sp3,inno3)*innop_prob';
                        
                        V(feasible,ife,educ,inno,in) = C(feasible,ife,educ,inno,in).^(1-gammac)/(1-gammac).*(C(feasible,ife,educ,inno,in)>=0) ...
                            + -(10^(5/gammac)).^(1-gammac)/(1-gammac).*(C(feasible,ife,educ,inno,in)<0) ...
                            + beta*vp(feasible);
                        
                        C(not_feasible,ife,educ,inno,in) = 0;
                        Ck(not_feasible,ife,educ,inno,in) = 0;
                        V(not_feasible,ife,educ,inno,in) = -(10^(5/gammac)).^(1-gammac)/(1-gammac);
                        Sp(not_feasible,ife,educ,inno,in)= Spgrid(1);
                end
            end
        end
    end
end


errors = sum(boundgrid(:))/length(C(:));
switch options.timer_on
    case {'Y'}
        if sum(boundgrid(:)) >= length(EDUC)*length(FE_pos)*length(INNO_pos)*length(N)
            fprintf('j : %i, Share of errors (increase grid) =%3.3f \n',par.age(j_pos),errors);
        end
end
end