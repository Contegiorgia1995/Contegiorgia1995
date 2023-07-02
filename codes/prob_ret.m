function [V,C,Sp] = prob_ret (par,options,Vp,Cp,j_pos)
% Solve household problem when retired
% keyboard
r_sav   = par.r_sav;
beta    = par.beta;
gammac  = par.gammac;
EDUC    = par.educ;
FE_pos  = par.inc.fe_pos;

if j_pos == par.Jd_pos-1 % consume all, next period is dead    
    S       = par.grids{1,j_pos};
    Sp      = zeros(length(S),length(FE_pos),length(EDUC)); % Last period, no savings
    C       = zeros(length(S),length(FE_pos),length(EDUC));
    V       = zeros(length(S),length(FE_pos),length(EDUC));
    
    for educ  = 1:length(EDUC)
        S   = par.grids{educ,j_pos};

        C(:,:,educ) =  (1+r_sav)* repmat(S',1,length(FE_pos)) + repmat(ret_rep(par,FE_pos,educ,options),length(S),1);
        V(:,:,educ) = C(:,:,educ).^(1-gammac)/(1-gammac);
    end
    
else % solve Euler equation
    S       = par.grids{1,j_pos};
    
    Sp      = zeros(length(S),length(FE_pos),length(EDUC)); 
    C       = zeros(length(S),length(FE_pos),length(EDUC));
    V       = zeros(length(S),length(FE_pos),length(EDUC));

    boundgrid = zeros(length(FE_pos),length(EDUC),1);

    for educ=1:length(EDUC)
        S       = par.grids{educ,j_pos};
        Spgrid  = par.grids{educ,j_pos+1};
        
        for ife=1:length(FE_pos)
            splVp = griddedInterpolant(Spgrid,Vp(:,ife,educ));
            Cpp   = Cp(:,ife,educ);

            ucp     = beta*(1+r_sav)* Cpp.^(-gammac);
            dispinc = ret_rep(par,ife,educ,options);

            [C(:,ife,educ),Sp(:,ife,educ),boundgrid(ife,educ)] = EGM(par,ucp,dispinc,Spgrid,S);
            V(:,ife,educ)  = C(:,ife,educ).^(1-gammac)/(1-gammac) + beta*splVp(Sp(:,ife,educ));
        end
    end
        
    errors = sum(boundgrid(:))/length(C(:));
    
    switch options.timer_on
        case {'Y'}
            if sum(boundgrid(:)) >= length(EDUC)*length(FE_pos);
                fprintf('j : %i, Share of errors (increase grid) =%3.3f \n',par.age(j_pos),errors);
            end
    end
end
end