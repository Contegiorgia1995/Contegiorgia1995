function [V,C,Sp] = prob_retT (par,options,Vp,Cp,j_pos)
% Solve household problem when retired and receiving OAS transfers

r_sav     = par.r_sav;
% r_debt     = par.r_sav;
beta      = par.beta;
gammac    = par.gammac;

S         = par.grids{1,j_pos};
H         = par.gridh{1,j_pos};
TAU       = par.gridtrans;
EDUC      = par.educ;

V         = zeros(length(S),length(H),length(EDUC),length(TAU));
C         = zeros(length(S),length(H),length(EDUC),length(TAU));
Sp        = zeros(length(S),length(H),length(EDUC),length(TAU));
boundgrid = zeros(length(H),length(EDUC),length(TAU));

for educ=1:length(EDUC)
    S         = par.grids{educ,j_pos};
    H         = par.gridh{educ,j_pos};
    Spgrid    = par.grids{educ,j_pos+1};
    for ih=1:length(H)
        splVp = griddedInterpolant(Spgrid,Vp(:,ih,educ));
        Cpp = Cp(:,ih,educ);
        ucp = beta*(1+r_sav)* Cpp.^(-gammac);
        ret = ret_rep(par,H(ih),educ,options);
        
        for it=1:length(TAU);
            dispinc = ret + TAU(it);
            
            [C(:,ih,educ,it),Sp(:,ih,educ,it),boundgrid(ih,educ,it)] = EGM(par,ucp,dispinc,Spgrid,S);
            
            V(:,ih,educ,it)  = C(:,ih,educ,it).^(1-gammac)/(1-gammac) + beta*splVp(Sp(:,ih,educ,it));
        end
    end
end

errors = sum(boundgrid(:))/length(C(:));

switch options.timer_on
    case {'Y'}
        if sum(boundgrid(:)) >= length(EDUC)*length(TAU)*length(H);
            fprintf('j : %i, Share of errors (increase grid) =%3.3f \n',par.age(j_pos),errors);
        end
end
end 
