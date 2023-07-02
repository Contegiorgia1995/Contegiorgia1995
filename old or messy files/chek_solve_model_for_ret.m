
V   = cell(par.Jd_pos,1);
C_   = cell(par.Jd_pos,1);
Ck  = cell(par.Jd_pos,1);
S   = cell(par.Jd_pos,1);
Tp  = cell(3,1);

%% Dead: V_Jd = 0
V{par.Jd_pos} = zeros(par.Ls(end),par.N_fe,length(par.educ));
C_{par.Jd_pos} = zeros(par.Ls(end),par.N_fe,length(par.educ));
S{par.Jd_pos} = zeros(par.Ls(end),par.N_fe,length(par.educ));

%for j_pos = par.Jd_pos-1:-1:par.Jr_pos
for j_pos = par.Jd_pos-1
Cp = C_{j_pos+1}

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
            %splVp = griddedInterpolant(Spgrid,Vp(:,ife,educ));
            Cpp   = Cp(:,ife,educ);

            ucp     = beta*(1+r_sav)* Cpp.^(-gammac);
            dispinc = ret_rep(par,ife,educ,options);

            %[C(:,ife,educ),Sp(:,ife,educ),boundgrid(ife,educ)] = EGM(par,ucp,dispinc,Spgrid,S);
            %V(:,ife,educ)  = C(:,ife,educ).^(1-gammac)/(1-gammac) + beta*splVp(Sp(:,ife,educ));
end
end
end
end



