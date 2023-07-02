function mup = trans_retT(par,mu0,S_pol,j_pos)
% keyboard
% Transition for retire with OAS as state variable
% Grids:
S           = par.grids{1,j_pos-1};
Sp          = par.grids{1,j_pos-1};
H           = par.gridh{1,j_pos-1};
Hp          = par.gridh{1,j_pos};
EDUC        = par.educ;
TAU         = par.gridtrans;

%% 1. Choice of savings 
mup                 = zeros(length(Sp),length(Hp),length(EDUC));
for itau = 1:length(TAU)    
    % Savings:
    Sp              = par.grids{1,j_pos}';
    QS              = sparse(length(S)*length(H)*length(EDUC),length(Sp));
    pos             = 1;
    npos            = length(S)*length(H);
    for educ = 1:length(EDUC)
        % vectorize policy
        Spol_vec                = S_pol(:,:,educ,itau);
        Spol_vec                = Spol_vec(:);
        Sp                      = par.grids{educ,j_pos}';
        % create transition matrix of savings: linear interpolant for policy
        fspace_s                = fundef({'spli',Sp,0,1});
        QS(pos:pos+npos-1,:)    = funbas(fspace_s,Spol_vec);
        pos                     = pos + npos ;
    end
    
    % Educ + H: Fixed
    QE                  = kron(speye(length(EDUC)*length(H)),ones(length(S),1));
    
    % Total Q
    Q                   = dprod(QE,QS);
    
    
    % new distribution:
    mu0int_vec          = mu0(:,:,:,itau);
    mu0int_vec          = mu0int_vec(:);
    mu_int_vec          = Q' * mu0int_vec;
    mup                 = mup + reshape(mu_int_vec,length(Sp),length(Hp),length(EDUC));
end
end