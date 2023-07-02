function mup = trans_ret(par,mu0,S_pol,j_pos)
% keyboard
% Transition for retire
% Grids:
S           = par.grids{1,j_pos-1};
Sp          = par.grids{1,j_pos-1};
FE_pos      = par.inc.fe_pos;
EDUC        = par.educ;

%% 1. Choice of savings
% Savings:
Sp              = par.grids{1,j_pos}';
QS              = sparse(length(S)*length(FE_pos)*length(EDUC),length(Sp));
pos             = 1;
npos            = length(S)*length(FE_pos);
for educ = 1:length(EDUC)
    % vectorize policy
    Spol_vec                = S_pol(:,:,educ);
    Spol_vec                = Spol_vec(:);
    Sp                      = par.grids{educ,j_pos}';
    % create transition matrix of savings: linear interpolant for policy
    fspace_s                = fundef({'spli',Sp,0,1});
    QS(pos:pos+npos-1,:)    = funbas(fspace_s,Spol_vec);
    pos                     = pos + npos ;
end

% Educ + H: Fixed
QE                  = kron(speye(length(EDUC)*length(FE_pos)),ones(length(S),1));

% Total Q
Q                   = dprod(QE,QS);

% new distribution:
mu0int_vec          = mu0(:);
mu_int_vec          = Q' * mu0int_vec;
mup                 = reshape(mu_int_vec,length(Sp),length(FE_pos),length(EDUC));

end