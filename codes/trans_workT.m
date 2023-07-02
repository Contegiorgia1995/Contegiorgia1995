function mup = trans_workT(par,mu0,S_pol,j_pos)
% keyboard
% Transition for workT
% Grids:
S         = par.grids{1,j_pos-1}';
Sp        = par.grids{1,j_pos}';
FE_pos      = par.inc.fe_pos;
INNO_pos    = par.inc.inno_pos;
EDUC      = par.educ;


% Savings:
Sp              = par.grids{1,j_pos}';
QS              = sparse(length(S)*length(FE_pos)*length(EDUC)*length(INNO_pos),length(Sp));
pos             = 1;
npos            = length(S)*length(FE_pos);
for inno = 1:length(INNO_pos)
    for educ = 1:length(EDUC)
        % vectorize policy
        Spol_vec                = S_pol(:,:,educ,inno);
        Spol_vec                = Spol_vec(:);
        Sp                      = par.grids{educ,j_pos}';
        % create transition matrix of savings: linear interpolant for policy
        fspace_s                = fundef({'spli',Sp,0,1});
        QS(pos:pos+npos-1,:)    = funbas(fspace_s,Spol_vec);
        pos                     = pos + npos ;
    end
end

% EDUC + FE : Fixed
QP   = kron(speye(length(FE_pos)*length(EDUC)),ones(length(S),1));
QP   = kron(ones(length(INNO_pos),1),QP);

% Total Q
Q                   = dprod(QP,QS);

% new distribution:
mu0int_vec          = mu0(:);
mu_int_vec          = Q' * mu0int_vec;
mup                 = reshape(mu_int_vec,length(Sp),length(FE_pos),length(EDUC));
end