function mup = trans_work_old(par,mu0,Spol,j_pos)
% keyboard
% Transition for work old
% Grids:
S         = par.grids{1,j_pos-1}';
Sp        = par.grids{1,j_pos}';
H         = par.gridh{1,j_pos-1}';
Hp        = par.gridh{1,j_pos}';
EDUC      = par.educ;
Ndeltas   = length(par.deltas{1,1}(j_pos,:));

mup       = zeros(length(Sp),length(Hp),length(EDUC));

for educ = 1:length(EDUC)
    % H policy function
    Hpol            = par.gridh{educ,j_pos-1}' * par.deltas{educ,1}(j_pos-1,:);
    Hpol_pr         = repmat(par.prob{educ,1}(j_pos-1,:),length(H),1);
    Hp              = par.gridh{1,j_pos}';
    QH              = sparse(length(H),length(Hp));
    for delta = 1:Ndeltas
        Hpol_vec                = Hpol(:,delta);
        Hpol_vec                = Hpol_vec(:);
        Hp                      = par.gridh{educ,j_pos}';
        fspace_h                = fundef({'spli',Hp,0,1});
        Qh_aux                  = funbas(fspace_h,Hpol_vec);
        Hpol_pr_vec             = squeeze(Hpol_pr(:,delta));
        Hpol_pr_vec             = Hpol_pr_vec(:);
        QH                      = QH + Qh_aux .* repmat(Hpol_pr_vec,1,size(Qh_aux,2));
    end
    QH                  = kron(QH,ones(length(S),1));
    
    % Savings:
    Sp                      = par.grids{educ,j_pos}';
    fspace_s                = fundef({'spli',Sp,0,1});
    Spol_vec                = Spol(:,:,educ);
    Spol_vec                = Spol_vec(:);
    QS                      = funbas(fspace_s,Spol_vec);
    
    % Total Q
    Q                   = dprod(QH,QS);
    
    mu0_aux                = mu0(:,:,educ);
    mu0_aux                = mu0_aux(:);
    mup_aux                = Q' * mu0_aux;
    mup(:,:,educ) = reshape(mup_aux,length(Sp),length(Hp));
end

end