function [mup,Q_ergo_out] = trans_fertility(par,mu0,Spol,fert_pol,j_pos,Q_ergo_in,options)
% keyboard
% Transition for work young
% Grids:
S           = par.grids{1,j_pos-1};
Sp          = par.grids{1,j_pos};
FE_pos      = par.inc.fe_pos;
INNO_pos    = par.inc.inno_pos;
EDUC        = par.educ;
N           = 1:length(par.N);
%% initial distribution in extended grid, position 1
mu0int_vec                  = mu0(:);

switch Q_ergo_in.created 
    case 'N'
      
        % Savings:
        Sp              = par.grids{1,j_pos}';
        QS              = sparse(length(S)*length(FE_pos)*length(EDUC)*length(INNO_pos),length(Sp));
        pos             = 1;
        npos            = length(S)*length(FE_pos);
        for inno = 1:length(INNO_pos)
            for educ = 1:length(EDUC)
                % vectorize policy
                Spol_vec                = Spol(:,:,educ,inno);
                Spol_vec                = Spol_vec(:);
                Sp                      = par.grids{educ,j_pos}';
                % create transition matrix of savings: linear interpolant for policy
                fspace_s                = fundef({'spli',Sp,0,1});
                QS(pos:pos+npos-1,:)    = funbas(fspace_s,Spol_vec);
                pos                     = pos + npos ;
            end
        end
%         QS                          = repmat(QS,length(N),1);

        % FE : Fixed
        QP   = kron(speye(length(FE_pos)),ones(length(S),1));
        QP   = kron(ones(length(EDUC)*length(INNO_pos),1),QP);
        
        % EDUC + INNO
        P_z        = zeros(length(EDUC),length(INNO_pos),length(INNO_pos));
        P_z(1,:,:) = par.inc.z_prob{1,1}(j_pos,:,:);
        P_z(2,:,:) = par.inc.z_prob{2,1}(j_pos,:,:);
        P_z(3,:,:) = par.inc.z_prob{3,1}(j_pos,:,:);
        
        P_z2                    = sparse(length(EDUC)*length(INNO_pos),length(EDUC)*length(INNO_pos));        
        for iz = 1:length(INNO_pos)
            row_1               = 1 + (iz-1) * length(EDUC);
            row_2               = iz * length(EDUC);
            P_z2(row_1:row_2,:) = dprod(squeeze(P_z(:,iz,:)),speye(length(EDUC)));
        end
        QZ                      = kron(P_z2,ones(length(S)*length(FE_pos),1));
        
        
        % Policy for N
        switch options.Fertility
            case 'Endo'
                fspace              = fundef({'spli',N,0,1});
                QN                  = funbas(fspace,fert_pol(:));
            case 'Exo'
                QN                  = ones(size(fert_pol(:)));
        end
        
        % Total Q
        Q                   = dprod(QP,QS);
        Q                   = dprod(QZ,Q);
        Q                   = dprod(QN,Q);

        
        Q_ergo_out.mat      = Q;
        Q_ergo_out.created = 'Y';
        
    case 'Y'
        Q   = Q_ergo_in.mat;
        Q_ergo_out = Q_ergo_in;
end

% new distribution: 
mu_int_vec = Q' * mu0int_vec;
mup       = reshape(mu_int_vec,length(Sp),length(FE_pos),length(EDUC),length(INNO_pos),length(N));

end