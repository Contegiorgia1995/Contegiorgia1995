function [mup, Q_ergo_out] = trans_educ(par,mu0,S_hsd,S_hsg,S_cg,j_pos,Q_ergo_in)
% Transition from age Je1 to Je2
% keyboard
% Grids:
S           = par.grids{1,j_pos-1};
Sp          = par.grids{1,j_pos};
FE_pos      = par.inc.fe_pos;
INNO_pos    = par.inc.inno_pos;
EDUC        = par.educ;
if j_pos-1 < par.Je1_pos + 1 % HS
    PSY         = par.psy_val_col;
else
    PSY         = 1;
end
%% 1. Choice of savings in education

% mu0int          = zeros(length(Sp),length(FE),length(phi_pos),length(PSY),length(INNO));

switch Q_ergo_in.created
    case 'N'
        Q = cell(length(FE_pos),length(PSY));
        for ife = 1:length(FE_pos)
            for ipsy = 1:length(PSY)
                % Saving policy function
                Spol          = zeros(length(S),length(EDUC),length(INNO_pos));
                if j_pos-1 < par.Je1_pos + 1 % HS
                    Spol(:,1,:) = reshape(squeeze(S_hsd(:,ife,:)),length(S),1,length(INNO_pos));
                    Spol(:,2,:) = repmat(squeeze(S_hsg(:,ife,ipsy)),1,1,length(INNO_pos));
                    Spol(:,3,:) = repmat(squeeze(S_cg(:,ife,ipsy)),1,1,length(INNO_pos));
                elseif j_pos-1 < par.Je2_pos + 2 % College
                    Spol(:,1,:) = reshape(squeeze(S_hsd(:,ife,:)),length(S),1,length(INNO_pos));
                    Spol(:,2,:) = reshape(squeeze(S_hsg(:,ife,:)),length(S),1,length(INNO_pos));
                    Spol(:,3,:) = repmat(squeeze(S_cg(:,ife)),1,1,length(INNO_pos));
                end
                    
                    
                    % 				Spol(:,:,3,:) = repmat(squeeze(S_st(:,ife,:)),1,1,1,length(INNO_pos));
                
                % H policy function
                % INNO function
                P_z        = zeros(length(EDUC),length(INNO_pos),length(INNO_pos));
                if j_pos-1 < par.Je1_pos + 1 % HS
                    P_z(1,:,:) = par.inc.z_prob{1,1}(j_pos,:,:);
                    P_z(2,:,:) = speye(length(INNO_pos));
                    P_z(3,:,:) = speye(length(INNO_pos));
                elseif j_pos-1 < par.Je2_pos + 2 % College
                    P_z(1,:,:) = par.inc.z_prob{1,1}(j_pos,:,:);
                    P_z(2,:,:) = par.inc.z_prob{2,1}(j_pos,:,:);
                    P_z(3,:,:) = speye(length(INNO_pos));
                end
                
                % Savings:
                Sp              = par.grids{1,j_pos}';
                QS              = sparse(length(S)*length(EDUC)*length(INNO_pos),length(Sp));
                pos             = 1;
                npos            = length(S);
                for inno = 1:length(INNO_pos)
                    for educ = 1:length(EDUC)
                        % vectorize policy
                        Spol_vec                = Spol(:,educ,inno);
                        Spol_vec                = Spol_vec(:);
                        Sp                      = par.grids{educ,j_pos}';
                        % create transition matrix of savings: linear interpolant for policy
                        fspace_s                = fundef({'spli',Sp,0,1});
                        QS(pos:pos+npos-1,:)    = funbas(fspace_s,Spol_vec);
                        pos                     = pos + npos ;
                    end
                end
                                
                % EDUC + INNO
                P_z2                    = sparse(length(EDUC)*length(INNO_pos),length(EDUC)*length(INNO_pos));
                for iz = 1:length(INNO_pos)
                    row_1               = 1 + (iz-1) * length(EDUC);
                    row_2               = iz * length(EDUC);
                    P_z2(row_1:row_2,:) = dprod(squeeze(P_z(:,iz,:)),speye(length(EDUC)));
                end
                QZ                      = kron(P_z2,ones(length(S),1));
                
                % Total Q
                %                 Q                   = dprod(QP,QS);
%                 Q{ife}          = dprod(QP2,QS);
                Q{ife,ipsy}       = dprod(QZ,QS);
            end
        end
        Q_ergo_out.mat   = Q;
        Q_ergo_out.created  = 'Y';
        
    case 'Y'
        Q   = Q_ergo_in.mat;
        Q_ergo_out = Q_ergo_in;
end

if j_pos-1 < par.Je1_pos + 1
    mup = zeros(length(Sp),length(FE_pos),length(PSY),length(EDUC),length(INNO_pos));
    for ife = 1:length(FE_pos)
        for ipsy = 1:length(PSY)
            % new distribution:
            mu0int_vec = mu0(:,ife,ipsy,:);
            mu_int_vec = Q{ife,ipsy}' * mu0int_vec(:);
            mup(:,ife,ipsy,:,:)    = reshape(mu_int_vec,length(Sp),1,1,length(EDUC),length(INNO_pos));
        end
    end
    mup        = squeeze(sum(mup,3));
    
else
    mup = zeros(length(Sp),length(FE_pos),length(EDUC),length(INNO_pos));
    for ife = 1:length(FE_pos)
            % new distribution:
            mu0int_vec = mu0(:,ife,:);
            mu_int_vec = Q{ife,1}' * mu0int_vec(:);
            mup(:,ife,:,:)    = reshape(mu_int_vec,length(Sp),1,length(EDUC),length(INNO_pos));
    end
end


end