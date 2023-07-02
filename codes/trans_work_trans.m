function [mup,muc,mu_ige,pop,Q_ergo_out] = trans_work_trans(par,mu0,Spol,phi_pol,j_pos,Q_ergo_in)
% keyboard
% Transition for fertility period
S           = par.grids{1,j_pos-1};
Sp          = par.grids{1,j_pos};
FE_pos      = par.inc.fe_pos;
INNO_pos    = par.inc.inno_pos;
EDUC        = par.educ;
N           = 1:length(par.N);
PHI         = 1:length(par.PHI);
%%

switch Q_ergo_in.created
    case 'N'
        % Saving policy function: Spol
        % Savings:
        Sp              = par.grids{1,j_pos}';
        QS              = sparse(length(S)*length(FE_pos)*length(EDUC)*length(INNO_pos)*length(N),length(Sp));
        pos             = 1;
        npos            = length(S)*length(FE_pos);
        for in = 1:length(N)
            for inno = 1:length(INNO_pos)
                for educ = 1:length(EDUC)
                    % vectorize policy
                    Spol_vec                = Spol(:,:,educ,inno,in);
                    Spol_vec                = Spol_vec(:);
                    Sp                      = par.grids{educ,j_pos}';
                    % create transition matrix of savings: linear interpolant for policy
                    fspace_s                = fundef({'spli',Sp,0,1});
                    QS(pos:pos+npos-1,:)    = funbas(fspace_s,Spol_vec);
                    pos                     = pos + npos ;
                end
            end
        end
        
        % FE : Fixed
        QP   = kron(speye(length(FE_pos)),ones(length(S),1));
        QP   = kron(ones(length(EDUC)*length(INNO_pos)*length(N),1),QP);
        
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
        QZ                  = kron(P_z2,ones(length(S)*length(FE_pos),1));
        QZ                  = kron(ones(length(N),1),QZ);
        
        % N: Fixed
        QN                  = kron(speye(length(N)),ones(length(S)*length(FE_pos)*length(EDUC)*length(INNO_pos),1));
        
        % Total Q
        Q                   = dprod(QP,QS);
        Q                   = dprod(QZ,Q);
        Q                   = dprod(QN,Q);
        
        % Outside main Q (for loop)
        % PHIp: linear interpolant for policy of PHI
        fspace          = fundef({'spli',PHI,0,1});
        QPhi            = funbas(fspace,phi_pol(:));
        QPhi            = full(QPhi);
        
        % Groups: INNO at time of choosing transfer
        QG              = kron(speye(length(INNO_pos)),ones(length(FE_pos)*length(EDUC)*length(S),1));
        QG              = kron(ones(length(N),1),QG);
        QG              = full(QG);
        
        Q_ergo_out.mat     = Q;
        Q_ergo_out.mat2    = QPhi;
        Q_ergo_out.mat3    = QG;
        Q_ergo_out.created = 'Y';
        
    case 'Y'
        Q      = Q_ergo_in.mat;
        QPhi   = Q_ergo_in.mat2;
        QG     = Q_ergo_in.mat3;
        Q_ergo_out = Q_ergo_in;
end

% new distribution:
mu0_vec     = mu0(:);
mup         = zeros(length(Sp),length(FE_pos),length(EDUC),length(INNO_pos),length(N),length(PHI),length(INNO_pos));


for ip = 1:length(PHI)
    for inno = 1:length(INNO_pos)
        mup_aux              = Q' * (mu0_vec.*QPhi(:,ip).*QG(:,inno));
        mup(:,:,:,:,:,ip,inno) = reshape(mup_aux,length(Sp),length(FE_pos),length(EDUC),length(INNO_pos),length(N));
    end
end

switch par.final_iter
    case 'Y'
        
        %% New cohort
        j_child_pos = par.Je1_pos;
        % S0c    = par.grids{j_child_pos};
        S0c      = par.grids{1,j_child_pos};
        PSY       = par.psy_val_hs;
        AGE_PROF  = par.inc.age_prof;
        
        muc      = zeros(length(S0c),length(FE_pos),length(PSY));
        mu_ige   = zeros(length(Sp),length(FE_pos),length(EDUC),length(INNO_pos),length(S0c),length(FE_pos),length(PSY));
        for educ_p = 1:length(EDUC)
            FE          = par.inc.fe{educ_p,1};
            INNO        = par.inc.z_val{educ_p,1}(j_pos-1,:);
            PSY_p    = par.psy_prob(educ_p,:)';
            psy_prob = repmat(reshape(PSY_p,1,1,length(PSY)),length(Sp),length(FE_pos),1);
            for fe_p0 = 1:length(FE_pos)
                for inno_p0 = 1:length(INNO_pos)
                    hp0     = exp(AGE_PROF{educ_p,1}(j_pos-1)) * exp(FE(fe_p0)) * exp(INNO(inno_p0));
                    prob_h0 = squeeze(par.PG(hp0));
                    
                    prob_n   = squeeze(sum(sum(sum(mup(:,fe_p0,educ_p,:,:,:,inno_p0),6),4),1));
                    ind_n    = (prob_n>0);
                    aux_n    = 1:length(N);
                    for in = aux_n(ind_n)
                        prob_p   = squeeze(sum(sum(mup(:,fe_p0,educ_p,:,in,:,inno_p0),1),4));
                        ind_p    = (prob_p>0);
                        aux_p    = 1:length(PHI);
                        for ip = aux_p(ind_p)
                            mu_aux                      = sum(sum(mup(:,fe_p0,educ_p,:,in,ip,inno_p0),1),4)*par.N(in)*par.fam_size/2;
                            muc(ip,:,:)                 = squeeze(muc(ip,:,:)) + (prob_h0 * mu_aux)*PSY_p';
                            
                            mu_ige_aux                   = sum(mup(:,fe_p0,educ_p,:,in,ip,inno_p0),4) * par.N(in) * par.fam_size/2;
                            mu_ige(:,fe_p0,educ_p,inno_p0,ip,:,:) = squeeze(mu_ige(:,fe_p0,educ_p,inno_p0,ip,:,:)) + repmat(kron(prob_h0',mu_ige_aux),1,1,length(PSY)).*psy_prob;
                            
                            %                 for ipsy = 1:length(PSY)
                            %                     mu_aux                      = sum(mup(:,h_p,educ_p,in,ip))*par.N(in)*par.fam_size/2;
                            %                     muc(ip,:,ipsy)              = squeeze(muc(ip,:,ipsy)) + prob_h0 * mu_aux;
                            %
                            %                     mu_ige_aux                  = mup(:,h_p,educ_p,in,ip) * par.N(in) * par.fam_size/2;
                            %                     mu_ige(:,h_p,educ_p,ip,:,1) = squeeze(mu_ige(:,h_p,educ_p,ip,:,1)) + kron(prob_h0,mu_ige_aux);
                            %                 end
                        end
                    end
                end
            end
        end
        pop_ige = sum(mu_ige(:));
        mu_ige = mu_ige./pop_ige;
        
    case 'N'
        %% New cohort
        j_child_pos = par.Je1_pos;
        % S0c    = par.grids{j_child_pos};
        S0c      = par.grids{1,j_child_pos};
        PSY       = par.psy_val_hs;
        AGE_PROF  = par.inc.age_prof;
        
        muc      = zeros(length(S0c),length(FE_pos),length(PSY));
        for educ_p = 1:length(EDUC)
            FE          = par.inc.fe{educ_p,1};
            INNO        = par.inc.z_val{educ_p,1}(j_pos-1,:);
            PSY_p    = par.psy_prob(educ_p,:)';
            for fe_p0 = 1:length(FE_pos)
                for inno_p0 = 1:length(INNO_pos)
                    hp0     = exp(AGE_PROF{educ_p,1}(j_pos-1)) * exp(FE(fe_p0)) * exp(INNO(inno_p0));
                    prob_h0 = squeeze(par.PG(hp0));
                    
                    prob_n   = squeeze(sum(sum(sum(mup(:,fe_p0,educ_p,:,:,:,inno_p0),6),4),1));
                    ind_n    = (prob_n>0);
                    aux_n    = 1:length(N);
                    for in = aux_n(ind_n)
                        prob_p   = squeeze(sum(sum(mup(:,fe_p0,educ_p,:,in,:,inno_p0),1),4));
                        ind_p    = (prob_p>0);
                        aux_p    = 1:length(PHI);
                        for ip = aux_p(ind_p)
                            mu_aux                      = sum(sum(mup(:,fe_p0,educ_p,:,in,ip,inno_p0),1),4)*par.N(in)*par.fam_size/2;
                            muc(ip,:,:)                 = squeeze(muc(ip,:,:)) + (prob_h0 * mu_aux)*PSY_p';
                            
                            
                        end
                    end
                end
            end
        end
        mu_ige = NaN;
end

pop = sum(muc(:));
muc = muc./pop;

mup    = sum(sum(sum(mup,5),6),7);

end


%                 for inno_p = 1:NZ
%                     mu_ige_aux = mup(:,fe_p,educ_p,inno_p,in,ip)*par.N(in)*par.fam_size/2;
%                     mu_ige(:,fe_p,educ_p,inno_p,ip,:) = squeeze(mu_ige(:,fe_p,educ_p,inno_p,ip,:,:)) + kron(prob_fe',mu_ige_aux);
%                 end