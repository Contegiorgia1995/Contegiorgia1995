function [mu_ige_out,mu_chetty,Inc_par_grid,mu_educ_ige,mu_par_ic] = trans_ige(par,mu_ige0,tau0,S_pol_c,Se_pol_c,options)
% Intergenerational distributions
% keyboard
switch options.timer_on
    case {'Y'}
        fprintf('\n Intergenerational distributions \n \n');
end


% Enter with mu_ige: (Spar,Hpar,Educpar,Sc,Hc,PSY);
% Child
j_pos_c   = par.Je1_pos;
Schild    = par.grids{1,j_pos_c};
FE_pos    = par.inc.fe_pos;
INNO_pos  = par.inc.inno_pos;
PSY       = par.psy_val_hs;
EDUC      = par.educ;
r_sav     = par.r_sav;
r_debt    = par.r_debt;

% Parent
j_pos_p   = find(par.age == 42,1,'first');
Spar      = par.grids{1,j_pos_p};

switch options.ComputeOtherMus
    case 'Y'
        mu_ige_out = cell(par.Jd_pos-1,1);
        %% 1. Choice of educ of child
        % Create Q
        tau0p          = ones(length(Schild),length(FE_pos),length(PSY),length(EDUC));
        tau0p(:,:,:,1) = tau0;
        tau0int_vec    = tau0p(:);
        % create transition matrix of educ: linear interpolant for policy of education
        fspaceergeduc   = fundef({'spli',EDUC,0,1});
        Qeduc           = funbas(fspaceergeduc,tau0int_vec);
        % create transition matrix of fixed states: Qfixed
        Qfixed          = kron(ones(length(EDUC),1),speye(length(Schild)*length(FE_pos)*length(PSY)));
        % create aggregate transition: Q
        Q               = dprod(Qeduc,Qfixed);
        
        clear tau0p tau0int_vec fspaceergeduc Qeduc Qfixed
        
        mu_int2  = zeros(length(Spar),length(FE_pos),length(EDUC),length(INNO_pos),length(Schild),length(FE_pos),length(PSY),length(EDUC));
        prob_e   = squeeze(sum(sum(sum(sum(sum(sum(mu_ige0,7),6),5),4),2),1));
        ind_e    = (prob_e>0);
        aux_e    = 1:length(EDUC);
        for educ_p = aux_e(ind_e)
            prob_s   = squeeze(sum(sum(sum(sum(sum(mu_ige0(:,:,educ_p,:,:,:,:),7),6),5),4),2));
            ind_s    = (prob_s>0);
            aux_s    = 1:length(Spar);
            for s_p = aux_s(ind_s)
                prob_fe   = squeeze(sum(sum(sum(sum(mu_ige0(s_p,:,educ_p,:,:,:,:),7),6),5),4));
                ind_fe    = (prob_fe>0);
                aux_fe    = 1:length(FE_pos);
                for fe_p = aux_fe(ind_fe)
                    for inno_p = 1:length(INNO_pos)
                   
                        mu0int          = zeros(length(Schild),length(FE_pos),length(PSY),length(EDUC));
                        mu0int(:,:,:,1) = squeeze(mu_ige0(s_p,fe_p,educ_p,inno_p,:,:,:));
                        mu0int_vec      = mu0int(:);

                        mu_int_vec      = Q' * mu0int_vec;
                        mup             = reshape(mu_int_vec',1,1,1,1,length(Schild),length(FE_pos),length(PSY),length(EDUC));
                        mu_int2(s_p,fe_p,educ_p,inno_p,:,:,:,:) = mup;
                    end
                end
            end
        end
        
        % load distribution of z0 conditional on education
        j_pos      = par.Je1_pos;
        S          = par.grids{1,j_pos};
        z_prob     = par.inc.z_prob;
        z0_prob    = [squeeze(z_prob{1,1}(par.Je1_pos,1,:))'; squeeze(z_prob{2,1}(par.Je2_pos,1,:))'; squeeze(z_prob{3,1}(par.Je2_pos+2,1,:))'];      % (educ,prob)
        INNO_pos   = par.inc.inno_pos;
        NZ0        = length(INNO_pos);       
        
        % create transition matrix of Z0 conditional on education
        Qz                = kron(z0_prob,ones(length(S)*length(FE_pos)*length(PSY),1));

        % create transition matrix of fixed states: Qfixed
        Qfixed             = speye(length(S)*length(FE_pos)*length(PSY)*length(EDUC)); 
        
        % create aggregate transition: Q
        Q   = dprod(Qz,Qfixed); 

        mu_aux       = mu_int2;
        mu_int2      = zeros(length(Spar),length(FE_pos),length(EDUC),length(INNO_pos),length(Schild),length(FE_pos),length(PSY),length(EDUC),length(INNO_pos));

        for educ_par = 1:length(EDUC)
            for ife_par = 1:length(FE_pos)
                for inno_p = 1:length(INNO_pos)
                    for is_par = 1:length(Spar)
                        mu0int_aux          = squeeze(mu_aux(is_par,ife_par,educ_par,inno_p,:,:,:,:));
                        mu0int_vec          = mu0int_aux(:);
                        clear mu0int_aux
                        
                        % new distribution:
                        mu_int_vec               = Q' * mu0int_vec;
                        clear mu0int_vec
                        
                        mu_int2(is_par,ife_par,educ_par,inno_p,:,:,:,:,:)  = reshape(mu_int_vec,1,1,1,1,length(S),length(FE_pos),length(PSY),length(EDUC),length(INNO_pos));
                        clear mu_int_vec
                    end
                end
            end
        end
        
        j_c_pos               = par.Je1_pos;
        mu_ige_out{j_c_pos}   = mu_int2;
        
        clear mu_int2 Q mu_int_vec mup mu0int
        
        %% 2. Remove Psychic Cost and Savings of parent and child from distribution
        mu_ige_0    = squeeze(sum(mu_ige_out{j_c_pos},7));
        mu_ige_0    = squeeze(sum(mu_ige_0,5));
        mu_ige_0    = squeeze(sum(mu_ige_0,1));
        
        
        %% 3. Moving h_child (parent fixed)
        % Parent (fixed)
        EDUC        = par.educ;

        
        j_pos_max   = find(par.age == 36,1,'first');
        
        for j_c_pos  = par.Je1_pos+1:j_pos_max
            % Child
            j_pos_c     =	j_c_pos;
            
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
            else
                P_z(1,:,:) = par.inc.z_prob{1,1}(j_pos,:,:);
                P_z(2,:,:) = par.inc.z_prob{2,1}(j_pos,:,:);
                P_z(3,:,:) = par.inc.z_prob{3,1}(j_pos,:,:);
            end
            
            % EDUC + INNO
            P_z2                    = sparse(length(EDUC)*length(INNO_pos),length(EDUC)*length(INNO_pos));
            for iz = 1:length(INNO_pos)
                row_1               = 1 + (iz-1) * length(EDUC);
                row_2               = iz * length(EDUC);
                P_z2(row_1:row_2,:) = dprod(squeeze(P_z(:,iz,:)),speye(length(EDUC)));
            end
            QZ                  = kron(P_z2,ones(length(FE_pos)*length(EDUC)*length(INNO_pos)*length(FE_pos),1));  %Parents Characteristics, Child's FE
            
            
            % FE : Fixed
            QP   = speye(length(FE_pos)*length(EDUC)*length(INNO_pos)*length(FE_pos));
            QP   = kron(ones(length(EDUC)*length(INNO_pos),1),QP);
            
            
            % Total Q
            Q                   = dprod(QZ,QP);
            
            clear QH QE
            
            % Distribution
            mu0int_vec = mu_ige_0(:);
            mu_int_vec = Q' * mu0int_vec;
            clear mu0int_vec Q
            mup        = reshape(mu_int_vec,length(FE_pos),length(EDUC),length(INNO_pos),length(FE_pos),length(EDUC),length(INNO_pos));
            clear mu_int_vec
            
            mu_ige_0              = mup;
            clear mup
            mu_ige_out{j_c_pos}   = mu_ige_0;
        end
        clear mu_ige_0
    case 'N'
        mu_ige_out = NaN;
end

%% Create ige dist of educ
% Enter with mu_ige: (Spar,Hpar,Educpar,Sc,Hc,PSY);
% Child
j_pos_c   = par.Je1_pos;
Schild    = par.grids{1,j_pos_c};
FE_pos    = par.inc.fe_pos;
PSY       = par.psy_val_hs;
EDUC      = par.educ;
% GROUP     = par.cutoffs;

% Parent
Spar      = par.grids{1,j_pos_p};

%% 1. Choice of educ of child
% Create Q
tau0p          = ones(length(Schild),length(FE_pos),length(PSY),length(EDUC));
tau0p(:,:,:,1) = tau0;
tau0int_vec    = tau0p(:);
% create transition matrix of educ: linear interpolant for policy of education
fspaceergeduc   = fundef({'spli',EDUC,0,1});
Qeduc           = funbas(fspaceergeduc,tau0int_vec);
% create transition matrix of fixed states: Qfixed
Qfixed          = kron(ones(length(EDUC),1),speye(length(Schild)*length(FE_pos)*length(PSY)));
% create aggregate transition: Q
Q               = dprod(Qeduc,Qfixed);

clear tau0p tau0int_vec fspaceergeduc Qeduc Qfixed

mu_int2  = zeros(length(Spar),length(FE_pos),length(EDUC),length(INNO_pos),length(Schild),length(FE_pos),length(PSY),length(EDUC));
prob_e   = squeeze(sum(sum(sum(sum(sum(sum(mu_ige0,7),6),5),4),2),1));
ind_e    = (prob_e>0);
aux_e    = 1:length(EDUC);
for educ_p = aux_e(ind_e)
    prob_s   = squeeze(sum(sum(sum(sum(sum(mu_ige0(:,:,educ_p,:,:,:,:),7),6),5),4),2));
    ind_s    = (prob_s>0);
    aux_s    = 1:length(Spar);
    for s_p = aux_s(ind_s)
        prob_fe   = squeeze(sum(sum(sum(sum(mu_ige0(s_p,:,educ_p,:,:,:,:),7),6),5),4));
        ind_fe    = (prob_fe>0);
        aux_fe    = 1:length(FE_pos);
        for fe_p = aux_fe(ind_fe)
            for inno_p = 1:length(INNO_pos)

                mu0int          = zeros(length(Schild),length(FE_pos),length(PSY),length(EDUC));
                mu0int(:,:,:,1) = squeeze(mu_ige0(s_p,fe_p,educ_p,inno_p,:,:,:));
                mu0int_vec      = mu0int(:);

                mu_int_vec      = Q' * mu0int_vec;
                mup             = reshape(mu_int_vec',1,1,1,1,length(Schild),length(FE_pos),length(PSY),length(EDUC));
                mu_int2(s_p,fe_p,educ_p,inno_p,:,:,:,:) = mup;
            end
        end
    end
end

mu_educ_ige = squeeze(sum(sum(sum(sum(sum(sum(mu_int2,1),2),4),5),6),7));

%% 4. For Chetty's calculation: Fix parent (evaluated at 40-44 yo) and move child (evaluated at 28-31) with savings
% mu_chetty    = mu_ige_out{par.Je1_pos}; % But including savings!
% Choice: education + draw of innovation
j_pos     = par.Je1_pos;
Schild    = par.grids{1,j_pos};
FE_pos    = par.inc.fe_pos;
EDUC      = par.educ;
w         = par.w;

% Income parent grid:
j_par_pos = find(par.age == 42,1,'first');
EDUC_PAR   = 1:length(EDUC);
curv = 3;
Inc_par_grid  = cell(length(EDUC),1);
mu_aux     = squeeze(sum(sum(sum(sum(sum(sum(mu_ige0,7),6),5),4),2),1));
ind_ed_par = logical((mu_aux>0));

for educ = EDUC_PAR(ind_ed_par)
    S_par       = par.grids{educ,j_par_pos};
    
    INNO       = repmat(reshape(par.inc.z_val{educ,1}(j_par_pos,:),1,1,length(INNO_pos)),length(S_par),length(FE_pos));
    FE         = repmat(reshape(par.inc.fe{educ,1},1,length(FE_pos)),length(S_par),1,length(INNO_pos));
    AGE_PROF   = par.inc.age_prof{educ,1}(j_par_pos);
    Lab        = par.w*exp(AGE_PROF) .* exp(FE) .* exp(INNO);
    
    Inc         = [];
    for iphi    = 1:length(Schild)
        S_tot       = S_par + Schild(iphi);
        r           = (r_sav .* (S_tot>=0) + r_debt .* (S_tot<0))';
        Sav        = repmat(r.*S_tot',1,length(FE_pos),length(INNO_pos));

        mu_aux     = squeeze(sum(sum(sum(mu_ige0(:,:,educ,:,iphi,:,:),7),6),5));
        ind        = logical((mu_aux(:)>0));
        
        Inc_aux    = Sav + Lab;
        Inc_aux    = Inc_aux(:);
        Inc_aux    = Inc_aux(ind);
        Inc        = [Inc; Inc_aux];
    end

    aux_min               = min(0.99*Inc(:));
    if aux_min < 0
        Inc_par_grid{educ,1}  = (linspace(-abs(aux_min)^(1/curv),max(1.01*Inc(:))^(1/curv),200)).^curv;
    else
        Inc_par_grid{educ,1}  = (linspace(aux_min^(1/curv),max(1.01*Inc(:))^(1/curv),200)).^curv;
    end
end


% 1. Move parents to income grid:
length_grid = max([length(Inc_par_grid{1,1}) length(Inc_par_grid{2,1}) length(Inc_par_grid{3,1})]);
mu_chetty_int   = zeros(length(EDUC),length_grid,length(Schild),length(FE_pos),length(PSY));
mu_aux          = mu_ige0;
clear mu_ige0


for educ_par = EDUC_PAR(ind_ed_par)
    S_par       = par.grids{educ_par,j_par_pos};
    r           = (r_sav .* (S_par>=0) + r_debt .* (S_par<0));
    S_par_int   = r.*S_par;
    fspace_inc  = fundef({'spli',Inc_par_grid{educ_par,1},0,1});
    
    AGE_PROF   = par.inc.age_prof{educ_par,1}(j_par_pos);
    INNO       = repmat(reshape(par.inc.z_val{educ_par,1}(j_par_pos,:),1,length(INNO_pos)),length(FE_pos),1,length(PSY));
    FE         = repmat(reshape(par.inc.fe{educ_par,1},length(FE_pos),1),1,length(INNO_pos),length(PSY));

    prob_s   = squeeze(sum(sum(sum(sum(sum(mu_aux(:,:,educ_par,:,:,:,:),7),6),5),4),2));
    ind_s    = (prob_s>0);
    aux_s    = 1:length(S_par_int);
    
    
    for s_par = aux_s(ind_s)
        % PSYc : Fixed
        QP        = kron(speye(length(PSY)),ones(length(FE_pos)*length(INNO_pos),1));
        for fe_ch = 1:length(FE_pos)
			prob_aux  = squeeze(mu_aux(s_par,:,educ_par,:,:,fe_ch,:));
        
			prob_sc   = squeeze(sum(sum(sum(prob_aux,4),2),1));
			ind_sc    = (prob_sc>0);
			aux_sc    = 1:length(Schild);
            for s_ch = aux_sc(ind_sc)
                sav    = S_par(s_par);
                r      = (r_sav .* (sav>=0) + r_debt .* (sav<0));
                sav    = r*sav;

                lab       = par.w*exp(AGE_PROF) .* exp(FE) .* exp(INNO);
                inc       = sav + lab;

                Qinc      = funbas(fspace_inc,inc(:));
                Q         = dprod(QP,Qinc);

                prob_aux2    = squeeze(prob_aux(:,:,s_ch,:));
                prob_aux2    = prob_aux2(:);

                mu_chetty_int(educ_par,:,s_ch,fe_ch,:) = mu_chetty_int(educ_par,:,s_ch,fe_ch,:) ...
                    + reshape(Q'*prob_aux2,1,length(Inc_par_grid{educ_par,1}),1,1,length(PSY)); 

            end
        end
        clear prob_aux prob_aux2 Qinc Q QP
    end
end
mu_par_ic       = mu_chetty_int;
clear INNO lab Lab mu_aux 

%% 2. Move child
% 2.1 Choice of education
j_pos     = par.Je1_pos;
S         = par.grids{1,j_pos};
FE_pos    = par.inc.fe_pos;
INNO_pos  = par.inc.inno_pos;
PSY       = par.psy_val_col;
EDUC      = par.educ;

% Create transition Q
tau0int_vec   = tau0(:);

% create transition matrix of educ: linear interpolant for policy of education
fspaceergeduc = fundef({'spli',EDUC,0,1});
Qeduc         = funbas(fspaceergeduc,tau0int_vec);

% create transition matrix of fixed states: Qfixed
Qfixed        = speye(length(S)*length(FE_pos)*length(PSY));

% create aggregate transition: Q
Q             = dprod(Qeduc,Qfixed);
clear Qeduc Qfixed

% Distributions
mu_int2  = zeros(length(EDUC),length_grid,length(S),length(FE_pos),length(PSY),length(EDUC));

for educ_par = EDUC_PAR(ind_ed_par)
    prob_inc = squeeze(sum(sum(sum(mu_chetty_int(educ_par,:,:,:,:),5),4),3));
    ind_inc  = logical((prob_inc>0));
    aux_inc  = 1:length(Inc_par_grid{educ_par,1});
    
%     for is_par = 1:200
    for is_par = aux_inc(ind_inc)
        mu0int_aux          = squeeze(mu_chetty_int(educ_par,is_par,:,:,:));
        mu0int_vec          = mu0int_aux(:);     
        clear mu0int_aux
        
        % new distribution:
        mu_int_vec               = Q' * mu0int_vec;
        clear mu0int_vec
        
        mu_int2(educ_par,is_par,:,:,:,:)  = reshape(mu_int_vec,1,1,length(S),length(FE_pos),length(PSY),length(EDUC));
        clear mu_int_vec
    end
end
mu_chetty           = mu_int2;
clear mu_int2 mu_chetty_int Q mu0int_vec mu0int_aux tau0int_vec

% 2.2 draw of innovation z0
% load distribution of z0 conditional on education
z_prob     = par.inc.z_prob;
% z0_prob    = [squeeze(z_prob{1,1}(par.Je2_pos,1,:))'; squeeze(z_prob{2,1}(par.Je2_pos+1,1,:))'; squeeze(z_prob{3,1}(par.Je2_pos+1,1,:))'];      % (educ,prob)
z0_prob    = [squeeze(z_prob{1,1}(par.Je1_pos,1,:))'; squeeze(z_prob{2,1}(par.Je2_pos,1,:))'; squeeze(z_prob{3,1}(par.Je2_pos+2,1,:))'];      % (educ,prob)
INNO_pos   = par.inc.inno_pos;

% create transition matrix of Z0 conditional on education
Qz                = kron(z0_prob,ones(length(S)*length(FE_pos)*length(PSY),1));

% create transition matrix of fixed states: Qfixed
Qfixed              = speye(length(S)*length(FE_pos)*length(PSY)*length(EDUC)); 

% create aggregate transition: Q
Q   = dprod(Qz,Qfixed); 

mu_int2  = zeros(length(EDUC),length_grid,length(S),length(FE_pos),length(PSY),length(EDUC),length(INNO_pos));

for educ_par = EDUC_PAR(ind_ed_par)
    prob_inc = squeeze(sum(sum(sum(sum(mu_chetty(educ_par,:,:,:,:,:),6),5),4),3));
    ind_inc  = (prob_inc>0);
    aux_inc  = 1:length(Inc_par_grid{educ_par,1});
    
    for is_par = aux_inc(ind_inc)
        mu0int_aux          = squeeze(mu_chetty(educ_par,is_par,:,:,:,:));
        mu0int_vec          = mu0int_aux(:);     
        clear mu0int_aux
       
        mup_vec = Q' * mu0int_vec;
        mu_int2(educ_par,is_par,:,:,:,:,:) = reshape(mup_vec,1,1,length(S),length(FE_pos),length(PSY),length(EDUC),length(INNO_pos));

    end
end
% Remove Psychic cost
mu_chetty   = squeeze(sum(mu_int2,5));
clear mu_int2 Q Qz Qfixed


% 2.2 Transition until Jc
for j_c_pos  = par.Je1_pos+1:par.Jc_pos
    % Child
    j_pos_c     =	j_c_pos;
    Schild      = par.grids{1,j_pos_c-1};
    if j_pos_c-1 <= par.Je1_pos + 1 % High School
        % Saving policy function
        Spol            = zeros(length(Schild),length(FE_pos),length(EDUC),length(INNO_pos));
        Spol(:,:,1,:)   = Se_pol_c{1,j_pos_c-par.Je1_pos};
        Spol(:,:,2,:)   = repmat(squeeze(Se_pol_c{2,j_pos_c-par.Je1_pos}(:,:,1)),1,1,1,length(INNO_pos));
        Spol(:,:,3,:)   = repmat(squeeze(Se_pol_c{3,j_pos_c-par.Je1_pos}(:,:,1)),1,1,1,length(INNO_pos));
        
        P_z        = zeros(length(EDUC),length(INNO_pos),length(INNO_pos));
        P_z(1,:,:) = par.inc.z_prob{1,1}(j_pos_c,:,:);
        P_z(2,:,:) = speye(length(INNO_pos));
        P_z(3,:,:) = speye(length(INNO_pos));
               
    elseif j_pos_c-1 <= par.Je2_pos + 1
        % Saving policy function
        Spol            = zeros(length(Schild),length(FE_pos),length(EDUC),length(INNO_pos));
        Spol(:,:,1,:)   = Se_pol_c{1,j_pos_c-par.Je1_pos};
        Spol(:,:,2,:)   = Se_pol_c{2,j_pos_c-par.Je1_pos};
        Spol(:,:,3,:)   = repmat(squeeze(Se_pol_c{3,j_pos_c-par.Je1_pos}(:,:)),1,1,1,length(INNO_pos));
        
        P_z(1,:,:) = par.inc.z_prob{1,1}(j_pos_c,:,:);
        P_z(2,:,:) = par.inc.z_prob{2,1}(j_pos_c,:,:);
        P_z(3,:,:) = speye(length(INNO_pos));
        
    else % Working
        % Saving policy function ---> Extend as previous ones, for below
        Spol = S_pol_c{j_c_pos-1};
        
        P_z(1,:,:) = par.inc.z_prob{1,1}(j_pos_c,:,:);
        P_z(2,:,:) = par.inc.z_prob{2,1}(j_pos_c,:,:);
        P_z(3,:,:) = par.inc.z_prob{3,1}(j_pos_c,:,:);
    end
    
    % new distribution:
    jp_pos_c      =	j_c_pos + 1;
    Spchild       = par.grids{1,jp_pos_c-1};
    mu_chetty1    = zeros(length(EDUC),length_grid,length(Spchild),length(FE_pos),length(EDUC),length(INNO_pos));
    
    % FE: Fixed
    QP                  = kron(speye(length(FE_pos)),ones(length(Schild),1));
    QP                  = kron(ones(length(EDUC)*length(INNO_pos),1),QP);

    
    % EDUC + INNO
    P_z2                    = sparse(length(EDUC)*length(INNO_pos),length(EDUC)*length(INNO_pos));
    for iz = 1:length(INNO_pos)
        row_1               = 1 + (iz-1) * length(EDUC);
        row_2               = iz * length(EDUC);
        P_z2(row_1:row_2,:) = dprod(squeeze(P_z(:,iz,:)),speye(length(EDUC)));
    end
    QZ                  = kron(P_z2,ones(length(Schild)*length(FE_pos),1));

    
    % Savings:
    Sp              = par.grids{1,j_c_pos}';
    QS              = sparse(length(Schild)*length(FE_pos)*length(EDUC)*length(INNO_pos),length(Sp));
    pos             = 1;
    npos            = length(Schild)*length(FE_pos);
    for inno = 1:length(INNO_pos)
        for educ = 1:length(EDUC)
            % vectorize policy

            Spol_vec                = Spol(:,:,educ,inno);
            Spol_vec                = Spol_vec(:);
            Sp                      = par.grids{educ,j_c_pos}';
            % create transition matrix of savings: linear interpolant for policy
            fspace_s                = fundef({'spli',Sp,0,1});
            QS(pos:pos+npos-1,:)    = funbas(fspace_s,Spol_vec);
            pos                     = pos + npos ;
        end
    end
    clear Spol_vec

    
    % Total Q
    Q                   = dprod(QP,QS);
    clear QP QS
    Q                   = dprod(QZ,Q);
    clear QZ
%     if min(Q(:))<0
%         fprintf('Stop');
%         keyboard
%     end

    for educ_par = EDUC_PAR(ind_ed_par)
        prob_inc = squeeze(sum(sum(sum(sum(mu_chetty(educ_par,:,:,:,:,:),6),5),4),3));
        ind_inc  = (prob_inc>0);
        aux_inc  = 1:length(Inc_par_grid{educ_par,1});

        for is_par = aux_inc(ind_inc)
            mu0int_aux          = squeeze(mu_chetty(educ_par,is_par,:,:,:,:));
            mu0int_vec          = mu0int_aux(:);     
            clear mu0int_aux

            mup_vec = Q' * mu0int_vec;
            mu_chetty1(educ_par,is_par,:,:,:,:) = reshape(mup_vec,1,1,length(Spchild),length(FE_pos),length(EDUC),length(INNO_pos));
        end
    end
    mu_chetty                       = mu_chetty1;
end
clear mu_chetty1


switch options.timer_on
    case {'Y'}
        fprintf('ige, time: %3.1f sec \n',toc)
end

end