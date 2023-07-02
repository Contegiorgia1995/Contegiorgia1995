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
Hchild    = par.gridh{2,j_pos_c};
PSY       = par.psy_val_hs;
EDUC      = par.educ;

% Parent
j_pos_p   = find(par.age == 40,1,'first');
Spar      = par.grids{1,j_pos_p};
Hpar      = par.gridh{2,j_pos_p};

switch options.ComputeOtherMus
    case 'Y'
        mu_ige_out = cell(par.Jd_pos-1,1);
        %% 1. Choice of educ of child
        % Create Q
        tau0p          = ones(length(Schild),length(Hchild),length(PSY),length(EDUC));
        tau0p(:,:,:,1) = tau0;
        tau0int_vec    = tau0p(:);
        % create transition matrix of educ: linear interpolant for policy of education
        fspaceergeduc   = fundef({'spli',EDUC,0,1});
        Qeduc           = funbas(fspaceergeduc,tau0int_vec);
        % create transition matrix of fixed states: Qfixed
        Qfixed          = kron(ones(length(EDUC),1),speye(length(Schild)*length(Hchild)*length(PSY)));
        % create aggregate transition: Q
        Q               = dprod(Qeduc,Qfixed);

        clear tau0p tau0int_vec fspaceergeduc Qeduc Qfixed
        
        mu_int2  = zeros(length(Spar),length(Hpar),length(EDUC),length(Schild),length(Hchild),length(PSY),length(EDUC));
        prob_e   = squeeze(sum(sum(sum(sum(sum(mu_ige0,6),5),4),2),1));
        ind_e    = (prob_e>0);
        aux_e    = 1:length(EDUC);
        for educ_p = aux_e(ind_e)
            prob_s   = squeeze(sum(sum(sum(sum(mu_ige0(:,:,educ_p,:,:,:),6),5),4),2));
            ind_s    = (prob_s>0);
            aux_s    = 1:length(Spar);
            for s_p = aux_s(ind_s);
                prob_h   = squeeze(sum(sum(sum(mu_ige0(s_p,:,educ_p,:,:,:),6),5),4));
                ind_h    = (prob_h>0);
                aux_h    = 1:length(Hpar);
                for h_p = aux_h(ind_h)
                
                    mu0int          = zeros(length(Schild),length(Hchild),length(PSY),length(EDUC));
                    mu0int(:,:,:,1) = squeeze(mu_ige0(s_p,h_p,educ_p,:,:,:));
        %             mu0int_vec      = reshape(mu0int,length(Schild)*length(Hchild)*length(PSY)*length(EDUC));
        %             mu0int_vec      = mu0int_vec';
                    mu0int_vec      = mu0int(:);

                    mu_int_vec      = Q' * mu0int_vec;
                    mup             = reshape(mu_int_vec',1,1,1,length(Schild),length(Hchild),length(PSY),length(EDUC));
                    mu_int2(s_p,h_p,educ_p,:,:,:,:) = mup;
                end
            end
        end


        j_c_pos               = par.Je1_pos;
        mu_ige_out{j_c_pos}   = mu_int2;
        
        clear mu_int2 Q mu_int_vec mup mu0int

        %% 2. Remove Psychic Cost and Savings of parent and child from distribution
        mu_ige_0    = squeeze(sum(mu_ige_out{j_c_pos},6));
        mu_ige_0    = squeeze(sum(mu_ige_0,4));
        mu_ige_0    = squeeze(sum(mu_ige_0,1));


        %% 3. Moving h_child (parent fixed)
        % Parent (fixed)
        Hpar        = par.gridh{1,j_pos_p};
        EDUC        = par.educ;
        Ndeltas     = length(par.deltas{1,1}(j_pos_p,:));
        
        j_pos_max   = find(par.age == 36,1,'first');

        for j_c_pos  = par.Je1_pos+1:j_pos_max
            % Child
            j_pos_c     =	j_c_pos;
            Hchild      = par.gridh{1,j_pos_c};
            H0child     = par.gridh{1,j_pos_c-1};
            if j_pos_c-1 <= par.Je1_pos+1 % First year of education
                % H policy function
                Hpol            = zeros(length(H0child),Ndeltas,length(EDUC));
                Hpol_pr         = zeros(length(H0child),Ndeltas,length(EDUC));

                Hpol(:,:,1)     = par.gridh{1,j_pos_c-1}' * par.deltas{1,1}(j_pos_c-1,:);
                Hpol_pr(:,:,1)  = repmat(par.prob{1,1}(j_pos_c-1,:),length(H0child),1);

                Haux            = repmat(par.gridh{2,j_pos_c-1}',1,Ndeltas);
                Hpol(:,:,2)     = Haux + par.Reducs(1,1) * Haux.^par.Reducs(1,2);
                Hpol_pr(:,:,2)  = repmat([1 0 0],length(H0child),1);

                Haux            = repmat(par.gridh{3,j_pos_c-1}',1,Ndeltas);                
                Hpol(:,:,3)     = Haux + par.Reducs(1,1) * Haux.^par.Reducs(1,2);
                Hpol_pr(:,:,3)  = repmat([1 0 0],length(H0child),1);

            elseif j_pos_c-1 <= par.Je2_pos + 1
                % H policy function
                Hpol            = zeros(length(H0child),Ndeltas,length(EDUC));
                Hpol_pr         = zeros(length(H0child),Ndeltas,length(EDUC));

                Hpol(:,:,1)     = par.gridh{1,j_pos_c-1}' * par.deltas{1,1}(j_pos_c-1,:);
                Hpol_pr(:,:,1)  = repmat(par.prob{1,1}(j_pos_c-1,:),length(H0child),1);

                Hpol(:,:,2)     = par.gridh{2,j_pos_c-1}' * par.deltas{2,1}(j_pos_c-1,:);
                Hpol_pr(:,:,2)  = repmat(par.prob{2,1}(j_pos_c-1,:),length(H0child),1);

                Haux            = repmat(par.gridh{3,j_pos_c-1}',1,Ndeltas);                
                Hpol(:,:,3)     = Haux + par.Reducs(2,1)* Haux.^par.Reducs(2,2);
                Hpol_pr(:,:,3)  = repmat([1 0 0],length(H0child),1);

            else % Working
                % H policy function
                Hpol            = zeros(length(H0child),Ndeltas,length(EDUC));
                Hpol_pr         = zeros(length(H0child),Ndeltas,length(EDUC));

                for educ = 1:length(EDUC)
                    Hpol(:,:,educ)     = par.gridh{educ,j_pos_c-1}' * par.deltas{educ,1}(j_pos_c-1,:);
                    Hpol_pr(:,:,educ)  = repmat(par.prob{educ,1}(j_pos_c-1,:),length(H0child),1);
                end
            end

            Hp              = par.gridh{1,j_pos_c}';
            QH              = sparse(length(H0child)*length(EDUC),length(Hp));
            for delta = 1:Ndeltas
                Qh_aux          = sparse(length(H0child)*length(EDUC),length(Hp));
                pos             = 1;
                npos            = length(H0child);
                for educ = 1:length(EDUC)
                    % vectorize policy
                    Hpol_vec                = Hpol(:,delta,educ);
                    Hpol_vec                = Hpol_vec(:);
                    Hp                      = par.gridh{educ,j_pos_c}';
                    % create transition matrix of savings: linear interpolant for policy
                    fspace_h                = fundef({'spli',Hp,0,1});
                    Qh_aux(pos:pos+npos-1,:)= funbas(fspace_h,Hpol_vec);
                    pos                     = pos + npos ;
                end
                Hpol_pr_vec     = squeeze(Hpol_pr(:,delta,:));
                Hpol_pr_vec     = Hpol_pr_vec(:);
                QH              = QH + Qh_aux .* repmat(Hpol_pr_vec,1,size(Qh_aux,2));
            end
            QH                  = kron(QH,ones(length(Hpar)*length(EDUC),1));  
            clear Hpol*
            
            % Fixed: EDUCchild, Hpar and EDUCpar
            QE                  = kron(ones(length(H0child),1),speye(length(EDUC)*length(Hpar)*length(EDUC)));

            % Total Q
            Q                   = dprod(QH,QE);

            clear QH QE
            
            % Distribution
            mu0int_vec = mu_ige_0(:);
            mu_int_vec = Q' * mu0int_vec;
            clear mu0int_vec Q
            mup        = reshape(mu_int_vec,length(Hpar),length(EDUC),length(Hchild),length(EDUC));
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
Hchild    = par.gridh{2,j_pos_c};
PSY       = par.psy_val_hs;
EDUC      = par.educ;
% GROUP     = par.cutoffs;

% Parent
Spar      = par.grids{1,j_pos_p};
Hpar      = par.gridh{2,j_pos_p};

% Choice of educ of child
% Create Q
tau0p          = ones(length(Schild),length(Hchild),length(PSY),length(EDUC));
tau0p(:,:,:,1) = tau0;
tau0int_vec    = tau0p(:);
% create transition matrix of educ: linear interpolant for policy of education
fspaceergeduc   = fundef({'spli',EDUC,0,1});
Qeduc           = funbas(fspaceergeduc,tau0int_vec);
% create transition matrix of fixed states: Qfixed
Qfixed          = kron(ones(length(EDUC),1),speye(length(Schild)*length(Hchild)*length(PSY)));
% create aggregate transition: Q
Q               = dprod(Qeduc,Qfixed);

clear tau0p tau0int_vec fspaceergeduc Qeduc Qfixed

mu_int2  = zeros(length(Spar),length(Hpar),length(EDUC),length(Schild),length(Hchild),length(PSY),length(EDUC));
prob_e   = squeeze(sum(sum(sum(sum(sum(mu_ige0,6),5),4),2),1));
ind_e    = (prob_e>0);
aux_e    = 1:length(EDUC);
for educ_p = aux_e(ind_e)
    prob_s   = squeeze(sum(sum(sum(sum(mu_ige0(:,:,educ_p,:,:,:),6),5),4),2));
    ind_s    = (prob_s>0);
    aux_s    = 1:length(Spar);
    for s_p = aux_s(ind_s);
        prob_h   = squeeze(sum(sum(sum(mu_ige0(s_p,:,educ_p,:,:,:),6),5),4));
        ind_h    = (prob_h>0);
        aux_h    = 1:length(Hpar);
        for h_p = aux_h(ind_h)

            mu0int          = zeros(length(Schild),length(Hchild),length(PSY),length(EDUC));
            mu0int(:,:,:,1) = squeeze(mu_ige0(s_p,h_p,educ_p,:,:,:));
%             mu0int_vec      = reshape(mu0int,length(Schild)*length(Hchild)*length(PSY)*length(EDUC));
%             mu0int_vec      = mu0int_vec';
            mu0int_vec      = mu0int(:);

            mu_int_vec      = Q' * mu0int_vec;
            mup             = reshape(mu_int_vec',1,1,1,length(Schild),length(Hchild),length(PSY),length(EDUC));
            mu_int2(s_p,h_p,educ_p,:,:,:,:) = mup;
        end
    end
end

mu_educ_ige = squeeze(sum(sum(sum(sum(sum(mu_int2,1),2),4),5),6));

%% 4. For Chetty's calculation: Fix parent (evaluated at 40-44 yo) and move child (evaluated at 28-31) with savings
% mu_chetty    = mu_ige_out{par.Je1_pos}; % But including savings!
% Choice: education + draw of innovation
j_pos     = par.Je1_pos;
Schild    = par.grids{1,j_pos};
Hchild    = par.gridh{1,j_pos};
EDUC      = par.educ;
w         = par.w;

% Income parent grid:
j_par_pos = find(par.age == 40,1,'first');
H_par     = par.gridh{1,j_par_pos};
min_inc   = w * min(H_par);
max_inc   = 0;
for educ = 1:length(EDUC)
    S_par       = par.grids{educ,j_par_pos};
    H_par       = par.gridh{educ,j_par_pos};
    min_inc     = min(min_inc, par.r_debt * min(S_par)+ w * min(H_par));
    max_inc     = max(max_inc, par.r_sav * max(S_par) + w * max(H_par));
end
curv = 2;
Inc_par_grid  = (linspace(-abs(min_inc)^(1/curv),max_inc^(1/curv),500)).^curv;
fspace_inc    = fundef({'spli',Inc_par_grid,0,1});

% 1. Move parents to income grid:
mu_chetty_int   = zeros(length(EDUC),length(Inc_par_grid),length(Schild),length(Hchild),length(PSY));
mu_aux          = mu_ige0;
clear mu_ige0

prob_e   = squeeze(sum(sum(sum(sum(sum(mu_aux,6),5),4),2),1));
ind_e    = (prob_e>0);
aux_e    = 1:length(EDUC);
for educ_par = aux_e(ind_e)
    prob_s   = squeeze(sum(sum(sum(sum(mu_aux(:,:,educ_par,:,:,:),6),5),4),2));
    ind_s    = (prob_s>0);
    aux_s    = 1:length(S_par);
    for s_par = aux_s(ind_s);
        S_par       = par.grids{educ_par,j_par_pos};
        H_par       = par.gridh{educ_par,j_par_pos};
        prob_aux    = squeeze(mu_aux(s_par,:,educ_par,:,:,:));
        
        prob_sc   = squeeze(sum(sum(sum(prob_aux,4),3),1));
        ind_sc    = (prob_sc>0);
        aux_sc    = 1:length(Schild);
        
        for s_ch = aux_sc(ind_sc)
            prob_aux2    = reshape(prob_aux(:,s_ch,:,:),length(Hpar),length(Hchild)*length(PSY));
            prob_aux2    = prob_aux2';

			r     		= (par.r_sav .* (S_par(s_par)>=0) + par.r_debt .* (S_par(s_par)<0))';
            inc         = r * S_par(s_par) + w * H_par;
            Qinc        = funbas(fspace_inc,inc');
    %         Qinc        = repmat(full(Qinc)',1,length(Schild),length(Hchild));


            mu_chetty_int(educ_par,:,s_ch,:,:) = squeeze(mu_chetty_int(educ_par,:,s_ch,:,:)) ...
                + reshape((prob_aux2 * Qinc)',length(Inc_par_grid),length(Hchild),length(PSY));        
        end
        clear prob_aux prob_aux2 Qinc
%         mu_chetty_int(educ_par,:,:,:) = squeeze(mu_chetty_int(educ_par,:,:,:)) ...
%             + reshape(full(Qinc)' * prob_aux',length(Inc_par_grid),length(Schild),length(Hchild));        
%             + Qinc' .* repmat(reshape(prob_aux',length(Hpar),length(Schild),length(Hchild))...
%             ,length(Inc_par_grid),1,1);
    end
end

mu_par_ic       = mu_chetty_int;
% mu_par_ic(educ_par,inc_par,S0_c,H0_c,psy_c);

mu_chetty_int  = squeeze(sum(mu_chetty_int,1));


% 2. Move child
% 2.1 Choice of education
j_pos     = par.Je1_pos;
S         = par.grids{1,j_pos};
H         = par.gridh{1,j_pos};
EDUC      = par.educ;

% Create transition Q
% extend policy for education tau
tau0p         = ones(length(S),length(H),length(PSY),length(EDUC));
tau0p(:,:,:,1)  = tau0;
tau0int_vec   = tau0p(:);
clear tau0p

% create transition matrix of educ: linear interpolant for policy of education
fspaceergeduc = fundef({'spli',EDUC,0,1});
Qeduc         = funbas(fspaceergeduc,tau0int_vec);

% create transition matrix of fixed states: Qfixed
Qfixed        = kron(ones(length(EDUC),1),speye(length(S)*length(H)*length(PSY)));

% create aggregate transition: Q
Q             = dprod(Qeduc,Qfixed);
clear Qeduc Qfixed

% Distributions
mu0int        = zeros(length(Inc_par_grid),length(S),length(H),length(PSY),length(EDUC));
mu_int2       = zeros(length(Inc_par_grid),length(S),length(H),length(PSY),length(EDUC));

prob_inc = squeeze(sum(sum(sum(mu_chetty_int,4),3),2));
ind_inc  = (prob_inc>0);
aux_inc  = 1:length(Inc_par_grid);
for is_par = aux_inc(ind_inc)
    mu0int_aux          = squeeze(mu0int(is_par,:,:,:,:));
    mu0int_aux(:,:,:,1) = mu_chetty_int(is_par,:,:,:);
    mu0int_vec          = mu0int_aux(:);     
    clear mu0int_aux
    % new distribution:
    mu_int_vec               = Q' * mu0int_vec;
    clear mu0int_vec
    mu_int2(is_par,:,:,:,:)  = reshape(mu_int_vec,1,length(S),length(H),length(PSY),length(EDUC));
    clear mu_int_vec
end

mu_chetty           = mu_int2;
clear mu_int2
j_c_pos  = par.Je1_pos;
Ndeltas  = length(par.deltas{1,1}(j_c_pos,:));

j_c_pos    = par.Je1_pos+1;
jp_pos_c   = j_c_pos + 1;
Spchild    = par.grids{1,jp_pos_c-1};
Hpchild    = par.gridh{1,jp_pos_c-1};
% mu_chetty1 = zeros(length(Inc_par_grid),length(Spchild),length(Hpchild),length(PSY),length(EDUC)); 

% 2.2 Transition until Jc
for j_c_pos  = par.Je1_pos+1:par.Jc_pos
    % Child
    j_pos_c     =	j_c_pos;
    Schild      = par.grids{1,j_pos_c-1};    
    Hchild      = par.gridh{1,j_pos_c-1};
    if j_pos_c-1 < par.Je1_pos + 1 % First year of education
        % Saving policy function
        Spol            = zeros(length(Schild),length(Hchild),length(PSY),length(EDUC));
        Spol(:,:,:,1)   = repmat(Se_pol_c{1,j_pos_c-par.Je1_pos},1,1,length(PSY));
        Spol(:,:,:,2)   = Se_pol_c{2,j_pos_c-par.Je1_pos};
        Spol(:,:,:,3)   = Se_pol_c{3,j_pos_c-par.Je1_pos};

        % H policy function
        Hpol            = zeros(length(Hchild),Ndeltas,length(EDUC));
        Hpol_pr         = zeros(length(Hchild),Ndeltas,length(EDUC));
        
        Hpol(:,:,1)     = par.gridh{1,j_pos_c-1}' * par.deltas{1,1}(j_pos_c-1,:);
        Hpol_pr(:,:,1)  = repmat(par.prob{1,1}(j_pos_c-1,:),length(Hchild),1);
        
        Haux            = repmat(par.gridh{2,j_pos_c-1}',1,Ndeltas);                
        Hpol(:,:,2)     = Haux + par.Reducs(1,1) * Haux.^par.Reducs(1,2);
        Hpol_pr(:,:,2)  = repmat([1 0 0],length(Hchild),1);
        
        Haux            = repmat(par.gridh{3,j_pos_c-1}',1,Ndeltas);                
        Hpol(:,:,3)     = Haux + par.Reducs(1,1)* Haux.^par.Reducs(1,2);
        Hpol_pr(:,:,3)  = repmat([1 0 0],length(Hchild),1);
        
    elseif j_pos_c-1 <= par.Je2_pos + 1
        % Saving policy function
        Spol            = zeros(length(Schild),length(Hchild),length(EDUC));
        Spol(:,:,1)   = Se_pol_c{1,j_pos_c-par.Je1_pos};
        Spol(:,:,2)   = Se_pol_c{2,j_pos_c-par.Je1_pos};
        Spol(:,:,3)   = Se_pol_c{3,j_pos_c-par.Je1_pos};

        % H policy function
        Hpol            = zeros(length(Hchild),Ndeltas,length(EDUC));
        Hpol_pr         = zeros(length(Hchild),Ndeltas,length(EDUC));
        
        Hpol(:,:,1)     = par.gridh{1,j_pos_c-1}' * par.deltas{1,1}(j_pos_c-1,:);
        Hpol_pr(:,:,1)  = repmat(par.prob{1,1}(j_pos_c-1,:),length(Hchild),1);
        
        Hpol(:,:,2)     = par.gridh{2,j_pos_c-1}' * par.deltas{2,1}(j_pos_c-1,:);
        Hpol_pr(:,:,2)  = repmat(par.prob{2,1}(j_pos_c-1,:),length(Hchild),1);
        
        Haux            = repmat(par.gridh{3,j_pos_c-1}',1,Ndeltas);                
        Hpol(:,:,3)     = Haux + par.Reducs(2,1) * Haux.^par.Reducs(2,2);
        Hpol_pr(:,:,3)  = repmat([1 0 0],length(Hchild),1);
        
    else % Working
        % Saving policy function ---> Extend as previous ones, for below
        Spol = S_pol_c{j_c_pos-1};

        % H policy function
        Hpol            = zeros(length(Hchild),Ndeltas,length(EDUC));
        Hpol_pr         = zeros(length(Hchild),Ndeltas,length(EDUC));
        
        for educ = 1:length(EDUC)
            Hpol(:,:,educ)     = par.gridh{educ,j_pos_c-1}' * par.deltas{educ,1}(j_pos_c-1,:);
            Hpol_pr(:,:,educ)  = repmat(par.prob{educ,1}(j_pos_c-1,:),length(Hchild),1);
        end
    end
    
    % new distribution:
    jp_pos_c      =	j_c_pos + 1;
    Spchild       = par.grids{1,jp_pos_c-1};
    Hpchild       = par.gridh{1,jp_pos_c-1};
    mu_chetty1    = zeros(length(Inc_par_grid),length(Spchild),length(Hpchild),length(PSY),length(EDUC));
    
    prob_psy = squeeze(sum(sum(sum(sum(mu_chetty,5),3),2),1));
    ind_psy  = (prob_psy>0);
    aux_psy  = 1:length(PSY);
    ipsy_one   = aux_psy(find(ind_psy == 1,1,'first'));
    for ipsy  = aux_psy(ind_psy)
%         if (j_pos_c-1 <= par.Je1_pos) ||  (j_pos_c-1 > par.Je1_pos && ipsy == ipsy_one )
        if ipsy == ipsy_one 
            % Savings:
            Sp              = par.grids{1,j_c_pos}';
            QS              = sparse(length(Schild)*length(Hchild)*length(EDUC),length(Sp));
            pos             = 1;
            npos            = length(Schild)*length(Hchild);
            for educ = 1:length(EDUC)
                % vectorize policy
                if j_pos_c-1 <= par.Je1_pos
                    Spol_vec                = Spol(:,:,ipsy,educ);
                else
                    Spol_vec                = Spol(:,:,educ);
                end
                Spol_vec                = Spol_vec(:);
                Sp                      = par.grids{educ,j_c_pos}';
                % create transition matrix of savings: linear interpolant for policy
                fspace_s                = fundef({'spli',Sp,0,1});
                QS(pos:pos+npos-1,:)    = funbas(fspace_s,Spol_vec);
                pos                     = pos + npos ;
            end
            clear Spol_vec

            % H:
            Hp              = par.gridh{1,j_pos_c}';
            QH              = sparse(length(Hchild)*length(EDUC),length(Hp)); 

            for delta = 1:Ndeltas
                Qh_aux          = sparse(length(Hchild)*length(EDUC),length(Hp));
                pos             = 1;
                npos            = length(Hchild);
                for educ = 1:length(EDUC)
                    % vectorize policy
                    Hpol_vec                = Hpol(:,delta,educ);
                    Hpol_vec                = Hpol_vec(:);
                    Hp                      = par.gridh{educ,j_pos_c}';
                    % create transition matrix of savings: linear interpolant for policy
                    fspace_h                = fundef({'spli',Hp,0,1});
                    Qh_aux(pos:pos+npos-1,:)= funbas(fspace_h,Hpol_vec);
                    pos                     = pos + npos ;
                end
                Hpol_pr_vec     = squeeze(Hpol_pr(:,delta,:));
                Hpol_pr_vec     = Hpol_pr_vec(:);
                QH              = QH + Qh_aux .* repmat(Hpol_pr_vec,1,size(Qh_aux,2));
            end
            QH                  = kron(QH,ones(length(Schild),1));
            clear Hpol_vec

            % Fixed: EDUCchild
            QE                  = kron(speye(length(EDUC)),ones(length(Schild)*length(Hchild),1));

            % Total Q
            Q                   = dprod(QH,QS);
            clear QH QS
            Q                   = dprod(QE,Q);
            clear QE
        end
   
        prob_inc = squeeze(sum(sum(sum(mu_chetty(:,:,:,ipsy,:),5),3),2));
        ind_inc  = (prob_inc>0);
        aux_inc  = 1:length(Inc_par_grid);
        for is_par = aux_inc(ind_inc)
            mu0int_aux                  = squeeze(mu_chetty(is_par,:,:,ipsy,:));
            mu0int_vec                  = mu0int_aux(:);
            clear mu0int_aux
            mu_int_vec                  = Q' * mu0int_vec;
            clear mu0int_vec
            mu_chetty1(is_par,:,:,ipsy,:)  = reshape(mu_int_vec,1,length(Spchild),length(Hpchild),1,length(EDUC));
        end      
    end
    mu_chetty                       = mu_chetty1;
end
clear mu_chetty1
% Remove Psychic cost
mu_chetty   = squeeze(sum(mu_chetty,4));

switch options.timer_on
    case {'Y'}
        fprintf('ige, time: %3.1f sec \n',toc)
end

end