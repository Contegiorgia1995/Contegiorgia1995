# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 13:49:56 2022

@author: Giorgia
"""

function [pol_dense,par,par2] = solve_ergodic_distribution1(par,pol,options)
% Interpolate policy functions into dense grid

%% Distribute policy functions
S           = pol.S;
Se          = pol.Se;
C           = pol.C;
Ce          = pol.Ce;
Ck          = pol.Ck;
Np          = pol.Np;
PHIp        = pol.PHIp;
Tp          = pol.Tp;
tau0        = pol.tau0;

clear pol


%% Interpolate policy functions in a dense grid
% Parameters
EDUC                            = par.educ;
N                               = par.N;
PSY                             = par.psy_val_hs;
FE_pos                          = par.inc.fe_pos;
INNO_pos                        = par.inc.inno_pos;

% Arrays for new policies
S_dense                         = cell(par.Jd_pos,1);
Se_dense                        = cell(3,3);
C_dense                         = cell(par.Jd_pos,1);
Ck_dense                        = cell(par.Jd_pos,1);
Ce_dense                        = cell(3,3);
Tp_dense                        = cell(par.Jd_pos,1);

%% 1. Create dense grid
curv                            = 3;
par.grids_dense                 = cell(3,par.NJ,1);
for educ = 1:length(EDUC)
    par.grids_dense{educ,par.Je1_pos} = par.grids{educ,par.Je1_pos};
end
n_2_s                           = 100;
for j_pos = par.Je1_pos+1:par.Jd_pos
    for educ = 1:length(EDUC)
        S0                          = par.grids{educ,j_pos};
        s_min                       = abs(min(S0));
        s_max                       = max(S0);
        par.grids_dense{educ,j_pos} = linspace(-s_min^(1/curv),s_max^(1/curv),n_2_s).^curv;
    end
end

%% 2. Interpolate policy function of savings for education
educ  = 1;
for t = 1:3
    S_dense_aux      = par.grids_dense{educ,par.Je1_pos+t-1};
    S_orig           = par.grids{educ,par.Je1_pos+t-1};
    Sp_max           = max(par.grids_dense{educ,par.Je1_pos+t});
    
    
    S_pol            = Se{educ,t};
    C_pol            = Ce{educ,t};
    if max(S_pol(:)) > Sp_max
        j_pos        = par.Je1_pos+t-1;
        fprintf('Attention: age = %i, educ = %i, extrapolation in ergodic distribution, max extrap = %3.2f pp  \n',par.age(j_pos),educ,log(max(S_pol(:))/Sp_max))
        S_pol        = min(S_pol,Sp_max);
    end
    
    Sp_min           = min(par.grids_dense{educ,par.Je1_pos+t});
    if min(S_pol(:)) < Sp_min
        j_pos        = par.Je1_pos+t-1;
        fprintf('Attention: age = %i, educ = %i, extrapolation in ergodic distribution, min extrap = %3.2f pp  \n',par.age(j_pos),educ,log(min(S_pol(:))/Sp_min))
        S_pol        = max(S_pol,Sp_min);
    end
    
    [S1,FE1,INNO1]   = ndgrid(S_orig,FE_pos,INNO_pos);
    [S2,FE2,INNO2]   = ndgrid(S_dense_aux,FE_pos,INNO_pos);
    
    pol                     = griddedInterpolant(S1,FE1,INNO1,squeeze(S_pol));
    Se_dense{educ,t}(:,:,:) = reshape(pol(S2,FE2,INNO2),length(S_dense_aux),length(FE_pos),length(INNO_pos));
    
    pol                     = griddedInterpolant(S1,FE1,INNO1,squeeze(C_pol));
    Ce_dense{educ,t}(:,:,:) = reshape(pol(S2,FE2,INNO2),length(S_dense_aux),length(FE_pos),length(INNO_pos));
end

educ  = 2;
for t = 1:3
    S_dense_aux      = par.grids_dense{educ,par.Je1_pos+t-1};
    S_orig           = par.grids{educ,par.Je1_pos+t-1};
    Sp_max           = max(par.grids_dense{educ,par.Je1_pos+t});
    
    S_pol            = Se{educ,t};
    C_pol            = Ce{educ,t};
    if max(S_pol(:)) > Sp_max
        j_pos        = par.Je1_pos+t-1;
        fprintf('Attention: age = %i, educ = %i, extrapolation in ergodic distribution, max extrap = %3.2f pp  \n',par.age(j_pos),educ,log(max(S_pol(:))/Sp_max))
        S_pol        = min(S_pol,Sp_max);
    end
    
    Sp_min           = min(par.grids_dense{educ,par.Je1_pos+t});
    if min(S_pol(:)) < Sp_min
        j_pos        = par.Je1_pos+t-1;
        fprintf('Attention: age = %i, educ = %i, extrapolation in ergodic distribution, min extrap = %3.2f pp  \n',par.age(j_pos),educ,log(min(S_pol(:))/Sp_min))
        S_pol        = max(S_pol,Sp_min);
    end
    
    if t <= 1
        psy_orig         = PSY;
        [S1,FE1,PSY1]    = ndgrid(S_orig,FE_pos,psy_orig);
        [S2,FE2,PSY2]    = ndgrid(S_dense_aux,FE_pos,psy_orig);
        
        pol                 = griddedInterpolant(S1,FE1,PSY1,squeeze(S_pol));
        Se_dense{educ,t}(:,:,:)    = reshape(pol(S2,FE2,PSY2), length(S_dense_aux),length(FE_pos),length(psy_orig));
        
        pol                 = griddedInterpolant(S1,FE1,PSY1,squeeze(C_pol));
        Ce_dense{educ,t}(:,:,:)    = reshape(pol(S2,FE2,PSY2), length(S_dense_aux),length(FE_pos),length(psy_orig));
        
    else
        [S1,FE1,INNO1]   = ndgrid(S_orig,FE_pos,INNO_pos);
        [S2,FE2,INNO2]   = ndgrid(S_dense_aux,FE_pos,INNO_pos);
        
        pol                     = griddedInterpolant(S1,FE1,INNO1,squeeze(S_pol));
        Se_dense{educ,t}(:,:,:) = reshape(pol(S2,FE2,INNO2),length(S_dense_aux),length(FE_pos),length(INNO_pos));
        
        pol                     = griddedInterpolant(S1,FE1,INNO1,squeeze(C_pol));
        Ce_dense{educ,t}(:,:,:) = reshape(pol(S2,FE2,INNO2),length(S_dense_aux),length(FE_pos),length(INNO_pos));
    end
end

educ  = 3;
for t = 1:3
    S_dense_aux          = par.grids_dense{educ,par.Je1_pos+t-1};
    S_orig               = par.grids{educ,par.Je1_pos+t-1};
    Sp_max               = max(par.grids_dense{educ,par.Je1_pos+t});
    
    psy_orig             = PSY;
    
    S_pol                = Se{educ,t};
    C_pol                = Ce{educ,t};
    if max(S_pol(:)) > Sp_max
        j_pos        = par.Je1_pos+t-1;
        fprintf('Attention: age = %i, educ = %i, extrapolation in ergodic distribution, max extrap = %3.2f pp  \n',par.age(j_pos),educ,log(max(S_pol(:))/Sp_max))
        S_pol            = min(S_pol,Sp_max);
    end
    
    Sp_min           = min(par.grids_dense{educ,par.Je1_pos+t});
    if min(S_pol(:)) < Sp_min
        j_pos        = par.Je1_pos+t-1;
        fprintf('Attention: age = %i, educ = %i, extrapolation in ergodic distribution, min extrap = %3.2f pp  \n',par.age(j_pos),educ,log(min(S_pol(:))/Sp_min))
        S_pol        = max(S_pol,Sp_min);
    end
    
    if t <= 1
        psy_orig         = PSY;
        
        [S1,FE1,PSY1]    = ndgrid(S_orig,FE_pos,psy_orig);
        [S2,FE2,PSY2]    = ndgrid(S_dense_aux,FE_pos,psy_orig);
        
        pol                 = griddedInterpolant(S1,FE1,PSY1,squeeze(S_pol));
        Se_dense{educ,t}(:,:,:)    = reshape(pol(S2,FE2,PSY2), length(S_dense_aux),length(FE_pos),length(psy_orig));
        
        pol                 = griddedInterpolant(S1,FE1,PSY1,squeeze(C_pol));
        Ce_dense{educ,t}(:,:,:)    = reshape(pol(S2,FE2,PSY2), length(S_dense_aux),length(FE_pos),length(psy_orig));
        
    else
        [S1,FE1]   = ndgrid(S_orig,FE_pos);
        [S2,FE2]   = ndgrid(S_dense_aux,FE_pos);
        
        pol                     = griddedInterpolant(S1,FE1,squeeze(S_pol));
        Se_dense{educ,t}(:,:) = reshape(pol(S2,FE2),length(S_dense_aux),length(FE_pos));
        
        pol                     = griddedInterpolant(S1,FE1,squeeze(C_pol));
        Ce_dense{educ,t}(:,:) = reshape(pol(S2,FE2),length(S_dense_aux),length(FE_pos));
    end
end
clear Se

%% 3. Interpolate policy function of savings for working young
for j_pos = par.Je2_pos+2:par.Jc_pos-1
    S_dense{j_pos}        = nan(length(par.grids_dense{1,j_pos}),length(FE_pos),length(EDUC),length(INNO_pos));
    C_dense{j_pos}        = nan(length(par.grids_dense{1,j_pos}),length(FE_pos),length(EDUC),length(INNO_pos));
    for educ = 1:length(EDUC)
        S_dense_aux      = par.grids_dense{educ,j_pos};
        S_orig           = par.grids{educ,j_pos}';
        
        Sp_max           = max(par.grids_dense{educ,j_pos+1});
        S_pol            = squeeze(S{j_pos}(:,:,educ,:));
        C_pol            = squeeze(C{j_pos}(:,:,educ,:));
        if max(S_pol(:)) > Sp_max
            fprintf('Attention: age = %i, educ = %i, extrapolation in ergodic distribution, max extrap = %3.2f pp  \n',par.age(j_pos),educ,log(max(S_pol(:))/Sp_max))
            S_pol            = min(S_pol,Sp_max);
        end
        
        Sp_min           = min(par.grids_dense{educ,j_pos+1});
        if min(S_pol(:)) < Sp_min
            fprintf('Attention: age = %i, educ = %i, extrapolation in ergodic distribution, max extrap = %3.2f pp  \n',par.age(j_pos),educ,log(min(S_pol(:))/Sp_min))
            S_pol            = max(S_pol,Sp_min);
        end
        
        [S1,FE1,INNO1]   = ndgrid(S_orig,FE_pos,INNO_pos);
        [S2,FE2,INNO2]   = ndgrid(S_dense_aux,FE_pos,INNO_pos);
        
        pol                         = griddedInterpolant(S1,FE1,INNO1,squeeze(S_pol));
        S_dense{j_pos}(:,:,educ,:)  = reshape(pol(S2,FE2,INNO2),length(S_dense_aux),length(FE_pos),1,length(INNO_pos));
        
        pol                         = griddedInterpolant(S1,FE1,INNO1,squeeze(C_pol));
        C_dense{j_pos}(:,:,educ,:)  = reshape(pol(S2,FE2,INNO2),length(S_dense_aux),length(FE_pos),1,length(INNO_pos));
        
    end
end

%% 4. Interpolate policy function of savings for fertility
j_pos                   = par.Jc_pos;
S_dense{j_pos}          = nan(length(par.grids_dense{1,j_pos}),length(FE_pos),length(EDUC),length(INNO_pos));
C_dense{j_pos}          = nan(length(par.grids_dense{1,j_pos}),length(FE_pos),length(EDUC),length(INNO_pos));
Ck_dense{j_pos}         = nan(length(par.grids_dense{1,j_pos}),length(FE_pos),length(EDUC),length(INNO_pos));
Np_dense                = nan(length(par.grids_dense{1,j_pos}),length(FE_pos),length(EDUC),length(INNO_pos));
Tp_dense{j_pos}         = nan(length(par.grids_dense{1,j_pos}),length(FE_pos),length(EDUC),length(INNO_pos));
for educ = 1:length(EDUC)
    S_dense_aux     = par.grids_dense{educ,j_pos};
    S_orig          = par.grids{educ,j_pos};
    Sp_max          = max(par.grids_dense{educ,j_pos+1});
    
    S_pol           = squeeze(S{j_pos}(:,:,educ,:));
    C_pol           = squeeze(C{j_pos}(:,:,educ,:));
    Ck_pol          = squeeze(Ck{j_pos}(:,:,educ,:));
    if max(S_pol(:)) > Sp_max
        fprintf('Attention: age = %i, educ = %i, extrapolation in ergodic distribution, max extrap = %3.2f pp  \n',par.age(j_pos),educ,log(max(S_pol(:))/Sp_max))
        S_pol       = min(S_pol,Sp_max);
    end
    
    Sp_min           = min(par.grids_dense{educ,j_pos+1});
    if min(S_pol(:)) < Sp_min
        fprintf('Attention: age = %i, educ = %i, extrapolation in ergodic distribution, max extrap = %3.2f pp  \n',par.age(j_pos),educ,log(min(S_pol(:))/Sp_min))
        S_pol            = max(S_pol,Sp_min);
    end
    
    [S1,FE1,INNO1]       = ndgrid(S_orig,FE_pos,INNO_pos);
    [S2,FE2,INNO2]       = ndgrid(S_dense_aux,FE_pos,INNO_pos);
    
    pol                         = griddedInterpolant(S1,FE1,INNO1,S_pol);
    S_dense{j_pos}(:,:,educ,:)  = reshape(pol(S2,FE2,INNO2),length(S_dense_aux),length(FE_pos),1,length(INNO_pos));
    
    pol                         = griddedInterpolant(S1,FE1,INNO1,C_pol);
    C_dense{j_pos}(:,:,educ,:)  = reshape(pol(S2,FE2,INNO2),length(S_dense_aux),length(FE_pos),1,length(INNO_pos));
    
    pol                         = griddedInterpolant(S1,FE1,INNO1,Ck_pol);
    Ck_dense{j_pos}(:,:,educ,:)  = reshape(pol(S2,FE2,INNO2),length(S_dense_aux),length(FE_pos),1,length(INNO_pos));
    
    Np_pol                      = squeeze(Np(:,:,educ,:));
    pol                         = griddedInterpolant(S1,FE1,INNO1,Np_pol,'nearest');
    Np_dense(:,:,educ,:)        = reshape(pol(S2,FE2,INNO2),length(S_dense_aux),length(FE_pos),1,length(INNO_pos));
    
    Tp_pol                      = squeeze(Tp{j_pos}(:,:,educ,:));
    pol                         = griddedInterpolant(S1,FE1,INNO1,Tp_pol,'nearest');
    Tp_dense{j_pos}(:,:,educ,:)  = reshape(pol(S2,FE2,INNO2),length(S_dense_aux),length(FE_pos),1,length(INNO_pos));
end


%% 5. Interpolate policy function of savings for work with child
for j_pos  = par.Jc_pos+1:par.Jc_pos+par.Je1_pos-3
    S_dense{j_pos}              = nan(length(par.grids_dense{1,j_pos}),length(FE_pos),length(EDUC),length(INNO_pos),length(N));
    C_dense{j_pos}              = nan(length(par.grids_dense{1,j_pos}),length(FE_pos),length(EDUC),length(INNO_pos),length(N));
    Ck_dense{j_pos}             = nan(length(par.grids_dense{1,j_pos}),length(FE_pos),length(EDUC),length(INNO_pos),length(N));
    Tp_dense{j_pos}             = nan(length(par.grids_dense{1,j_pos}),length(FE_pos),length(EDUC),length(INNO_pos),length(N));
    for educ = 1:length(EDUC)
        S_dense_aux       = par.grids_dense{educ,j_pos};
        S_orig            = par.grids{educ,j_pos};
        Sp_max            = max(par.grids_dense{educ,j_pos+1});
        
        S_pol             = squeeze(S{j_pos}(:,:,educ,:,:));
        C_pol             = squeeze(C{j_pos}(:,:,educ,:,:));
        Ck_pol            = squeeze(Ck{j_pos}(:,:,educ,:,:));
        Tp_pol            = squeeze(Tp{j_pos}(:,:,educ,:,:));
        
        if max(S_pol(:)) > Sp_max
            fprintf('Attention: age = %i, educ = %i, extrapolation in ergodic distribution, max extrap = %3.2f pp  \n',par.age(j_pos),educ,log(max(S_pol(:))/Sp_max))
            S_pol         = min(S_pol,Sp_max);
        end
        
        Sp_min           = min(par.grids_dense{educ,j_pos+1});
        if min(S_pol(:)) < Sp_min
            fprintf('Attention: age = %i, educ = %i, extrapolation in ergodic distribution, min extrap = %3.2f pp  \n',par.age(j_pos),educ,log(min(S_pol(:))/Sp_min))
            S_pol            = max(S_pol,Sp_min);
        end
        
        switch options.Fertility
            case 'Endo'
                [S1,FE1,INNO1,N1]           = ndgrid(S_orig,FE_pos,INNO_pos,N);
                [S2,FE2,INNO2,N2]           = ndgrid(S_dense_aux,FE_pos,INNO_pos,N);
                
                pol                           = griddedInterpolant(S1,FE1,INNO1,N1,S_pol);
                S_dense{j_pos}(:,:,educ,:,:)  = reshape(pol(S2,FE2,INNO2,N2),length(S_dense_aux),length(FE_pos),1,length(INNO_pos),length(N));
                
                pol                         = griddedInterpolant(S1,FE1,INNO1,N1,C_pol);
                C_dense{j_pos}(:,:,educ,:,:)  = reshape(pol(S2,FE2,INNO2,N2),length(S_dense_aux),length(FE_pos),1,length(INNO_pos),length(N));
                
                pol                         = griddedInterpolant(S1,FE1,INNO1,N1,Ck_pol);
                Ck_dense{j_pos}(:,:,educ,:,:) = reshape(pol(S2,FE2,INNO2,N2),length(S_dense_aux),length(FE_pos),1,length(INNO_pos),length(N));
                
                pol                         = griddedInterpolant(S1,FE1,INNO1,N1,Tp_pol,'nearest');
                Tp_dense{j_pos}(:,:,educ,:,:) = reshape(pol(S2,FE2,INNO2,N2),length(S_dense_aux),length(FE_pos),1,length(INNO_pos),length(N));
            case 'Exo'
                
                [S1,FE1,INNO1]       = ndgrid(S_orig,FE_pos,INNO_pos);
                [S2,FE2,INNO2]       = ndgrid(S_dense_aux,FE_pos,INNO_pos);
                
                pol                         = griddedInterpolant(S1,FE1,INNO1,S_pol);
                S_dense{j_pos}(:,:,educ,:,1)  = reshape(pol(S2,FE2,INNO2),length(S_dense_aux),length(FE_pos),1,length(INNO_pos));
                
                pol                         = griddedInterpolant(S1,FE1,INNO1,C_pol);
                C_dense{j_pos}(:,:,educ,:,1)  = reshape(pol(S2,FE2,INNO2),length(S_dense_aux),length(FE_pos),1,length(INNO_pos));
                
                pol                         = griddedInterpolant(S1,FE1,INNO1,Ck_pol);
                Ck_dense{j_pos}(:,:,educ,:,1)  = reshape(pol(S2,FE2,INNO2),length(S_dense_aux),length(FE_pos),1,length(INNO_pos));
                
                Tp_pol                      = squeeze(Tp{j_pos}(:,:,educ,:));
                pol                         = griddedInterpolant(S1,FE1,INNO1,Tp_pol,'nearest');
                Tp_dense{j_pos}(:,:,educ,:,1)  = reshape(pol(S2,FE2,INNO2),length(S_dense_aux),length(FE_pos),1,length(INNO_pos));
                
        end
    end
end

%% 6. Interpolate policy function of savings for work + transfers to child
j_pos                       = par.Jc_pos+par.Je1_pos-2;

S_dense{j_pos}              = nan(length(par.grids_dense{1,j_pos}),length(FE_pos),length(EDUC),length(INNO_pos),length(N));
C_dense{j_pos}              = nan(length(par.grids_dense{1,j_pos}),length(FE_pos),length(EDUC),length(INNO_pos),length(N));
Ck_dense{j_pos}             = nan(length(par.grids_dense{1,j_pos}),length(FE_pos),length(EDUC),length(INNO_pos),length(N));
Tp_dense{j_pos}             = nan(length(par.grids_dense{1,j_pos}),length(FE_pos),length(EDUC),length(INNO_pos),length(N));
PHIp_dense                  = nan(length(par.grids_dense{1,j_pos}),length(FE_pos),length(EDUC),length(INNO_pos),length(N));
for educ = 1:length(EDUC)
    S_dense_aux       = par.grids_dense{educ,j_pos};
    S_orig            = par.grids{educ,j_pos};
    Sp_max            = max(par.grids_dense{educ,j_pos+1});
    
    Tp_pol            = squeeze(Tp{j_pos}(:,:,educ,:,:));
    PHIp_pol          = squeeze(PHIp(:,:,educ,:,:));
    S_pol             = squeeze(S{j_pos}(:,:,educ,:,:));
    C_pol             = squeeze(C{j_pos}(:,:,educ,:,:));
    Ck_pol            = squeeze(Ck{j_pos}(:,:,educ,:,:));
    
    %     aux               = reshape(1*(H_dense_aux<= CUTOFFS(educ,1)) + ...
    %         2.*(H_dense_aux> CUTOFFS(educ,1)).*(H_dense_aux<= CUTOFFS(educ,2)) +...
    %         3.*(H_dense_aux> CUTOFFS(educ,2)),1,length(H_dense_aux));
    %     Grs_dense{j_pos}(:,:,educ,:)  = repmat(aux,length(par.grids_dense{1,j_pos}),1,1,length(N));
    
    if max(S_pol(:)) > Sp_max
        fprintf('Attention: age = %i, educ = %i, extrapolation in ergodic distribution, max extrap = %3.2f pp  \n',par.age(j_pos),educ,log(max(S_pol(:))/Sp_max))
        S_pol         = min(S_pol,Sp_max);
    end
    
    Sp_min           = min(par.grids_dense{educ,j_pos+1});
    if min(S_pol(:)) < Sp_min
        fprintf('Attention: age = %i, educ = %i, extrapolation in ergodic distribution, max extrap = %3.2f pp  \n',par.age(j_pos),educ,log(min(S_pol(:))/Sp_min))
        S_pol            = max(S_pol,Sp_min);
    end
    
    switch options.Fertility
        case 'Endo'
            [S1,FE1,INNO1,N1]           = ndgrid(S_orig,FE_pos,INNO_pos,N);
            [S2,FE2,INNO2,N2]           = ndgrid(S_dense_aux,FE_pos,INNO_pos,N);
            
            pol                           = griddedInterpolant(S1,FE1,INNO1,N1,S_pol);
            S_dense{j_pos}(:,:,educ,:,:)  = reshape(pol(S2,FE2,INNO2,N2),length(S_dense_aux),length(FE_pos),1,length(INNO_pos),length(N));
            
            pol                         = griddedInterpolant(S1,FE1,INNO1,N1,C_pol);
            C_dense{j_pos}(:,:,educ,:,:)  = reshape(pol(S2,FE2,INNO2,N2),length(S_dense_aux),length(FE_pos),1,length(INNO_pos),length(N));
            
            pol                         = griddedInterpolant(S1,FE1,INNO1,N1,Ck_pol);
            Ck_dense{j_pos}(:,:,educ,:,:) = reshape(pol(S2,FE2,INNO2,N2),length(S_dense_aux),length(FE_pos),1,length(INNO_pos),length(N));
            
            pol                         = griddedInterpolant(S1,FE1,INNO1,N1,Tp_pol,'nearest');
            Tp_dense{j_pos}(:,:,educ,:,:) = reshape(pol(S2,FE2,INNO2,N2),length(S_dense_aux),length(FE_pos),1,length(INNO_pos),length(N));
            
            pol                         = griddedInterpolant(S1,FE1,INNO1,N1,PHIp_pol,'nearest');
            PHIp_dense(:,:,educ,:,:)      = reshape(pol(S2,FE2,INNO2,N2),length(S_dense_aux),length(FE_pos),1,length(INNO_pos),length(N));
            
        case 'Exo'
            
            [S1,FE1,INNO1]       = ndgrid(S_orig,FE_pos,INNO_pos);
            [S2,FE2,INNO2]       = ndgrid(S_dense_aux,FE_pos,INNO_pos);
            
            pol                         = griddedInterpolant(S1,FE1,INNO1,S_pol);
            S_dense{j_pos}(:,:,educ,:,1)  = reshape(pol(S2,FE2,INNO2),length(S_dense_aux),length(FE_pos),1,length(INNO_pos));
            
            pol                         = griddedInterpolant(S1,FE1,INNO1,C_pol);
            C_dense{j_pos}(:,:,educ,:,1)  = reshape(pol(S2,FE2,INNO2),length(S_dense_aux),length(FE_pos),1,length(INNO_pos));
            
            pol                         = griddedInterpolant(S1,FE1,INNO1,Ck_pol);
            Ck_dense{j_pos}(:,:,educ,:,1)  = reshape(pol(S2,FE2,INNO2),length(S_dense_aux),length(FE_pos),1,length(INNO_pos));
            
            Tp_pol                      = squeeze(Tp{j_pos}(:,:,educ,:));
            pol                         = griddedInterpolant(S1,FE1,INNO1,Tp_pol,'nearest');
            Tp_dense{j_pos}(:,:,educ,:,1)  = reshape(pol(S2,FE2,INNO2),length(S_dense_aux),length(FE_pos),1,length(INNO_pos));
            
            pol                         = griddedInterpolant(S1,FE1,INNO1,PHIp_pol,'nearest');
            PHIp_dense(:,:,educ,:,1)        = reshape(pol(S2,FE2,INNO2),length(S_dense_aux),length(FE_pos),1,length(INNO_pos));
            
    end
end

%% 6. Interpolate policy function of savings for work old
for j_pos       = par.Jc_pos+par.Je1_pos-1: par.Jr_pos-1
    S_dense{j_pos}                  = nan(length(par.grids_dense{1,j_pos}),length(FE_pos),length(EDUC),length(INNO_pos));
    C_dense{j_pos}                  = nan(length(par.grids_dense{1,j_pos}),length(FE_pos),length(EDUC),length(INNO_pos));
    for educ = 1:length(EDUC)
        S_dense_aux = par.grids_dense{educ,j_pos};
        S_orig      = par.grids{educ,j_pos};
        Sp_max      = max(par.grids_dense{educ,j_pos+1});

        S_pol       = squeeze(S{j_pos}(:,:,educ,:));
        C_pol       = squeeze(C{j_pos}(:,:,educ,:));
        if max(S_pol(:)) > Sp_max
            fprintf('Attention: age = %i, educ = %i, extrapolation in ergodic distribution, max extrap = %3.2f pp  \n',par.age(j_pos),educ,log(max(S_pol(:))/Sp_max))
            S_pol       = min(S_pol,Sp_max);
        end
        
        Sp_min           = min(par.grids_dense{educ,j_pos+1});
        Sp_min           = min(par.grids_dense{educ,j_pos+1});
        if min(S_pol(:)) < Sp_min
            fprintf('Attention: age = %i, educ = %i, extrapolation in ergodic distribution, max extrap = %3.2f pp  \n',par.age(j_pos),educ,log(min(S_pol(:))/Sp_min))
            S_pol            = max(S_pol,Sp_min);
        end
        
        [S1,FE1,INNO1]       = ndgrid(S_orig,FE_pos,INNO_pos);
        [S2,FE2,INNO2]       = ndgrid(S_dense_aux,FE_pos,INNO_pos);
        
        pol                         = griddedInterpolant(S1,FE1,INNO1,S_pol);
        S_dense{j_pos}(:,:,educ,:)  = reshape(pol(S2,FE2,INNO2),length(S_dense_aux),length(FE_pos),1,length(INNO_pos));

        pol                         = griddedInterpolant(S1,FE1,INNO1,C_pol);
        C_dense{j_pos}(:,:,educ,:)  = reshape(pol(S2,FE2,INNO2),length(S_dense_aux),length(FE_pos),1,length(INNO_pos));
    end
end

%% 8. Retirement
for j_pos       = par.Jr_pos:par.Jd_pos
    S_dense{j_pos}          = nan(length(par.grids_dense{1,j_pos}),length(FE_pos),length(EDUC));
    for educ = 1:length(EDUC)
        S_dense_aux = par.grids_dense{educ,j_pos};
        S_orig      = par.grids{educ,j_pos};
        j_pos_1     = min(j_pos+1,par.Jd_pos);
        Sp_max      = max(par.grids_dense{educ,j_pos_1});
        S_pol       = squeeze(S{j_pos}(:,:,educ));
        C_pol       = squeeze(C{j_pos}(:,:,educ));
        if max(S_pol(:)) > Sp_max
            fprintf('Attention: age = %i, educ = %i, extrapolation in ergodic distribution, max extrap = %3.2f pp  \n',par.age(j_pos),educ,log(max(S_pol(:))/Sp_max))
            S_pol       = min(S_pol,Sp_max);
        end
        
        Sp_min           = min(par.grids_dense{educ,j_pos_1});
        if min(S_pol(:)) < Sp_min
            fprintf('Attention: age = %i, educ = %i, extrapolation in ergodic distribution, max extrap = %3.2f pp  \n',par.age(j_pos),educ,log(min(S_pol(:))/Sp_min))
            S_pol            = max(S_pol,Sp_min);
        end
        
        [S1,FE1]       = ndgrid(S_orig,FE_pos);
        [S2,FE2]       = ndgrid(S_dense_aux,FE_pos);
        
         pol                         = griddedInterpolant(S1,FE1,S_pol);
        S_dense{j_pos}(:,:,educ)  = reshape(pol(S2,FE2),length(S_dense_aux),length(FE_pos),1);

        pol                         = griddedInterpolant(S1,FE1,C_pol);
        C_dense{j_pos}(:,:,educ)  = reshape(pol(S2,FE2),length(S_dense_aux),length(FE_pos),1);

    end
end

clear S
%% Save interpolated policy functions

% 9. Rename grid vectors
par2                        = par;
par2                        = rmfield(par2,'grids');
par2                        = rmfield(par2,'grids_dense');
par2.grids                  = par.grids_dense;

% Policy functions
Se                          = Se_dense;
S                           = S_dense;
Ce                          = Ce_dense;
C                           = C_dense;
Ck                          = Ck_dense;
Np                          = Np_dense;
PHIp                        = PHIp_dense;
% Grs                         = Grs_dense;
pol_dense.Se                = Se;
pol_dense.S                 = S;
pol_dense.Ce                = Ce;
pol_dense.C                 = C;
pol_dense.Ck                = Ck;
pol_dense.Np                = Np;
pol_dense.PHIp              = PHIp;
pol_dense.Tp                = Tp_dense;
% pol_dense.Grs               = Grs;
pol_dense.tau0              = tau0;

end