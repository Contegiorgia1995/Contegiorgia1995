function [IGE] = gen_ige_sample(par,pol,muc,options)
% Intergenerational distributions

%Enter with mu_ige_aux: (Spar,FEpar,Educpar,INNOpar,FEc,Sc);
%Enter with mu_ige_aux: (Spar,FEpar,Psypar);
Nsample   = options.Nsample;
rng(666);   % Set seed for shocks

FE_pos    = par.inc.fe_pos;
INNO_pos  = par.inc.inno_pos;
EDUC      = par.educ;
PSY       = 1:length(par.psy_val_col);
r_sav     = par.r_sav;
r_debt    = par.r_debt;

j_child_pos      = par.Je1_pos;

%% Create States and Income of Parents
IGE.States     = cell(par.Jd_pos-1,2);
IGE.lab        = cell(par.Jd_pos-1,2);
IGE.sav        = cell(par.Jd_pos-1,2);

%% Initialize Sample
siz        = size(muc);
draws      = rand(Nsample,1);
probs      = cumsum(muc(:));
init_ind   = nan(Nsample,1);
for in = 1:Nsample
    init_ind(in) = find((draws(in)>=[0; probs(1:end-1)]).*(draws(in)<probs));
end

[Spos,FEpos,Psypos]         = ind2sub(siz,init_ind);
IGE.States{j_child_pos-1,1} = [par.PHI(Spos)' FEpos Psypos pol.tau0(init_ind)];

for igen = 1:2
    %% Education Stage
    % Add draw of z
    j_pos       = par.Je1_pos;
    z_prob      = par.inc.z_prob;
    z0_prob    = [squeeze(z_prob{1,1}(par.Je1_pos,1,:))'; squeeze(z_prob{2,1}(par.Je2_pos,1,:))'; squeeze(z_prob{3,1}(par.Je2_pos+2,1,:))'];      % (educ,prob)      % (educ,prob)
    IGE.States{j_pos,igen} = nan(Nsample,5);
    draws      = rand(Nsample,1);
    for educ = 1:length(EDUC)
        p_min = cumsum([0 z0_prob(educ,1:end-1)]);
        p_max = cumsum(z0_prob(educ,:));
        for inno = 1:length(INNO_pos)
            ind = logical((IGE.States{j_pos-1,igen}(:,4) == educ).*(draws>=p_min(inno)).*(draws<p_max(inno)));
            IGE.States{j_pos,igen}(ind,:) =  [IGE.States{j_pos-1,igen}(ind,:) repmat(inno,sum(ind),1)];
        end
    end
    
    if igen == 1
       IGE.States{j_pos-1,igen}(:,5) = Spos; 
    elseif igen == 2
        fspace          = fundef({'spli',par.PHI,0,1});
        QPhi            = funbas(fspace,IGE.States{j_pos-1,igen}(:,1));
        [~,ip]          = max(QPhi,[],2);
        IGE.States{j_child_pos-1,igen}(:,5) = ip;
    end
    
    % Create labor and savings income
    IGE.lab{j_pos,igen} = nan(Nsample,1);
    IGE.sav{j_pos,igen} = nan(Nsample,1);
    IGE.States{j_pos+1,igen} = nan(Nsample,4);
    for educ = 1:length(EDUC)
        ind = logical((IGE.States{j_pos,igen}(:,4) == educ));
        
        S_pol            = pol.Se{educ,1};
        S_orig           = par.grids{educ,j_pos};
        Ssim            = IGE.States{j_pos,igen}(ind,1);
        FEsim           = IGE.States{j_pos,igen}(ind,2);
        
        rsim            = r_sav.*(Ssim>=0) + r_debt.*(Ssim<0);
        
        FE              = par.inc.fe{educ,1}(FEsim);
        AGE_PROF        = par.inc.age_prof{educ,1}(j_pos);
        
        if educ == 1
            [S1,FE1,INNO1]   = ndgrid(S_orig,FE_pos,INNO_pos);
            
            INNOsim         = IGE.States{j_pos,igen}(ind,5);
            INNO            = par.inc.z_val{educ,1}(j_pos,INNOsim)';
            
            IGE.lab{j_pos,igen}(ind,1) = par.w*exp(AGE_PROF) .* exp(FE) .* exp(INNO);
            
            polin             = griddedInterpolant(S1,FE1,INNO1,squeeze(S_pol));
            Spsim             = polin(Ssim,FEsim,INNOsim);
            IGE.States{j_pos+1,igen}(ind,1) = Spsim;
        else
            [S1,FE1,PSY1]   = ndgrid(S_orig,FE_pos,PSY);            
            PSYsim         = IGE.States{j_pos,igen}(ind,3);
            
            IGE.lab{j_pos,igen}(ind,1) = 0;
            
            polin             = griddedInterpolant(S1,FE1,PSY1,squeeze(S_pol));
            Spsim             = polin(Ssim,FEsim,PSYsim);
            IGE.States{j_pos+1,igen}(ind,1) = Spsim;
        end
        IGE.sav{j_pos,igen}(ind,1) = rsim.*Ssim;
    end
    
    % Education Stage: Future States
    j_pos      = j_pos + 1;
    draws      = rand(Nsample,1);
    for educ = 1:length(EDUC)
        if educ == 1
            for inno = 1:length(INNO_pos)
                z_prob = squeeze(par.inc.z_prob{educ,1}(j_pos,inno,:));
                p_min = cumsum([0; z_prob(1:end-1)]);
                p_max = cumsum(z_prob(:));
                for inno_p = 1:length(INNO_pos)
                    ind = logical((IGE.States{j_pos-1,igen}(:,5) == inno).*(IGE.States{j_pos-1,igen}(:,4) == educ).*(draws>=p_min(inno_p)).*(draws<p_max(inno_p)));
                    IGE.States{j_pos,igen}(ind,2)      =  IGE.States{j_pos-1,igen}(ind,2);
                    IGE.States{j_pos,igen}(ind,3)      =  IGE.States{j_pos-1,igen}(ind,4);
                    IGE.States{j_pos,igen}(ind,4)      =  repmat(inno_p,sum(ind),1);
                end
            end
        else
            ind = logical((IGE.States{j_pos-1,igen}(:,4) == educ));
            IGE.States{j_pos,igen}(ind,2) =  IGE.States{j_pos-1,igen}(ind,2);
            IGE.States{j_pos,igen}(ind,3) =  IGE.States{j_pos-1,igen}(ind,4);
            IGE.States{j_pos,igen}(ind,4) =  IGE.States{j_pos-1,igen}(ind,5);
        end
    end
    
    while j_pos <= par.Je2_pos+1    
        % Create labor and savings income
        IGE.lab{j_pos,igen} = nan(Nsample,1);
        IGE.sav{j_pos,igen} = nan(Nsample,1);
        IGE.States{j_pos+1,igen} = nan(Nsample,4);
        for educ = 1:length(EDUC)
            ind = logical((IGE.States{j_pos,igen}(:,3) == educ));

            S_pol            = pol.Se{educ,j_pos-(par.Je1_pos-1)};
            if educ == 3
                S_pol        = repmat(S_pol,1,1,length(INNO_pos));
            end
            S_orig           = par.grids{educ,j_pos};
            Ssim            = IGE.States{j_pos,igen}(ind,1);
            FEsim           = IGE.States{j_pos,igen}(ind,2);
            INNOsim         = IGE.States{j_pos,igen}(ind,4);
            rsim            = r_sav.*(Ssim>=0) + r_debt.*(Ssim<0);

            FE              = par.inc.fe{educ,1}(FEsim);
            AGE_PROF        = par.inc.age_prof{educ,1}(j_pos);

            [S1,FE1,INNO1]   = ndgrid(S_orig,FE_pos,INNO_pos);
            INNO            = par.inc.z_val{educ,1}(j_pos,INNOsim)';

            polin             = griddedInterpolant(S1,FE1,INNO1,squeeze(S_pol));
            Spsim             = polin(Ssim,FEsim,INNOsim);
            IGE.States{j_pos+1,igen}(ind,1) = Spsim;
            
            if educ < 3            
                IGE.lab{j_pos,igen}(ind,1) = par.w*exp(AGE_PROF) .* exp(FE) .* exp(INNO);
            else
                IGE.lab{j_pos,igen}(ind,1) = par.w_college*exp(AGE_PROF) .* exp(FE);
            end
            IGE.sav{j_pos,igen}(ind,1) = rsim.*Ssim;
        end

        % Education Stage: Future States
        j_pos      = j_pos + 1;
        draws      = rand(Nsample,1);
        for educ = 1:length(EDUC)
            if educ <= 2
                for inno = 1:length(INNO_pos)
                    z_prob = squeeze(par.inc.z_prob{educ,1}(j_pos,inno,:));
                    p_min = cumsum([0; z_prob(1:end-1)]);
                    p_max = cumsum(z_prob(:));
                    for inno_p = 1:length(INNO_pos)
                        ind = logical((IGE.States{j_pos-1,igen}(:,4) == inno).*(IGE.States{j_pos-1,igen}(:,3) == educ).*(draws>=p_min(inno_p)).*(draws<p_max(inno_p)));
                        IGE.States{j_pos,igen}(ind,2)      =  IGE.States{j_pos-1,igen}(ind,2);
                        IGE.States{j_pos,igen}(ind,3)      =  IGE.States{j_pos-1,igen}(ind,3);
                        IGE.States{j_pos,igen}(ind,4)      =  repmat(inno_p,sum(ind),1);
                    end
                end
            else
                ind = logical((IGE.States{j_pos-1,igen}(:,3) == educ));
                IGE.States{j_pos,igen}(ind,2) =  IGE.States{j_pos-1,igen}(ind,2);
                IGE.States{j_pos,igen}(ind,3) =  IGE.States{j_pos-1,igen}(ind,3);
                IGE.States{j_pos,igen}(ind,4) =  IGE.States{j_pos-1,igen}(ind,4);
            end
        end
    end
    
    %% Iterate until Fertility - 1
    while j_pos <= par.Jc_pos-1
        
        % Create labor and savings income
        IGE.lab{j_pos,igen} = nan(Nsample,1);
        IGE.sav{j_pos,igen} = nan(Nsample,1);
        IGE.States{j_pos+1,igen} = nan(Nsample,4);
        for educ = 1:length(EDUC)
            ind = logical((IGE.States{j_pos,igen}(:,3) == educ));
            
            S_pol            = squeeze(pol.S{j_pos,1}(:,:,educ,:));
            S_orig           = par.grids{educ,j_pos};
            Ssim             = IGE.States{j_pos,igen}(ind,1);
            
            [S1,FE1,INNO1]   = ndgrid(S_orig,FE_pos,INNO_pos);
            
            INNOsim         = IGE.States{j_pos,igen}(ind,4);
            INNO            = par.inc.z_val{educ,1}(j_pos,INNOsim)';
            FEsim           = IGE.States{j_pos,igen}(ind,2);
            FE              = par.inc.fe{educ,1}(FEsim);
            AGE_PROF        = par.inc.age_prof{educ,1}(j_pos);
            
            IGE.lab{j_pos,igen}(ind,1) = par.w*exp(AGE_PROF) .* exp(FE) .* exp(INNO);
           
            
            rsim            = r_sav.*(Ssim>=0) + r_debt.*(Ssim<0);
            polin             = griddedInterpolant(S1,FE1,INNO1,squeeze(S_pol));
            Spsim             = polin(Ssim,FEsim,INNOsim);
            IGE.States{j_pos+1,igen}(ind,1) = Spsim;
            
            IGE.sav{j_pos,igen}(ind,1) = rsim.*Ssim;
        end
        
        % Future States
        j_pos      = j_pos + 1;
        draws      = rand(Nsample,1);
        for educ = 1:length(EDUC)
            for inno = 1:length(INNO_pos)
                z_prob = squeeze(par.inc.z_prob{educ,1}(j_pos,inno,:));
                p_min = cumsum([0; z_prob(1:end-1)]);
                p_max = cumsum(z_prob(:));
                for inno_p = 1:length(INNO_pos)
                    ind = logical((IGE.States{j_pos-1,igen}(:,4) == inno).*(IGE.States{j_pos-1,igen}(:,3) == educ).*(draws>=p_min(inno_p)).*(draws<p_max(inno_p)));
                    IGE.States{j_pos,igen}(ind,2) =  IGE.States{j_pos-1,igen}(ind,2);
                    IGE.States{j_pos,igen}(ind,3) =  IGE.States{j_pos-1,igen}(ind,3);
                    IGE.States{j_pos,igen}(ind,4) =  repmat(inno_p,sum(ind),1);
                end
            end
        end
    end
    
    
    %% Fertility Period
    % Create labor and savings income
    IGE.lab{j_pos,igen} = nan(Nsample,1);
    IGE.sav{j_pos,igen} = nan(Nsample,1);
    IGE.States{j_pos+1,igen} = nan(Nsample,5);
    for educ = 1:length(EDUC)
        ind = logical((IGE.States{j_pos,igen}(:,3) == educ));
        
        S_pol            = squeeze(pol.S{j_pos,1}(:,:,educ,:));
        Np_pol           = squeeze(pol.Np(:,:,educ,:));
        S_orig           = par.grids{educ,j_pos};
        Ssim             = IGE.States{j_pos,igen}(ind,1);
        
        [S1,FE1,INNO1]      = ndgrid(S_orig,FE_pos,INNO_pos);
        
        INNOsim         = IGE.States{j_pos,igen}(ind,4);
        INNO            = par.inc.z_val{educ,1}(j_pos,INNOsim)';
        FEsim           = IGE.States{j_pos,igen}(ind,2);
        FE              = par.inc.fe{educ,1}(FEsim);
        
        AGE_PROF        = par.inc.age_prof{educ,1}(j_pos);
        
        IGE.lab{j_pos,igen}(ind,1) = par.w*exp(AGE_PROF) .* exp(FE) .* exp(INNO);
        
        
        rsim            = r_sav.*(Ssim>=0) + r_debt.*(Ssim<0);
        polin             = griddedInterpolant(S1,FE1,INNO1,squeeze(S_pol));
        Spsim             = polin(Ssim,FEsim,INNOsim);
        IGE.States{j_pos+1,igen}(ind,1) = Spsim;
        
        polin             = griddedInterpolant(S1,FE1,INNO1,squeeze(Np_pol),'nearest');
        Npsim             = polin(Ssim,FEsim,INNOsim);
        IGE.States{j_pos+1,igen}(ind,5) = Npsim;
        
        IGE.sav{j_pos,igen}(ind,1) = rsim.*Ssim;
    end
    
    % Future States
    j_pos      = j_pos + 1;
    draws      = rand(Nsample,1);
    for educ = 1:length(EDUC)
        for inno = 1:length(INNO_pos)
            z_prob = squeeze(par.inc.z_prob{educ,1}(j_pos,inno,:));
            p_min = cumsum([0; z_prob(1:end-1)]);
            p_max = cumsum(z_prob(:));
            for inno_p = 1:length(INNO_pos)
                ind = logical((IGE.States{j_pos-1,igen}(:,4) == inno).*(IGE.States{j_pos-1,igen}(:,3) == educ).*(draws>=p_min(inno_p)).*(draws<p_max(inno_p)));
                IGE.States{j_pos,igen}(ind,2) =  IGE.States{j_pos-1,igen}(ind,2);
                IGE.States{j_pos,igen}(ind,3) =  IGE.States{j_pos-1,igen}(ind,3);
                IGE.States{j_pos,igen}(ind,4) =  repmat(inno_p,sum(ind),1);
            end
        end
    end

    
    %% Iterate until Transfers
    N           = 1:length(par.N);
    while j_pos <= par.Jc_pos+par.Je1_pos-3
        % Create labor and savings income
        IGE.lab{j_pos,igen} = nan(Nsample,1);
        IGE.sav{j_pos,igen} = nan(Nsample,1);
        IGE.States{j_pos+1,igen} = nan(Nsample,5);
        for educ = 1:length(EDUC)
            ind = logical((IGE.States{j_pos,igen}(:,3) == educ));
            
            INNOsim         = IGE.States{j_pos,igen}(ind,4);
            INNO            = par.inc.z_val{educ,1}(j_pos,INNOsim)';
            FEsim           = IGE.States{j_pos,igen}(ind,2);
            FE              = par.inc.fe{educ,1}(FEsim);
            Nsim            = IGE.States{j_pos,igen}(ind,5);
            AGE_PROF        = par.inc.age_prof{educ,1}(j_pos);
            
            IGE.lab{j_pos,igen}(ind,1) = par.w*exp(AGE_PROF) .* exp(FE) .* exp(INNO);
            
            S_pol            = squeeze(pol.S{j_pos,1}(:,:,educ,:,:));
            S_orig           = par.grids{educ,j_pos};
            Ssim            = IGE.States{j_pos,igen}(ind,1);
            
            switch options.Fertility
                case 'Endo'
                    [S1,FE1,INNO1,N1]      = ndgrid(S_orig,FE_pos,INNO_pos,N);
                    
                    polin             = griddedInterpolant(S1,FE1,INNO1,N1,squeeze(S_pol));
                    Spsim             = polin(Ssim,FEsim,INNOsim,Nsim);
                    
                case 'Exo'
                    [S1,FE1,INNO1]      = ndgrid(S_orig,FE_pos,INNO_pos);
                    
                    polin             = griddedInterpolant(S1,FE1,INNO1,squeeze(S_pol));
                    Spsim             = polin(Ssim,FEsim,INNOsim);
            end
            rsim            = r_sav.*(Ssim>=0) + r_debt.*(Ssim<0);      
            IGE.States{j_pos+1,igen}(ind,1) = Spsim;
            IGE.sav{j_pos,igen}(ind,1) = rsim.*Ssim;
        end
        
        % Future States
        j_pos      = j_pos + 1;
        draws      = rand(Nsample,1);
        for educ = 1:length(EDUC)
            for inno = 1:length(INNO_pos)
                z_prob = squeeze(par.inc.z_prob{educ,1}(j_pos,inno,:));
                p_min = cumsum([0; z_prob(1:end-1)]);
                p_max = cumsum(z_prob(:));
                for inno_p = 1:length(INNO_pos)
                    ind = logical((IGE.States{j_pos-1,igen}(:,4) == inno).*(IGE.States{j_pos-1,igen}(:,3) == educ).*(draws>=p_min(inno_p)).*(draws<p_max(inno_p)));
                    IGE.States{j_pos,igen}(ind,2) =  IGE.States{j_pos-1,igen}(ind,2);
                    IGE.States{j_pos,igen}(ind,3) =  IGE.States{j_pos-1,igen}(ind,3);
                    IGE.States{j_pos,igen}(ind,4) =  repmat(inno_p,sum(ind),1);
                    IGE.States{j_pos,igen}(ind,5) =  IGE.States{j_pos-1,igen}(ind,5);
                end
            end
        end
    end
    
    
    %% Transfers to Children
    N           = 1:length(par.N);    
    % Create labor and savings income
    IGE.lab{j_pos,igen} = nan(Nsample,1);
    IGE.sav{j_pos,igen} = nan(Nsample,1);
    IGE.States{j_pos+1,igen} = nan(Nsample,4);
    for educ = 1:length(EDUC)
        ind = logical((IGE.States{j_pos,igen}(:,3) == educ));

        S_pol            = squeeze(pol.S{j_pos,1}(:,:,educ,:,:));
        S_orig           = par.grids{educ,j_pos};
        Ssim            = IGE.States{j_pos,igen}(ind,1);
        
        INNOsim         = IGE.States{j_pos,igen}(ind,4);
        INNO            = par.inc.z_val{educ,1}(j_pos,INNOsim)';
        FEsim           = IGE.States{j_pos,igen}(ind,2);
        FE              = par.inc.fe{educ,1}(FEsim);
        Nsim            = IGE.States{j_pos,igen}(ind,5);
        AGE_PROF        = par.inc.age_prof{educ,1}(j_pos);
        
        switch options.Fertility
            case 'Endo'
                [S1,FE1,INNO1,N1]      = ndgrid(S_orig,FE_pos,INNO_pos,N);

                polin             = griddedInterpolant(S1,FE1,INNO1,N1,squeeze(S_pol));
                Spsim             = polin(Ssim,FEsim,INNOsim,Nsim);

            case 'Exo'
                [S1,FE1,INNO1]      = ndgrid(S_orig,FE_pos,INNO_pos);

                polin             = griddedInterpolant(S1,FE1,INNO1,squeeze(S_pol));
                Spsim             = polin(Ssim,FEsim,INNOsim);
        end

        

        IGE.lab{j_pos,igen}(ind,1) = par.w*exp(AGE_PROF) .* exp(FE) .* exp(INNO);
        
        rsim            = r_sav.*(Ssim>=0) + r_debt.*(Ssim<0);
        IGE.States{j_pos+1,igen}(ind,1) = Spsim;

        IGE.sav{j_pos,igen}(ind,1) = rsim.*Ssim;
        
    end

    % Future States
    j_pos      = j_pos + 1;
    draws      = rand(Nsample,1);
    for educ = 1:length(EDUC)
        for inno = 1:length(INNO_pos)
            z_prob = squeeze(par.inc.z_prob{educ,1}(j_pos,inno,:));
            p_min = cumsum([0; z_prob(1:end-1)]);
            p_max = cumsum(z_prob(:));
            for inno_p = 1:length(INNO_pos)
                ind = logical((IGE.States{j_pos-1,igen}(:,4) == inno).*(IGE.States{j_pos-1,igen}(:,3) == educ).*(draws>=p_min(inno_p)).*(draws<p_max(inno_p)));
                IGE.States{j_pos,igen}(ind,2) =  IGE.States{j_pos-1,igen}(ind,2);
                IGE.States{j_pos,igen}(ind,3) =  IGE.States{j_pos-1,igen}(ind,3);
                IGE.States{j_pos,igen}(ind,4) =  repmat(inno_p,sum(ind),1);
            end
        end
    end
    
    
    %% Create next generation
    if igen == 1
        if igen == 1
            IGE.States{j_child_pos-1,igen+1} = nan(Nsample,4);
        end
        j_pos      = j_pos - 1;
        for educ = 1:length(EDUC)
            ind       = logical((IGE.States{j_pos,igen}(:,3) == educ));
            S_orig    = par.grids{educ,j_pos};
            Ssim      = IGE.States{j_pos,igen}(ind,1);
            INNOsim         = IGE.States{j_pos,igen}(ind,4);
            FEsim           = IGE.States{j_pos,igen}(ind,2);
            PHI_pol          = squeeze(pol.PHIp(:,:,educ,:,:));
            Nsim = IGE.States{j_pos,igen}(ind,5);


            switch options.Fertility
                case 'Endo'
                    [S1,FE1,INNO1,N1]      = ndgrid(S_orig,FE_pos,INNO_pos,N);
                    
                    polin             = griddedInterpolant(S1,FE1,INNO1,N1,squeeze(PHI_pol),'nearest');
                    PhiSim            = polin(Ssim,FEsim,INNOsim,Nsim);
                    
                case 'Exo'
                    [S1,FE1,INNO1]      = ndgrid(S_orig,FE_pos,INNO_pos);
                    
                    polin             = griddedInterpolant(S1,FE1,INNO1,squeeze(PHI_pol),'nearest');
                    PhiSim            = polin(Ssim,FEsim,INNOsim);
            end
                        
            IGE.States{j_child_pos-1,igen+1}(ind,1) = par.PHI(PhiSim)';
        end
        
        % Fix with fertility rate
        switch options.Fertility
            case 'Endo'
                Nsim = IGE.States{j_pos,igen}(:,5);
                for in = N
                    ind         = logical((Nsim == in));
                    if in == 1
                        IGE.States{j_child_pos-1,igen+1}(ind,1) = nan;
                    end
                    if in > 2
                        count_ind   = sum(ind);
                        in0         = in;

                        while in0>2
                            count_tot   = length(IGE.States{j_pos,igen});

                            % Fix child states
                            IGE.States{j_child_pos-1,igen+1}(count_tot+1:count_tot+count_ind,1) = IGE.States{j_child_pos-1,igen+1}(ind,1);

                            % Fix parent states
                            for j_pos_fix = par.Je1_pos-1:j_pos + 1
                                IGE.States{j_pos_fix,igen}(count_tot+1:count_tot+count_ind,:) = IGE.States{j_pos_fix,igen}(ind,:);
                                if j_pos_fix >= par.Je1_pos && j_pos_fix <= j_pos
                                    IGE.lab{j_pos_fix,igen}(count_tot+1:count_tot+count_ind,:)    = IGE.lab{j_pos_fix,igen}(ind,:);
                                    IGE.sav{j_pos_fix,igen}(count_tot+1:count_tot+count_ind,:)    = IGE.sav{j_pos_fix,igen}(ind,:);          
                                end
                            end
                            in0         = in0 - 1;
                        end
                    end
                end
        end
        
        Nsampletot = length(IGE.States{j_pos,igen});
        draws      = rand(Nsampletot,1);
        drawsh0    = rand(Nsampletot,1);
        prob_h0    = nan(Nsampletot,length(FE_pos));
        for educ = 1:length(EDUC)
            ind      = logical((IGE.States{j_pos,igen}(:,3) == educ));

            % Add FE
            INNOsim         = IGE.States{j_pos,igen}(ind,4);
            INNO            = par.inc.z_val{educ,1}(j_pos,INNOsim)';
            FEsim           = IGE.States{j_pos,igen}(ind,2);
            FE              = par.inc.fe{educ,1}(FEsim);
            AGE_PROF        = par.inc.age_prof{educ,1}(j_pos);

            hp0             = exp(AGE_PROF) .* exp(FE) .* exp(INNO);
            prob_h0(ind,:)    = squeeze(par.PG(hp0))';
            
            % Add Psy
            PSY_p    = par.psy_prob(educ,:);
            
            Sample   = 1:Nsampletot;
            for in = Sample(ind)
                IGE.States{j_child_pos-1,igen+1}(in,2) = find((drawsh0(in)>=cumsum([0 prob_h0(in,1:end-1)])).*(drawsh0(in)<cumsum(prob_h0(in,:))));
                IGE.States{j_child_pos-1,igen+1}(in,3) = find((draws(in)>=cumsum([0 PSY_p(1,1:end-1)])).*(draws(in)<cumsum(PSY_p)));
            end
        end
        
         % Add education
        S_orig                  = par.PHI;
        [S1,FE1,PSY1]           = ndgrid(S_orig,FE_pos,PSY);
        polin                   = griddedInterpolant(S1,FE1,PSY1,squeeze(pol.tau0),'nearest');
        
        Ssim            = IGE.States{j_child_pos-1,igen+1}(:,1);
        FEsim           = IGE.States{j_child_pos-1,igen+1}(:,2);
        PSYsim          = IGE.States{j_child_pos-1,igen+1}(:,3);
        
        tau_sim                                 = polin(Ssim,FEsim,PSYsim);
        IGE.States{j_child_pos-1,igen+1}(:,4)   = tau_sim;
        
        Nsample    = Nsampletot;
        j_pos      = j_pos + 1;
    end
    
    %% After Children
    while j_pos <= par.Jr_pos-1      
        % Create labor and savings income
        IGE.lab{j_pos,igen} = nan(Nsample,1);
        IGE.sav{j_pos,igen} = nan(Nsample,1);
        IGE.States{j_pos+1,igen} = nan(Nsample,4);
        for educ = 1:length(EDUC)
            ind = logical((IGE.States{j_pos,igen}(:,3) == educ));
            
            S_pol            = squeeze(pol.S{j_pos,1}(:,:,educ,:));
            S_orig           = par.grids{educ,j_pos};
            Ssim            = IGE.States{j_pos,igen}(ind,1);
            
            [S1,FE1,INNO1]   = ndgrid(S_orig,FE_pos,INNO_pos);
            
            INNOsim         = IGE.States{j_pos,igen}(ind,4);
            INNO            = par.inc.z_val{educ,1}(j_pos,INNOsim)';
            FEsim           = IGE.States{j_pos,igen}(ind,2);
            FE              = par.inc.fe{educ,1}(FEsim);
            AGE_PROF        = par.inc.age_prof{educ,1}(j_pos);
            
            IGE.lab{j_pos,igen}(ind,1) = par.w*exp(AGE_PROF) .* exp(FE) .* exp(INNO);
            
            rsim            = r_sav.*(Ssim>=0) + r_debt.*(Ssim<0);
            polin             = griddedInterpolant(S1,FE1,INNO1,squeeze(S_pol));
            Spsim             = polin(Ssim,FEsim,INNOsim);
            IGE.States{j_pos+1,igen}(ind,1) = Spsim;
            
            IGE.sav{j_pos,igen}(ind,1) = rsim.*Ssim;
        end
        
        % Future States
        j_pos      = j_pos + 1;
        draws      = rand(Nsample,1);
        for educ = 1:length(EDUC)
            for inno = 1:length(INNO_pos)
                z_prob = squeeze(par.inc.z_prob{educ,1}(j_pos,inno,:));
                p_min = cumsum([0; z_prob(1:end-1)]);
                p_max = cumsum(z_prob(:));
                for inno_p = 1:length(INNO_pos)
                    ind = logical((IGE.States{j_pos-1,igen}(:,4) == inno).*(IGE.States{j_pos-1,igen}(:,3) == educ).*(draws>=p_min(inno_p)).*(draws<p_max(inno_p)));
                    IGE.States{j_pos,igen}(ind,2) =  IGE.States{j_pos-1,igen}(ind,2);
                    IGE.States{j_pos,igen}(ind,3) =  IGE.States{j_pos-1,igen}(ind,3);
                    IGE.States{j_pos,igen}(ind,4) =  repmat(inno_p,sum(ind),1);
                end
            end
        end
    end
    while j_pos <= par.Jd_pos-1      
        % Create labor and savings income
        IGE.lab{j_pos,igen} = nan(Nsample,1);
        IGE.sav{j_pos,igen} = nan(Nsample,1);
        IGE.States{j_pos+1,igen} = nan(Nsample,3);
        for educ = 1:length(EDUC)
            ind = logical((IGE.States{j_pos,igen}(:,3) == educ));
            
            S_pol            = squeeze(pol.S{j_pos,1}(:,:,educ));
            S_orig           = par.grids{educ,j_pos};
            Ssim            = IGE.States{j_pos,igen}(ind,1);
            
            [S1,FE1]   = ndgrid(S_orig,FE_pos);
            
            FEsim           = IGE.States{j_pos,igen}(ind,2);
            
            IGE.lab{j_pos,igen}(ind,1) = ret_rep(par,FEsim,educ,options);
            
            rsim              = r_sav.*(Ssim>=0) + r_debt.*(Ssim<0);
            polin             = griddedInterpolant(S1,FE1,squeeze(S_pol));
            Spsim             = polin(Ssim,FEsim);
            IGE.States{j_pos+1,igen}(ind,1) = Spsim;
            
            IGE.sav{j_pos,igen}(ind,1) = rsim.*Ssim;
        end
        
        % Future States
        j_pos      = j_pos + 1;
        for educ = 1:length(EDUC)
            ind = logical((IGE.States{j_pos-1,igen}(:,3) == educ));
            IGE.States{j_pos,igen}(ind,2) =  IGE.States{j_pos-1,igen}(ind,2);
            IGE.States{j_pos,igen}(ind,3) =  IGE.States{j_pos-1,igen}(ind,3);
        end
    end
end

% %% Create additional distributions for moments
% j_par_pos    = par.Jc_pos + 3;
% j_ch_pos     = par.Je2_pos;
% 
% trans_fe_ige = nan(length(FE_pos),length(FE_pos));
% for ife_par = 1:length(FE_pos)
%     ind_par = logical((IGE.States{j_par_pos,1}(:,2) == ife_par));
%     tot_par = sum(ind_par);
%     for ife_ch = 1:length(FE_pos)
%         ind_ch = logical(ind_par.*(IGE.States{j_ch_pos,2}(:,2) == ife_ch));
%         tot_ch = sum(ind_ch);
%         trans_fe_ige(ife_par,ife_ch) = tot_ch/tot_par;
%     end
% end
% 
% trans_educ_fe_ige = nan(length(EDUC),length(FE_pos),length(FE_pos));
% for educ_par = 1:length(EDUC)
%     for ife_par = 1:length(FE_pos)
%         ind_par = logical((IGE.States{j_par_pos,1}(:,3) == educ_par).*(IGE.States{j_par_pos,1}(:,2) == ife_par));
%         tot_par = sum(ind_par);
%         for ife_ch = 1:length(FE_pos)
%             ind_ch = logical(ind_par.*(IGE.States{j_ch_pos,2}(:,2) == ife_ch));
%             tot_ch = sum(ind_ch);
%             trans_educ_fe_ige(educ_par,ife_par,ife_ch) = tot_ch/tot_par;
%         end
%     end
% end
% 
% trans_dist = nan(length(EDUC),length(FE_pos));
% for educ_par = 1:length(EDUC)
%     for ife_par = 1:length(FE_pos)
%         ind_par = logical((IGE.States{j_par_pos,1}(:,3) == educ_par).*(IGE.States{j_par_pos,1}(:,2) == ife_par));
%         tot_par = sum(ind_par);
%         
%         tot_trans = sum(IGE.States{j_child_pos-1,2}(ind_par,1));
%         trans_dist(educ_par,ife_par) = tot_trans/tot_par;
%     end
% end


end