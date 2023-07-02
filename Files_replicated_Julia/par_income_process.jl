function par_income_process(par::par_,options)
    
    EDUC         = par.educ
    age          = par.age[1:par.Jr_pos]
    
    N_markov     = par.N_Markov
    N_fe         = par.N_fe
    inc.inno_pos = 1:N_markov
    inc.fe_pos   = 1:N_fe
    
    age_start       = [par.Je1_pos par.Je2_pos par.Je2_pos+2]
    
    age_prof     = @(e,a) (e==1) .* ( 0.0507705  .* a  -0.0005522    .*a.^2)...
                        + (e==2) .* ( 0.0668012  .* a -0.0007312  .*a.^2)...
                        + (e==3) .* ( 0.1221627  .* a -0.0013147  .*a.^2);
    
    inc.age_prof = cell(length(EDUC),1);
    
    for ie = 1:length(EDUC)
        inc.age_prof{ie,1}      = zeros(par.Jr_pos,1);
        for ia = par.Je1_pos:par.Jr_pos
            inc.age_prof{ie,1}(ia,1) = age_prof(ie,age(ia));
        end
        inc.age_prof{ie,1}(par.Je1_pos:par.Jr_pos,1)      = inc.age_prof{ie,1}(par.Je1_pos:par.Jr_pos,1) - inc.age_prof{ie,1}(par.Je1_pos,1);
    end
    
                     
    inc.inno_rho    = [0.875153     0.961399        0.966375];    % autocorrelation of innovation for education groups
    inc.z0_var      = [0.22832      0.10134         0.0660616];  % Variance of iid shocks to innovation for education groups
    inc.inno_var    = [0.0624742    0.02341         0.0291241];   % Variance of iid shocks to innovation for education groups
    switch options.AdultRisk
        case 'N'
            inc.inno_var = 1e-6*inc.inno_var;
            inc.z0_var   = 1e-6*inc.z0_var;
    end
    inc.z_var       = zeros(length(EDUC),length(age));
    for ie = 1:length(EDUC)
        rho      = inc.inno_rho(ie);
        z0_var   = inc.z0_var(ie);
        inno_var = inc.inno_var(ie);
        inc.z_var(ie,age_start(ie)) = z0_var;
        for ia = age_start(ie)+1:par.Jr_pos
            inc.z_var(ie,ia) = rho^2 * inc.z_var(ie,ia-1) + inno_var;
        end
    end
    
    inc.z_val         = cell(length(EDUC),1);
    inc.z_prob        = cell(length(EDUC),1);
    
    #Z0 PROBABILITIES (AND VALUES) 
    for ie = 1:length(EDUC) 
        inc.z_val{ie,1}         = zeros(par.Jr_pos,N_markov);
        inc.z_prob{ie,1}        = zeros(par.Jr_pos,N_markov,N_markov);
        
        % Initial draw of z
        ia                      = age_start(ie);
        nu                      = sqrt(N_markov-1) * inc.z_var(ie,ia)^.5;
        etagrid                 = linspace(-nu, nu, N_markov)';
        inc.z_val{ie,1}(ia,:)   = etagrid;
        
        grid_med                 = .5*(etagrid(2:end)+etagrid(1:end-1));
        aux_cdf                  = normcdf(grid_med,0,inc.z_var(ie,ia)^.5);
        aux_pdf                  = [aux_cdf(1)
                                    aux_cdf(2:end)-aux_cdf(1:end-1)];
        aux_pdf                  = [aux_pdf
                                    1-sum(aux_pdf)];
    
        inc.z_prob{ie,1}(ia,:,:) = repmat(aux_pdf',N_markov,1);
        for ia = age_start(ie)+1:par.Jr_pos
            [etagrid, P]             = discretize_nonstationary_ar_rouwen(inc.inno_rho(ie),inc.z_var(ie,ia-1)^.5,inc.z_var(ie,ia)^.5, N_markov);
            inc.z_val{ie,1}(ia,:)    = etagrid;
            inc.z_prob{ie,1}(ia,:,:) = P;
        end
    end
    
    end
    return inc