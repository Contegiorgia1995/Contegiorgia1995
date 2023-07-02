function [mu,mu_ige0,muc,pop_final,Q_ergo] = solve_ergodic_distribution2(par2,pol_dense,options)
% Find ergodic distribution across cohorts
switch options.timer_on
    case {'Y'}
        fprintf('\n Find ergodic distribution across cohorts \n \n');
end

%% Compute demographic transition
% Find ergodic distribution across cohorts
%%%%%%%%%%%%%%%%%%%%%%%%%%%% OUTPUT: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mu{j}(a) is the pdf of state a at age j (i.e. at the begining of age j)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% distribute policy functions
pol_mat = fieldnames(pol_dense);
for i = 1 : length(pol_mat)
    eval([cell2mat(pol_mat(i)) '= pol_dense.(cell2mat(pol_mat(i)));']);
end
clear pol_dense


% Empty cell for measures at each age
mu        =  cell(par2.Jd_pos-1,1);

Q_ergo =  cell(par2.Jd_pos-1,1);
for jj = 1:par2.Jd_pos-1
    Q_ergo{jj}.created = 'N';
end


%% mu{1}: guess

% Load initial guess for tau and Vc0
filename = strcat('dist_ss');
if  exist(sprintf('%s.mat',filename),'file') == 2
    guess = load(sprintf('%s.mat',filename),'mu');
    j_pos = par2.Je1_pos-1;
    muc   = guess.mu{j_pos};
else
    muc = zeros(2,2);
end
j_pos     = par2.Je1_pos;
S_grid    = par2.grids{1,j_pos};
FE_pos    = par2.inc.fe_pos;
PSY       = par2.psy_val_hs;

% If guess is of different size of current grid start new guess
muc_aux = zeros(length(S_grid),length(FE_pos),length(PSY));
if ndims(muc) ~= ndims(muc_aux)
    % start with a guess for distribution over (s,educ)
    % Guess 1: 'symmetric': one obs for each grid point (s,educ)
    muc    = zeros(length(S_grid),length(FE_pos),length(PSY));
    nn     = length(S_grid)*length(FE_pos)*length(PSY); %size of grid to get probability measure
    muc(:) = 1/nn;
elseif sum((size(muc) -size(muc_aux))) ~= 0
    % start with a guess for distribution over (s,educ)
    % Guess 1: 'symmetric': one obs for each grid point (s,educ)
    muc    = zeros(length(S_grid),length(FE_pos),length(PSY));
    nn     = length(S_grid)*length(FE_pos)*length(PSY); %size of grid to get probability measure
    muc(:) = 1/nn;
end

pop = zeros(options.maxiterD,1);
cohort = 0; dif = 1;

switch options.timer_on
    case {'Y'}
        fprintf('cohort: %3.0f, dif:%3.7f, pop:%1.2f, time =%3.2f \n',cohort,dif,sum(muc(:)),toc);
end
while ((dif>options.tolD ) && cohort <=options.maxiterD)
    cohort = cohort + 1;
    j_pos      = par2.Je1_pos-1;
    mu{j_pos} = muc; % update initial distribution for cohort
    
    %% Age j=Je1: new states are education and innovation
    j_pos = par2.Je1_pos;
    [mu{j_pos},Q_ergo{j_pos-1}] = trans_draw(par2,mu{j_pos-1},tau0,Q_ergo{j_pos-1});
    
    switch options.timer_on
        case {'Y'}
            tot = sum(mu{j_pos}(:));
            fprintf('age:%i, pop:%1.2f, time: %3.1f sec \n',par2.age(j_pos),tot,toc)
    end
    
    %% Age j=Je1+1:Je2+3
    for j_pos = par2.Je1_pos + 1:par2.Je2_pos + 2
        [mu{j_pos},Q_ergo{j_pos-1}] = trans_educ(par2,mu{j_pos-1},Se{1,j_pos - par2.Je1_pos},Se{2,j_pos - par2.Je1_pos},Se{3,j_pos - par2.Je1_pos},j_pos,Q_ergo{j_pos-1});
        % Only need to pass policy functions given education, iid shock does
        % not matter once we condition on education choice
        switch options.timer_on
            case {'Y'}
                tot = sum(mu{j_pos}(:));
                fprintf('age:%i, pop:%1.2f, time: %3.1f sec \n',par2.age(j_pos),tot,toc)
        end
    end
%     mu{j_pos}                = squeeze(sum(mu{j_pos},3)); % Remove Psychic Cost
 
    
    %% Age j=Je2 + 2: Jc
    for j_pos = par2.Je2_pos+3:par2.Jc_pos
        [mu{j_pos},Q_ergo{j_pos-1}] = trans_work(par2,mu{j_pos-1},S{j_pos-1},j_pos,Q_ergo{j_pos-1});
        
        switch options.timer_on
            case {'Y'}
                tot = sum(mu{j_pos}(:));
                fprintf('age:%i, pop:%1.2f, time: %3.1f sec \n',par2.age(j_pos),tot,toc)
        end
    end
    
    %% Age Jc
    j_pos = par2.Jc_pos + 1;
    fert_pol = Np;
    [mu{j_pos},Q_ergo{j_pos-1}] = trans_fertility(par2,mu{j_pos-1},S{j_pos-1},fert_pol,j_pos,Q_ergo{j_pos-1},options);
    switch options.timer_on
        case {'Y'}
            tot = sum(mu{j_pos}(:));
            fprintf('age:%i, pop:%1.2f, time: %3.1f sec \n',par2.age(j_pos),tot,toc)
    end
    
    %% Age Jc+1
    for j_pos = par2.Jc_pos+2:par2.Jc_pos+par2.Je1_pos-2
        [mu{j_pos},Q_ergo{j_pos-1}] = trans_work_with_child(par2,mu{j_pos-1},S{j_pos-1},j_pos,Q_ergo{j_pos-1});
        switch options.timer_on
            case {'Y'}
                tot = sum(mu{j_pos}(:));
                fprintf('age:%i, pop:%1.2f, time: %3.1f sec \n',par2.age(j_pos),tot,toc)
        end
    end
    
    %% Age Jc+2: phi transfer + new cohort
    j_pos       = par2.Jc_pos+par2.Je1_pos-1;
    %     j_pos_child = par2.Je1_pos;
    phi_pol = PHIp;
    [mu{j_pos},muc,mu_ige0,pop(cohort),Q_ergo{j_pos-1}] = trans_work_trans(par2,mu{j_pos-1},S{j_pos-1},phi_pol,j_pos,Q_ergo{j_pos-1});
    switch options.timer_on
        case {'Y'}
            tot = sum(mu{j_pos}(:));
            fprintf('age:%i, pop:%1.2f, time: %3.1f sec \n',par2.age(j_pos),tot,toc)
    end
    
    pop_final  = pop(cohort);
    j_0_pos = par2.Je1_pos-1;
    dif = max(abs(mu{j_0_pos}(:)-muc(:)));
    % Describe muc
    mu{1}      = muc;     mu{2}      = muc;
    
    switch options.timer_on
        case {'Y'}
            fprintf('cohort: %3.0f, dif:%3.7f, pop:%1.2f, time =%3.2f \n',cohort,dif,pop(cohort),toc)
    end
end

end