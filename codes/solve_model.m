function [par,pol,guess,flags] = solve_model(par,options)
%% Solve model - April 2015

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes:
% Time length: 2 years
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch options.start_Calibrate 
	case 'N'
		%% Load initial guess for tau and Vc0
		filename = strcat('guess');
		if  exist(sprintf('%s.mat',filename),'file') == 2
			load(sprintf('%s.mat',filename));
			Vc0   = guess.Vc0;
		else
			Vc0   = zeros(1,1,1);
		end
	case 'Y'
		%% Load initial guess for tau and Vc0 from HPC
		Vc0   = options.guess.Vc0;
end

% If guess is of different size of current grid start new guess
Vc0_aux  = zeros(length(par.PHI),par.N_fe,length(par.psy_val_hs));
if ndims(Vc0)~= ndims(Vc0_aux)
    Vc0  = zeros(length(par.PHI),par.N_fe,length(par.psy_val_hs));
elseif sum(size(Vc0) - size(Vc0_aux)) ~= 0
    Vc0  = zeros(length(par.PHI),par.N_fe,length(par.psy_val_hs));
end

%% Empty matrix
V   = cell(par.Jd_pos,1);
C   = cell(par.Jd_pos,1);
Ck  = cell(par.Jd_pos,1);
S   = cell(par.Jd_pos,1);
Tp  = cell(3,1);

%% Dead: V_Jd = 0
V{par.Jd_pos} = zeros(par.Ls(end),par.N_fe,length(par.educ));
C{par.Jd_pos} = zeros(par.Ls(end),par.N_fe,length(par.educ));
S{par.Jd_pos} = zeros(par.Ls(end),par.N_fe,length(par.educ));

%% Retirement: for ages j=Jr:Jd-1
for j = par.Jd_pos-1:-1:par.Jr_pos
    [V{j},C{j},S{j}] = prob_ret(par,options,V{j+1},C{j+1},j);
    switch options.timer_on
        case {'Y'}
            fprintf('j: %i, time: %3.0f sec\n',par.age(j),toc);
    end
end

%% Working: for age j=Jr-1: 
j = par.Jr_pos-1;
[V{j},C{j},S{j}] = prob_workT(par,options,V{j+1},C{j+1},j);
switch options.timer_on
    case {'Y'}
        fprintf('j: %i, time: %3.0f sec\n',par.age(j),toc);
end

%% Working old: for ages j=Jc+Je1 : Jr -2
for j = par.Jr_pos-2:-1:par.Jc_pos+par.Je1_pos-1
    [V{j},C{j},S{j}] = prob_work_old(par,options,V{j+1},C{j+1},j);
    switch options.timer_on
        case {'Y'}
            fprintf('j: %i, time: %3.0f sec\n',par.age(j),toc);
    end
end

%% Iterate over years of educ and Vchild0
iter = 0; distV = 1;
while (distV>options.tolV && iter <=options.maxiter)
    ticiter = tic;
    %% Transfer money to children, age j = Jc+Je1
    j = par.Jc_pos+par.Je1_pos-2;
    [V{j},C{j},Ck{j},S{j},PHIp,Tp{j}] = prob_work_trans(par,options,V{j+1},C{j+1},Vc0,j);
    switch options.timer_on
        case {'Y'}
            fprintf('j: %i, time: %3.0f sec\n',par.age(j),toc);
    end
    
    %% Working with Children at home: for ages j=Jc+1:Jc+Je1-3: age  = 32
    for j = par.Jc_pos+par.Je1_pos-3:-1:par.Jc_pos+1
        [V{j},C{j},Ck{j},S{j},Tp{j}] = prob_work_with_child(par,options,V{j+1},C{j+1},j);
        switch options.timer_on
            case {'Y'}
                fprintf('j: %i, time: %3.0f sec\n',par.age(j),toc);
        end
    end
    
    %% Fertility: for age j=Jc
    j=par.Jc_pos;
    [V{j},C{j},Ck{j},S{j},Np,Tp{j}] = prob_fertility(par,options,V{j+1},C{j+1},j);
    switch options.timer_on
        case {'Y'}
            fprintf('j: %i, time: %3.0f sec\n',par.age(j),toc);
    end
    
    %% Work young: for ages j=Je2+1:Jc-1
    for j = par.Jc_pos-1:-1:par.Je2_pos+2
        [V{j},C{j},S{j}] = prob_work_young(par,options,V{j+1},C{j+1},j);
        switch options.timer_on
            case {'Y'}
                fprintf('j: %i, time: %3.0f sec\n',par.age(j),toc);
        end
    end
    
    %% Solve Education choice
    j = par.Je2_pos+1;
    [Ve,Ce,Se,Vc1,tau1] = prob_educ(par,options,V{j+1},C{j+1});
    
    
    difV  = abs(Vc1-Vc0);
    distV = max(difV(:));
    tau0  = tau1;

%     pace    = 0.25;
    pace    = 0.9;
    Vc0 = pace*Vc1+(1-pace)*Vc0;
    

    iter = iter+1;
    switch options.model_iter_on
            case {'Y'}
                fprintf('HH iter %03i, norm Vchild : %9.4f, tfs = %03.0f sec\n', iter, distV,toc(ticiter));
    end
    
%     fprintf(' \n')
end

if (distV>options.tolV) == 1
    flags.converge = 'N';
else
    flags.converge = 'Y';
end

% keyboard

%% Output
% Policy functions:
pol.C         = C;
pol.Ce        = Ce;
pol.Ck        = Ck;
pol.S         = S;
pol.Se     	  = Se;
pol.V         = V;
pol.Ve        = Ve;
pol.Vc0       = Vc0;
pol.Np        = Np;
pol.PHIp      = PHIp;
pol.Tp        = Tp;
pol.tau0      = tau1;
pol.tau       = tau1;

% Initial guess
guess.Vc0    = Vc0;

% flags
flags.iter = iter;
end