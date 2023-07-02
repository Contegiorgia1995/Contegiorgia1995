%% Start
datestr(now)
pwd
fprintf('sobolmin = %i\n',sobolmin);
warning('off','MATLAB:nearlySingularMatrix')

%% Options
options.print_moms      = 'Y';          % display moments
options.timer_on    	= 'N';          % display timer
options.model_iter_on   = 'Y';          % display model solve iterations
options.calib_steps     = 'Y';          % display calibration steps
options.save_mat    	= 'N';          % save policy functions + LE_dist
options.save_guess	    = 'Y';          % save guess
options.save_dist 	    = 'Y';          % save distribution
options.start_Calibrate = 'N';          % Default  = N: use standard initial guess; Y: load specific initial guess
options.maxiter         = 25;           % max iter for policy functions
options.tolV            = 1e-3;         % tolerance for policy functions
options.maxiterD        = 100;           % max iter for demographics
options.tolD            = 1e-6;         % tolerance for demographics
options.Nsim            = 1e5;          % Number of simulations
options.exerciseSS.On   = 'N';          % Default = N; if = Y: solve SS exercise (different wage and savings options)
options.ComputeMoments  = 'Y';          % Y: compute moments
options.transfers       = 'On';         % 'On': with OAS; 'Reduced'; 'Off'
options.ChildCost       = 'OppC'; 		% 'OC': opportunity cost or 'Constant'
options.Fertility       = 'Endo';       % 'Endo': endogenous fertility; 'Exo': constant fertility
options.Ck              = 'Yes';        % 'Yes': children consume at home; 'No'
options.ParentTrans     = 'Endo';       % Endo: Benchmark case with parents transfers; 'Exo': exogenously given
options.ParentTransVal  = .05;          % Value for exogenous trasnfers
options.AdultRisk       = 'Y';          % Y: benchmark case; N: no adult risk
options.ExPolRetirement = 'N';          % N: benchmark case; Y: Exercise Retirement Policy
options.PolRetMult      = 1;            % multiplier of retirment benefit in first bracket. Benchmark = 1;
options.ComputeOtherMus = 'N';
options.ComputeLE_age   = 'N';
options.equilibrium     = 'partial';    % 'partial' or 'general'
options.MeanRevEx       = 'Off';
options.GridH0orig      = 'N';
options.psy_avg         = 'N';
options.fert_transfer   = 'N';
options.init_transfer   = 'N';
options.patient         = 'N';
options.noborrowingwedge = 'N';
options.solve_h0        = 'Y';



%% MODEL
%1. Initiate storage
res_x                   = [];
res_MOM                 = [];
error_x                 = [];
% 2. Create 'result'.mat files
parentpath              = cd(cd('..'));
filename                = strcat(parentpath,'/RESULTS/','results_',num2str(options.n_fold,'%i'));
save(filename,'res_x','res_MOM');
clear res*

parentpath              = cd(cd('..'));
filename                = strcat(parentpath,'/RESULTS/','errors_',num2str(options.n_fold,'%i'));
save(filename,'error_x');
clear error_x

parentpath              = cd(cd('..'));
filename_diary          = strcat(parentpath,'/RESULTS/','diary_',num2str(options.n_fold,'%i'),'.txt');
diary(filename_diary);

% 3. Set bounds on sobol search
options.xnames          = {'gamman', 'ln_level', 'sigmah0', 'psy_col' , 'cor_psy', 'psy_hs', 'R_hs_curv', 'R_col_curv'  , 'R_hs_level'  ,'R_col_level' , 'Child_C_Tot', 'child_cost_inc_curv'   , 'beta_inc'};
%xlb                     = [ 0.140  , 0.470     , 0.180    , 18.000     , 1.400    , 1.800  , 0.10      , 0.10           , 0.020         , 0.01         , 0.400        , 0.500                   , 0.12      ];
%xub                     = [ 0.250  , 0.600     , 0.350    , 30.000    , 2.200    , 5.000   , 0.6       , 0.7            , 0.200          , 0.20         , 0.700        , 0.800                   , 0.21      ];  
  
xlb                     = [ 0.170  , 0.700     , 0.200    , 20.000     , 1.400    , 2.000  , 0.10      , 0.10           , 0.080         , 0.01         , 0.450        , 0.500                   , 0.14      ];
xub                    = [ 0.300  , 0.820     , 0.350    , 35.000    , 2.900    , 4.500   , 0.60      , 0.70           , 0.200         , 0.20         , 0.750        , 0.800                   , 0.25      ];  

switch options.patient
    case 'Y'
        xlb                    = [ 0.100  , 0.450     , 0.200    , 20.000     , 1.400    , 3.000  , 0.10      , 0.10           , 0.080         , 0.01         , 0.550        , 0.500                   , 0.12      ];
        xub                    = [ 0.250  , 0.650     , 0.350    , 50.000    , 2.900    , 8.000   , 0.60      , 0.70           , 0.200         , 0.20         , 0.850        , 0.800                   , 0.25      ];  
end

p                       = sobolset(numel(xlb));
sobolmax                = sobolmin+sobol_size-1;

%% Initialise sobol sequence and then search
for sss = sobolmin:sobolmax
    diary off; diary on
    sobshock        = p(sss,1:numel(xlb));
    options.soboln  = sss;
    x               = xlb+sobshock.*(xub-xlb);
    fprintf('start sobol: %i, time:%3.2f sec \n',sss,toc);
    fprintf('gamman = %3.4f, ln_level = %3.4f, mu = 1.00, sigmah0 = %3.4f, psy max = %3.4f, psy cor = %1.3f, psy hs = %1.3f, R hs curv = %3.4f, R col curv = %3.4f, R hs level = %3.4f, R col level = %3.4f,  Child_C_Tot= %1.3f,  child_cost_inc_curv= %1.3f,  BetaInc= %1.3f   \n',x(1),x(2),x(3),x(4),x(5),x(6),x(7),x(8),x(9),x(10),x(11),x(12),x(13));
    try
        calibrate(x,options);
    catch
        fprintf('ERROR HERE: Try-Catch at start_model \n')
        error_type  = 1;
        parentpath  = cd(cd('..'));
        filename    = strcat(parentpath,'/RESULTS/','errors_',num2str(options.n_fold,'%i'));
        Y           = load(sprintf('%s.mat',filename));
        aux         = [reshape(x,1,numel(x)) options.n_fold options.soboln error_type];
        error_x     = [ Y.error_x;        reshape(aux,1,numel(aux))];
        save(filename,'error_x');
    end
    diary off; diary on
end

filedist            = strcat('dist_ss.mat');
fileguess           = strcat('guess.mat');
delete(filedist,fileguess);

fprintf('CONGRATULATIONS!!\n');
