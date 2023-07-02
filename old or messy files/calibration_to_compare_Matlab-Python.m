par.time_period = 2 %%% added myself from somewhere
par.age             = 0:par.time_period:80;
% par.NJ              = length(par.age);
% par.Je1             = 16;                                       % HS age
% par.Je1_pos         = find(par.age == par.Je1,1,'first');
% par.Je2             = 18;                                       % College age
% par.Je2_pos         = find(par.age == par.Je2,1,'first');
% par.Jc              = 28;                                       % Fertility
% par.Jc_pos          = find(par.age == par.Jc,1,'first');
% % par.Jt              = 40;                                       % OAS
% % par.Jt_pos          = find(par.age == par.Jt,1,'first');
% par.Jr              = 66;                                       % Retirement
% par.Jr_pos          = find(par.age == par.Jr,1,'first');
par.Jd              = 80;                                       % Die, source: WB - life expectancy = 79
par.Jd_pos          = find(par.age == par.Jd,1,'first');

par.educ     = [1 2 3];

par.gammac          = 0.5;

%%%%%%%%%%%%%%%%%%%%ù
prime        = 3/100; % Smets and Wouters AER 2007
par.r        = (1+prime)^par.time_period-1;
par.r_sav     = (1+prime)^par.time_period-1;            % Roys, Seshadri - 2014 XXX CONSIDER ADDING TAXES
par_est.int_iota = 0.1; % XXX Add to estimation? XXX
par.r_debt   =  ((1+prime+par_est.int_iota)^par.time_period - 1);
%%%%%%%%%%%%%%%%%%

par.grids            = cell(length(par.educ),par.Jd_pos);


N_fe         = par.N_fe;
inc.fe_pos   = 1:N_fe;
par.inc.fe           = cell(length(par.educ),1);
par.inc.fe_pos       = 
FE_pos               = par.inc.fe_pos;


