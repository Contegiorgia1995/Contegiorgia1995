#using DataFrames


using Parameters############we are here

@with_kw mutable struct options_{S<:String, I<: Int64, F<: Float64}
    print_moms ::S     = "Y"        # display moments
    timer_on ::S   	   = "N"            # display timer
    model_iter_on ::S  = "Y"         # display model solve iterations
    calib_steps ::S    = "Y"         # display calibration steps
    save_mat   ::S 	   = "N"            # save policy functions + LE_dist
    save_guess	::S    = "Y"         # save guess
    save_dist ::S	   = "Y"         # save distribution
    start_Calibrate :: S = "N"       # Default  = N: use standard initial guess; Y: load specific initial guess
    maxiter :: I       = 25          # max iter for policy functions
    tolV  ::F          = 1e-3        # tolerance for policy functions
    maxiterD  ::I      = 100        # max iter for demographics
    tolD  ::F          = 1e-6       # tolerance for demographics
    Nsim   ::F         = 1e5        # Number of simulations
    exerciseSS_On ::S  = "N"        # Default = N; if = Y: solve SS exercise (different wage and savings options)
    ComputeMoments :: S = "Y"       # Y: compute moments
    transfers ::S      = "On"       # 'On': with OAS; 'Reduced'; 'Off'
    ChildCost  ::S     = "OppC"		# 'OC': opportunity cost or 'Constant'
    Fertility  :: S    = "Endo"     # 'Endo': endogenous fertility; 'Exo': constant fertility
    Ck  :: S           = "Yes"      # 'Yes': children consume at home; 'No'
    ParentTrans :: S   = "Endo"     # Endo: Benchmark case with parents transfers; 'Exo': exogenously given
    ParentTransVal:: F = .05;       # Value for exogenous trasnfers
    AdultRisk ::S      = "Y"        # Y: benchmark case; N: no adult risk
    ExPolRetirement::S = "N"         # N: benchmark case; Y: Exercise Retirement Policy
    PolRetMult ::I     = 1            # multiplier of retirment benefit in first bracket. Benchmark = 1;
    ComputeOtherMus::S = "N"
    ComputeLE_age ::S  = "N"
    equilibrium ::S    = "partial"    # 'partial' or 'general'
    MeanRevEx   ::S    = "Off"
    GridH0orig  ::S    = "N"
    psy_avg ::S        = "N"
    fert_transfer ::S  = "N"
    init_transfer ::S  = "N"
    patient ::S        = "N"
    noborrowingwedge::S= "N"
    solve_h0 ::S       = "Y"
    
end

options = options_()
fieldnames(typeof(options))


@with_kw mutable struct par_est_{F<:Float64}
    sigmah0::F  =  0.24776039834302887588
    psy_max::F =  0.75658886515041290366
    psy_cor::F   =  0.25397448399310357248
    psy_hs::F    = 25.71397930469017012456
    R_hs_curv::F  = 1.7323110070687675055
    R_col_curv::F  = 3.9825431046569619298
    R_hs_level::F = 0.55750684703128172703
    R_col_level::F = 0.23634828161700616178
    Nanny_oppC ::F =  0.16301011354749064819
    # par_est.Child_Ccurv 
    Child_C_Tot ::F = 0.10918570913971439862 
    Child_Ccurv ::F = 0.66823638808333596373
    beta_inc ::F = 0.6646077270329715514
    child_cost_inc_curv ::F = 0.18544353931303594885 
    #if options.solve_h0 
    meanh0::F = 0.47056974917650218337

end

par_est = par_est_()
fieldnames(typeof(par_est))


@with_kw mutable struct par_{F<:Float64}
    inc::F
    gamman::F
    ln_level::F
    mu::F
    sigmah0::F
    psy_max::F
    psy_cor::F
    psy_hs::F
    R_hs_curv::F
    R_col_curv::F
    R_hs_level::F
    R_col_level::F
    Nanny_oppC::F
    Child_Ccurv::F
    Child_C_Tot::F
    child_cost_inc_curv::F
    beta_inc::F
    meanh0::F
    mean_inc::F
    w::F
    time_period::F
    age::F
    NJ::F
    Je1::F
    Je1_pos::F
    Je2::F
    Je2_pos::F
    Jc::F
    Jc_pos::F
    Jt::F
    Jt_pos::F
    Jr::F
    Jr_pos::F
    Jd::F
    Jd_pos::F
    educ::F
    dist_educ::F
    N_markov::F
    N_fe::F
    Reducs::F
    H0::F
    p_model_data::F
    beta::F
    gammac::F
    lambdan::F
    pe2::F
    pe1::F
    psy_val_hs::F
    psy_val_col::F
    psy_prob::F
    r::F
    r_sav::F
    int_iota::F
    r_debt::F
    r_col::F
    col_fact::F
    Lambda::F
    pn::F
    w_college::F
    gov_debt::F
    prod_alpha::F
    prod_delta::F
    prod_s::F
    prod_rho::F
    fam_size::F
    N::F
    Ls_0::F
    Ls_1::F
    Ls::F
    grids::F
    PHI::F
    fert_trans::F
    init_trans ::F
end