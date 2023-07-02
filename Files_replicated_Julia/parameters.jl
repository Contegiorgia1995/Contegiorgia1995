using Parameters
using Distributions
using UnPack
using Infinity
using Statistics
@with_kw mutable struct party5{F <: Float64, I<:Int64, MI<: Matrix{Int64}, MF<:Matrix{Float64}, SR<: StepRange{Int64, Int64}, FU <:Function, Struct<:Struct, V<:Vector}

    mean_inc :: F           = 70179.45
                           # mean income at age 40-43 (to compute OAS transfers)
    w :: I                  = w ################defining in function parameters()
    time_period::I          = 2
    age::SR                 = 0:time_period:80
    NJ::I                   = length(age) 
    Je1::I                  = 16
    Je1_pos::I              = findfirst(item -> item == Je1,age)
    Je2::I                  = 18                                       # College age
    Je2_pos::I              = findfirst(item -> item == Je2,age)
    Jc::I                   = 28
    Jc_pos::I               = findfirst(item -> item == Jc,age)
    Jr::I                   = 66
    Jr_pos::I               = findfirst(item -> item == Jr,age)
    Jd::I                   = 80
    Jd_pos::I               = findfirst(item -> item == Jd,age)
    educ :: MI              = educ
    dist_educ:: MF          = [0.094 0.601 (1-.601-.094)]
    N_markov:: I            = 5
    N_fe:: I                = 5
    sigmah0 ::F             = par_est.sigmah0
    inc :: Struct           = inc
    Reducs:: MF             = Reducs ################defining in function parameters()
    beta_inc ::F            = par_est.beta_inc
    sigmah0     ::F         = par_est.sigmah0
    meanh0  ::F             = par_est.meanh0
    H0 ::  MF               = H0 ################defining in function parameters()
    PG :: Function                  = PG ################defining in function parameters()
    p_model_data            = p_model_data
    beta                    = beta
    gammac                  = 0.5
    gamman                  = par_est.gamman
    lambdan::F              = par_est.ln_level 
    pe2::F                  = pe2
    pe1::F                  = pe1
    psy_val_hs::F           = psy_val_hs
    psy_val_col::F          = psy_val_col      
    psy_prob::F             = psy_prob
    r::F
    Nanny_oppC :: I         = par_est.Nanny_oppC
    Child_Ccurv :: F        = par_est.Child_Ccurv
    Child_C_Tot ::F         = par_est.Child_C_Tot
    child_cost_inc_curv  ::F= par_est.child_cost_inc_curv
    r_sav::F                = r_sav 
    int_iota::F             = int_iota
    r_debt::F               = r_debt
    r_col::F                = r_col
    col_fact::F             = col_fact
    lambda::F               = lambda
    pn::F                   = pn
    w_college::F            = w_college
    gov_debt::F             = gov_debt
    prod_alpha::F           = prod_alpha
    prod_delta::F           = prod_delta
    prod_s::F               = prod_s
    prod_rho::F             = prod_rho
    fam_size::F             = 2
    N::V                    = N
    Ls_0::F
    Ls_1::F
    Ls::F
    grids::F
    PHI::F
    fert_trans::F
    init_trans ::F
    #
end

function parameters(par_est, options)
    educ = [1 2 3]
    mean_inc        = 70179.45
    mean_inc_age42  = 75630 
    time_period = 2
    N_fe = 5
    meanh0 = par_est.meanh0

    

    if options.exerciseSS_On == "N"
        w       = 1
    end
    
    Reducs = zeros(2,2)
    Reducs[1,:] = [par_est.R_hs_level par_est.R_hs_curv]
    Reducs[2,:] = [par_est.R_col_level par_est.R_col_curv]

    pp                  = LinRange(0.025,0.975,N_fe)

    if options.GridH0orig == "N"
        sigma_grid          = 1.5*par_est.sigmah0 # Larger grid of H0 so that it includes more points if beta_inc is large.
    else
        sigma_grid          = 1.5*options.SigmaH0orig# Larger grid of H0 so that it includes more points if beta_inc is large.
    end
    H0                      = quantile(LogNormal(log(meanh0),sigma_grid), pp) 
    H0                      = reshape(H0, (1,5))

    function PG(x::Vector{Float64})

        v                       = H0[1:end-1]
        v                       = append!(v, ∞)
        v                       = reshape(v, (1,5))
    
    
        logn_x1 = zeros(14,5)
        for i = 1:length(x)
            for j = 1:length(v)
                logn_x1[i,j] = cdf(LogNormal(log(meanh0) + par_est.beta_inc * log(x[i]), par_est.sigmah0), v[j])
            end
        end
        

        w = H0[1:end-1]
        
        w = pushfirst!(w, 0)
        w = reshape(w, (1,5))
        
        
        logn_x2 = zeros(14,5)
        for i = 1:length(x)
            for j = 1:length(w)
                logn_x2[i,j] = cdf(LogNormal(log(meanh0) + par_est.beta_inc * log(x[i]), par_est.sigmah0), w[j])
            end
        end

        PG_values = transpose(logn_x1 - logn_x2)

        return PG_values
    end

    PG(x)

    H0 = transpose(H0)
    
    ######################define inc.fe
    inc_fe           = Dict()
    for i=1:length(educ)
        inc_fe[i] = []
    end
    log_H0 = zeros(length(H0))
    for i=1:length(H0)
        log_H0[i] = log(H0[i])
    end
    log_H0_1 = zeros(length(H0))
    for i=1:length(H0)
        log_H0_1[i] = log(H0[i] + Reducs[1,1] * H0[i]^Reducs[1,2])
    end
    log_H0_2 = zeros(length(H0))
    for i=1:length(H0)
        log_H0_2[i] = log(H0[i] + Reducs[2,1] * H0[i]^Reducs[2,2])
    end

    inc_fe[1]      = log_H0 
    inc_fe[2]      =  log_H0_1
    inc_fe[3]      =  log_H0_1 
    ###############################

    mean_inc_data           = time_period*mean_inc_age42
    mean_inc_model          = time_period

    p_model_data            = mean_inc_model / mean_inc_data

    beta_annual         = .975
    if options.patient == "Y"
        beta_annual         = .99   
    end
    beta                = beta_annual^time_period
    #gammac          = 0.5
    #gamman          = par_est.gamman                           
    #lambdan         = par_est.ln_level


    pe2                 = time_period * 2 * 6588 * p_model_data

    pe1                 = pe2*0.09

    if options.exerciseSS_On == "Y"
        if options.exerciseSS_change_pe == "Y"

            pe2             = par.pe2 * options.exerciseSS_p_mult
            pe1             = par.pe1 * options.exerciseSS_p_mult
        end
    end

    pe2_total           = (4/par.time_period)*pe2
    pe1_total           = (4/par.time_period)*pe1
    
    gridpsy             = 100

    if  options.exerciseSS.On == "Y"
        if options.exerciseSS.change_psy == "Y"
            par_est.psy_hs  = par_est.psy_hs * options.exerciseSS.w_ss^options.psy_adj;
            par_est.psy_max = par_est.psy_max * options.exerciseSS.w_ss^options.psy_adj;
        end
    
    psy_val_hs      = LinRange(0,par_est.psy_hs,gridpsy)
    psy_val_col     = LinRange(0,par_est.psy_max,gridpsy)  
    
    psy_prob        = zeros(3,gridpsy)

    Ngrid               = gridpsy+1

    curv_psy            = par_est.psy_cor
    grid_cdf            = LinRange(0^(1/curv_psy),1^(1/curv_psy),Ngrid).^curv_psy
    psy_prob[1,:]       = grid_cdf[2:end]-grid_cdf[1:end-1]
    # HS grad
    curv_psy            = 1
    grid_cdf            = LinRange(0^(1/curv_psy),1^(1/curv_psy),Ngrid).^curv_psy
    psy_prob[2,:]       = grid_cdf[2:end]-grid_cdf[1:end-1]
    #CO grad: more prob in lower values
    curv_psy            = par_est.psy_cor
    grid_cdf            = LinRange(0^(1/curv_psy),1^(1/curv_psy),Ngrid).^curv_psy
    grid_cdf_inv        = grid_cdf[2:end]-grid_cdf[1:end-1]
    psy_prob[3,:]       = grid_cdf_inv[gridpsy:-1:1]

    prob_educ_data      = [9.4 60.1 30.5]./100
    prob_psy            = prob_educ_data * psy_prob
    avg_psy             = [prob_psy*psy_val_hs; prob_psy*psy_val_col]

    if options.psy_avg  == "Y"
        gridpsy         = 2
        psy_val_hs      = [avg_psy[1]*0.995; avg_psy[1]]
        psy_val_col     = [avg_psy[2]*0.995; avg_psy[2]]
        
        psy_prob        = repeat(0.5,(3,gridpsy))

    # 7. Partial Equilibrium prices and taxes
    # 7.1 Wages: internally, to match mean income and education shares
    # 7.2 Interest rate
    # prime        = 1/beta_annual - 1;
    prime               = 3/100 # Smets and Wouters AER 2007
    r                   = (1+prime)^time_period-1
    r_sav               = (1+prime)^time_period-1            # Roys, Seshadri - 2014 XXX CONSIDER ADDING TAXES
    par_est.int_iota        = 0.1 # XXX Add to estimation? XXX
    if options.noborrowingwedge == 'Y'
            par_est.int_iota = 0.01
    end
    r_debt              =  ((1+prime+par_est.int_iota)^time_period - 1)
    r_col               = ((1+prime+0.00925)^time_period - 1); #See /DK/data/Data_US/Education/Loans/College Loans.xls
    col_fact            = r_col/(1-(1+r_col)^(-10))*(1-(1+r_debt)^(-10))/r_debt # Assume debt is repaid in 20 years, i.e. 10 periods.
    # switch options.noborrowingwedge
    #    case 'Y'
    #        par.r_col    = ((1+prime+0.0)^par.time_period - 1); %See /DK/data/Data_US/Education/Loans/College Loans.xls
    #        par.col_fact = 1;
    # end
    # Borrowing contraint: By education
    borr                = w * [10000 24000 34000].* p_model_data #self-reported limits on unsecured credit by family type from the SCF. Based on Abbot et al (2016).
    borr_sch            = w * [0 0 2*col_fact*23000].* p_model_data #Less borrowing if not finished college

    # 7.3 Government social security taxes
    if options.ExPolRetirement == "N"
            lambda      = .124                                     # Social Security Ppayroll Tax, Krueger, Ludwig - 2015
    elseif options.ExPolRetirement == "Y"
            lambda      = options.PolRetTaxes
    end     
    pn                  = 0.54 * time_period * mean_inc * p_model_data   # price of childcare, source: Folbre 2008, page 129. Wage of child care in 2000 = $7.43, mean wage = 13.74, 7.43/13.74=.54
    w_college           = 0.56 * w                             # wage at college, eg inc = w_college * w * h; source: IPUMS
    if options.exerciseSS.On =="Y"
    #         par.pn				= par.pn  * options.exerciseSS.p_mult ;  
    end

    gdp_pc              = 44308;
    gov_debt            = 0.2 * gdp_pc * p_model_data
                    
    # Source: AGMV 2013
    prod_alpha          = 0.350;
    prod_delta          = 1-(1-0.0650)^par.time_period;
    prod_s              = [0.160 0.390 0.450];
    prod_rho            = 0.680;
    
    ## 9. Grids for savings and discrete choices
    # 8.1 Number of children
    fam_size            = 2                                     #each child of the hh represents par.fam_size children
    if options.Fertility == "Endo"
            N       = [0:3]
    elseif options.Fertility == "ENdo"
            N       = options.Fertility_Exo_N/fam_size
    end
    # N = 0: no child
    # N = 1: fam_size children
    # N = 2: 2*fam_size children
    # etc.

    # 8.4 Savings

end




