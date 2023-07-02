

import time
import itertools as it
import numpy as np
from scipy import optimize
from scipy import stats
from __future__ import print_function
from numpy import *
from statistics import NormalDist
import scipy.interpolate as spi
from scipy.io import loadmat
import os
from types import SimpleNamespace
import numpy as np

import types, copy


model = SimpleNamespace() #model
model.par = SimpleNamespace() #param
model.sol = SimpleNamespace() 
model.options = SimpleNamespace()
model.inc = SimpleNamespace()
model.par.inc = SimpleNamespace()
model.par_temp = SimpleNamespace()
model.guess = SimpleNamespace()
model.options.guess = SimpleNamespace()
model.pol = SimpleNamespace()
model.flags = SimpleNamespace()
model.par2 = SimpleNamespace()
model.pol_dense = SimpleNamespace()

options = model.options
par = model.par
par_temp = model.par_temp
par_est = model.par
par2 = model.par2
inc = model.inc
par.inc = model.par.inc
guess = model.guess
options.guess = model.options.guess 
pol = model.pol
flags = model.flags
pol_dense = model.pol_dense 
## Options

options.exerciseSS_On   = 'N'         # Default = N; if = Y: solve SS exercise (different wage and savings options)
options.Fertility       = 'Endo'       # 'Endo': endogenous fertility; 'Exo': constant fertility
options.Fertility_Exo_N = 2.15       # 'Exo' case: number of children
options.ParentTrans     = 'Endo'      # Endo: Benchmark case with parents transfers; 'Exo': exogenously given
options.ParentTransVal  = 0.05         # Value for exogenous trasnfers
options.ExPolRetirement = 'N'         # N: benchmark case; Y: Exercise Retirement  icy
options.GridH0orig      = 'N'        # Force initial grid of H0 -- To use when doing sensitivity analysis of parameters
options.patient         = 'N'         # 'N': Baseline Discount Factor; 'Y': Higher beta (lower discounting)
options.noborrowingwedge= 'N'         # 'N': Baseline Wedge between borrowing and lending; 'Y': Smaller wedge
options.fert_transfer   = 'N'
options.init_transfer   = 'N'
options.psy_avg         = 'N'
options.AdultRisk       = 'Y' #from start
options.timer_on    	= 'N';          # display timer
options.start_Calibrate = 'N'
options.AdultRisk       = 'Y'
options.Ck              = 'Yes'
options.maxiter         = 25
options.tolV            = 1e-3
options.solve_h0        = 'N'
options.maxiterD        = 100;  # max iter for demographics
options.AdultRisk       = 'Y'         
#inc.age_prof = cell(length(EDUC),1); #from file par_income_process
#inc.fe_pos   = np.arange(1,N_fe,N_fe)





x = np.matrix([[0, 0.24776039834302887588],
              [1, 0.75658886515041290366],  
              [2,  0.25397448399310357248],
              [3, 25.71397930469017012456],
              [4, 1.7323110070687675055],
              [5, 3.9825431046569619298],
              [6, 0.55750684703128172703],
              [7, 0.23634828161700616178],
              [8, 0.16301011354749064819],
              [9, 0.10918570913971439862],
              [10, 0.66823638808333596373],
              [11, 0.6646077270329715514],
              [12, 0.18544353931303594885],
              [13, 0.47056974917650218337]])
                   

x_val = np.array([ 0.24776039834302887588,
             0.75658886515041290366,
             0.25397448399310357248,
             25.71397930469017012456,
             1.7323110070687675055,
             3.9825431046569619298,
             0.55750684703128172703,
             0.23634828161700616178,
             0.16301011354749064819,
             0.10918570913971439862,
             0.66823638808333596373,
             0.6646077270329715514,
             0.18544353931303594885,
             0.47056974917650218337])


par_est.gamman      = x[0,1]
par_est.ln_level    = x[1,1]
par_est.mu          = 0.5
par_est.sigmah0     = x[2,1]
par_est.psy_max     = x[3,1]
par_est.psy_cor     = x[4,1]
par_est.psy_hs      = x[5,1]
par_est.R_hs_curv   = x[6,1]
par_est.R_col_curv  = x[7,1]
par_est.R_hs_level  = x[8,1]
par_est.R_col_level = x[9,1]
par_est.Nanny_oppC  = 1
par_est.Child_Ccurv = 0.6445
# par_est.Child_Ccurv = 1.000
par_est.Child_C_Tot = x[10,1] #Multiply by (1-tax) in child cost function to obtain parameter from paper
par_est.child_cost_inc_curv = x[11,1]
par_est.beta_inc    = x[12,1]
#if options.solve_h0 = 'N'
par_est.meanh0      = x[13,1]
#par.meanh0 = 2
#par.sigmah0 = 0.5

##############  PROBNLEMS:
    #parameters:  par.PG becasue of lognormcdf
                # options.init_trans_size
    #prob_workT:  C[inno,educ,not_feasible,ife] = 0
                # V[inno,educ,not_feasible,ife] = (-(10**(5/gammac))**(1-gammac))/(1-gammac),
                # S[inno,educ,not_feasible,ife] seems to work but would also work if Sp[inno,educ,not_feasible,ife]= Spgrid[0] was not commented
