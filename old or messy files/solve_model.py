# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 13:35:40 2022

@author: Giorgia
"""

from types import SimpleNamespace
import time
import itertools as it
import numpy as np
from scipy import optimize
from scipy.interpolate import RegularGridInterpolator
from scipy.io import loadmat
import os



model = SimpleNamespace() #model
model.par = SimpleNamespace() #param
model.sol = SimpleNamespace() 
model.options = SimpleNamespace()
options = model.options
par = model.par
par_est = model.par
model.inc = SimpleNamespace()
par = model.par
par.inc = model.par.inc
inc = model.inc
os.chdir('/Users/Giorgia/Dropbox/ReplicationProject/DairuchPaper/DairuchFiles/run_generic')

guess = loadmat('/Users/Giorgia/Dropbox/ReplicationProject/DairuchPaper/DairuchFiles/run_generic/guess.mat')
guess.keys()
#type(guess['guess']),guess['guess'].shape
#type(annots['guess'][0][0]),annots['guess'][0][0].shape
tau_guess = guess['guess'][0][0]['tau']
Vc0_guess = guess['guess'][0][0]['Vc0']

def tic():
    #Homemade version of matlab tic and toc functions
    import time
    global startTime_for_tictoc
    startTime_for_tictoc = time.time()

def toc():
    import time
    if 'startTime_for_tictoc' in globals():
        print("Elapsed time is " + str(time.time() - startTime_for_tictoc) + " seconds.")
    else:
        print("Toc: start time not set")
tic()
toc()

options.start_Calibrate = 'N'

function [par,pol,guess,flags] = solve_model(par,options)

def solve_model(par,options):
##Solve model - April 2015

# Notes:
# Time length: 2 years

    if options.start_Calibrate == 'N':
	# Load initial guess for tau and Vc0
        filename = str('guess')
        file = str(filename + '.mat')
        file_path = r'/Users/Giorgia/Dropbox/ReplicationProject/DairuchPaper/DairuchFiles/run_generic/'+file
        exist = os.path.isfile(file_path)

        if exist:
            guess = loadmat(file)
            #guess.keys()
            tau = guess['guess'][0][0]['tau']
            Vc0 = guess['guess'][0][0]['Vc0']
            Vc0 = np.reshape(Vc0,[10,20,20])

        else:
            Vc0   = np.repeat(np.zeros([1,1]),1)

    elif options.start_Calibrate == 'Y':
		## Load initial guess for tau and Vc0 from HPC
        Vc0   = options.guess.Vc0


# If guess is of different size of current grid start new guess
    Vc0_aux  = np.zeros([len(par.psy_val_hs),len(par.PHI),par.N_fe])
    if Vc0.ndim != Vc0_aux.ndim:
        Vc0  = np.zeros([len(par.psy_val_hs),len(par.PHI),par.N_fe])
    elif sum(len(Vc0[:,0,0]) - len(Vc0_aux[:,0,0]) + len(Vc0[0,:,0]) - len(Vc0_aux[0,:,0]) + len(Vc0[0,0,:]) - len(Vc0_aux[0,0,:])) != 0:########### Vc0 is (20,20,10), Vc0_aux is (100,5,20). I want to get -75=  
        Vc0  = np.zeros([len(par.psy_val_hs),len(par.PHI),par.N_fe])
############################################################################

## Empty matrix
    V   = {}
    for i in range (0, par.Jd_pos):
        V[i] = [[]]
    
    C   = {}
    for i in range (0, par.Jd_pos):
        C[i] = [[]]
    Ck   = {}
    for i in range (0, par.Jd_pos):
        Ck[i] = [[]]
        
    S   = {}
    for i in range (0, par.Jd_pos):
        S[i] = [[]]
        
    Tp   = {}
    for i in range (0, par.Jd_pos):
        Tp[i] = [[]]
    
    ## Dead: V_Jd = 0
    V[par.Jd_pos-1] = np.zeros([len(par.educ),int(par.Ls[-1]),par.N_fe])
    C[par.Jd_pos-1] = np.zeros([len(par.educ),int(par.Ls[-1]),par.N_fe])
    S[par.Jd_pos-1] = np.zeros([len(par.educ),int(par.Ls[-1]),par.N_fe])

    ## Retirement: for ages j=Jr:Jd-1
    for j in range(par.Jd_pos-2,par.Jr_pos-2,-1):
        V[j],C[j],S[j] = prob_ret(par,options,V[j+1],C[j+1],j)
        if options.timer_on == 'Y':
                print(f'j: {par.age[j]}, time: {toc()} sec\n')
    
    
