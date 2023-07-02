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
import types, copy





def solve_model(par,options):
##Solve model - April 2015
################################premises
    os.chdir('/Users/Giorgia/Dropbox/ReplicationProject/DairuchPaper/DairuchFiles/run_generic')
    
    guess_ = loadmat('/Users/Giorgia/Dropbox/ReplicationProject/DairuchPaper/DairuchFiles/run_generic/guess.mat')
    
    ##type(guess['guess']),guess['guess'].shape
    ##type(annots['guess'][0][0]),annots['guess'][0][0].shape
    model.guess = SimpleNamespace()
    guess       = model.guess
    guess.tau = guess_['guess'][0][0]['tau']
    guess.Vc0 = guess_['guess'][0][0]['Vc0']
    #guess.keys()
    
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
  #############################################################actual function starts here
    r_sav = par.r_sav
    w = par.w
# Notes:
# Time length: 2 years

    if options.start_Calibrate == 'N':
	# Load initial guess for tau and Vc0
        filename = str('guess')
        file = str(filename + '.mat')
        file_path = r'/Users/Giorgia/Dropbox/ReplicationProject/DairuchPaper/DairuchFiles/run_generic/'+file
        exist = os.path.isfile(file_path)

        if exist:
            guess_ = loadmat(file)
            #guess.keys()
            #tau = guess['guess'][0][0]['tau']
            #Vc0 = guess['guess'][0][0]['Vc0']
            #Vc0 = np.reshape(Vc0,[10,20,20])###################wrong!!!!!
            Vc0 = guess.Vc0

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
    
    ## Working: for age j=Jr-1: 
    j = par.Jr_pos-2
    V[j],C[j],S[j]  = prob_workT(par,options,V[j+1],C[j+1],j);
    if options.timer_on == 'Y':
            print(f'j: {par.age[j]}, time: {toc()} sec\n')

   ## Working old: for ages j=Jc+Je1 : Jr -2
    for j in range (par.Jr_pos-3,(par.Jc_pos + par.Je1_pos-3),-1):
        V[j],C[j],S[j] = prob_work_old(par,options,V[j+1],C[j+1],j)
        if options.timer_on == 'Y':
                print(f'j: {par.age[j]}, time: {toc()} sec\n')
   
    ## Iterate over years of educ and Vchild0
    iter_ = 0
    distV = 1
    while(distV > options.tolV and iter_ <= options.maxiter):
        ticiter = tic
        ## Transfer money to children, age j = Jc+Je1
        j = par.Jc_pos+par.Je1_pos-3
        V[j],C[j],Ck[j],S[j],PHIp,Tp[j] = prob_work_trans(par,options,V[j+1],C[j+1],Vc0,j)
        if options.timer_on == 'Y':
            print(f'j: {par.age[j]}, time: {toc()} sec\n')  
            
    ## Working with Children at home: for ages j=Jc+1:Jc+Je1-3: age  = 32
        for j in range(par.Jc_pos+par.Je1_pos-4,par.Jc_pos-1,-1):
            V[j],C[j],Ck[j],S[j],Tp[j] = prob_work_with_child(par,options,V[j+1],C[j+1],j)
            if options.timer_on == 'Y':
                print(f'j: {par.age[j]}, time: {toc()} sec\n')  
##################################################################
##############################################################THINGS ARE LOOKING GOOD BUT FOR VERY NEGATIVE NUMBERS IN FIRST ROWS BOTH C AND V ARE WORNG. THE APROXIMATION MESSES THINGS UP. gegm USES splVp AS INPUT SO ALSO C AND S (OUTPUT OF GEGM) ARE WORNG
    ## Fertility: for age j=Jc
        j=par.Jc_pos-1
        V[j],C[j],Ck[j],S[j],Np,Tp[j] = prob_fertility(par,options,V[j+1],C[j+1],j)
        if options.timer_on == 'Y':
            print(f'j: {par.age[j]}, time: {toc()} sec\n')  
        
        ## Work young: for ages j=Je2+1:Jc-1
        for j in range(par.Jc_pos-2,par.Je2_pos,-1):
            V[j],C[j],S[j] = prob_work_young(par,options,V[j+1],C[j+1],j)
            if options.timer_on == 'Y':
                print(f'j: {par.age[j]}, time: {toc()} sec\n')
        
        ## Solve Education choice
        j = par.Je2_pos
        Ve,Ce,Se,Vc1,tau1 = prob_educ(par,options,V[j+1],C[j+1])
        if options.timer_on == 'Y':
            print(f'j: {par.age[j]}, time: {toc()} sec\n, iter:{iter_}, distV: {distV}')
        
        
        difV  = np.abs(Vc1-Vc0)
        distV = np.max(difV[:])
        tau0  = tau1
    
    #    pace    = 0.25;
        pace    = 0.9
        Vc0 = pace*Vc1+(1-pace)*Vc0
        
    
        iter_ = iter_+1
        #if options.model_iter_on == 'Y':
           # fprintf('HH iter %03i, norm Vchild : %9.4f, tfs = %03.0f sec\n', iter, distV,toc(ticiter));

        
    #     fprintf(' \n')

    
    if (distV>options.tolV) == 1:
        flags.converge = 'N'
    else:
        flags.converge = 'Y'

    
    # keyboard
    model.pol = SimpleNamespace()
    pol = model.pol
    ## Output
    # Policy functions:
    pol.C         = copy.deepcopy(C)
    pol.Ce        = copy.deepcopy(Ce)
    pol.Ck        = copy.deepcopy(Ck)
    pol.S         = copy.deepcopy(S)
    pol.Se     	  = copy.deepcopy(Se)
    pol.V         = copy.deepcopy(V)
    pol.Ve        = copy.deepcopy(Ve)
    pol.Vc0       = copy.deepcopy(Vc0)
    pol.Np        = copy.deepcopy(Np)
    pol.PHIp      = copy.deepcopy(PHIp)
    pol.Tp        = copy.deepcopy(Tp)
    pol.tau0      = copy.deepcopy(tau1)
    pol.tau       = copy.deepcopy(tau1)
    
    # Initial guess
    guess.Vc0    = Vc0
    
    # flags
    flags.iter_ = iter_
    
    return pol, par, guess, flags

            

parentpath1    = str('/Users/Giorgia/Dropbox/ReplicationProject/DairuchPaper/Replication/Files_replicated_Python')

filename1    = str(parentpath1 + '/RESULTS/' +'pol')

dic = {"pol1":pol1, "par":par,"iter": iter_ ,"distV":distV}

scipy.io.savemat(str(filename1 + '.mat'), {str(filename1):dic})
Y = loadmat(str(filename1 + '.mat'))####not working



print(new_dict)   
import numpy as np
np.save('/Users/Giorgia/Dropbox/ReplicationProject/DairuchPaper/Replication/Files_replicated_Python/dic.npy', dic) 

new_dict = np.load('dic.npy', allow_pickle='TRUE')
print(new_dict.item())


