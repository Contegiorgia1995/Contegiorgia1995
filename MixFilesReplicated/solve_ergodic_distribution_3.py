# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 15:43:01 2022

@author: Giorgia
"""
def solve_ergodic_distribution3(par2,pol_dense,mu,pop_final,Q_ergo,options):
    
# Find distribution for a given cohort
    if options.timer_on == 'Y':
        print(f'\n Find ergodic distribution for cohort \n \n')


## distribute policy functions
    pol_mat = copy.deepcopy(pol_dense.__dict__)##### this will be useful!!! converse SImpleNameSPce to dictionary
    pol_mat = pol_mat.keys()
    for i in range (0,len(pol_mat)):
        eval([literal_eval(pol_mat[i]) '= pol_dense.(cell2mat(pol_mat(i)));']);

    clear pol_dense
    
    dist_age  = np.zeros([par2.Jd_pos-1,1])
    mu_cs     =  {}

    for i in range (0, par2.Jd_pos-1):
        mu_cs[i] = [[]]
    
    ## Age j=Jc + 3: Jr-1: 
    for j_pos in range( par2.Jc_pos+par2.Je1_pos-1,par2.Jr_pos-1):
        mu[j_pos],Q_ergo[j_pos-1] = trans_work(par2,mu[j_pos-1],S[j_pos-1],j_pos,Q_ergo[j_pos-1])
        if options.timer_on == 'Y':
            tot = np.sum(np.ravel(mu[j_pos][:]))
                print(f'age:{par2.age[j_pos]}, pop:{tot}, time: {toc} sec \n')

    
    ## age j = Jr: draw of transfers
    j_pos = par2.Jr_pos-1
    mu[j_pos] = trans_workT(par2,mu[j_pos-1],S[j_pos-1],j_pos)
        if options.timer_on == 'Y':
            tot = np.sum(np.ravel(mu[j_pos][:]))
                print(f'age:{par2.age[j_pos]}, pop:{tot}, time: {toc} sec \n')
    
    
    ## Age j=Jr + 1: state are transfers
    for j_pos in range (par2.Jr_pos,par2.Jd_pos):
        mu[j_pos] = trans_workT(par2,mu[j_pos-1],S[j_pos-1],j_pos)
        if options.timer_on == 'Y':
            tot = np.sum(np.ravel(mu[j_pos][:]))
                print(f'age:{par2.age[j_pos]}, pop:{tot}, time: {toc} sec \n')

    
    pop_g   = pop_final**(1/par2.Jc)
    sum_mu  = 0
    for j_pos in range (0,par2.Jd_pos-1):
        yy        = par2.age[j_pos]
        mu_cs[j_pos] = mu[j_pos]/(pop_g**yy)
        sum_mu       = sum_mu + np.sum(np.ravel(mu_cs[j_pos][:])
    
    for j_pos in range (0,par2.Jd_pos-1):
        mu_cs[j_pos] = mu[j_pos]/sum_mu
        dist_age[j_pos] = np.sum(np.ravel(mu_cs[j_pos][:])

    
    return function [mu,mu_cs,dist_age] 
