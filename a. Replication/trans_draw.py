# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 19:37:04 2022

@author: Giorgia
"""
def trans_draw(par,mu0,tau0,Q_ergo_in):
# From mu0 to mup
# keyboard
# Choice: education
    j_pos     = par.Je1_pos-1
    S         = par.grids[0][j_pos]
    FE_pos    = par.inc.fe_pos
    EDUC      = par.educ
    PSY       = par.psy_val_hs
    
    
    ## 1. Choice of educ
    mu0int          = np.zeros([len(EDUC),len(PSY),len(S),len(FE_pos)])
    mu0int[0,:,:,:] = mu0    # initial distribution in extended grid, position 1
    mu0int_vec      = mu0int[:]
    
    if Q_ergo_in.created == 'N':

            # extend policy for education tau
            tau0p         = np.ones([len(EDUC),len(PSY),len(S),len(FE_pos)])
            tau0p[0,:,:,:]  = tau0
            tau0int_vec   = tau0p[:]
    
            # create transition matrix of educ: linear interpolant for policy of education
            fspaceergeduc   = fundef({'spli',EDUC,0,1})
            Qeduc           = funbas(fspaceergeduc,tau0int_vec)
    
            # create transition matrix of fixed states: Qfixed
            Qfixed              = np.kron(np.ones([len(EDUC),1]),scipy.sparse.identity(len(S)*len(FE_pos)*len(PSY)))###?????????
    
            # create aggregate transition: Q
            Q   = dprod(Qeduc,Qfixed)
            Q_ergo_out.mat1 = Q
            Q_ergo_out.created = 'Y'
    if Q_ergo_in.created == 'Y':
            Q   = Q_ergo_in.mat1
            Q_ergo_out = Q_ergo_in
    # new distribution: 
    mu_int_vec      = Q.T * mu0int_vec
    muint             = reshape(mu_int_vec,length(S),length(FE_pos),length(PSY),length(EDUC));
    
    ## 2. draw of innovation z0
    # load distribution of z0 conditional on education
    z_prob     = par.inc.z_prob
    z0_prob    = [z_prob[0][0][par.Je1_pos-1].T, z_prob[1][0][par.Je2_pos-1].T, z_prob[2][0][par.Je2_pos+1].T]      # (educ,prob)
    INNO_pos   = par.inc.inno_pos
    NZ0        = len(INNO_pos)         
    
    # vectorize distribution and policy
    mu0int_vec    = muint[:]
    
    if Q_ergo_in.created == 'N':
        # create transition matrix of Z0 conditional on education
        Qz                = np.kron(z0_prob,np.ones([len(EDUC),len(PSY),len(S),len(FE_pos)]))

        # create transition matrix of fixed states: Qfixed
        Qfixed              = scipy.sparse.identity(len(S)*len(FE_pos)*len(PSY)*len(EDUC))

        # create aggregate transition: Q
        Q   = dprod(Qz,Qfixed) 
        Q_ergo_out.mat2 = Q
        Q_ergo_out.created = 'Y'
    
   if Q_ergo_in.created == 'Y':
        Q   = Q_ergo_in.mat2
        Q_ergo_out = Q_ergo_in

    
    # new distribution: 
   mup_vec = Q.T * mu0int_vec
   mup     = np.reshape(mup_vec,[NZ0,len(EDUC),len(PSY),len(S),len(FE_pos)])
    

   return mup, Q_ergo_out 