# -*- coding: utf-8 -*-
"""
Created on Mon Oct 24 10:28:34 2022

@author: Giorgia
"""
#https://github.com/QuantEcon/BasisMatrices.jl







from ast import literal_eval
def solve_ergodic_distribution2(par2,pol_dense,options):
# Find ergodic distribution across cohorts
    if options.timer_on == 'Y':
        print(f'\n Find ergodic distribution across cohorts \n \n')


# Compute demographic transition
# Find ergodic distribution across cohorts
###################################### OUTPUT:#####################################
#% mu{j}(a) is the pdf of state a at age j (i.e. at the begining of age j)
#################################################################################

## distribute policy functions
# from varname import varname
# def function():
#     return varname()

# func = function()
    pol_mat = copy.deepcopy(pol_dense.__dict__)##### this will be useful!!! converse SImpleNameSPce to dictionary
    pol_mat = pol_mat.keys()
    for i in range (0,len(pol_mat)):
        eval([literal_eval(pol_mat[i]) '= pol_dense.(cell2mat(pol_mat(i)));']);
    #######################################need to finish but eval is ugly
    #del pol_dense
    
    
    # Empty cell for measures at each age
    mu        =  {}
    
    for i in range (0, par2.Jd_pos-1):
        mu[i] = [[]]
    
    Q_ergo =  {}
    for i in range (0, par2.Jd_pos-1):
        Q_ergo[i] = [[]]

    for jj in range (0,par2.Jd_pos-1):#structure array 
        Q_ergo[jj].created = 'N'

    
    
    ## mu{1}: guess
    
    # Load initial guess for tau and Vc0

    os.chdir('/Users/Giorgia/Dropbox/ReplicationProject/DairuchPaper/DairuchFiles/run_generic')        
    filename = str('dist_ss')
    mu = str(filename + '.mat')
    file_path = r'/Users/Giorgia/Dropbox/ReplicationProject/DairuchPaper/DairuchFiles/run_generic/'+mu
    exist = os.path.isfile(file_path)

    if exist:
        guess_ = loadmat(mu)
        j_pos = par2.Je1_pos-2
        guess.mu = copy.deepcopy(guess_["mu"])
        for j in range (0,3):
            guess.mu[j][0] = guess.mu[j][0].T
            for i in range(0,10):
                guess.mu[j][0][i] = guess.mu[j][0][i].T
        j = 3######################################################careful here that dimensions 3 and 10 are inverted here (in matlab (10,3))
        guess.mu[j][0] = guess.mu[j][0].T
        for i in range (0,3):
            guess.mu[j][0][i] = np.reshape(guess.mu[j][0][i],[10,20,20])
            for k in range (0,10):
                guess.mu[j][0][i][k] = guess.mu[j][0][i][k].T
        j = 4# (3, 10, 100, 100) vs (100   100    10     3)
        guess.mu[j][0] = guess.mu[j][0].T
        for k in range (0,10):
            guess.mu[j][0][i][k] = guess.mu[j][0][i][k].T
                
        for j in range (5,8):# (3, 100, 100)  vs 100   100     3
            guess.mu[j][0] = guess.mu[j][0].T
            for i in range(0,3):
                guess.mu[j][0][i] = guess.mu[j][0][i].T
        for j in range (8,10):# (4, 3, 100, 100)  vs 100   100     3     4
            guess.mu[j][0] = guess.mu[j][0].T
            # for i in range(0,3):
            #     guess.mu[j][0][i] = guess.mu[j][0][i].T
        j = 10 # (3, 20, 4, 3, 100, 100) vs 100   100     3     4    20     3 ######################gotta fix this
        guess.mu[j][0] = guess.mu[j][0].T
        
        muc   = guess.mu[j_pos][0]
    
    else:
        muc = np.zeros([2,2])

    j_pos     = par2.Je1_pos-1
    S_grid    = par2.grids[0][j_pos]
    FE_pos    = par2.inc.fe_pos
    PSY       = par2.psy_val_hs
    
    # If guess is of different size of current grid start new guess
    muc_aux = np.zeros([len(PSY),len(S_grid),len(FE_pos)])
    if muc.ndim != muc_aux.ndim:
        # start with a guess for distribution over (s,educ)
        # Guess 1: 'symmetric': one obs for each grid point (s,educ)
        muc    = np.zeros([len(PSY),len(S_grid),len(FE_pos)])
        nn     = len(PSY)*len(S_grid)*len(FE_pos) #size of grid to get probability measure
        muc[:] = 1/nn
    elif np.sum(np.array(np.shape(muc)) - np.array(np.shape(muc_aux))) != 0:
        # start with a guess for distribution over (s,educ)
        # Guess 1: 'symmetric': one obs for each grid point (s,educ)
        muc    = np.zeros([len(PSY),len(S_grid),len(FE_pos)])
        nn     = len(PSY)*len(S_grid)*len(FE_pos)#size of grid to get probability measure
        muc[:] = 1/nn

    
    pop = np.zeros([options.maxiterD,1])
    cohort = 0
    dif = 1
    
    if options.timer_on == 'Y':
        print(f'cohort: {cohort}, dif:{dif}, pop:{np.sum(muc[:])}, time ={toc}')

    while ((dif>options.tolD ) and cohort <=options.maxiterD):
        cohort = cohort + 1
        j_pos      = par2.Je1_pos-2
        mu[j_pos] = muc # update initial distribution for cohort
        
        ## Age j=Je1: new states are education and innovation
        j_pos = par2.Je1_pos-1
        [mu[j_pos],Q_ergo[j_pos-1]] = trans_draw(par2,mu[j_pos-1],tau0,Q_ergo[j_pos-1])#### here the shit starts
        
        if options.timer_on == 'Y':
            tot = sum(mu{j_pos}(:))
            print(f'age:{par2.age(j_pos)}, pop:{tot}, time:{toc} \n')

        
        ## Age j=Je1+1:Je2+3
        for j_pos in range(par2.Je1_pos,par2.Je2_pos + 2):
            mu[j_pos],Q_ergo[j_pos-1] = trans_educ(par2,mu[j_pos-1],Se[0][j_pos - par2.Je1_pos],Se[1][j_pos - par2.Je1_pos],Se[2][j_pos - par2.Je1_pos],j_pos,Q_ergo[j_pos-1])
            # Only need to pass policy functions given education, iid shock does
            # not matter once we condition on education choice
            if options.timer_on == 'Y':
                tot = np.sum(mu[j_pos][:])
                print(f'age:{par2.age(j_pos)}, pop:{toc}, time: {toc} \n')

    #     mu{j_pos}                = squeeze(sum(mu{j_pos},3)); % Remove Psychic Cost
     
        
       ## Age j=Je2 + 2: Jc
        for j_pos in range(par2.Je2_pos+2,par2.Jc_pos):
            mu[j_pos],Q_ergo[j_pos-1] = trans_work(par2,mu[j_pos-1],S[j_pos-1],j_pos,Q_ergo[j_pos-1])
            
            if options.timer_on == 'Y':
                tot = np.sum(mu[j_pos][:])
                print(f'age:{par2.age(j_pos)}, pop:{toc}, time: {toc} \n')

        
        ## Age Jc
        j_pos = par2.Jc_pos
        fert_pol = Np
        mu[j_pos],Q_ergo[j_pos-1] = trans_fertility(par2,mu[j_pos-1],S[j_pos-1],fert_pol,j_pos,Q_ergo[j_pos-1],options)
        switch options.timer_on
            case {'Y'}
                tot = sum(mu{j_pos}(:));
                fprintf('age:%i, pop:%1.2f, time: %3.1f sec \n',par2.age(j_pos),tot,toc)
        
        ## Age Jc+1
        for j_pos in range(par2.Jc_pos+1,par2.Jc_pos+par2.Je1_pos-2):
            mu[j_pos],Q_ergo[j_pos-1] = trans_work_with_child(par2,mu[j_pos-1],S[j_pos-1],j_pos,Q_ergo[j_pos-1])
            if options.timer_on == 'Y':
                tot = np.sum(mu[j_pos][:])
                print(f'age:{par2.age(j_pos)}, pop:{toc}, time: {toc} \n')
        
        ## Age Jc+2: phi transfer + new cohort
        j_pos       = par2.Jc_pos+par2.Je1_pos-2
        #     j_pos_child = par2.Je1_pos;
        phi_pol = PHIp
        mu[j_pos],Q_ergo[j_pos-1] = trans_work_trans(par2,mu[j_pos-1],S[j_pos-1],j_pos,Q_ergo[j_pos-1])
        if options.timer_on == 'Y':
            tot = np.sum(mu[j_pos][:])
            print(f'age:{par2.age(j_pos)}, pop:{toc}, time: {toc} \n')
        
        pop_final  = pop[cohort]
        j_0_pos = par2.Je1_pos-2
        dif = np.max(abs(mu[j_0_pos][:]-muc[:]))
        # Describe muc
        mu[0]      = muc
        mu[1]      = muc
        
        if options.timer_on == 'Y':
            tot = np.sum(mu[j_pos][:])
            print(f'cohort:{cohort}, pop:{pop(cohort)}, time: {toc} \n')

    return mu,mu_ige0,muc,pop_final,Q_ergo