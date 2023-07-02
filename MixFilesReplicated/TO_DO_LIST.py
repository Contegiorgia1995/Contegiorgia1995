# -*- coding: utf-8 -*-
"""
Created on Sun Sep 18 10:12:36 2022

@author: Giorgia
"""

###########
#1.##fix interpolation
#############################UNDERSTAND HOW EXP_PROB WORKS IN prob_educ_col_grad as that determines whetehr we use spi.intep1d of interp2d and it messes up GEM_col >>>>>>>>>>>>>> how does len innop_prob change?????
#kinda got that but Cpp           = np.squeeze(Cpaux[:,:,ife]).T

#####IndexError: too many indices for array: array is 2-dimensional, but 3 were indexed


#2.#check GEG, GEM_with child; GEG_college if they work once :
    # for i in range(int(imin+1),int(imax-1)):
    #     if i == imin + 1:
    #         if len(innop_prob)>1:
    #             inno3,Sp3  = np.meshgrid(INNO_pos,Spgrid)
    #            #Vpaux          = splVp(Sp3,inno3)*innop_prob.T
    #             Vpaux          = np.dot(splVp(Sp3[:,0],inno3[0,:]).T,innop_prob)
    #         else:
    #             Vpaux          = splVp(Spgrid).T######################problem
    #is out of comment
#3. fix repmat at the end of prob_educ_co_grad and prob_educ_hs_grad
#4. adjust solve_model
#5. adjust text
    
#########################
#fixed C-V in Vaux (prob_work_trans and prob_work_fertility)
#fixed GEM---- maybe need to check GEGM with child and college too: ucp redenined before EGM and Spgrid redefined inside function "c" of GEGM or it becomes messy and SplVp etc also would be affected
##