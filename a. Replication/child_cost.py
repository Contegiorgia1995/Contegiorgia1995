# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 17:27:21 2022

@author: Giorgia
"""

def ChildCost(par,labor_inc,n,T,options): ##ok<INUSD>
    ## Child cost function - Source: Angrist and Evans 1998. Assume 3rd child cost equal to any child + linear cost (see table 10)
    # 1/3 is loss of hours when raising kids oneself, 1/8 when hiring nannys. Difference must be covered by nannies.
    # C  = @(h,n,T) T * n^0.7 *(par.w*h*(1-par.lambda) * 1/8  + (1/3-1/8) * pn ) ...
    # 		  + (1-T) * n^0.7 *  (par.w*h*(1-par.lambda) * 1/3);
    # Source: Angrist and Evans 1998. Assume 3rd child cost equal to any child + linear cost (see table 10)
    
    n                   = n*2/par.fam_size
    h                   = labor_inc/(par.w*(1-par.Lambda))
    
    
    #lost_wnanny         = par.Child_C_Tot*par.Nanny_oppC
    l_ret               = par.Child_Ccurv
    child_cost_inc_curv = par.child_cost_inc_curv
    # lost_total          = par.Child_C_Tot * (1-par.lambda)**(child_cost_inc_curv-1)
    lost_total          = par.Child_C_Tot
    C                   = n**l_ret * par.w*(1-par.Lambda) *  h**child_cost_inc_curv * lost_total
    
    return C
# lost_total          = par.Child_C_Tot
    # C                   = n**l_ret * labor_inc**child_cost_inc_curv * lost_total
    
    # l_ret               = 0.4757
    # l_ret               = 0.7;
    # if options.ChildCost == 'OppC_Const':
    #     
    ##         C           = T *       (n*par.fam_size)**l_ret * (labor_inc* lost_wnanny   + (lost_total-lost_wnanny) * par.pn ) + (1-T) * (n*par.fam_size)^l_ret * (labor_inc* lost_total);
    ##             If estimating lost_total, lost_wnanny, l_ret: use below
    #           C           = T *       n**l_ret * (labor_inc* lost_wnanny   + (lost_total-lost_wnanny) * par.pn ) + (1-T) * n**l_ret * (labor_inc* lost_total);
    # elif options.ChildCost == 'Constant':
    ##        C           = (n*par.fam_size)**l_ret * lost_total * par.pn
    #           C           = n**l_ret * lost_total * par.pn
    # elif options.ChildCost == 'OppC':
    ##        C           = (n*par.fam_size)**l_ret * (labor_inc* lost_total);
   ##        C           = min(n**l_ret * (labor_inc* lost_total),labor_inc)
    #         C           = n^l_ret .* labor_inc.^child_cost_inc_curv .* lost_total;

    

    
    
    ### Economies of scale:
    ## Source: Folbre 2008, table 6.4, also used in Cordoba, Liu, Ripoll 2015
    # Y                   = [150 2*100 3*85]';
    # l_ret               = 0.4757;
    #% Y                   = [55.2 82.6 113]';
    # N                   = (1:3)';
    #lY                  = log(Y);
    # lN                  = log(N);
    # X                   = [ones(3,1) lN];
    # b                   = (X.T*X)\X.T*lY;
    # l_ret               = b(2)
    # l_ret               = 0.6445;
