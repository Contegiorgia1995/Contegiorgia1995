# -*- coding: utf-8 -*-
"""
Created on Sun Oct  9 18:05:49 2022

@author: Giorgia
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 17:45:09 2022

@author: Giorgia
"""
#https://github.com/pig2015/mathpy/blob/master/polation/globalspline.py
        
    r_sav     = par.r_sav
    r_debt    = par.r_debt
    beta      = par.beta
    gammac    = par.gammac
    w         = par.w
    Lambda    = par.Lambda
    lambdan   = par.lambdan
    gamman    = par.gamman
    
    S         = par.grids[0][j_pos]
    EDUC      = par.educ
    
    AGE_PROF   = par.inc.age_prof
    INNO_pos   = par.inc.inno_pos
    INNO_posT = np.reshape(INNO_pos,[1,5])
    FE_pos     = par.inc.fe_pos
    
    N         = np.arange(0,4,1)
    
    Sp = np.zeros([len(INNO_pos),len(EDUC),len(N),len(S),len(FE_pos)])
    C  = np.zeros([len(INNO_pos),len(EDUC),len(N),len(S),len(FE_pos)])
    V  = np.zeros([len(INNO_pos),len(EDUC),len(N),len(S),len(FE_pos)])
    Ck = np.zeros([len(INNO_pos),len(EDUC),len(N),len(S),len(FE_pos)])
    Tp = np.zeros([len(INNO_pos),len(EDUC),len(N),len(S),len(FE_pos)])
    boundgrid = np.zeros([len(INNO_pos),len(N),len(FE_pos),len(EDUC)])
    
    ## Without children and phi = 1
    if options.Fertility == 'Endo': # Only need to solve case without children in benchmark case where we allow them to choose endogenous fertility
            i_n   = 0
            
            for educ in range (0,len(EDUC)):
                #parameters(par_est, options)
                S          = par.grids[educ][j_pos]
                Spgrid  = par.grids[educ][j_pos+1]
                
                FE        = par.inc.fe[educ]
                INNO      = par.inc.z_val[educ][j_pos,:]
                
                
                r         = (r_sav * (Spgrid>=0) + r_debt * (Spgrid<0)).T
                rs        = (r_sav * (S>=0) + r_debt * (S<0))
                
                inno3 = np.repeat(INNO_posT,len(S),axis=0)
                INNOp_prob = par.inc.z_prob[educ][:,j_pos+1,:].T # Note I need j_pos + 1 here!

                INNO2,S2  = np.meshgrid(INNO_pos,Spgrid)
                
                
                for ife in range (0,len(FE_pos)):
                    splVp = interp2d(Spgrid, INNO_pos,Vp[:,educ,i_n,:,ife])
                    #splVp1 = extrap1d(splVp)
                    #splVp = LinearNDInterpolator(np.reshape(list(zip(S2,INNO2)),[80,5,2]),Vp[:,educ,i_n,:,ife].T )
        
                    for inno  in range (0,len(INNO_pos)):
                        
                   
                        innop_prob = INNOp_prob[inno,:]
                        
                        
                        h   = exp(AGE_PROF[educ][0][j_pos]) * exp(FE[ife])* exp(INNO[inno])
                        
                        Cpp = Cp[:,educ,i_n,:,ife].T
                        
                        if min(np.ravel(Cpp.any())<0 and abs(min(np.ravel(Cpp.any()))))>1e-6:
                            Cpp = max(Cpp,0)
                        ucp = beta*np.reshape(((1+r)),[80,])*np.dot(Cpp**(-gammac), np.reshape(innop_prob,[5,]))

                        dispinc = w*h*(1-Lambda)
                        feasible = ((dispinc + (1+rs)*S - Spgrid[0])>0)
                        
                        not_feasible = (dispinc + (1+rs)*S - Spgrid[0]<=0)
                        
                        C[inno,educ,i_n,feasible,ife],Sp[inno,educ,i_n,feasible,ife],boundgrid[inno,i_n,ife,educ] = EGM(par,ucp.T,dispinc,Spgrid,S[feasible])
                        Ck[inno,educ,i_n,feasible,ife] = np.zeros(len(C[inno,educ,i_n,feasible,ife]))
                        
                        
                        Sp3 = np.repeat(np.reshape(Sp[inno,educ,i_n,feasible,ife],[80,1]),len(INNO_pos),axis = 1)
                        
                        vp = np.dot(splVp(Sp3[:,0],inno3[0,:]).T,innop_prob)
                        
                        V[inno,educ,i_n,feasible,ife] = C[inno,educ,i_n,feasible,ife]**(1-gammac)/(1-gammac) + beta*vp[feasible]
                        
                        C[inno,educ,i_n,not_feasible,ife] = 0;
                        V[inno,educ,i_n,not_feasible,ife] = -(10**(5/gammac))**(1-gammac)/(1-gammac)
                        Sp[inno,educ,i_n,not_feasible,ife]= Spgrid[0]
                        
                        
                        Tp[inno,educ,i_n,:,ife] = np.zeros(len(C[inno,educ,i_n,:,ife]))

    
    ## Case with Children
    if options.Fertility == 'Endo': #Only need to solve case without children in benchmark case where we allow them to choose endogenous fertility
            in_1    = 1
    elif options.Fertility == 'Exo':
            in_1    = 0

    
    for educ in range (0,len(EDUC)):
        
        INNOp_prob  = par.inc.z_prob[educ][:,j_pos+1,:].T  # Note I need j_pos + 1 here!
        
        S          = par.grids[educ][j_pos]
        Spgrid     = par.grids[educ][j_pos+1]
    
        r         = (r_sav * (Spgrid>=0) + r_debt * (Spgrid<0)).T
        rs        = (r_sav * (S>=0) + r_debt * (S<0))
        
        inno3 = np.repeat(INNO_posT,len(S),axis=0)
        
        FE        = par.inc.fe[educ]
        INNO      = par.inc.z_val[educ][j_pos,:]
        
        INNO2,S2  = np.meshgrid(INNO_pos,Spgrid)
        
        for i_n in range(in_1,len(N)):
            n             = par.N[i_n]
            for ife in range (0,len(FE_pos)):
                splVp =  GlobalSpline2D(Spgrid, INNO_pos,Vp[:,educ,i_n,:,ife])
                for inno in range (0,len(INNO_pos)):
                    innop_prob      = INNOp_prob[inno,:]
                    h   = exp(AGE_PROF[educ][0][j_pos]) * exp(FE[ife])* exp(INNO[inno])
                    
                    Cpp = np.squeeze(Cp[:,educ,i_n,:,ife].T)
                    if min(np.ravel(Cpp.any())<0 and abs(min(np.ravel(Cpp.any()))))>1e-6:
                        Cpp = max(Cpp,0)
                    ucp = beta*np.reshape(((1+r)),[80,1])*np.dot(Cpp**(-gammac), np.reshape(innop_prob,[5,1]))
                    
                    labor_inc = w*h*(1-Lambda)
                    
                    n_final = n* (par.fam_size/2)
                    
                    ChildCost_0 = ChildCost(par,labor_inc,n_final,0,options)
                    ChildCost_1 = ChildCost(par,labor_inc,n_final,1,options)
                    ChildCost_opt = min(ChildCost_0,ChildCost_1)
                    Tp[inno,educ,i_n,:,ife]   = (ChildCost_1>ChildCost_0)
                    
                    dispinc = labor_inc- ChildCost_opt
                    feasible = (dispinc + (1+rs)*S - Spgrid[0]>0)
                    not_feasible = (dispinc + (1+rs)*S - Spgrid[0]<=0)
                    
                    
                    if options.Ck == 'Yes':
                        
                        C[inno,educ,i_n,feasible,ife],Sp[inno,educ,i_n,feasible,ife],boundgrid[inno,i_n,ife,educ] = GEGM_withchild(par,ucp.T,dispinc,Spgrid,S[feasible],splVp,innop_prob,n_final) 
                       # parameters(par_est, options)
                        Spgrid     = par.grids[educ][j_pos+1]#                        EGM_withchild(par,ucp,dispinc,Spgrid,S,n_final); %%%% To do: Generalized EGM with child
                        
                        f_n     = (lambdan/n_final**(1-gamman))**(1/gammac)
                        Ck[inno,educ,i_n,feasible,ife]     = f_n * C[inno,educ,i_n,feasible,ife]
                        
                        
                        Sp3 = np.repeat(np.reshape(Sp[inno,educ,i_n,:,ife],[len(Spgrid),1]),len(INNO_pos),axis = 1)   #######wrong!! but Sp[] ok                 
                        vp = np.dot(splVp(Sp3[:,0],inno3[0,:]).T,innop_prob)
                        
                        
                        V[inno,educ,i_n,feasible,ife] = (C[inno,educ,i_n,feasible,ife]**(1-gammac))/(1-gammac)*(C[inno,educ,i_n,feasible,ife]>=0) + lambdan* n_final**(gamman) * (Ck[inno,educ,i_n,feasible,ife]**(1-gammac))/(1-gammac)*(C[inno,educ,i_n,feasible,ife]>=0) -((10**(5/gammac))**(1-gammac))/(1-gammac)*(C[inno,educ,i_n,feasible,ife]<0) + beta*vp[feasible]
                        # Fix bad extrapolation:
                        Sp[inno,educ,i_n,feasible,ife]  = Sp[inno,educ,i_n,feasible,ife]*(C[inno,educ,i_n,feasible,ife] >=0) + Spgrid[0]*(C[inno,educ,i_n,feasible,ife]<0)
                        C[inno,educ,i_n,feasible,ife]  = C[inno,educ,i_n,feasible,ife]*(C[inno,educ,i_n,feasible,ife] >= 0) + 0*(C[inno,educ,i_n,feasible,ife] <0)
                        Ck[inno,educ,i_n,feasible,ife] = Ck[inno,educ,i_n,feasible,ife]*(Ck[inno,educ,i_n,feasible,ife]>= 0) + 0*(Ck[inno,educ,i_n,feasible,ife]<0)

                        C[inno,educ,i_n,not_feasible,ife] = 0
                        Ck[inno,educ,i_n,not_feasible,ife] = 0
                        V[inno,educ,i_n,not_feasible,ife] = -(10**(5/gammac))**(1-gammac)/(1-gammac)
                        Sp[inno,educ,i_n,not_feasible,ife]= Spgrid[0]


                        
                    elif options.Ck == 'No':
                        C[inno,educ,i_n,feasible,ife],Sp[inno,educ,i_n,feasible,ife],boundgrid[inno,i_n,ife,educ]= GEGM(par,ucp.T,dispinc,Spgrid,S[feasible],splVp,innop_prob)
                        
                        Ck[inno,educ,i_n,feasible,ife] = np.zeros(len(C[inno,educ,i_n,feasible,ife]))
                        
                        Sp3         = np.repeat(np.reshape(Sp[inno,educ,i_n,:,ife],[len(feasible),1]),len(INNO_pos),axis = 1)    
                        vp          = np.dot(GlobalSpline2D(Sp3[:,0],inno3[0,:]),innop_prob)### wrong interpolation
                        

                        V[inno,educ,i_n,feasible,ife] = C[inno,educ,i_n,feasible,ife]**(1-gammac)/(1-gammac)*(C[inno,educ,i_n,feasible,ife]>=0) -(10**(5/gammac))**(1-gammac)/(1-gammac)*(C[inno,educ,i_n,feasible,ife]<0) + beta*vp[feasible]

                        
                        C[inno,educ,i_n,not_feasible,ife] = 0
                        Ck[inno,educ,i_n,not_feasible,ife] = 0
                        V[inno,educ,i_n,not_feasible,ife] = -(10**(5/gammac))**(1-gammac)/(1-gammac)
                        Sp[inno,educ,i_n,not_feasible,ife]= Spgrid[0]
                        
def extrap2d(interpolator):
    xs = interpolator.x
    ys = interpolator.y
    zs = interpolator.z
    zs = np.reshape(zs, (-1, len(xs)))
def pointwise(x, y):
    if x < xs[0] or y < ys[0]:
        x1_index = np.argmin(np.abs(xs - x))
        x2_index = x1_index + 1
        y1_index = np.argmin(np.abs(ys - y))
        y2_index = y1_index + 1
        x1 = xs[x1_index]
        x2 = xs[x2_index]
        y1 = ys[y1_index]
        y2 = ys[y2_index]
        z11 = zs[x1_index, y1_index]
        z12 = zs[x1_index, y2_index]
        z21 = zs[x2_index, y1_index]
        z22 = zs[x2_index, y2_index]

        return (z11 * (x2 - x) * (y2 - y) +
        z21 * (x - x1) * (y2 - y) +
        z12 * (x2 - x) * (y - y1) +
        z22 * (x - x1) * (y - y1)
       ) / ((x2 - x1) * (y2 - y1) + 0.0)


    elif x > xs[-1] or y > ys[-1]:
        x1_index = np.argmin(np.abs(xs - x))
        x2_index = x1_index - 1
        y1_index = np.argmin(np.abs(ys - y))
        y2_index = y1_index - 1
        x1 = xs[x1_index]
        x2 = xs[x2_index]
        y1 = ys[y1_index]
        y2 = ys[y2_index]
        z11 = zs[x1_index, y1_index]
        z12 = zs[x1_index, y2_index]
        z21 = zs[x2_index, y1_index]
        z22 = zs[x2_index, y2_index]#

        return (z11 * (x2 - x) * (y2 - y) +
        z21 * (x - x1) * (y2 - y) +
        z12 * (x2 - x) * (y - y1) +
        z22 * (x - x1) * (y - y1)
        ) / ((x2 - x1) * (y2 - y1) + 0.0)
    else:
        return interpolator(x, y)
    
    
from scipy import interpolate
import numpy as np


class GlobalSpline2D(interpolate.interp2d):
    def __init__(self, x , y , z , kind='linear'):
        if kind == 'linear':
            if len(x) < 2 or len(y) < 2:
                raise self.get_size_error(2, kind)
        elif kind == 'cubic':
            if len(x) < 4 or len(y) < 4:
                raise self.get_size_error(4, kind)
        elif kind == 'quintic':
            if len(x) < 6 or len(y) < 6:
                raise self.get_size_error(6, kind)
        else:
            raise ValueError('unidentifiable kind of spline')

        super().__init__(x, y, z, kind=kind)
        self.extrap_fd_based_xs = self._linspace_10(self.x_min, self.x_max, -4)
        self.extrap_bd_based_xs = self._linspace_10(self.x_min, self.x_max, 4)
        self.extrap_fd_based_ys = self._linspace_10(self.y_min, self.y_max, -4)
        self.extrap_bd_based_ys = self._linspace_10(self.y_min, self.y_max, 4)

    @staticmethod
    def get_size_error(size, spline_kind):
        return ValueError('length of x and y must be larger or at least equal '
                          'to {} when applying {} spline, assign arrays with '
                          'length no less than '
                          '{}'.format(size, spline_kind, size))

    @staticmethod
    def _extrap1d(xs, ys, tar_x):
        if isinstance(xs, np.ndarray):
            xs = np.ndarray.flatten(xs)
        if isinstance(ys, np.ndarray):
            ys = np.ndarray.flatten(ys)
        assert len(xs) >= 4
        assert len(xs) == len(ys)
        f = interpolate.InterpolatedUnivariateSpline(xs, ys)
        return f(tar_x)

    @staticmethod
    def _linspace_10(p1, p2, cut=None):
        ls = list(np.linspace(p1, p2, 10))
        if cut is None:
            return ls
        assert cut <= 10
        return ls[-cut:] if cut < 0 else ls[:cut]

    def _get_extrap_based_points(self, axis, extrap_p):
        if axis == 'x':
            return (self.extrap_fd_based_xs if extrap_p > self.x_max else
                    self.extrap_bd_based_xs if extrap_p < self.x_min else [])
        elif axis == 'y':
            return (self.extrap_fd_based_ys if extrap_p > self.y_max else
                    self.extrap_bd_based_ys if extrap_p < self.y_min else [])
        assert False, 'axis unknown'
        
    def __call__(self, x_, y_, **kwargs):
        xs = np.atleast_1d(x_)
        ys = np.atleast_1d(y_)

        if xs.ndim != 1 or ys.ndim != 1:
            raise ValueError("x and y should both be 1-D arrays")

        pz_yqueue = []
        for y in ys:
            extrap_based_ys = self._get_extrap_based_points('y', y)

            pz_xqueue = []
            for x in xs:
                extrap_based_xs = self._get_extrap_based_points('x', x)

                if not extrap_based_xs and not extrap_based_ys:
                    # inbounds
                    pz = super().__call__(x, y, **kwargs)[0]

                elif extrap_based_xs and extrap_based_ys:
                    # both x, y atr outbounds
                    # allocate based_z from x, based_ys
                    extrap_based_zs = self.__call__(x,
                                                    extrap_based_ys,
                                                    **kwargs)
                    # allocate z of x, y from based_ys, based_zs
                    pz = self._extrap1d(extrap_based_ys, extrap_based_zs, y)

                elif extrap_based_xs:
                    # only x outbounds
                    extrap_based_zs = super().__call__(extrap_based_xs,
                                                       y,
                                                       **kwargs)
                    pz = self._extrap1d(extrap_based_xs, extrap_based_zs, x)

                else:
                    # only y outbounds
                    extrap_based_zs = super().__call__(x,
                                                       extrap_based_ys,
                                                       **kwargs)
                    pz = self._extrap1d(extrap_based_ys, extrap_based_zs, y)

                pz_xqueue.append(pz)

            pz_yqueue.append(pz_xqueue)

        zss = pz_yqueue
        if len(zss) == 1:
            zss = zss[0]
        return np.array(zss)    