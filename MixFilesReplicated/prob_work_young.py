import numpy as np
import scipy.interpolate as spi
from scipy.interpolate import griddata
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import LinearNDInterpolator
from scipy.interpolate import interpn
from scipy.interpolate import interp2d

def prob_work_young(par,options,Vp,Cp,j_pos):


    # Solve household problem when young
    r_sav     = par.r_sav
    r_debt    = par.r_debt
    beta      = par.beta
    gammac    = par.gammac
    w         = par.w
    Lambda    = par.Lambda
    
    S         = par.grids[0][j_pos]
    EDUC      = par.educ
    
    AGE_PROF   = par.inc.age_prof
    INNO_pos   = par.inc.inno_pos
    INNO_posT = np.reshape(INNO_pos,[1,5])

    FE_pos     = par.inc.fe_pos
    
    Sp = np.zeros([len(INNO_pos),len(EDUC),len(S),len(FE_pos)])
    C  = np.zeros([len(INNO_pos),len(EDUC),len(S),len(FE_pos)])
    V  = np.zeros([len(INNO_pos),len(EDUC),len(S),len(FE_pos)])
    boundgrid = np.zeros([len(INNO_pos),len(FE_pos),len(EDUC)])

    
    for educ in range (0,len(EDUC)):
        
        S          = par.grids[educ][j_pos]
        Spgrid  = par.grids[educ][j_pos+1]
        
        
        FE        = par.inc.fe[educ]
        INNO      = par.inc.z_val[educ][j_pos,:]
        
        inno3 = np.repeat(INNO_posT,len(S),axis=0)
        INNOp_prob = par.inc.z_prob[educ][:,j_pos+1,:].T # Note I need j_pos + 1 here, because it is in continuation value!
        
        r         = (r_sav * (Spgrid>=0) + r_debt * (Spgrid<0)).T
        rs        = (r_sav * (S>=0) + r_debt * (S<0))
        
        INNO2,S2  = np.meshgrid(INNO_pos,Spgrid)
        
        for ife in range (0,len(FE_pos)):
            splVp = GlobalSpline2D(Spgrid, INNO_pos,Vp[:,educ,:,ife] )
            for inno  in range (0,len(INNO_pos)):
            
                innop_prob = INNOp_prob[inno,:]
                
                
                h   = exp(AGE_PROF[educ][0][j_pos]) * exp(FE[ife])* exp(INNO[inno])
                
                Cpp = Cp[:,educ,:,ife].T
                
                if min(np.ravel(Cpp.any())<0 and abs(min(np.ravel(Cpp.any()))))>1e-6:
                    Cpp = max(Cpp,0)
                ucp = beta*np.reshape(((1+r)),[80,1])*np.dot(Cpp**(-gammac), np.reshape(innop_prob,[size(innop_prob),1]))
                
                dispinc = w*h*(1-Lambda)
                
                feasible = (dispinc + (1+rs)*S - Spgrid[0]>0)
                
                not_feasible = (dispinc + (1+rs)*S - Spgrid[0]<=0)
        
                
                C[inno,educ,feasible,ife],Sp[inno,educ,feasible,ife],boundgrid[inno,ife,educ] = GEGM(par,ucp.T,dispinc,Spgrid,S[feasible],splVp,innop_prob)

                #             EGM(par,ucp,dispinc,Spgrid,S);
                
                Sp3 = np.repeat(np.reshape(Sp[inno,educ,feasible,ife],[80,1]),len(INNO_pos),axis = 1)
                vp = np.dot(splVp(Sp3[:,0],inno3[0,:]).T,innop_prob)
                
                V[inno,educ,feasible,ife] = C[inno,educ,feasible,ife]**(1-gammac)/(1-gammac) + beta*vp[feasible]

                
                C[inno,educ,not_feasible,ife] = 0;
                V[inno,educ,not_feasible,ife] = -(10**(5/gammac))**(1-gammac)/(1-gammac)
                Sp[inno,educ,not_feasible,ife]= Spgrid[0]

    
    errors = sum(boundgrid[:])/len(C[:])
    if options.timer_on == 'Y':
            if sum(boundgrid[:]) >= len(EDUC)*len(FE_pos)*len(INNO_pos):
                print('j : {par.age(j_pos)}, Share of errors (increase grid) ={errors} \n')

    return V,C,Sp
