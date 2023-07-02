# -*- coding: utf-8 -*-
"""
Created on Sun Jul 24 09:59:08 2022

@author: Giorgia
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 11:58:53 2022

@author: Giorgia
"""


#approx_2d with arrays
from types import SimpleNamespace
import time
import itertools as it
import numpy as np
from scipy import optimize

#approx_2d
def approx_2d(X,Y,x):
    Nx = len(X)
    #Mx = round(Nx/2)
    j = 1
    y = np.empty(len(x))
    y[:] = np.nan
    for i in range(0,len(x)):
        #while (j < Nx-1 and x[i] > X[j]):
        while (j<Nx-1 and x[i]>X[j]):
            j = j+1

    #First order aproximation
        dy = (Y[j]-Y[j-1])/(X[j]-X[j-1])############### problem, in matlab index j is different, each element of a matrix has a single number
        dx = x[i]-X[j-1]
        y[i] = Y[j-1] + dy*dx

    #approximation outside bounds: linear extrapolation (OLS)
        Nextrap = 10
        Nextrap = min(Nextrap,Nx)
        if x[i] < X[0]:#CHAGED SIGN FOR ATTEMPT
            X_reg = np.array([np.ones(Nextrap), np.reshape(X[0:Nextrap],Nextrap)]).T #I transpose it to get same result as Matlab, but I X is a marix then it does not work well
            Y_reg = np.reshape(Y[0:Nextrap],Nextrap).T
            beta_nom = np.dot(X_reg.T,Y_reg)
            beta_denom = np.dot(X_reg.T, X_reg)
            beta = np.linalg.inv(beta_denom).dot(beta_nom)##### equivalent of blacklash operator to solve linear system in matlab
            new_array = np.array([1,x[i]])
            y[i] = np.dot(new_array,beta) #####I get reversed order from if: (instead of y = [2.71, 2.71, 5] I get y = [5,5,2.71])
        elif x[i] > X[Nx-1]:
            X_reg = np.array([np.ones(Nextrap), np.reshape(X[Nx-Nextrap:], Nextrap)]).T #I transpose it to get same result as Matlab
            Y_reg = np.reshape(Y[Nx - Nextrap:],Nextrap).T
            beta_nom = np.dot(X_reg.T,Y_reg)
            beta_denom = np.dot(X_reg.T, X_reg)
            beta = np.linalg.inv(beta_denom).dot(beta_nom)
            y[i] = np.dot(np.array([1,x[i]]),beta)########bad things happening PROBLEM WITH IF STATEMENT, BOTH RESULTS ARE THE SAME BECAUSE IT OPERATES ONLY ON x[i]>X[i]
    return y




