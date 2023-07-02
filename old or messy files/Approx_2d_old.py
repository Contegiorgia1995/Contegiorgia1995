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
        if x[i] < X[i]:#CHAGED SIGN FOR ATTEMPT
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

X = [6,5,1]
Y = [1,2,3]
x = [2,2,2]

approx_2d(np.append(S_bc,a_endo),np.append(c_bc,c),S)

################## For matrixes__ still something to fix, but should only need the one above

def approx_2d_mat(X,Y,x):
    #X = np.array(X)
    #Y = np.array(Y)
    #x = np.array(x)
        
    
    Nx = len(X)
    #Mx = round(Nx/2)
    j = 1
    y = np.empty(len(x))
    ind_Y = np.arange(np.size(Y))
    ind_X = np.arange(np.size(X))
    ind_x = np.arange(np.size(x))
    
    Y_res = Y[np.unravel_index(ind_Y, Y.shape, 'F')]
    X_res = X[np.unravel_index(ind_X, X.shape, 'F')]
    x_res = x[np.unravel_index(ind_x, x.shape, 'F')]
    for i in range(0,len(x)):
        #while (j < Nx-1 and x[i] > X[j]):
        while (j<Nx-1 and x_res[i]>X_res[j]):
            j = j+1
       
    #First order aproximation: For dy to work X and Y must have the same lenght!!!!!!!
        dy = np.divide((Y_res[j]-Y_res[j-1]),(X_res[j]-X_res[j-1]))############### problem, in matlab index j is different, each element of a matrix has a single number
        dx = x_res[i]-X_res[j-1]
        y[i] = Y_res[j-1] + dy*dx

    #approximation outside bounds: linear extrapolation (OLS)
        Nextrap = 10
        Nextrap = min(Nextrap,Nx)
        if x_res[i] < X_res[i]:#CHAGED SIGN FOR ATTEMPT
            X_reg = np.array([np.ones(Nextrap), np.reshape(X_res[0:Nextrap],(Nextrap))]).T #I transpose it to get same result as Matlab, but I X is a marix then it does not work well
            Y_reg = np.reshape(Y_res[0:Nextrap],Nextrap).T
            beta_nom = np.dot(X_reg.T,Y_reg)
            beta_denom = np.dot(X_reg.T, X_reg)
            beta = np.linalg.inv(beta_denom).dot(beta_nom)
            y[i] = np.dot(np.array([1,x_res[i]]),beta) #####I get reversed order from if: (instead of y = [2.71, 2.71, 5] I get y = [5,5,2.71])
        elif x_res[i] > X_res[Nx-1]:
            # X_reg = np.array([np.ones(Nextrap), np.reshape(X[Nx-Nextrap:], (Nextrap))]).T #Here I have a huge problem with shapes, also in Matlab it struggles to understand how can convert (NX-Nextrap = 0 to the end) to an array of 3 entires only
            # Y_reg = np.reshape(Y[Nx - Nextrap:],Nextrap).T
            X_reg = np.array([np.ones(Nextrap), np.reshape(X[0,Nx-Nextrap:], (Nextrap))]).T #Here I have a huge problem with shapes, also in Matlab it struggles to understand how can convert (NX-Nextrap = 0 to the end) to an array of 3 entires only
            Y_reg = np.reshape(Y[0,Nx - Nextrap:],Nextrap).T ##just an idea of picking the first row only of X and Y and rshape it so that it fits the dimension of Nextrap
            beta_nom = np.dot(X_reg.T,Y_reg)
            beta_denom = np.dot(X_reg.T, X_reg)
            beta = np.linalg.inv(beta_denom).dot(beta_nom)
            y[i] = np.dot(np.array([1,x_res[i]]),beta)
    return y#
    
   
    

 
X = [6,5,1]
Y = [1,2,3]
x = [2,2,2]
# X = [7,1,5,4]
# Y = [1,2,3,4]
# x = [8,2,1,2]
approx_2d(X,Y,x)

X = np.array([[6,5,1],[4,3,2]])

Y = np.array([[1,2,3],[3,2,1]])

x= np.array([[8,8,2],[2,2,2]])