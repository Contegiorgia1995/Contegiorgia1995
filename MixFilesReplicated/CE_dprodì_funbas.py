# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 12:30:03 2022

@author: Giorgia
"""

import numpy as np
from numba import jit
from scipy.interpolate import fitpack as spl # For splines - temporary
import functools

@jit
def dprod(A,B):
    """ The direct sum of two matrices with the same number
	of rows is equivalent to computing row-wise tensor
	(Kronecker) products.
    
    Parameters
    ----------
    A : array_like
        In multidimensional approximation this is the 
        n - 1 dimensional basis matrix
    B: array_like
        In multidimensional approximation this is the 
        1 dimensional basis matrix
    Returns
    -------
    Res :  ndarray
        Matrix of shape (m,n), where m = no. rows in A and
        n = no. columns in A * no. columns in B
    Notes
    -----
    """

    nobsa , na = A.shape
    nobsb , nb = B.shape
    Res = np.empty((nobsa,nb*na))
    if nobsa != nobsb:
        return 'A and B must have same number of rows'
    for t in range(nobsa):
        for ia in range(na):
            for ib in range(nb):            
                Res[t,nb*(ia)+ib] = A[t,ia] * B[t, ib]
    return Res

@jit
def cheb_nodes(p, nodetype=0):

    """Chebyshev nodes - for 1 dimension.
    Returns Chebyshev nodes
    Parameters
    ----------
    p : array_like
        Parameter array containing:
         - the number of nodes
         - the lower bound of the approximation
         - the upper bound of the approximation
    nodetype : int
         - if 0 (default value) then use the usual nodes
         - if 1 then extend it to the endpoints
         - if 2 then us Lobatto nodes
    Returns
    -------
    x :  an array containing the Chebyshev nodes
    Notes
    -----
    """
    n , a , b = p[0] , p[1] , p[2]
    s = (b-a) / 2 
    m = (b+a) / 2  
    if (nodetype < 2):  # usual nodes
        k = np.pi*np.linspace(0.5,n-0.5,n)  
        x = m - np.cos(k[0:n]/n) * s  
        if (nodetype == 1):  # Extend nodes to endpoints
            aa = x[0]  
            bb = x[-1]  
            x = (bb*a - aa*b)/(bb-aa) + (b-a)/(bb-aa)*x
    else: # Lobatto nodes
        k = np.pi*np.linspace(0,n-1,n)
        x = m - np.cos(k[0:n]/(n-1)) * s
    return x
@jit
def wealth_knot(p,degree = 1,curv  = 0.15,lower_bound = 1e-7):
    """Knots for wealth distributions - for 1 dimension.
    Returns commonly used knots that can be used with splines
    Parameters
    ----------
    p : array_like
        Parameter array containing:
         - the number of nodes
         - the lower bound of the approximation
         - the upper bound of the approximation
    deg : int
            Degree of the splines.
            Default: create nodes instead of knots
    curve : float
            Weight on the lower end
            Default: 0.15
    lower_bound : float
            the correction for the lower bound
    Returns
    -------
    x :  an array containing the knots
    Notes
    -----
    Use these knots if you have boundary problems.
    """
    n , a , b = p[0] , p[1] , p[2]
    knots = np.linspace((a + lower_bound)**curv, (b+ lower_bound)**curv,n + 1 - degree) **(1.0/curv) -lower_bound
    return knots
@jit
def cheb_basex(p, x): 
    """Cheb basis matrix - for 1 dimension.
    Returns a matrix whose columns are the values of the (first kind) Chebyshev 
    polynomial of maximum degree n-1 evaluated at the points `x`. Degree 0 is 
    the constant 0.
    Parameters
    ----------
    p : array_like
        Parameter array containing:
         - the order of approximation - the highest degree polynomial is n-1
         - the lower bound of the approximation
         - the upper bound of the approximation
    x : array_like
        Points at which to evaluate the b-splines.
    Returns
    -------
    bas : ndarray
        Matrix of shape (m,n), where ``m = len(x)`` and
        ``n - 1 = order(polynomial)``
    Notes
    -----
    Orthogonal polynomial
    """
    n , a , b = p[0] , p[1] , p[2]
    z = (2/(b-a)) * (x-(a+b)/2)
    m = z.shape[0]
    bas = np.empty((m, n));
    bas[:, 0] = 1.0
    bas[:, 1] = z[:]
    z = 2 * z
    for i in range(m):
        for j in range(2,n):
            bas[i, j] = z[i] * bas[i, j-1] - bas[i, j-2]
    return bas
@jit
def mono_basex(p, x): 
    """Monomials basis matrix- for 1 dimension.
    Returns a matrix whose columns are the values of the monomials of maximum 
    order n - 1 evaluated at the points `x`. Degree 0 is the constant 0.
    Parameters
    ----------
    p : array_like
        Parameter array containing:
         - the order of approximation - the highest degree polynomial is n-1
         - the lower bound of the approximation
         - the upper bound of the approximation
    x : array_like
        Points at which to evaluate the b-splines.
    Returns
    -------
    bas : ndarray
        Matrix of shape (m,n), where ``m = len(x)`` and
        ``n - 1 = order(polynomial)``
    Notes
    -----
    Also known as the Vandermonde matrix
    """
    n , a , b = p[0] , p[1] , p[2]
    z = (2/(b-a)) * (x-(a+b)/2)
    m = z.shape[0]
    bas = np.empty((m, n));
    bas[:, 0] = 1.0
    for i in range(m):
        for j in range(1,n):
            bas[i, j] = z[i] * bas[i, j-1]
    return bas
@jit
def spli_basex(p, x ,knots=None , deg = 3 , order = 0 ):
    """Vandermonde type matrix for splines.
    Returns a matrix whose columns are the values of the b-splines of deg
    `deg` associated with the knot sequence `knots` evaluated at the points
    `x`.
    Parameters
    ----------
    p : array_like
        Parameter array containing:
         - the number of knots
         - the lower bound of the approximation
         - the upper bound of the approximation
    x : array_like
        Points at which to evaluate the b-splines.
    deg : int
        Degree of the splines.
        Default: cubic splines
    knots : array_like
        List of knots. The convention here is that the interior knots have
        been extended at both ends by ``deg + 1`` extra knots - see augbreaks.
        If not given the default is equidistant grid
    order : int
        Evaluate the derivative of the spline
    Returns
    -------
    vander : ndarray
        Vandermonde like matrix of shape (m,n), where ``m = len(x)`` and
        ``n = len(augbreaks) - deg - 1``
    Notes
    -----
    The knots exending the interior points are usually taken to be the same
    as the endpoints of the interval on which the spline will be evaluated.
    """
    n , a , b = p[0] , p[1] , p[2]
    if knots is None:
        knots = np.linspace(a , b , n + 1 - deg)
    augbreaks = np.concatenate(( a * np.ones((deg)),knots, b * np.ones((deg))))
    m = len(augbreaks) - deg - 1
    v = np.empty((m, len(x)))
    d = np.eye(m, len(augbreaks))
    for i in range(m):
        v[i] = spl.splev(x, (augbreaks, d[i], deg),order)
    return v.T
@jit
def cheb_diff(p):
    """Differentiating matrix for Chebyshev polynomials
    Returns a matrix which multiplied from the right with the coefficients
    of the Chebyshev polynomial returns the derivative of the respective 
    Chebyshev polynomial. Can be used instead to evaluate the basis matrix of
    the derivative of a Chebyshev polynomial.
    Parameters
    ----------
    p : array_like
        Parameter array containing:
         - the number of knots = degree + 1 of the polynomial
         - the lower bound of the approximation
         - the upper bound of the approximation
    Returns
    -------
    D : ndarray
       Returns an upper triangular derivative operator matrix 
    Notes
    -----
    See usage in funbas
    """
    n , a , b = p[0] , p[1] , p[2]  
    D = np.zeros((n,n))
    for j in range(n):
        for i in range(int((n-j)/2)):
            D[j,j+1+2*i] = 4*((2*i+j+1))/(b-a)
    D[0,:] = D[0,:]/2
    return D
@jit
def mono_diff(p):  
    """Differentiating matrix for monomials
    Returns a matrix which multiplied from the right with the coefficients
    of the monomial returns the derivative of the respective 
    monomial. Can be used instead to evaluate the basis matrix of
    the derivative of a monomial.
    Parameters
    ----------
    p : array_like
        Parameter array containing:
         - the number of knots = degree + 1 of the polynomial
         - the lower bound of the approximation
         - the upper bound of the approximation
    Returns
    -------
    D : ndarray
       Returns an upper triangular derivative operator matrix 
    Notes
    -----
    See usage in funbas
    """
    n , a , b = p[0] , p[1] , p[2] 
    D = np.zeros((n,n))
    for j in range(n-1):
        D[j,j+1] = (j+1)/(b-a)*2
    return D
def funbas(p,x,order = None, polynomial = None):
    """Creating a multidimensional approximation basis matrix
    Returns a matrix which is a tensor product of the basis matrices 
    Parameters
    ----------
    p : array_like
        Parameter matrix where each row contains:
         - the number of knots(for splines) or  degree + 1
         of the polynomial (Chebyshev or Monomials)
         - the lower bound of the approximation
         - the upper bound of the approximation
    x : array_like
        Matrix of evaluation points, one column for each dimension
        (note the difference with p)
    order : array_like
        Specifies which derivatives should be evaluated - default is zero
        Can only contain natural numbers (integration is not done yet)
    polynomial : array_like
        Specifies the type of approximation to be used, one for each dimension:
        default is Chebyshev for all dimensions. 
        Can only contain the following strings:
        - 'cheb' for Chebyshev
        - 'mono' for Monomials
        - 'spli' for Cubic Splines - cannot pre-specify knots
        
    Returns
    -------
    Phi0 : ndarray
       Returns the n dimensional (n = no. of rows in p) basis matrix.
    Notes
    -----
    To perform a multidimensional approximation, first create a product of the evaluation
    points for each dimension - see the example
    """
    
        
    if x.ndim==1:
        if order is None:
            order = []
            order[0] = 0
        if polynomial is None:
            polynomial = []
            polynomial[0] = 'cheb'
        if polynomial=='cheb':
            Phi0=cheb_basex(p, x) @ np.linalg.matrix_power(cheb_diff(p),order[0])
        elif polynomial=='mono':
            Phi0=mono_basex(p, x) @ np.linalg.matrix_power(mono_diff(p),order[0])
        elif polynomial=='spli':
            Phi0=spli_basex(p, x, order = order[0]) 
    else:
        Phi0=np.ones((x.shape[0],1))  
        if order is None:
            order = np.zeros(x.shape[1], dtype=np.int)
        if polynomial is None:
            polynomial = ["cheb" for j in range(x.shape[1])]
        for j in range(x.shape[1]):
            if polynomial[j]=='cheb':
                Phi1=cheb_basex(p[j,:], x[:,j]) @ np.linalg.matrix_power(cheb_diff(p[j,:]),order[j])
            elif polynomial[j]=='mono':
                Phi1=mono_basex(p[j,:], x[:,j]) @ np.linalg.matrix_power(mono_diff(p[j,:]),order[j])
            elif polynomial[j]=='spli':
                Phi1=spli_basex(p[j,:], x[:,j], order = order[j]) 
            Phi0=dprod(Phi1,Phi0)
    return Phi0
@jit
def goldenx(f,a,b,tol,*arg):
    """Vectorized golden section search to maximize univariate functions simultaneously
    Returns the maximum and the maximal value of f
    Parameters
    ----------
    f : function
        Function that maps from R^n to R^n such that it is an augmented univariate function
    a : array_like
        The lower bound for the maximization, for each dimension
    b : array_like
        The upper bound for the maximization, for each dimension
    tol : float
        Specifies the default tolerance value (for the argument of f) - default is 10**(-10)        
    Returns
    -------
    x1 : ndarray
       Returns the n dimensional solution of the maximization.
    f1 : ndarray
       Returns the n dimensional maximum values.
    Notes
    -----
    What is missing: additional arguments can be passed to f but not as a tuple (difference to scipy optimize).
    As a consequence, tolerance has to be defined all the time. The sig function is not nice. But pretty fast
    """ 
    
    alpha1 = (3.0 - np.sqrt(5)) / 2.0
    alpha2 = 1.0 - alpha1
    d  = b - a
    x1 = a + alpha1 * d
    x2 = a + alpha2 * d
    s  = np.ones(x1.shape)
    f1 = f(x1,*arg)
    f2 = f(x2,*arg)
    d = alpha1 * alpha2 * d
#    if((b-a).min() <0):
#        print('Fix the bounds yo')
#        conv = 0.0
#    else:
#       conv = 2.0   
    conv = 2.0
    while conv > tol:        
        i = np.greater(f2,f1)
        not_i = np.logical_not(i)
        sig = np.ones(x1.shape)
        sig[not_i] = -1.0
        x1[i] = x2[i]
        f1[i] = f2[i]
        d = alpha2 * d
        x2 = x1 + np.multiply(np.multiply(np.multiply(s,(i^(not_i))),d),sig)
        s = np.sign(x2-x1)
        f2 = f(x2,*arg)
        conv = np.max(d)          
    return x1, f1
def equidistant_nonlin_grid(orig_grid,newgrid_size):
    """Create a fine 1D grid,equidistant between the points of the original grid 
    Returns the new grid
    Parameters
    ----------
    orig_grid : (n,0) array_like
    
    newgrid_size : int
    Specifies the size of the new grid
     
    Returns
    -------
    fine_grid : ndarray
       Returns the finer grid containing the original gridpoints.
    Notes
    -----
    """ 
    orig_grid_size = orig_grid.shape[0]
    no_equidistant = np.int((newgrid_size)/(orig_grid_size-1))
    remainder = newgrid_size - no_equidistant *(orig_grid_size-1) 
    fine_grid = np.zeros((newgrid_size))
    prev_grid_end = 0
    for i in range(orig_grid_size-1):
        endpoint_included = False
        if i< remainder:
            points_between = no_equidistant + 1
        elif(i ==(orig_grid_size-2)):
            points_between = no_equidistant
            endpoint_included = True
        else:
            points_between = no_equidistant
        fine_grid[prev_grid_end:(prev_grid_end + points_between) ] = np.linspace(orig_grid[i],orig_grid[i+1],points_between,endpoint = endpoint_included)
        prev_grid_end = prev_grid_end + points_between
    return fine_grid