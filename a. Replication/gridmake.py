# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 15:55:00 2022

@author: Giorgia
"""
#https://github.com/maitlahcen/CompEcon-python/blob/master/compecon/tools.py
    
from functools import reduce
import numpy as np
from scipy.linalg import qz
import time

from scipy.sparse import identity


def gridmake(*arrays):
    """
    Forms grid points
    USAGE
    X = gridmake(x1,x2,...,xn)
    X1,..., Xn = gridmake(x1,x2,...,xn)
    Expands matrices into the associated grid points.
    If N is the 2xd array that indexes the size of the inputs, GRIDMAKE returns a sum(N[0]) by prod(N[1]) array.
    The output can also be returned as either
      d matrices or
      sum(N(:,2)) matrices
    Note: the grid is expanded so the last variable change most quickly.
    Example:
    X = gridmake([1, 2, 3], [4, 5])
    array([[1, 1, 2, 2, 3, 3],
           [4, 5, 4, 5, 4, 5]])
    Also the inputs need not be vectors.
    Y = gridmake(X, [10, 20])
    array([[ 1,  1,  1,  1,  2,  2,  2,  2,  3,  3,  3,  3],
           [ 4,  4,  5,  5,  4,  4,  5,  5,  4,  4,  5,  5],
           [10, 20, 10, 20, 10, 20, 10, 20, 10, 20, 10, 20]])
    """
    if len(arrays) == 1:
        return arrays[0]

    arrays = np.atleast_2d(*arrays)
    n = len(arrays)
    idx = np.indices([a.shape[1] for a in arrays]).reshape([n, -1])
    return np.vstack(arrays[k][:, idx[k]] for k in range(n))
