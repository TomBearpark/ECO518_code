import numpy as np
from sympy import Matrix
mat =  np.array(
   [[-2,0,-1,0], 
    [0,0,0,0], 
    [1,0,0,0],
    [0,1,0,0]])

m = Matrix(mat)
m
P, J = m.jordan_form()
P
J
