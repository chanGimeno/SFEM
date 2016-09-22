import numpy as np
import scipy as sp
from scipy.special import factorial

def h_k(k,x):
    if k == 0:
        return np.ones_like(x)
    elif k == 1:    
        return x
    else:        
        return x*h_k(k-1,x)-(k-1)*h_k(k-2,x)

def h_k_norm(k,x):
    """
    Returns the normalized Hermite polynomial of degree k evaluated at x 
    """
    return h_k(k,x)/np.sqrt(factorial(k))

def hermPol(p,x):
    """
    Returns the normalized Hermite polynomials of maximum degree p 
    """

    h = np.empty(p, dtype=object)

    for k in range(p):
        h[k] = h_k_norm(k,x)

    return h

# def hermPol_lamda(p):
#     """
#     Returns the normalized Hermite polynomials of maximum degree p 
#     """

#     h = np.empty(p, dtype=object)
#     for k in range(p):
#         h[k] = lamda x: h_k(k,x)
        
#     return h
#     """
#     Returns the normalized Hermite polynomial of degree k evaluated at x 
#     """
#     return h_k(k,x)/np.sqrt(factorial(k))

# def hermPol(p):
#     """
#     Returns the normalized Hermite polynomials of maximum degree p 
#     """

#     h = np.empty(p, dtype=object)
#     H = np.empty(p, dtype=object)

#     for k in range(p):
        
#         if k == 0:
#             h[k] = h_0(k)
#         elif k == 1:    
#             h[k] = h_1(k)
#         else:        
#             h[k] = h_k(k)
        
#     # normalize Hermite polynomials
#     for i in range(p):
#         H[i] = h_norm(i)
#         # h[k] = H

#     return H

if __name__=='__main__':

    from pylab import *

    p = 4

    x = np.linspace(-2,2,101)

    h = hermPol(p, x)

    figure()

    for i in range(p):
        # plot(x,h_k_norm(i,x))
        plot(x,h[i])

    show()
