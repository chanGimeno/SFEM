import numpy as np
import scipy as sp
from scipy.special import factorial

# In Python only parameters are held in the lambda scope.
# L, l, alpha, w are taken from the surrounding scope.
# ==> An intermediate scope is needed
def h_0(k):
    return lambda x: np.ones_like(x)

def h_1(k):
    return lambda x: x

def h_k(k):
    return lambda x: x*h[k-1](x)-(k-1)*h[k-2](x)

def h_norm(kk):
    return lambda xx: h[kk](xx)/np.sqrt(factorial(kk))

def hermPol(p):
    """
    Returns the normalized Hermite polynomials of maximum degree p 
    """

    h = np.empty(p, dtype=object)
    H = np.empty(p, dtype=object)

    for k in range(p):
        
        if k == 0:
            h[k] = h_0(k)
        elif k == 1:    
            h[k] = h_1(k)
        else:        
            h[k] = h_k(k)
        
    # normalize Hermite polynomials
    for i in range(p):
        H[i] = h_norm(i)
        # h[k] = H

    return H

if __name__=='__main__':

    from pylab import *

    p = 4

    h = hermPol(p)

    x = np.linspace(-2,2,101)

    figure()

    for i in range(p):
        # plot(x,h[i](x))
        plot(x,h[i](x)/np.sqrt(factorial(i)))

    H = np.empty(p, dtype=object)
    # normalize Hermite polynomials
    for k in range(p):
        H[k] = h_norm(k)
        # h[k] = H
        
    figure()

    for i in range(p):
        plot(x,H[i](x))
        
    show()
