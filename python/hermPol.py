import numpy as np
import scipy as sp
import scipy.special

# In Python only parameters are held in the lambda scope.
# L, l, alpha, w are taken from the surrounding scope.
# ==> An intermediate scope is needed
def h_k(k):
    return lambda x: sp.special.eval_hermitenorm(k,x)/np.sqrt(sp.special.factorial(k))

def hermPol(p):
    """
    Returns the normalized Hermite polynomials of maximum degree p 
    """

    h = np.empty(p, dtype=object)
    for k in range(p):
        # h[k] = lambda x,k: sp.special.eval_hermitenorm(k,x)/np.sqrt(sp.special.factorial(k))
        h[k] = h_k(k)
        
    return h

if __name__=='__main__':

    from pylab import *

    p = 4

    h = hermPol(p)

    x = np.linspace(-2,2,101)

    figure()

    for i in range(p):
        plot(x,h[i](x))

    show()
