import numpy as np
from roots import roots
import scipy.optimize
import scipy as sc

# In Python only parameters are held in the lambda scope.
# L, l, alpha, w are taken from the surrounding scope.
# ==> An intermediate scope is needed
def phi_even(L, l, alpha, w):
    return lambda z: alpha * np.sin(w*(z - L/2.))

def phi_odd(L, l, alpha, w):
    return lambda z: alpha * np.cos(w*(z - L/2.))

def EigFcnKL(m , L , l):
    """
    Returns the eigenvalues and eigenfunctions of the KL expansion with
    exponential kernel

    Parameters
    ----------
    m: number of random variables
    L: Length of the domain
    l: correlation length

    Returns
    -------
    lambda: vector of eigenvalues
    phi: cell array of eigenfunctions
    w: vector 
    alpha: vector 
    """

    Lambda = np.zeros(m);
    w      = np.zeros(m);
    alpha  = np.zeros(m);
    phi    = np.empty(m, dtype='object');  # roughly like the matlab "cell"
    
    for j in range(m):

        i = j + 1

        if np.mod(i,2) == 0: # even

            lowerBound = (i-0.999)*np.pi/L
            upperBound = 0.999*i*np.pi/L
            f = lambda ww: 1./l*np.tan(ww*L/2.)+ww
            w[j] = sc.optimize.brentq(f, lowerBound, upperBound, xtol=1e-6, rtol=1.e-10)
            Lambda[j] = 2.*l / (1. + w[j]**2 * l**2)
            alpha[j] = 1. / np.sqrt(L/2.-np.sin(w[j]*L)/(2.*w[j]))
            Alpha = alpha[j]
            W = w[j]
            phi[j] = phi_even(L, l, Alpha, W)
            
        else: # odd
            lowerBound = (i-0.999)*np.pi/L
            upperBound = 0.999*i*np.pi/L
            f = lambda ww: 1./l-ww*np.tan(ww*L/2.)
            w[j] = sc.optimize.brentq(f, lowerBound, upperBound, xtol=1e-6, rtol=1.e-10)
            Lambda[j] = 2.*l / (1. + w[j]**2 * l**2)
            alpha[j] = 1. / np.sqrt(L/2. + np.sin(w[j]*L)/(2.*w[j]))
            Alpha = alpha[j]
            W = w[j]
            phi[j] = phi_odd(L, l, Alpha, W)

    return Lambda, phi, w, alpha


