import numpy as np
from roots import roots
import scipy.optimize
import scipy as sc

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
    """

    # BEWARE: The lambda is just looking up the global value of the variables.
    
    Lambda = np.zeros(m);
    w      = np.zeros(m);
    alpha  = np.zeros(m);
    phi    = np.empty(m, dtype='object');
    
    for j in range(m):

        i = j + 1

        if np.mod(i,2) == 0: # even

            lowerBound = (i-0.999)*np.pi/L
            upperBound = 0.999*i*np.pi/L
            # midPoint = (lowerBound + upperBound)/2.
            f = lambda ww: 1./l*np.tan(ww*L/2.)+ww
            # w[j] = roots(f, lowerBound, upperBound, eps=1.e-4)
            w[j] = sc.optimize.brentq(f, lowerBound, upperBound, xtol=1e-6, rtol=1.e-10)
            print 'i*np.pi/L, w', i*np.pi/L, w[j],
            Lambda[j] = 2.*l / (1. + w[j]**2 * l**2)
            
            alpha[j] = 1. / np.sqrt(L/2. - np.sin(w[j]*L)/(2.*w[j]))
            print 'alpha', alpha[j]
            aalpha = alpha[j]
            W = w[j]
            phi[j] = lambda z, aalpha, W: aalpha * np.sin(W*(z - L/2.))

        else: # odd
            lowerBound = (i-0.999)*np.pi/L
            upperBound = 0.999*i*np.pi/L
            # midPoint = (lowerBound + upperBound)/2.
            f = lambda ww: 1./l-ww*np.tan(ww*L/2.)
            # w[j] = roots(f, lowerBound, upperBound, eps=1.e-4)
            w[j] = sc.optimize.brentq(f, lowerBound, upperBound, xtol=1e-6, rtol=1.e-10)
            print 'i*np.pi/L, w', i*np.pi/L, w[j],
        
            Lambda[j] = 2.*l / (1. + w[j]**2 * l**2)
        
            alpha[j] = 1. / np.sqrt(L/2. + np.sin(w[j]*L)/(2.*w[j]))
            print 'alpha', alpha[j]
            aalpha = alpha[j]
            W = w[j]
            phi[j] = lambda z, aalpha, W: aalpha * np.cos(W*(z - L/2.))

    return Lambda, phi, w, alpha


