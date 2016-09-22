import numpy as np
from scipy.optimize import fsolve


def EigFcnKL(m , L , l):
    """
    Returns the eigenvalues and eigenfunctions of the KL expansion with
    exponential kernel

    m: number of random variables
    L: Length of the domain
    l: correlation length
    lambda: vector of eigenvalues
    phi: cell array of eigenfunctions
    """

    Lambda = np.zeros(m);
    phi = np.empty(m, dtype='object');

    for j in range(m):

        i = j + 1

        if np.mod(i,2) == 0: # even

            lowerBound = (i-0.999)*np.pi/L
            upperBound = 0.999*i*np.pi/L
            midPoint = (lowerBound + upperBound)/2.
            f = lambda w: 1./l*np.tan(w*L/2.)+w
            w = fsolve(f, midPoint)
            print 'i*np.pi/L, w', i*np.pi/L, w
            Lambda[j] = 2.*l / (1. + w**2 * l**2)
            
            alpha = 1. / np.sqrt(L/2. - np.sin(w*L)/(2.*w))
            phi[j] = lambda z: alpha * np.sin(w*(z - L/2.))
    
        else: # odd
            lowerBound = (i-0.999)*np.pi/L
            upperBound = 0.999*i*np.pi/L
            midPoint = (lowerBound + upperBound)/2.
            f = lambda w: 1./l-w*np.tan(w*L/2.)
            w = fsolve(f, midPoint)
            print 'i*np.pi/L, w', i*np.pi/L, w
        
            Lambda[j] = 2.*l / (1. + w**2 * l**2)
        
            alpha = 1. / np.sqrt(L/2. + np.sin(w*L)/(2.*w))
            phi[i] = lambda z: alpha * np.cos(w*(z - L/2.))

    return Lambda, phi


