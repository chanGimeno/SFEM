import numpy as np

def relVarianceErrorKL(m, L, l, Lambda, phi, z):
    """
    Returns the relative variance error of the KL expansion with
    exponential kernel
    
    Parameters
    ----------
    m: number of random variables
    L: Length of the domain
    l: correlation length
    lambda: vector of eigenvalues
    phi: cell array of eigenfunctions
    z: longitudinal coordinate

    Returns
    -------
    errVar : 2-rank array 
         relative variance error
    """
    
    errVar = np.ones((len(z),m))
    errVar[:,0] = errVar[:,0] - Lambda[0]*phi[0](z)**2
    for i in range(1,m):
        errVar[:,i] = errVar[:,i-1] - Lambda[i]*phi[i](z)**2

    return errVar

