import numpy as np

def avgVarianceErrorKL(m, L, Lambda):
    """
    Returns the average variance error of the KL expansion with
    exponential kernel

    Parameters
    ----------
    m: number of random variables
    L: Length of the domain
    l: correlation length
    lambda: vector of eigenvalues
    phi: cell array of eigenfunctions

    Returns
    -------
    errVar : 1-rank array 
         average variance error
    """

    errVar = np.ones(m);
    errVar0 = 1.
    errVar[0] = errVar0 - 1./L * Lambda[0]
    for i in range(1,m):
        errVar[i] = errVar[i-1] - 1./L * Lambda[i]

    return errVar
