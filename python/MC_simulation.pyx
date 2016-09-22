import numpy as np
import scipy.stats

ctypedef void (*func_type)()

def monte_carlo_simulation(double L,
                           double l,
                           double q,
                           double mu_U,
                           double sigma_U,
                           double[:] Lambda,
                           func_type[:] phi,
                           int N_terms,
                           double[:] z_mid,
                           int N_repetitions,
                           int N_realizations):


    cdef int N = X.shape[0]
    cdef int D = X.shape[1]
    cdef double[:] Y = np.zeros(N)
    cdef int i, j, ii
    cdef double r = 0

    cdef int i, j, ii
    cdef int N_elements = len(z_mid)
    cdef double[:] mean_U  = np.zeros(N_repetitions)
    cdef double[:] std_U   = np.zeros(N_repetitions)
    cdef double[:] delta_U = np.zeros(N_repetitions)
    cdef double[:,:] u_tip = np.zeros((N_repetitions, N_realizations))
    cdef double[:,:] hlds  = np.zeros((N_realizations,N_terms))
    cdef double[:,:] lhd   = np.zeros((N_realizations,N_terms))
    cdef double[:,:] UU    = np.zeros((N_realizations,N_terms))
    cdef double[:,:] U     = np.zeros(N_terms)
    cdef double[:,:] U     = np.zeros(N_terms)
    cdef int[:] pseed
    cdef double[:] U_RF = np.zeros(N_elements)
    cdef double[:] D_RF = np.zeros(N_elements)
    cdef double[:] u_RR = np.zeros(N_elements+1)
    cdef double[:,:] u_tip = np.zeros((N_repetitions, N_realizations))

    # Experiment repetitions
    for j in range(N_repetitions):

        # "random" sequences
        if case=='pseudo':
            UU = randn(N_realizations,N_terms)
        elif case=='halton':
            import ghalton
            # sequencer = ghalton.Halton(N_terms)
            # permutation seed | 100 permutation rules avaible in ghalton
            pseed = np.arange(100); np.random.shuffle(pseed)
            sequencer = ghalton.GeneralizedHalton(N_terms, pseed[j])
            hlds = np.array(sequencer.get(N_realizations))
            # transform to normal
            UU = scipy.stats.norm.ppf(hlds, loc=0., scale=1.)
        elif case=='lhs':
            from pyDOE import lhs
            lhd = lhs(N_terms, samples=N_realizations)
            # transform to normal
            UU = scipy.stats.norm.ppf(lhd, loc=0., scale=1.)

        # loop over realizations
        for i in range(N_realizations):

            # Gaussian RVs
            U = UU[i,:]

            # Initialize Gaussian RF 
            U_RF = np.zeros_like(z_mid)

            # Calculate Gaussian RF 
            for ii in range(N_terms):
                U_RF = U_RF + np.sqrt(Lambda[ii])*phi[ii](z_mid)*U[ii]

            U_RF *= sigma_U
            U_RF += mu_U

            # Transformation to distribution of axial resistance (lognormal)
            D_RF = np.exp(U_RF)

            # finite element solution (random realization)
            u_RR = febar(D_RF, q, L, N_elements)

            # store tip displacements
            u_tip[j,i] = u_RR[-1]


    return u_tip
