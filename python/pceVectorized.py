import numpy as np
import scipy.stats
from scipy.special import factorial
from sklearn.neighbors import KernelDensity
from pylab import *
import string

from EigFcnKL import EigFcnKL
from relVarianceErrorKL import relVarianceErrorKL
from avgVarianceErrorKL import avgVarianceErrorKL
from febar import febar
from assembleK import assembleK
from hermPol import  hermPol
from seqPCE import seqPCE

###############################################################
# Vectorized PCE
###############################################################

# Linear elestic bar of length L=2m
L = 2.

# The axial resistance is described by a random field with lognormal marginal distribution
# mean: mu_D = 100kN
# coefficient of variation: delta_D = std_D / mu_D = 0.2
mu_D    = 100.
delta_D = 0.2
std_D   = delta_D * mu_D

# The auto-correlation coefficient function of the underlying Gaussian random
# field is given by the following exponential model, with correlation length=1m
l = 0.2 # [m]
# l = 1. # [m]

# distributed load
q    = 1.  # [kN/m]


# eigenvalues and eigenfunctions of the KL expansion with exponential kernel
m_max = 200
Lambda, phi, w, alpha = EigFcnKL(m_max , L , l)

# Make a grid to evaluate the eigenfunctions and the errors
N_nodes    = 21
N_elements = N_nodes - 1
z     = np.linspace(0, L, N_nodes)
z_mid = (z[:-1] + z[1:]) / 2.

# relative variance error at midpoints
errVar = relVarianceErrorKL(m_max, L, l, Lambda, phi, z_mid)

# average variance error at midpoints
avgErr = avgVarianceErrorKL(m_max, L, Lambda)

# D(z,theta) is described by a homogeneous RF with lognormal marginal distribution
# Transform the parameters of the lognormal distribution to a normal distribution
sigma_U =  np.sqrt(np.log(delta_D**2 + 1.))
mu_U =  np.log(mu_D) - sigma_U**2/2.

# find the first m-terms that ensure an average variance error < 0.05 
avgErr_max = 0.05

m = 0
while m < m_max:
    m = m + 1
    if avgErr[m] < avgErr_max:
        break
print 'N_terms, avgErr[m]', m+1, avgErr[m]

# we will work here with d instead of m
d = m+1;

# We should investigate order 1-3 of the Hermite polynomials
p_max = 3

# Calculate Hermite Polynomials
hermPolNorm = hermPol(p_max+1)

# PCE sequences
p = p_max
nump = int(np.round(factorial(p+d)/(factorial(p)*factorial(d))))
print 'nump', nump
Seq = seqPCE( d,p )

# Monte Carlo parameters
N_realizations = 1
# realization
i = 0


# "random" sequences
np.random.seed(seed=12345)
UU = np.random.randn(N_realizations,d)

# Evaluate PCEs
Psi = np.ones(nump)
for kk in range(nump):
    for jj in range(d):
        Psi[kk]= Psi[kk] * hermPolNorm[Seq[kk,jj]](UU[i,jj])

# print Psi
