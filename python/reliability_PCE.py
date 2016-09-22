from pylab import *
from scipy.special import factorial
from scipy.stats import norm
from scipy.optimize import fsolve
from pce1D import pce1D

# plot text defaults
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('text', usetex=True)

# clear all
# clc

# Reliability analysis with PCE
# g = R - S

# S is deterministic
S = 3

# R is lognormal
# mean of R
mu_R = 10
# CV of R
delta_R = 0.3
# mean and standard deviation of underlying normal
sigma = sqrt(log(delta_R**2+1));
mu = log(mu_R)-0.5*sigma**2;

# Exact probabability of failure
Pf_exact = norm.cdf(log(S), mu, sigma)

# Coefficients of polynomial chaos
# maximum order 
p = 11
a = zeros(p+1)
var_PCE = zeros(p)
Pf_PCE  = zeros(p)

for i in range(p+1):  
    a[i] = exp(mu+0.5*sigma**2)*sigma**i/sqrt(factorial(i))

# Evaluate the variance of the PCE approximation
for i in range(p):
    var_PCE[i] = sum(a[1:i+2]**2)


# Plot var_PCE/var_exact in terms of p
var_exact = (delta_R*mu_R)**2

figure(1)
plot(range(p),var_PCE/var_exact,'o-')
xlabel(r'p')
ylabel(r"Var(R)$_{PCE}$/Var(R)$_{exact}$", fontsize=16)


# Evaluate the probability of failure for each polynomial degree
for i in range(p):
    def f(x): return pce1D(a,i,x)-S
    u0 = fsolve(f,0)
    Pf_PCE[i] = norm.cdf(u0)


# Plot logPf_PCE/logPf_exact in terms of p
figure(2)
plot(range(p),log(Pf_PCE)/log(Pf_exact),'o-')
xlabel(r'p')
ylabel(r'log(Pf$_{PCE}$)/log(Pf$_{exact}$)')


# plot PCE approximation
figure(3)
xp=arange(-5,5.1,0.1)
for i in range(p):
    plot(xp,pce1D(a,i,xp));
xlabel(r'U')
ylabel(r'R$_{PCE}$(U)')
legend(('p = 1','p = 2','p = 3','p = 4','p = 5','p = 6'),loc='best')


