from pylab import *

# clear all;
# clc;

# uniformly distributed pseudo-random samples
rand

# set seed
#rng(0);

# get n  pseudo random samples from the uniform distribution
n = 100;
v = rand(n,1);

# get Latin Hypercube samples
vlhs = lhsdesign(n,1);

# get Quasi random samples (Halton sequence)
p = haltonset(1,'Skip',1e3,'Leap',1e2);
vqrs = p(1:n,:);

# transform to samples of the normal distribution
mu = 10;
sigma = 5;
x = norminv(v,mu,sigma);
#x = normrnd(mu,sigma,n,1);

xlhs = norminv(vlhs,mu,sigma);
xqrs = norminv(vqrs,mu,sigma);

# histogram
figure(1);
subplot(3,1,1)       
hist(x)
title('Pseudo random')

subplot(3,1,2)       
hist(xlhs);
title('Latin Hypercube')

subplot(3,1,3)       
hist(xqrs);
title('Quasi random')


[f,xi1] = ksdensity(x);
[flhs,xi2] = ksdensity(xlhs);
[fqrs,xi3] = ksdensity(xqrs);

#------------------------------------------------------------------------------------------
# http://stackoverflow.com/questions/28572731/matlab-ksdensity-equivalent-in-python

# kernel density estimation -- SKLEARN
import numpy as np
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity

np.random.seed(12345)
# similar to MATLAB ksdensity example x = [randn(30,1); 5+randn(30,1)];
Vecvalues=np.concatenate((np.random.normal(0,1,30), np.random.normal(5,1,30)))[:,None]
Vecpoints=np.linspace(-8,12,100)[:,None]
kde = KernelDensity(kernel='gaussian', bandwidth=0.5).fit(Vecvalues)
logkde = kde.score_samples(Vecpoints)
plt.plot(Vecpoints,np.exp(logkde))
plt.show()

# kernel density estimation -- SCIPY

# That form of the ksdensity call automatically generates an arbitrary
# x. scipy.stats.gaussian_kde() returns a callable function that can be
# evaluated with any x of your choosing. The equivalent x would be
# np.linspace(data.min(), data.max(), 100).

import numpy as np
from scipy import stats

data = ...
kde = stats.gaussian_kde(data)
x = np.linspace(data.min(), data.max(), 100)
p = kde(x)
#------------------------------------------------------------------------------------------



figure(2);

plot(xi1,f,'LineWidth',2);
hold on;
plot(xi2,flhs,'LineWidth',2,'Color','green');
plot(xi3,fqrs,'LineWidth',2,'Color','red');
plot(xi1,normpdf(xi1,mu,sigma),'LineWidth',2,'Color','black');
legend('PRS','LHS','QRS','Exact')
title('Probability density function')
hold off;



# lognormal random samples
u = [];
ulhs = [];
uqrs = [];

# LHS and QRS samples from the uniform distribution
vlhs = lhsdesign(1000,1);
p = haltonset(1,'Skip',1e3,'Leap',1e2);
vqrs = p(1:1000,:);

for i=1:100

n = i*10;
u = [u ; normrnd(0,1,10,1)];
ulhs = norminv(vlhs(1:n));
uqrs = norminv(vqrs(1:n));

y = exp(u);
ylhs = exp(ulhs);
yqrs = exp(uqrs);

# estimation of mean
mu_y(i) = 1/n*sum(y);
# mu_y(i) = mean(y);

mu_y_lhs(i) = mean(ylhs);
mu_y_qrs(i) = mean(yqrs);


# estimation of standard deviation
sigma_y(i) = sqrt(1/(n-1)*sum((y-mu_y(i)).^2));
# sigma_y(i) = std(y);

sigma_y_lhs(i) = std(ylhs);
sigma_y_qrs(i) = std(yqrs);


end

# exact mean
mu_y_exact = sqrt(exp(1));

figure(3);
plot(10:10:1000,mu_y,'LineWidth',2);
hold on;
plot(10:10:1000,mu_y_lhs,'LineWidth',2,'Color','green');
plot(10:10:1000,mu_y_qrs,'LineWidth',2,'Color','red');
plot([0 1000],[mu_y_exact mu_y_exact],'LineWidth',2,'Color','black');
legend('PRS','LHS','QRS','Exact')
title('Mean')
hold off;

# exact standard deviation
sigma_y_exact = sqrt(exp(1)*(exp(1)-1));

figure(4);
plot(10:10:1000,sigma_y,'LineWidth',2);
hold on;
plot(10:10:1000,sigma_y_lhs,'LineWidth',2,'Color','green');
plot(10:10:1000,sigma_y_qrs,'LineWidth',2,'Color','red');
plot([0 1000],[sigma_y_exact sigma_y_exact],'LineWidth',2,'Color','black');
legend('PRS','LHS','QRS','Exact')
title('Standard deviation')
hold off;




