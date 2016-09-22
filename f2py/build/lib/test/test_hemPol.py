import numpy as np
import sys
sys.path.insert(0,'../')
from PCE import hermPol, hermPol2

if __name__=='__main__':

    from pylab import *

    p = 4

    x = np.linspace(-2,2,101)

    h = hermPol(p, x)

    figure()

    for i in range(p):
        plot(x,h[i,:])


    h = hermPol2(p, x)

    figure()

    for i in range(p):
        plot(x,h[i,:])

    show()
