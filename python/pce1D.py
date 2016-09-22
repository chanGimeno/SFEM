from hermPol import hermPol

def pce1D( a,p,x ):
    """
    %   Returns the 1D Hermite pce expansion
    %   a: vector of coefficients
    %   p: polynomial degree
    %   x: inpute of expansion
    """
    h = hermPol(p+1)
    val = 0

    for i in range(p+1):
   
        val = val+a[i]*h[i](x)

    return val

if __name__=='__main__':

    from pylab import *

    p = 4
    a = random(p+1)
    h = hermPol(p+1)

    x = random()

    val = pce1D( a,p,x )

    print val
