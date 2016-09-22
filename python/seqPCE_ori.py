import numpy as np
from scipy.special import factorial

def seqPCE( d,p ):
    """
    % returns the sequences of multidimensional polynomial is the PC expansion
    %  d: dimension of the input space
    %  q: maximum polynomial degree
    %  Seq: matrix of dimensions P x d, where P = (q+d-1)!/((d-1)!q!)
    """

    n = 0;
    Seq = [];
    seq = np.empty(d);

    for q in range(0,p+1):

        First = 0; 
    
        # Box positions of each ball
        BoxPos = np.zeros(d-1);
    
        # Number of polynomials for current q
        nump = int(factorial(q+d-1)/(factorial(d-1)*factorial(q)))

        for jj in range(nump):

            if First == 0:
                # initialize: fill the first d - 1 boxes
                for i in range(d-1):
                    ii = i + 1
                    BoxPos[i] = ii;

                First = 1;
            else:
                # find rightmost ball
                # [pp,ii] = np.max(BoxPos);
                ii = np.argmax(BoxPos);
                pp = BoxPos[ii]

                # if rightmost ball in rightmost box
                if pp == q+d-1: # ### q+d-2 ???
                    # find rightmost ball that can be shifted
                    kk = 1;
                    while(kk == 1):
                        kk = BoxPos[ii]-BoxPos[ii-1];
                        ii = ii-1;

                    # shift rightmost ball
                    BoxPos[ii] = BoxPos[ii]+1;
                    # return balls to immediate right
                    for i in range(ii+1,d-1):
                        BoxPos[i] = BoxPos[i-1]+1;

                else:
                    # normal increment
                    BoxPos[ii] = pp+1;


                    
            # convert balls to polynomial sequences
            n = n+1;
            
            seq[0] = BoxPos[0]-1;
            for i in range(1,d-1):
                seq[i] = BoxPos[i]-BoxPos[i-1]-1; 

            seq[d-1] = q+d-1-BoxPos[d-2];
            Seq.append(seq)

    return Seq

            
if __name__=='__main__':

    from pylab import *

    d = 4;
    p = 2;

    for q in range(0,p+1):
        nump = factorial(q+d-1)/(factorial(d-1)*factorial(q));
        print nump
    
    Seq = seqPCE( d,p )

