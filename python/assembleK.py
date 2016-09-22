import numpy as np
from scipy import sparse

def assembleK(Ke, nele, nnode):
    """
    assembles the global stiffness matrix of bar
    Ke: vector of element stiffness
    nele: number of elements
    nnode: number of nodes
    """

    # initialize global stiffness matrix of bar
    K = np.zeros((nnode,nnode))
    K = sparse.csr_matrix(K, shape=(nnode,nnode))

    # assembles the global stiffness matrix
    idx  = range(0,nele)
    idx2 = range(1,nele+1)
    K = K + sparse.csr_matrix(( Ke, (idx, idx )), shape=(nnode,nnode))
    K = K + sparse.csr_matrix((-Ke, (idx, idx2)), shape=(nnode,nnode))
    K = K + sparse.csr_matrix((-Ke, (idx2,idx )), shape=(nnode,nnode))
    K = K + sparse.csr_matrix(( Ke, (idx2,idx2)), shape=(nnode,nnode))

    return K


if __name__=='__main__':

    from pylab import *

    nele = 3
    q    = 1.
    L    = 10.
    # D    = np.random.random(nele)
    D = np.array([0.4, 0.7, 0.25])

    # number of nodes
    nnode = nele+1

    # element length
    h = 1.*L/nele

    # assemble stiffness matrix
    Ke = 1.*D/h

    K = assembleK(Ke,nele,nnode)





