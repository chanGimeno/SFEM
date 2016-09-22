import numpy as np
import scipy.sparse.linalg
import scipy as sp 
from assembleK import assembleK

def febar( D, q, L, nele ):
    """
    returns the nodal displacements of a bar subjected to uniform axial loading
    D: vector with axial resistances at each element
    q: uniformly distributed load
    L: length of bar
    nele: number of elements
    """

    # number of nodes
    nnode = nele+1

    # element length
    h = 1.*L/nele

    # consistent nodal forces
    F = np.zeros((nnode,1))

    # assemble force vector
    F[:]  = h*q
    F[0]  = F[0]/2 # boundary conditions
    F[-1] = F[-1]/2 # boundary conditions

    # assemble stiffness matrix
    Ke = 1.*D/h

    K = assembleK(Ke,nele,nnode)

    # reduce stiffness matrix and force vector to apply fixed end boundary
    # condition

    Kred = K[1:nnode,1:nnode]
    Fred = F[1:nnode]

    # solution for nodal displacements
    # u = np.linalg.solve(Kred, Fred)
    u = sp.sparse.linalg.spsolve(Kred, Fred)

    # add fixed end boundary condition
    u = np.hstack((0.,u))

    return u

