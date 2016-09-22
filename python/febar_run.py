import numpy as np
import string
from febar import febar

nele = 3
q    = 1.
L    = 10.
# D    = np.random.random(nele)
D = np.array([0.4, 0.7, 0.25])

u = febar( D, q, L, nele )

print u
print string.join(["\n%12.5f"%uu for uu in u])
