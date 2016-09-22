import numpy as np

# exponential auto-correlation function
def auto_correlation_KL(x, y, Lambda, phi, m):
    """
    exponential auto-correlation. KL decomposition
    """
    distance = np.linalg.norm(x-y)
    AC = 0.
    for i in range(m):
        AC = AC + Lambda[i]*phi[i](x)*phi[i](y)
    return AC
