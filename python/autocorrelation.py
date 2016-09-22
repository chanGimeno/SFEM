import numpy as np

# exponential auto-correlation function
def auto_correlation(x,y,l_corr,sigma=1.):
    """
    exponential auto-correlation
    """
    distance = np.linalg.norm(x-y)
    return sigma**2 * np.exp(-distance/l_corr)
