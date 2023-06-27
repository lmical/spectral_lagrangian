import numpy as np

def compute_integral_error(c,c_exact):  # c is dim x times
    times=np.shape(c)[1]
    error=0.
    for t in range(times):
        error = error + (c[0,t]-c_exact[0,t])**2.
    error = np.sqrt(error/times) 
    return error

def compute_integral_error_approx(tt,c,tt_exact,c_exact):  # c is dim x times
    error=0.
    iex=0
    for it, t in enumerate(tt):
        while (t> tt_exact[iex]+1e-5):
            iex+=1
        error = error + np.linalg.norm(c[:,it]-c_exact[:,iex],2)**2.
    error = np.sqrt(error/len(tt)) 
    return error