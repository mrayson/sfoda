"""
Numpy-like array operations
"""

import numpy as np

def accum1d(xin, yin, xout, method='mean'):
    """
    Accumulates values onto a vector by averaging multiple entries

    Inputs:
    ---
        xin: input vector to map with output
        yin: input vector of values to average
        xout: output vector
        method: 'mean' or 'sum' [default: 'mean']

    returns:
    ---
        yout: vector with average values
    """
    
    # Find the index of each input value in the output array
    idx = np.searchsorted(xout, xin)

    # This averages
    if method == 'mean':
        # Find the number of values in each bin
        N = np.bincount(idx)

        tmp = np.bincount(idx, weights = yin / N[idx])

    elif method == 'sum':

        tmp = np.bincount(idx, weights = yin)

    # Ensure that the output is on the same axis
    #iny = np.zeros(Z.shape, np.bool)
    #iny[idx] = True
    #outidx = np.unique(idx)

    nx = xout.shape[0]
    # need to check if the vector is longer than the inputs
    #  This will happen if the values in xin are beyond the range of xout
    if idx.max() == nx:
        iout = idx.max()
        tmp = tmp[:-1]
    else:
        iout = idx.max()+1

    yout = np.zeros_like(Z)

    yout[0:iout] = tmp

    return yout
 
