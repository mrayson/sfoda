"""
Numpy-like array operations
"""

import numpy as np
import pdb

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

    yout = np.zeros((nx,))

    yout[0:iout] = tmp

    return yout
 
def depthint(y, z, ztop=None, zbed=None, cumulative=False, axis=0):
    """
    Integrate the variable "y" along its first dimension

    Inputs (optional):
        
        ztop : (scalar) set to top of water column if different from z[-1]
        zbed : (scalar) set to bottom of water column if different from z[0]
    """
    Nz = z.shape[0]

    # Reshape the y variable
    y = y.swapaxes(0, axis)

    assert y.shape[0] == Nz

    # Calculate the vertical grid spacing
    zmid = np.zeros((Nz+1,))
    zmid[1:-1] = 0.5*(z[1:]+z[0:-1])
    if ztop is None:
        zmid[0] = z[0]
    else:
        zmid[0] = ztop

    if zbed is None:
        zmid[-1] = z[-1]
    else:
        zmid[-1] = zbed

    dz = np.abs(zmid[1:] - zmid[0:-1])

    # Perform the integration
    if cumulative:
        return np.cumsum(y*dz[:,np.newaxis], axis=0).swapaxes(axis,0), dz
    else:
        return np.sum(y*dz[:,np.newaxis], axis=0).swapaxes(axis,0), dz

def depthavg(y, z, ztop=None, zbed=None, axis=0):
    """
    Depth average the variable along the first dimension

    Inputs (optional):
        ztop : (scalar) set to top of water column if different from z[-1]
        zbed : (scalar) set to bottom of water column if different from z[0]
    """

    y_dz, dz = depthint(y, z, ztop=ztop, zbed=zbed, axis=axis)
    H = dz.sum()

    return y_dz/H

def grad_z(y, z, axis=0):
    """
    Compute the vertical gradient

    "z" can be an array same size as y, or vector along the first axis of "y"

    Takes the derivative along the dimension specified by axis(=0)
    """
    Nz = z.shape[0]

    # Reshape the y variable
    y = y.swapaxes(0, axis)
    assert y.shape[0] == Nz

    z = z.swapaxes(0, axis)
    assert z.shape == (Nz,) or z.shape == y.shape

    dy_dz = np.zeros_like(y)
    
    # Second-order accurate for mid-points
    ymid = 0.5*(y[1:,...]+y[0:-1,...])

    zmid = 0.5*(z[1:,...]+z[0:-1,...])

    dzmid  = zmid[1:,...] - zmid[0:-1,...] 

    dy_dz[1:-1, ...] = (ymid[1:,...] - ymid[0:-1,...])/\
            dzmid[:,...]

    # First-order accurate for top and bottom cells
    dy_dz[0,...] = (y[1,...] - y[0,...])/dzmid[0,...]
    dy_dz[-1,...] = (y[-1,...] - y[-2,...])/dzmid[-1,...]

    return dy_dz.swapaxes(axis, 0)

