"""
Isopycnal slicing 

Borrowed from Rob Hetland's Octant package:

        github.com/hetland/octant
"""

import numpy as np

import pdb

def isoslice(var, prop, isoval=0, axis=0, masking=True):
    """
    result = isoslice(variable,property[,isoval=0])
    
    result is a a projection of variable at property == isoval in the first
    nonsingleton dimension.  In the case when there is more than one zero
    crossing, the results are averaged.
    
    EXAMPLE:
    Assume two three dimensional variable, s (salt), z (depth), and
    u (velicity), all on the same 3D grid.  x and y are the horizontal 
    positions, with the same horizontal dimensions as z (the 3D depth 
    field).  Here, assume the vertical dimension, that will be projected,
    is the first.  
    
    s_at_m5  = isoslice(s,z,-5);        # s at z == -5
    h_at_s30 = isoslice(z,s,30);       # z at s == 30
    u_at_s30 = isoslice(u,s,30);       # u at s == 30
    """
    if var.ndim<2:
        raise ValueError('variable must have at least two dimensions')
    if not prop.shape == var.shape:
        raise ValueError('dimension of var and prop must be identical')

    var = var.swapaxes(0, axis)
    prop = prop.swapaxes(0, axis)

    prop= prop - isoval
    sz = np.shape(var)
    var = var.reshape(sz[0],-1)
    prop = prop.reshape(sz[0],-1)
    #find zero-crossings (zc == 1)
    zc =  np.where( (prop[:-1,:]*prop[1:,:])<=0.0 ,1.0, 0.0)
    varl = var[:-1,:]*zc
    varh = var[1:,:]*zc
    propl = prop[:-1,:]*zc
    proph = prop[1:,:]*zc
    result = varl - propl*(varh-varl)/(proph-propl)
    result = np.where(zc==1., result, 0.)
    szc = zc.sum(axis=0)
    szc = np.where(szc==0., 1, szc)
    result = result.sum(axis=0)/szc
    if masking:
        result = np.ma.masked_where(zc.sum(axis=0)==0, result)
        if np.all(result.mask):
            raise Warning('property==%f out of range (%f, %f)' % \
                           (isoval, (prop+isoval).min(), (prop+isoval).max()))
    result = result.reshape(sz[1:])
    return(result)


