"""
Barycentric interpolation on scattered data

"""

from scipy.spatial import Delaunay
import numpy as np

class BarycentricInterp(object):
    """
    Triangular barycentric interpolation class
    """
    
    tri = None
	
    def __init__(self, xyin, xyout, **kwargs):
        self.__dict__.update(kwargs)
        
        # Triangulate the input 
        if self.tri is None:
            self.tri = Delaunay(xyin)

        self.weights, self.verts, self.simplex = self._calc_weights(xyout)
        
   
    def __call__(self, z, xyout=None, mask=False):
        """
        Perform the interpolation
        """
        if xyout is not None:
            self.weights, self.verts = self._calc_weights(xyout)
            
        zout = (z[self.verts]*self.weights).sum(axis=1)

        if mask:
            zout[self.simplex==-1] = np.nan

        return zout

    def _calc_weights(self, xyout):
        # Find the simplex (i.e. the cell or polygon number) where each point lies
        s = self.tri.find_simplex(xyout)
        
        # Compute the barycentric coordinates (these are the weights)
        X = self.tri.transform[s,:2]
        Y = xyout - self.tri.transform[s,2]
        b = np.einsum('ijk,ik->ij', X, Y)
        weights = np.c_[b, 1 - b.sum(axis=1)]
        
        # These are the vertices of the output points
        verts = self.tri.simplices[s]
        
        return weights, verts, s
 
