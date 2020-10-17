"""
Barycentric interpolation on scattered data

"""

from scipy.spatial import Delaunay
import numpy as np

class BarycentricInterp(object):
    """
    Triangular barycentric interpolation class
    """
    def __init__(self, xyin, xyout, **kwargs):
        self.__dict__.update(kwargs)
        
        # Triangulate the input 
        tri = Delaunay(xyin)
        
        # Find the simplex (i.e. the cell or polygon number) where each point lies
        s = tri.find_simplex(xyout)
        
        # Compute the barycentric coordinates (these are the weights)
        X = tri.transform[s,:2]
        Y = xyout - tri.transform[s,2]
        b = np.einsum('ijk,ik->ij', X, Y)
        self.weights = np.c_[b, 1 - b.sum(axis=1)]
        
        # These are the vertices of the output points
        self.verts = tri.simplices[s]
    
    def __call__(self, z):
        
        return (z[self.verts]*self.weights).sum(axis=1)
