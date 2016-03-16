"""
Cython utility functions for geometric searching
"""
import numpy as np
cimport numpy as np
cimport cython

from cython.parallel import prange


#cdef class Point:
#
#    cdef double x
#    cdef double y
#
#    def __init__(self, double x, double y):
#        self.x = x
#        self.y = y
cdef struct _Point:
    double x
    double y

ctypedef _Point Point

cdef inline int ccw(Point A, Point B, Point C) nogil:
    return (C.y-A.y)*(B.x-A.x) > (B.y-A.y)*(C.x-A.x)

cdef inline int intersect(Point A, Point B, Point C, Point D) nogil:
    return (ccw(A,C,D) != ccw(B,C,D)) & (ccw(A,B,C) != ccw(A,B,D))

cdef inline int check_edge_crossing(\
        double xnew, double ynew,\
        double xold, double yold,\
        double xp1, double yp1,\
        double xp2, double yp2 ) nogil:
    """
    Check to see if a particle has crossed an edge of a cell
    """
    cdef Point p1, p2, A, B

    #p1 = Point(xold,yold)
    #p2 = Point(xnew,ynew)
    p1.x = xold
    p1.y = yold
    p2.x = xnew
    p2.y = ynew
    
    #A = Point(xp, yp)
    #B = Point(xp, yp)
    A.x = xp1
    A.y = yp1
    B.x = xp2
    B.y = yp2
    
    if intersect(p1,p2,A,B):
        return 1
    else:
        return 0
 
@cython.boundscheck(False)
cpdef check_cell_crossing(\
        np.ndarray[np.int32_t,ndim=1] cells_i,\
        np.ndarray[np.double_t,ndim=1] xnew,\
        np.ndarray[np.double_t,ndim=1] ynew,\
        np.ndarray[np.double_t,ndim=1] xold,\
        np.ndarray[np.double_t,ndim=1] yold,\
        np.ndarray[np.double_t,ndim=1] xp,\
        np.ndarray[np.double_t,ndim=1] yp,\
        np.ndarray[np.int64_t,ndim=1] nfaces,\
        np.ndarray[np.int32_t,ndim=2] cells):
    """
    Cythonized cell crossing code
    """
    cdef:
        int nc = xnew.shape[0]
        int nf, ii, nn, cell_i
        int pt1, pt2

        np.ndarray[np.int16_t, ndim=1] changedcell = np.zeros( (nc,), np.int16)
        np.ndarray[np.int64_t, ndim=1] neigh = np.zeros( (nc,), np.int)

    neigh[:] = -1

    # Loop through all cells
    #for ii in range(nc):
    for ii in prange(nc, nogil=True):
        cell_i = cells_i[ii]
        nf = nfaces[cell_i]
        for nn in range(nf):
            pt1 = nn
            #pt2 = mod(nn+1,nf-1)
            #pt2 = (nn+1)%(nf-1) # hopefully this calls the c modulus function
            pt2 = pt1+1
            if pt2 == nf:
                pt2 = 0

            if check_edge_crossing(\
                xnew[ii], ynew[ii], xold[ii], yold[ii],\
                xp[cells[cell_i, pt1]], yp[cells[cell_i, pt1]],\
                xp[cells[cell_i, pt2]], yp[cells[cell_i, pt2]],):

                changedcell[ii] = 1
                neigh[ii] = nn

                break

    return changedcell, neigh
             


    '''
    # Loop through each face on the edge
    for ii in range(self.maxfaces):
        # find the index of the first and second edge
        ind = ii>=self.nfaces[cell_i]-1
        pt1 = ii*np.ones(xold.shape,np.int)
        pt1[ind] = self.nfaces[cell_i][ind]-1

        pt2 = pt1+1
        pt2[ind]=0

        # create two points along the edge
        A = Point(self.xp[self.cells[cell_i,pt1]],self.yp[self.cells[cell_i,pt1]])
        B = Point(self.xp[self.cells[cell_i,pt2]],self.yp[self.cells[cell_i,pt2]])
    
        # Check if the line between the two particles crosss the edge
        hasleftcell = intersectvec(p1,p2,A,B)
        changedcell[hasleftcell] = True
        neigh[hasleftcell] = pt1[hasleftcell]
    '''
    
    
 

    
