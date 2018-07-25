# cython: profile=True

# Cython utilities to speed up certain tasks

import numpy as np
from numpy cimport ndarray, int32_t, int64_t, double_t
cimport numpy as np
cimport cython


@cython.boundscheck(False)
cpdef create_pnt2cells(ndarray[int32_t,ndim=2] cells,
    ndarray[int64_t,ndim=1] nfaces):
    """
    build hash table for point->cell lookup
    """

    cdef int ii, jj, cc
    cdef int nc = nfaces.shape[0]

    cdef dict _pnt2cells = {}

    for ii in range(nc):
        for jj in range(nfaces[ii]):
            cc = cells[ii,jj]
            if not _pnt2cells.has_key(cc):
                _pnt2cells[cc] = set()
            _pnt2cells[cc].add(ii)
    
    return _pnt2cells

@cython.boundscheck(False)
cpdef create_pnt2edges(ndarray[int64_t, ndim=2] edges,
	ndarray[int64_t, ndim=1] mark,
	int DELETED_EDGE):
    """
    Creates the node->edge lookup table
    """
    cdef int e, p
    cdef int ne = edges.shape[0]
    cdef dict p2e = {}

    for e in range(ne):
	# skip deleted edg
        if edges[e,2] == DELETED_EDGE:
            continue
        if mark[e] == DELETED_EDGE:
            continue

        for p in edges[e,:2]:
            if not p2e.has_key(p):
                p2e[p] = []
            p2e[p].append(e)
	
    #for p in p2e.keys():
    #    p2e[p] = array(p2e[p])

    return p2e

@cython.boundscheck(False)
cpdef make_neigh_from_cells(ndarray[int32_t,ndim=2] cells,
    ndarray[int64_t,ndim=1] nfaces,
	dict pnt2cells):
    """
    Find the neighbouring cells
    """

    cdef int i, j
    cdef int nc = cells.shape[0]
    cdef int maxfaces = cells.shape[1]
    cdef ndarray[int64_t, ndim=2] neigh = np.zeros((nc, maxfaces),np.int)
    pnt2cells = create_pnt2cells(cells, nfaces)

    for i in range(nc):
        # find the neighbors:
        # the first neighbor: need another cell that has
        # both self.cells[i,0] and self.cells[i,1] in its
        # list.
        my_set = set([i])
        n = nfaces[i] * [-1]
        for j in range(nfaces[i]):
            adj1 = pnt2cells[cells[i,j]]
            adj2 = pnt2cells[cells[i,(j+1)%nfaces[i]]]
            neighbor = adj1.intersection(adj2).difference(my_set)
            if len(neighbor) == 1:
                n[j] = neighbor.pop()
        
        neigh[i,0:nfaces[i]] = n
    
    return neigh

@cython.boundscheck(False)
@cython.wraparound(False)
def make_edges_from_cells(ndarray[int32_t,ndim=2] cells,
	ndarray[int64_t,ndim=1] nfaces,
	dict pnt2cells):


    cdef int i, j, n, pnt_a, pnt_b
    cdef int nc = nfaces.shape[0]
    cdef int maxfaces = nfaces.shape[1]


    # iterate over cells, and for each cell, if it's index
    # is smaller than a neighbor or if no neighbor exists,
    # write an edge record
    cdef list edges = list()
    cdef set my_set = set()
    cdef set adj1 = set()
    cdef set adj2 = set()
    cdef set neighbor = set()
    default_marker = 0

    # this will get built on demand later.
    #self._pnt2edges = None
    
    for i in range(nc):
        # find the neighbors:
        # the first neighbor: need another cell that has
        # both self.cells[i,0] and self.cells[i,1] in its
        # list.
        if i%1000==0:
            print(i,nc)

        my_set = set([i])
        #n = [-1,-1,-1]
        #n = nfaces[i] * [-1]
        #for j in 0,1,2:
        #for j in range(self.MAXFACES):
        for j in range(nfaces[i]):
            pnt_a = cells[i,j]
            pnt_b = cells[i,(j+1)%nfaces[i]]
            
                
            adj1 = pnt2cells[pnt_a] # cells that use pnt_a
            adj2 = pnt2cells[pnt_b] # cells that use pnt_b

            # the intersection is us and our neighbor
            #  so difference out ourselves...
            #neighbor = adj1.intersection(adj2).difference(my_set)
            neighbor = (adj1 & adj2) - my_set
            ##adj3 = set.intersection(adj1, adj2)
            ##neighbor = set.difference(adj3, my_set)

            # and maybe we ge a neighbor, maybe not (we're a boundary)
            if len(neighbor) == 1:
                n = neighbor.pop()
            else:
                n = -1

            # Use a loop instead
            #n = -1
            #for a1 in adj1:
            #    for a2 in adj2:
            #        if a1 == a2: # intersection
            #            for m1 in my_set:
            #                if a1 != m1:
            #                    n = a1 # difference
                
            if n==-1 or i<n:
                # we get to add the edge:
                edges.append((pnt_a,
                              pnt_b,
                              default_marker,
                              i,n))

    #self.edge = np.array(edges,np.int32)
    Ne = len(edges)
    #edges = np.asarray(edges)
    #alledges = np.array([edges[ii,0:2] for ii in range(Ne)])
    #mark = np.array([edges[ii,2] for ii in range(Ne)])
    #grad = np.array([edges[ii,3:5] for ii in range(Ne)])

    alledges = np.array([edge[0:2] for edge in edges])
    mark = np.array([edge[2] for edge in edges])
    grad = np.array([edge[3:5] for edge in edges])

        
    #return edges[:,0:2], edges[:,2], edges[3:5]
    return alledges, mark, grad

@cython.boundscheck(False)
cpdef ensure_ccw(ndarray[int32_t,ndim=2] cells,
	ndarray[int64_t,ndim=1] nfaces,
	ndarray[double_t,ndim=1] Ac):
    """
    Ensure that the nodes are rotated counter-clockwise

    Modifies the cells and Ac arrays
    """
    cdef int i
    cdef int nc = nfaces.shape[0]

    for i in range(nc):
        if Ac[i] < 0:
            cells[i,0:nfaces[i]] = cells[i,0:nfaces[i]][::-1] # reverse order
            Ac[i] *= -1.
        else:
	    # Do nothing
            continue	
            #cells[i,0:nfaces[i]] =  cells[i,0:nfaces[i]]

@cython.boundscheck(False)
@cython.wraparound(False)
cdef cell2edges(int cell_i,
	ndarray[int32_t,ndim=2] cells,
	ndarray[int64_t,ndim=1] nfaces,
	dict pnt2edges):
    """
    Finds the edge indices for individual grid cells
    """    
    cdef:
        int i
        int nf = nfaces[cell_i]

        ndarray[int32_t,ndim=1] edges = np.zeros( (nf,), np.int32)
        np.ndarray[np.int32_t,ndim=1] pnts = np.zeros( (nf,), np.int32)

    #if cells[cell_i,0] == -1:
    #    raise "cell %i has been deleted"%cell_i
    
    # return indices to the edges for this cell:
    pnts = cells[cell_i,0:nf] # the vertices

    for i in range(nf):
        edges[i] = find_edge( pnts[(i)%nf], pnts[(i+1)%nf], pnt2edges )

    #edges = [ find_edge( pnts[(i)%nf], pnts[(i+1)%nf], pnt2edges ) \
    #   for i in range(nf) ]
    #
    #edges = [ self.find_edge( (pnts[(i+1)%nf], pnts[(i+2)%nf]) ) for i in range(nf) ]
    # This ordering ensures that the 'face' and 'neigh' arrays correspond
    #edges = [ self.find_edge( (pnts[(i+1)%nf], pnts[(i+2)%nf]) ) for i in range(-1,nf-1) ]

    return edges

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef cell_edge_map(ndarray[int32_t,ndim=2] cells,
	ndarray[int64_t,ndim=1] nfaces,
	dict pnt2edges):
    """ cell2edges for the whole grid
    return an integer valued [Nc,3] array, where [i,k] is the edge index
    opposite point self.cells[i,k]

    """
    cdef:
        int nc = cells.shape[0]
        int maxfaces = cells.shape[1]
        int i
    
    cdef ndarray[int32_t, ndim=2] cem = 999999*np.ones( (nc, maxfaces), np.int32)

    for i in range(nc):
        cem[i,0:nfaces[i]] = cell2edges(i, cells, nfaces, pnt2edges)

    return cem


@cython.boundscheck(False)
@cython.wraparound(False)
cdef find_edge(int32_t node0, int32_t node1, dict pnt2edges):

    #cdef np.int32_t el0, el1, e
    cdef:
        #int ne0 = pnt2edges[node0].shape[0]
        #int ne1 = pnt2edges[node1].shape[0]
        int e
        #np.ndarray[np.int64_t,ndim=1] el0 = np.zeros( (ne0,), np.int64)
        #np.ndarray[np.int64_t,ndim=1] el1 = np.zeros( (ne1,), np.int64)
	

    #print node0, type(node0)
    #print pnt2edges[node0]

    el0 = pnt2edges[node0]
    el1 = pnt2edges[node1]

    for e in el0:
        if e in el1:
            return e
    return -1


