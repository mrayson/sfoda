# -*- coding: utf-8 -*-
"""
    Kriging interpolation library

    Fast version - uses numba
"""
from scipy import spatial
import numpy as np
from numba import jit
from soda.utils.linsolve import linsolve

import pdb

@jit(nopython=False)
def semivariogram(D, vrange, nugget, sill):
    """ Semivariogram functions"""
    if D > vrange:
        F = sill
    else:
        tmp = D/vrange
        F = nugget + (sill-nugget)*(1.5*tmp - 0.5*tmp**3)
    return F


@jit(nopython=False)
def get_weights(Ns, C, gamma, dist, xin, yin, vrange, nugget, sill):
    """ 
    Calculates the kriging weights point by point
    """
    eps = 1e-10
    for i in range(0,Ns):
        C[i,i] = semivariogram(eps, vrange, nugget, sill)
        for j in range(i+1,Ns):
            D = np.sqrt((xin[i]-xin[j])**2+(yin[i]-yin[j])**2)
            C[i,j] = semivariogram(D+eps, vrange, nugget, sill)
            C[j,i] = C[i,j]

    C[Ns,Ns]=0
   
    # Loop through each model point and calculate the vector D
    for j in range(0,Ns):
        gamma[j]= semivariogram(dist[j]+eps, vrange, nugget, sill)

    # Solve the matrix to get the weights
    linsolve(C,gamma,Ns+1)
    #gamma = np.linalg.solve(C,gamma)

    return gamma
    #return gamma[:-1]
 

@jit(nopython=False)
def get_all_weights(Nc, Ns, W, C, gamma, \
         dist, ind, xyin, vrange, nugget, sill):

    for ii in range(Nc):
        gamma = get_weights(Ns, C, gamma,\
            dist[ii,:], xyin[ind[ii,:],0], xyin[ind[ii,:],1],\
                vrange, nugget, sill)

        W[:,ii] = gamma[:-1]

        # Reset these arrays
        gamma[:] = 0
        gamma +=1
        C[:] = 0
        C += 1

    return W

######
# Distance calculation functions
######
def get_distance_kdtree(xyin, xyout, nnear, anisofac=1.):
        # Compute the spatial tree
        kd = spatial.cKDTree(xyin)
        
        # Perform query on all of the points in the grid
        dist_old, idx = kd.query(xyout, k=nnear)

        # Compute the distance with anisotropy
        dx = xyout[:,0, np.newaxis] - xyin[idx,0]
        dy = xyout[:,1, np.newaxis] - xyin[idx,1]
    
        dist = np.abs(anisofac*dx+1j*dy)

        return dist, idx
 
def nearest(trid, xpt, ypt, nnear):

    # Find the cell idx
    cidx = trid.find_simplex([xpt,ypt])
    cidx = [cidx]
    nodes = set(trid.simplices[cidx])
    flag=True
    ii = 0
    while flag:
        ii+=1
        if ii > 20:
            
            flag=False
        # Loop through the neighbours
        for cc in cidx:
            neighs = trid.neighbors[cc]
            for ff in neighs:
                if ff == -1:
                    break
                newnodes = trid.simplices[ff]
                for nn in newnodes:
                    if nn not in nodes and nn != -1:
                        nodes.add(nn)

                        if len(nodes) == nnear:
                            flag=False
                            
            cidx = neighs

    return list(nodes)
    
def find_tri_nearest(xyin, xyout, nnear, anisofac=1.):
    trid = spatial.Delaunay(xyin)
    nx = xyout.shape[0]
    idx = np.zeros((nx,nnear), np.int32)
    for ii in range(nx):
        if ii%10000 == 0:
            print(ii, nx)
        xpt, ypt = xyout[ii,0], xyout[ii,1]
        mynodes = nearest(trid, xpt, ypt, nnear)

        idx[ii,:] = mynodes[0:nnear]
    
    # Compute the distance
    dx = xyout[:,0, np.newaxis] - xyin[idx,0]
    dy = xyout[:,1, np.newaxis] - xyin[idx,1]
    
    dist = np.abs(anisofac*dx+1j*dy)
    
    return dist, idx       

#######

class kriging(object):
    
    """ Class for kriging interpolation"""
    
    ### Properties ###
    maxdist = 1000
    NNear = 12
    
    # Variogram paramters
    varmodel = 'spherical'
    nugget = 0.1
    sill = 0.8
    vrange = 250.0

    # Anisotropy factor for distances (scalar or vector)
    anisofac = 1.
    
    verbose = True
    
    def __init__(self,XYin, XYout, distance_func=get_distance_kdtree, **kwargs):
        self.__dict__.update(kwargs)
        self.distance_func = distance_func
        
        self.XYin = XYin
        self.XYout = XYout
        
        self._build_weights()
        
    def __call__(self,Zin):
        """
        Calls the interpolation function with the scalar in Zin
        """
        self.Z = np.zeros((self.Nc,))
        for ii in range(0,self.Nc):
            self.Z[ii] = np.dot(self.W[:,ii],Zin[self.ind[ii,:]])
            
        return self.Z
                
    def _build_weights(self):
        """ Calculates the kriging weights for all of the points in the grid"""
        
        dist, self.ind = self.distance_func(self.XYin, self.XYout,\
                self.NNear, anisofac=self.anisofac)       
        self.Nc = np.size(self.ind,axis=0)
        print('%d interpolation points.'%self.Nc)

        # Initialise some array to calculate the weights
        Ns = self.NNear
    
        # Construct the LHS matrix C
        C=np.ones((Ns+1,Ns+1))
        gamma = np.ones((Ns+1,))
        self.W = np.zeros((Ns,self.Nc))

        get_all_weights(self.Nc, Ns, self.W, C, gamma, \
            dist, self.ind, self.XYin,\
            self.vrange, self.nugget, self.sill)

        ##get_weights = np.vectorize(self.get_weights)
        #W = [self.get_weights(dist[ii,:], self.XYin[self.ind[ii,:],0],\
        #        self.XYin[self.ind[ii,:],1]) for ii in range(self.Nc)]
        #self.W = np.array(W).squeeze().T

        # Print percentages

        # Now loop through and get the weights for each point
        #self.W = np.zeros((self.NNear,self.Nc))
        #p0=0
        #pstep=5
        #for ii in range(0,self.Nc):
        #    
        #    if self.verbose:
        #        pfinish = float(ii)/float(self.Nc)*100.0
        #        if  pfinish> p0:
        #            print '%3.1f %% complete...'%pfinish
        #            p0+=pstep
        #                        
        #    W = self.get_weights(dist[ii,:],\
        #        self.XYin[self.ind[ii,:],0],\
        #        self.XYin[self.ind[ii,:],1])
        #    
        #    self.W[:,ii] = W.T 
                
        
    #def get_weights(self,dist,xin,yin):
    #    
    #    """ Calculates the kriging weights point by point"""
    #    
    #    eps = 1e-10
    #    Ns = len(dist)
    #    #Ns = dist.shape[0]
    #    
    #    # Construct the LHS matrix C
    #    C=np.ones((Ns+1,Ns+1))
    #    for i in range(0,Ns):
    #        C[i,i]=0
    #        for j in range(i+1,Ns):
    #            D = np.sqrt((xin[i]-xin[j])**2+(yin[i]-yin[j])**2)
    #            C[i,j] = self.semivariogram(D+eps)
    #            C[j,i] = C[i,j]

    #    C[Ns,Ns]=0

    #    ###
    #    # Old method
    #    ###
    #    # Calculate the inverse of C 
    #    #Cinv = np.linalg.inv(C)
    #    
    #    # Loop through each model point and calculate the vector D
    #    gamma = np.ones((Ns+1,1))
    #    
    #    for j in range(0,Ns):
    #        gamma[j,0]= self.semivariogram(dist[j]+eps)

    #    # Solve the matrix to get the weights
    #    #W = np.dot(Cinv,gamma)
    #    #W = W[:-1,:]

    #    W = np.linalg.solve(C,gamma)

    #    #print np.size(gamma,axis=0),np.size(gamma,axis=1)   
    #    return W[:-1]
    #    #
    #    #return 1.0/float(Ns)*np.ones((Ns,1))
    #    
    #def semivariogram(self,D):
    #    """ Semivariogram functions"""
    #    if self.varmodel == 'spherical':
    #        if D > self.vrange:
    #            F = self.sill
    #        else:
    #            tmp = D/self.vrange
    #            F = self.nugget + (self.sill-self.nugget)*(1.5*tmp - 0.5*tmp**3)
    #    return F
  
        
        
