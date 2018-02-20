# -*- coding: utf-8 -*-
"""
Tools for filtering unstructured grid data

Created on Tue Dec 18 13:19:39 2012

@author: mrayson
"""

import numpy as np
from scipy import spatial, sparse
#from numba import jit

import pdb

class ufilter(object):
    """
    Unstructured grid filter class
    """
    
    c = 4. # uses point c x p for filter
    filtertype = 'gaussian' # 'gaussian' or 'lanczos'
    kmax = 50 # Maximum number of points to use in filter matrix
    vectorized=True
    
    def __init__(self,X,delta_f,**kwargs):
        
        self.X = X
        self.delta_f = delta_f
        
        self.__dict__.update(kwargs)
        
        s = np.shape(self.X)
        self.n = s[0] # Number of points
        if len(s)>1:
            self.m = s[1] # Number of dimensions
        else:
            self.m = 1
        
        self.GetP()
        
        #if self.vectorized:
        #    self.BuildFilterMatrix2()
        #else:
        #    self.BuildFilterMatrix()
        self.build_filter_matrix()
        
    def __call__(self,y):
        """
        Performs the filtering operation on data in vector y
        """
        #print np.argwhere(self.G.sum(axis=1) == 0.0)
        return self.G.dot(y)
        
    #@jit
    def build_filter_matrix(self):
        """
        Builds a sparse matrix, G, used to filter data in vector, y, via:
            y_filt = G x y
            
        """
        # Compute the spatial tree
        print 'Building KD-Tree...'
        kd = spatial.cKDTree(self.X)
        eps=1e-6

        # Initialise the sparse matrix
        print 'Creating an %d x %d sparse matrix...'%(self.n, self.n)
        #self.G = sparse.lil_matrix((self.n,self.n))
        weights = np.zeros((self.n, self.kmax))
        rows = np.zeros((self.n, self.kmax), dtype=np.int)
        cols = np.zeros((self.n, self.kmax), dtype=np.int)
        mask = np.zeros((self.n, self.kmax), dtype=np.bool)
        
        printstep = 5 
        printstep0 = 0
        ii=0
        for nn in range(self.n):
            ii+=1
            perccomplete = float(nn)/float(self.n)*100.0
            if perccomplete > printstep0:
                print '%d %% complete...'%(int(perccomplete))
                printstep0+=printstep
            #print nn
            # Find all of the points within c * p distance from point
            #dx, i = kd.query(self.X[nn,:]+eps,k=self.n+1,distance_upper_bound=self.c*self.p)
            dx, i = kd.query(self.X[nn,:]+eps, k=self.kmax, distance_upper_bound=self.c*self.p)
            
            # Calculate the filter weights on these points
            ind = dx != np.inf
            Gtmp = self.gaussian(dx[ind])
            
            # Insert the points into the sparse matrix
            I = i[ind]
            #self.G[nn,I] = Gtmp

            nnear = ind.sum()
            weights[nn,0:nnear] = Gtmp
            cols[nn,0:nnear] = I
            rows[nn,0:nnear] = nn
            mask[nn,0:nnear] = True

        # build a sparse matrix here
        print 'Building the CSR matrix...'
        self.G = sparse.csr_matrix((weights[mask], (rows[mask],cols[mask])),\
                shape=(self.n, self.n))
        print 'Done'


    def BuildFilterMatrix(self):
        """
        Builds a sparse matrix, G, used to filter data in vector, y, via:
            y_filt = G x y
            
        """
        # Compute the spatial tree
        print 'Building KD-Tree...'
        kd = spatial.cKDTree(self.X)
        eps=1e-6
        # Initialise the sparse matrix
        print 'Creating an %d x %d sparse matrix...'%(self.n, self.n)
        self.G = sparse.lil_matrix((self.n,self.n))
        
        printstep = 5 
        printstep0 = 0
        ii=0
        for nn in range(self.n):
            ii+=1
            perccomplete = float(nn)/float(self.n)*100.0
            if perccomplete > printstep0:
                print '%d %% complete...'%(int(perccomplete))
                printstep0+=printstep
            #print nn
            # Find all of the points within c * p distance from point
            dx, i = kd.query(self.X[nn,:]+eps,k=self.n+1,distance_upper_bound=self.c*self.p)
            
            # Calculate the filter weights on these points
            ind = dx != np.inf
            Gtmp = self.gaussian(dx[ind])
            
            # Insert the points into the sparse matrix
            I = i[ind]
            self.G[nn,I] = Gtmp
            #self.G[nn,I] = dx[ind] # testing

    def BuildFilterMatrix2(self):
        """
        Builds a sparse matrix, G, used to filter data in vector, y, via:
            y_filt = G x y
        
        Vectorized version of the above
        """
        # Compute the spatial tree
        print 'Building KD-Tree...'
        kd = spatial.cKDTree(self.X)
        eps=1e-6

        # Initialise the sparse matrix
        print 'Creating an %d x %d sparse matrix...'%(self.n, self.n)
        self.G = sparse.lil_matrix((self.n,self.n))
        
        # Find all of the points within c * p distance from point

        print 'Querying matrix...'
        dx, i = kd.query(self.X+eps,k=self.kmax,distance_upper_bound=self.c*self.p)
        
        ind = np.isinf(dx)
        # Calculate the filter weights
        if self.filtertype=='gaussian':
            Gtmp = self.gaussian(dx)
        elif self.filtertype=='lanczos':
            Gtmp = self.Lanczos(dx)
        
        # Set the weighting to zero for values outside of the range
        Gtmp[ind]=0
        i[ind]=0
        
        # Normalise the filter weights
        sumG = np.sum(Gtmp,axis=1)
        Gout = [Gtmp[ii,:]/G for ii, G in enumerate(sumG)]

        
        for nn,gg in enumerate(Gout):
            self.G[nn,i[nn,:]]=gg

        # sets the first value to 1.
        self.G[0,0] = 1.
            
#        ind = [nn*self.n + i[nn,:] for nn,G in enumerate(G)]
#        ind = np.array(ind)
#        G=np.array(G)
#        pdb.set_trace()
#        self.G[ind.ravel()] = G.ravel()

#        N = [nn + i[nn,:]*0 for nn,G in enumerate(G)]
#        N = np.array(N)
#        G=np.array(G)
#        self.G[N,i]=G
         
        #[self.G[nn,i[nn,:]] for nn,G in enumerate(tmp)]
            
#        # Calculate the filter weights on these points
#        ind = dx != np.inf
#        Gtmp = self.Gaussian(dx[ind])
#        
#        # Insert the points into the sparse matrix
#        I = i[ind]
#        self.G[nn,I] = Gtmp
#        #self.G[nn,I] = dx[ind] # testing
 
    def GetP(self):
        """
        Calculate the 'p' parameter
        """
        #self.p = self.delta_f**2/40.0
        self.p = self.delta_f
        
    #def kmax(self):
    #    """
    #    Estimate the maximum number of points in the search radius
    #    """        
    #    return np.round(self.c*self.p/self.dxmin)
        
    #@jit
    def gaussian(self,dx):
        """
        Calculate the Gaussian filter weights
        """      
        #Gtmp = 1.0 / (4.0*np.pi*self.p) * np.exp(-dx/np.sqrt(self.p))/dx
        #Gtmp = Gtmp[1:] # discard closest point (self)
        if self.m == 1:
            # 1D gaussian
            coef = 1.0/(np.sqrt(2.0*np.pi)*self.p)
        elif self.m == 2:
            # 2D gaussian
            coef = 1.0/(2.0*np.pi*self.p**2)
            
        Gtmp = coef * np.exp(- dx**2 / (2*self.p**2))
        
        return Gtmp / np.sum(Gtmp)
        
    def Lanczos(self,dx):
        """
        Lanczos filter weights
        
        !!!Need to check this!!!
        """        
        a = self.p
        
        Gtmp = np.sinc(dx) * np.sinc(dx/a)
        
        return Gtmp / np.sum(Gtmp)
