"""
Unstructured grid plotting module
"""
from .hybridgrid import HybridGrid, circumcenter
from .gridsearch import GridSearch
from sfoda.utils.barycentric import BarycentricInterp

import numpy as np
from matplotlib.collections import PolyCollection, LineCollection
import matplotlib.pyplot as plt
from matplotlib import tri
from scipy import spatial

import pdb


class Plot(HybridGrid):
    """
    Plotting module

    Methods:
        - plotmesh: 2D mesh outline plot
        - plotcelldata: Face plot of cell (face) centered data
        - plotedgedata: 
    """
    _xlims = None
    _ylims = None
    _xy = None
    clim = None
    def __init__(self, xp, yp, cells, **kwargs):

        HybridGrid.__init__(self, xp, yp, cells, lightmode=True, **kwargs)
        # Calculate the area
        self.Ac = self.calc_area()


    ############
    # Main user routines
    ############
    def plotmesh(self,ax=None,facecolors='none',linewidths=0.2,**kwargs):
        """
        Plots the outline of the grid mesh
        """
        fig = plt.gcf()
        if ax==None:
            ax = fig.gca()
    
        xlim=self.get_xlims()
        ylim=self.get_ylims()
        collection = PolyCollection(self.xy(), facecolors=facecolors,\
            linewidths=linewidths, **kwargs)
        
        ax.add_collection(collection)
    
        ax.set_aspect('equal')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        return ax, collection

    def plotcelldata(self, z, xlims=None, ylims=None, colorbar=True,\
            vmin=None, vmax=None, **kwargs):
        """
        Plot cell centered data
        """
        ax=plt.gca()
        fig = plt.gcf()
        # Find the colorbar limits if unspecified
        #if self.clim is None:
        #    self.clim = [z.min(),z.max()]
        # Set the xy limits
        if xlims is None or ylims is None:
            xlims=self.get_xlims()
            ylims=self.get_ylims()
        
        collection = PolyCollection(self.xy(),closed=False,**kwargs)
        collection.set_array(z)
        collection.set_clim(vmin=vmin,vmax=vmax)
        collection.set_edgecolors(collection.to_rgba(z))    
        
        ax.add_collection(collection)

        ax.set_aspect('equal')
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)

        axcb=None
        if colorbar:
            axcb = fig.colorbar(collection)

    
        return fig, ax, collection, axcb

    def plotedgedata(self, z, xlims=None,ylims=None, colorbar=True, **kwargs):
        """
          Plot the unstructured grid edge data
        """
        ax=plt.gca()
        fig = plt.gcf()

        assert(z.shape[0] == self.Ne),\
            ' length of z scalar vector not equal to number of edges, Ne.'
        
        # Find the colorbar limits if unspecified
        if self.clim is None:
            self.clim = [z.min(),z.max()]
        # Set the xy limits
        if xlims is None or ylims is None:
            xlims=self.xlims()
            ylims=self.ylims()
        
        xylines = [self.xp[self.edges],self.yp[self.edges]]

        # Create the inputs needed by line collection
        Ne = xylines[0].shape[0]

        # Put this into the format needed by LineCollection
        linesc = [list(zip(xylines[0][ii,:],xylines[1][ii,:])) for ii in range(Ne)]

        collection = LineCollection(linesc,array=z,**kwargs)

        collection.set_clim(vmin=self.clim[0],vmax=self.clim[1])
        
        ax.add_collection(collection)

        ax.set_aspect('equal')
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)

        axcb=None
        if colorbar:
            axcb = fig.colorbar(collection)
        
        return fig, ax, collection, axcb

    def contourf(self, zdata, clevs=20, \
        xlims=None, ylims=None, colorbar=True, filled=True,\
        **kwargs):
        """
        Filled contour plot
        """
        self.build_tri_voronoi()
       
        fig = plt.gcf()
        ax = fig.gca()

        # data must be on nodes for tricontourf
        #zdata = self.cell2node(z)
        #zdata[np.isnan(zdata)] = 0.

        if isinstance(clevs,int): # is integer
            V = np.linspace(self.clim[0],self.clim[1],clevs)
        else:
            V = clevs
            
        if filled:
            camp = ax.tricontourf(self._triv, zdata, V, **kwargs)
        else:
            camp = ax.tricontour(self._triv, zdata, V, **kwargs)
                
        ax.set_aspect('equal')
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)
        
        if filled and colorbar:
            axcb = fig.colorbar(camp)
        else:
            axcb = None
        
        return fig, ax, camp, axcb

    def calc_grad(self, phi, k=0):
        """
        Computes the horizontal gradient of variable, phi, along layer, k.
        
        Uses divergence (Green's) theorem to compute the gradient. This can be noisy 
        for triangular control volumes.
        
        Returns: d(phi)/dx, d(phi)/dy
        
        Based on MATLAB code sungradient.m
        """
        if self.face is None:
            self.calc_all_properties()

        def _GradientAtFace(phi,jj,k):
            
            grad1 = self.grad[:,0]
            grad2 = self.grad[:,1]
            nc1 = grad1[jj]
            nc2 = grad2[jj]
                    
            # check for edges (use logical indexing)
            ind1 = nc1==-1
            nc1[ind1]=nc2[ind1]
            ind2 = nc2==-1
            nc2[ind2]=nc1[ind2]
            
            # check depths (walls)
            indk = (k>=self.Nk[nc1]) | (k>=self.Nk[nc2])
            ind3 = (indk) &  (self.Nk[nc2]>self.Nk[nc1])
            nc1[ind3]=nc2[ind3]
            ind4 = (indk) & (self.Nk[nc1]>self.Nk[nc2])
            nc2[ind4]=nc1[ind4]
            
            # Calculate gradient across face            
            return (phi[nc1]-phi[nc2]) / self.dg[jj]
            
        ne = self.face #edge-indices
        mask = ne == self._FillValue
        ne[mask]=0
        
        Gn_phi = _GradientAtFace(phi,ne,k)
        Gn_phi[mask]=0
        
        Gx_phi = Gn_phi * self.n1[ne] * self.DEF * self.df[ne]
        Gy_phi = Gn_phi * self.n2[ne] * self.DEF * self.df[ne]
        dX = np.sum(Gx_phi,axis=1)/self.Ac
        dY = np.sum(Gy_phi,axis=1)/self.Ac

        return dX, dY
 
    ##########
    # Private routines
    ##########
    def get_xlims(self):
        if self._xlims is None:
            self._xlims = [self.xp.min(),self.xp.max()]
        return self._xlims

    def get_ylims(self):
        if self._ylims is None:
            self._ylims = [self.yp.min(),self.yp.max()]
        return self._ylims

    def xy(self):
        """ 
        Returns a list of Nx2 vectors containing the grid cell node coordinates
            
        Used by spatial ploting routines 
        """
        if self._xy is None:
            eptr, eidx = self.to_metismesh() 
            
            self._xy = [ np.asarray([ self.xp[ eidx[eptr[ii]:eptr[ii+1]] ],\
                self.yp[eidx[eptr[ii]:eptr[ii+1] ] ] ]).T\
                        for ii in range(eptr.shape[0]-1)]


        return self._xy

    def build_tri(self):
        """
        Create a matplotlib triangulation object from the grid
        
        Used primarily to contour the data (but also for interpolation)      
        """
        if '_tri' not in self.__dict__:
            maxfaces = self.nfaces.max()
            if maxfaces==3:
                self._tri =tri.Triangulation(self.xp,self.yp,self.cells)

            else:

                # Need to compute a delaunay triangulation for mixed meshes
                pts = np.vstack([self.xp,self.yp]).T
                D = spatial.Delaunay(pts)
                self._tri =tri.Triangulation(self.xp,self.yp,D.simplices)

                # Compute a mask by removing triangles outside of the polygon
                xy = D.points
                cells=D.simplices
                xv,yv = circumcenter(xy[cells[:,0],0],\
                    xy[cells[:,0],1],\
                    xy[cells[:,1],0],\
                    xy[cells[:,1],1],\
                    xy[cells[:,2],0],\
                    xy[cells[:,2],1])
                mask = self.find_cell(xv,yv)
                self._tri.set_mask(mask==-1)

    def build_tri_voronoi(self):
        """
        Create a matplotlib triangulation object from the grid voronoi points
        
        Used primarily to for interpolation      
        """
        if '_triv' not in self.__dict__:
            # Need to compute a delaunay triangulation for mixed meshes
            if self.VERBOSE:
                print('Computing delaunay triangulation and computing mask...')
            pts = np.vstack([self.xv,self.yv]).T
            D = spatial.Delaunay(pts)
            self._triv =tri.Triangulation(self.xv,self.yv,D.simplices)

            # Compute a mask by removing triangles outside of the polygon
            xy = D.points
            cells=D.simplices
            xv,yv = circumcenter(xy[cells[:,0],0],xy[cells[:,0],1],\
                xy[cells[:,1],0],xy[cells[:,1],1],xy[cells[:,2],0],xy[cells[:,2],1])
            mask = self.find_cell(xv,yv)
            self._triv.set_mask(mask==-1)

    def find_cell(self,x,y):
        """
        Return the cell index that x and y lie inside of 

        return -1 for out of bounds
        """
        if '_tsearch' not in self.__dict__:
            self._tsearch=GridSearch(self.xp,self.yp,self.cells,nfaces=self.nfaces,\
                edges=self.edges,mark=self.mark,grad=self.grad,neigh=self.neigh,\
                xv=self.xv,yv=self.yv)
        
        return self._tsearch(x,y)

    def interpolate(self, z, x, y,  kind='linear'):
        """
        Interpolate data in 'z' on current grid onto points 'x' and 'y'

        kind = 'linear', 'cubic', 'barycentric' (default), 'nearest'
        """
        if kind == 'linear':
            self.build_tri_voronoi()
            F = tri.LinearTriInterpolator(self._triv, z)
            return F(x, y)

        elif kind == 'cubic':
            self.build_tri_voronoi()
            F = tri.CubicTriInterpolator(self._triv, z)
            return F(x, y)
        
        elif kind == 'barycentric':
            self.build_tri_voronoi()
            xyin = np.vstack([self.xv,self.yv]).T
            xyout = np.vstack([x,y]).T
            F = BarycentricInterp(xyin, xyout)
            return F(z)

        elif kind == 'nearest':
            if isinstance(x, float):
                dist, idx = self.find_nearest([x,y])
                return z[idx]
            elif isinstance(x, np.ndarray):
                sz = x.shape
                xy = np.array([x.ravel(),y.ravel()]).T
                dist, idx = self.find_nearest(xy)
                return z[idx].reshape(sz)

    def cell2node(self,cell_scalar):
        """
        Map a cell-based scalar onto a node
        
        This is calculated via a mean of the cells connected to a node(point)
        """
        # Area weighted interpolation
        node_scalar = [np.sum(cell_scalar[list(self.pnt2cells(ii))] * \
                self.Ac[list(self.pnt2cells(ii))]) / \
                np.sum( self.Ac[list(self.pnt2cells(ii))])\
                    for ii in range(self.Np)]
        return np.array(node_scalar)

    def find_nearest(self, xy, NNear=1):
        """
        Returns the grid indices of the closest points to the nx2 array xy
        
        Uses the scipy KDTree routine
        """
        
        if '_kd' not in self.__dict__:
            self._kd = spatial.cKDTree(np.vstack((self.xv,self.yv)).T)
    
        # Perform query on all of the points in the grid
        dist,ind=self._kd.query(xy,k=NNear)
        
        return dist, ind


