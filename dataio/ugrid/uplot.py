"""
Unstructured grid plotting module
"""
from .hybridgrid import HybridGrid

import numpy as np
from matplotlib.collections import PolyCollection, LineCollection
import matplotlib.pyplot as plt

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
    
        xlim=self.xlims()
        ylim=self.ylims()
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
            xlims=self.xlims()
            ylims=self.ylims()
        
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

 
    ##########
    # Private routines
    ##########
    def xlims(self):
        if self._xlims is None:
            self._xlims = [self.xp.min(),self.xp.max()]
        return self._xlims

    def ylims(self):
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


