# -*- coding: utf-8 -*-
"""

Tools for handling SUNTANS output data

Created on Mon Sep 24 16:55:45 2012

@author: mrayson
"""

from netCDF4 import MFDataset, Dataset 
from cftime import num2pydate
import numpy as np
from datetime import datetime
import os, time, getopt, sys
from scipy import spatial
from scipy import sparse
import operator

import xarray as xray

from matplotlib import tri
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, LineCollection
import matplotlib.animation as animation

from sfoda.utils import othertime
from sfoda.utils.mynumpy import grad_z
from sfoda.utils.timeseries import timeseries
from sfoda.utils.ufilter import ufilter
from sfoda.utils.myproj import MyProj
from sfoda.ugrid.hybridgrid import HybridGrid, circumcenter
from sfoda.ugrid.gridsearch import GridSearch
from .suntans_ugrid import ugrid

import pdb

# Constants
GRAV=9.81
FILLVALUE=-999999

###############################################################        
# Dictionary with lookup table between object variable name and netcdf file
# variable name
suntans_gridvars = {'xp':'xp',\
                    'yp':'yp',\
                    'xv':'xv',\
                    'yv':'yv',\
                    'xe':'xe',\
                    'ye':'ye',\
                    'lonv':'lonv',\
                    'latv':'latv',\
                    'lonp':'lonp',\
                    'latp':'latp',\
                    'cells':'cells',\
                    'face':'face',\
                    'nfaces':'nfaces', \
                    'edges':'edges',\
                    'neigh':'neigh',\
                    'grad':'grad',\
                    #'gradf':'gradf',\
                    'mark':'mark',\
                    'normal':'normal',\
                    'mnptr':'mnptr',\
                    'eptr':'eptr',\
                    'n1':'n1',\
                    'n2':'n2',\
                    'df':'df',\
                    'dg':'dg',\
                    'def':'def',\
                    'Ac':'Ac',\
                    'dv':'dv',\
                    'dz':'dz',\
                    'z_r':'z_r',\
                    'z_w':'z_w',\
                    'Nk':'Nk',\
                    'Nke':'Nke',\
                    'time':'time'\
                    }
suntans_dimvars = {'Np':'Np',\
                   'Ne':'Ne',\
                   'Nc':'Nc',\
                   'Nkmax':'Nk',\
                   'Nk':'Nk',\
                   'maxfaces':'numsides'\
                   }

########
# Utility functions
########
def calc_z(dz):
    z_bot = np.cumsum(dz)
    z_top = np.hstack((0.0,z_bot[:-1]))
    return 0.5*(z_bot+z_top)

##########
# Classes
##########
class Grid(object):
    """ Class for handling SUNTANS grid data"""

    MAXFACES=3 # Default number of faces
    _FillValue=FILLVALUE
    gridvars = suntans_gridvars
    griddims = suntans_dimvars

    # Some grid properties
    DEF = None

    VERBOSE=True

    ###
    # Grid projection details
    projstr = None
    utmzone = None
    isnorth = None

    def __init__(self,infile ,**kwargs):
               
        self.__dict__.update(kwargs)

        if isinstance(infile,list):
            infile2=infile[0]
        else:
            infile2=infile
        
        if os.path.isdir(infile2):
            # Load ascii grid file
            self.infile = infile
            self.__loadascii()            
            
        else:
            # Load the grid fromm a netcdf file
            self.infile = infile
            self.ncfile = infile
            self.__loadGrdNc()
        

        # Find the grid limits
        self.xlims = [self.xp.min(),self.xp.max()]
        self.ylims = [self.yp.min(),self.yp.max()]
        
        # Cell polygon attribute for plotting  
        #self.cells[self.cells.mask]=0
        #self.grad[self.grad.mask]=0
        #self.face[self.face.mask]=0
        self.xy = self.cellxy(self.xp, self.yp)

        ###
        # Find the lat/lon of the cell-centres
        self.lonv, self.latv = self.to_latlon(self.xv, self.yv)
        self.lonp, self.latp = self.to_latlon(self.xp, self.yp)

        
    
    def __loadascii(self):
        """
        Load the grid variables from the ascii files: points.dat, edges.dat, cells.dat
        """
        pointdata = readTXT(self.infile+'/points.dat')
        celldata = readTXT(self.infile+'/cells.dat')
        edgedata = readTXT(self.infile+'/edges.dat')
        
        self.xp = pointdata[:,0]
        self.yp = pointdata[:,1]
        #self.dv = pointdata[:,2] # zero to start
        self.Np = len(self.xp)
        
        # Work out if cells.dat is in the quad grid format based on number of
        # columns
        self.Nc = celldata.shape[0]
        if celldata.ndim==2:
            ncols = celldata.shape[1]
            #print '!!!cells.dat has %d columns!!!'%ncols
            if ncols==8: # Old format
                self.xv = celldata[:,0]
                self.yv = celldata[:,1]
                self.cells = np.asarray(celldata[:,2:5],int)
                self.neigh = np.asarray(celldata[:,5:8])
                self.nfaces = 3*np.ones((self.Nc,),np.int)
                self.maxfaces=self.MAXFACES
            elif ncols==9: # New format
                nfaces = celldata[:,0]
                self.nfaces = np.zeros((self.Nc,),np.int)
                self.nfaces[:] = nfaces
                self.xv = celldata[:,1]
                self.yv = celldata[:,2]
                self.cells = np.asarray(celldata[:,3:6],int)
                self.neigh = np.asarray(celldata[:,6:9])
                self.nfaces = 3*np.ones((self.Nc,),np.int)
                self.maxfaces=self.MAXFACES

            elif ncols==11: # Quad grid format
                nfaces = celldata[:,0]
                self.nfaces = np.zeros((self.Nc,),np.int)
                self.nfaces[:] = nfaces
                self.xv = celldata[:,1]
                self.yv = celldata[:,2]
                self.cells = np.asarray(celldata[:,3:7],int)
                self.neigh = np.asarray(celldata[:,7:11])
                self.maxfaces=4
        else: # Uneven number of cells
            celldata = celldata.tolist()
            nfaces = [ll[0] for ll in celldata]
            self.nfaces=np.array(nfaces,np.int)
            self.maxfaces=self.nfaces.max()
            self.cells = self._FillValue*np.ones((self.Nc,self.maxfaces),int)
            self.neigh = self._FillValue*np.ones((self.Nc,self.maxfaces),int)
            self.xv = np.zeros((self.Nc,))
            self.yv = np.zeros((self.Nc,))
            for ii in range(self.Nc):
                nf = self.nfaces[ii]
                self.xv[ii] = celldata[ii][1]
                self.yv[ii] = celldata[ii][2]
                self.cells[ii,0:nf] = celldata[ii][3:3+nf]
                self.neigh[ii,0:nf] = celldata[ii][3+nf:3+2*nf]

        
        self.edges = np.asarray(edgedata[:,0:2],int)
        self.Ne = self.edges.shape[0]
        self.mark = np.asarray(edgedata[:,2],int)
        self.grad = np.asarray(edgedata[:,3:5],int)
        if np.size(edgedata,1)==6:
            self.edge_id = np.asarray(edgedata[:,5],int)
        else:
            self.edge_id = np.zeros((self.Ne,),int)
        
        # Load the vertical grid info from vertspace.dat if it exists
        try:
            vertspace=readTXT(self.infile+'/vertspace.dat')
        except:
            if self.VERBOSE:
                print('Warning could not find vertspace.dat in folder, setting Nkmax=1')
            vertspace=0.0
            
        self.setDepth(vertspace)

        self.maskgrid()

    def __loadGrdNc(self):
        
        """
        Load the grid variables into the object from a netcdf file
        
        Try statements are for backward compatibility  
        
        Variables loaded are presently:
        'xp','yp','xv','yv','xe','ye','cells','face','nfaces','edges','neigh','grad',
        'gradf','mark','normal','n1','n2','df','dg','def','Ac','dv','dz','z_r','z_w','Nk','Nke'
        """
        
        self.__openNc()
        nc=self.nc

        # Get the dimension sizes
        for vv in list(self.griddims.keys()):
           try:
               setattr(self,vv,nc.dimensions[self.griddims[vv]].__len__())
           except:
               if self.VERBOSE:
                   print('Cannot find dimension: %s'%self.griddims[vv])
       
       
        for vv in list(self.gridvars.keys()):
            try:
                if vv=='def': # Cannot have this attribute name in python!
                    setattr(self,'DEF',nc.variables[self.gridvars[vv]][:])
                else:
                    setattr(self,vv,nc.variables[self.gridvars[vv]][:])
            except:
                if self.VERBOSE:
                    print('Cannot find variable: %s'%self.gridvars[vv])
         
        if 'Nk' in self.__dict__:
            self.Nk-=1 #These need to be zero based
        
        if 'nfaces' not in self.__dict__:
            self.MAXFACES = self.cells.shape[1]
            self.nfaces = self.MAXFACES*np.ones((self.Nc,),np.int)
            self.maxfaces = self.MAXFACES

        # If edges, grad or neigh have not been stored then calculate them
        if 'edges' not in self.__dict__:
            self.reCalcGrid()
        elif 'grad' not in self.__dict__:
            self.reCalcGrid()
        #elif not self.__dict__.has_key('neigh'):
        #    self.reCalcGrid()

        # Set the mark equal zero if doesn't exist
        if 'mark' not in self.__dict__:
            self.mark = np.zeros((self.Ne))

        # Check the _FillValue attribute is consistent with the grid
        # If the maximum cells value exceeds the number of points the fillvalue
        # is likely to be  999999
        if self.cells.max()>self.Np:
            if self.VERBOSE:
                print('Changing the _FillValue from {} to {}'.format(\
                    self._FillValue, self.cells.max()))
            self._FillValue=self.cells.max()

        if type(self.cells) != type(np.ma.MaskedArray()):
            self.maskgrid()
        else:
            self.cellmask=self.cells.mask

        #if type(self.DEF) == type(np.ma.MaskedArray()):
        #    if np.all(self.DEF.mask):
        #        self.calc_def()
        try:
            self.calc_def()
        except:
            print('No def array...')
       
    def maskgrid(self):
        """
        Mask the cells, face and neigh arrays
        """
        self.cellmask = self.cells==int(self._FillValue)
        #for ii in range(self.Nc):
        #    self.cellmask[ii,self.nfaces[ii]::]=True
            
        self.cells[self.cellmask]=0
        self.cells =\
            np.ma.masked_array(self.cells,mask=self.cellmask,fill_value=0)
        
        if 'face' in self.__dict__:
            self.face[self.cellmask]=0
            self.face =\
                np.ma.masked_array(self.face,mask=self.cellmask,fill_value=0)
            
        if 'neigh' in self.__dict__:
            self.neigh =\
                np.ma.masked_array(self.neigh,mask=self.cellmask,fill_value=0)
        
             
    def convert2hybrid(self):
        """
        Converts the suntans grid to a HybridGrid type
        """
        return HybridGrid(self.xp, self.yp, self.cells, nfaces=self.nfaces,\
            xv=self.xv, yv=self.yv,\
            mark=self.mark,\
            edges=self.edges,\
            grad=self.grad,\
            neigh=self.neigh,\
            _FillValue=self._FillValue)


    def reCalcGrid(self):
        """
        Recalculate some of the important grid properties manually

        Method for finding the following arrays: grad,edges,neigh,mark...
        """
        if self.VERBOSE:
            print('Re-calculating the grid variables...')

        grd = self.convert2hybrid()

        self.edges=grd.edges
        self.grad=grd.grad
        self.neigh=grd.neigh
        self.face=grd.face
        self.mark=grd.mark

        #self.calc_def()
        self.DEF = grd.DEF

    def calc_def(self):
        """
        Recalculate the edge to face distance
        """
        ne = np.array(self.face)

        #try:
        #    mask = ne.mask.copy()
        #except:
        #    mask = np.zeros(self.DEF.shape,np.bool)
        ne[self.cellmask]=0

        def dist(x0,x1,y0,y1):
            return np.sqrt( (x0-x1)**2. + (y0-y1)**2.)

        self.DEF = dist(self.xv,self.xe[ne].T,self.yv,self.ye[ne].T).T

        self.DEF = np.ma.masked_array(self.DEF,mask=self.cellmask)


    def plot(self,**kwargs):
        """
          Plot the unstructured grid data
        """
        if 'clim' in self.__dict__:
            clim = self.clim
        else:
            clim = [self.dv.min(), self.dv.max()]
       
        self.fig,self.ax,self.patches,self.cb=unsurf(self.xy,self.dv,xlim=self.xlims,ylim=self.ylims,\
            clim=clim,**kwargs)
        
        plt.title('SUNTANS Grid Bathymetry [m]')
        
    def plotvtk(self):
        """
          Plot the unstructured grid data using vtk libraries
        """
        if 'clim' in self.__dict__:
            clim = self.clim
        else:
            clim = [self.dv.min(), self.dv.max()]
        points = np.column_stack((self.xp,self.yp,0.0*self.xp))
        self.fig, h, ug,d, title=unsurfm(points,self.cells,self.dv,clim=clim,title='SUNTANS Grid Bathymetry [m]',\
            colormap='gist_earth')
        
    
    def plotBC(self):
        """
        Plot the boundary markers and the grid nodes
        """
        # Find the edge points
        xe = np.mean(self.xp[self.edges],axis=1)
        ye = np.mean(self.yp[self.edges],axis=1)
        plt.plot(self.xp,self.yp,'.')
        plt.plot(xe,ye,'k+')
        plt.plot(xe[self.mark==1],ye[self.mark==1],'ro')
        plt.plot(xe[self.mark==2],ye[self.mark==2],'yo')
        plt.plot(xe[self.mark==3],ye[self.mark==3],'go')
        plt.plot(xe[self.mark==4],ye[self.mark==4],'co')
        plt.legend(('Node','Edge','Marker=1','Marker=2','Marker=3','Marker=4'))
        plt.axis('equal')
    
    def plotmesh(self,ax=None,facecolors='none',linewidths=0.2,**kwargs):
        """
        Plots the outline of the grid mesh
        """
        fig = plt.gcf()
        if ax is None:
            ax = fig.gca()
    
        xlim=self.xlims
        ylim=self.ylims
        collection = PolyCollection(self.xy,facecolors=facecolors,\
            linewidths=linewidths,**kwargs)
        
        ax.add_collection(collection)
    
        ax.set_aspect('equal')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        return ax, collection
    
    def plothist(self):
        """
        Histogram plot of distances between cell voronoi points
        """
        if 'dg' not in self.__dict__:
            self.calc_dg()
            
        dg = self.dg
        ind = np.argwhere(dg!=0) # boundary edges have dg = 0
            
        fig = plt.gcf()
        ax = fig.gca()
        
        plt.hist(self.dg[ind],bins=100,log=False)
        textstr='Min = %3.1f m\nMean = %3.1f m\nMedian = %3.1f m\nMax = %3.1f m'%\
           (np.min(self.dg[ind]),np.mean(self.dg[ind]),np.median(self.dg[ind]),np.max(self.dg[ind]))
        textstr += '\nNc = %d\nNp = %d\nNe = %d\nmaxfaces = %d'\
            %(self.Nc,self.Np,self.Ne,self.maxfaces)
        plt.text(0.7,0.6,textstr,transform=ax.transAxes)
        plt.xlabel('dg [m]')
        plt.ylabel('Edge Count')

    def plotedgedata(self,z,xlims=None,ylims=None,**kwargs):
        """
          Plot the unstructured grid edge data
        """

        assert(z.shape[0] == self.Ne),\
            ' length of z scalar vector not equal to number of edges, Ne.'
        
        # Find the colorbar limits if unspecified
        if self.clim is None:
            self.clim=[]
            self.clim.append(np.min(z))
            self.clim.append(np.max(z))
        # Set the xy limits
        if xlims is None or ylims is None:
            xlims=self.xlims 
            ylims=self.ylims
        
        xylines = [self.xp[self.edges],self.yp[self.edges]]
        self.fig,self.ax,self.collection,self.cb=edgeplot(xylines,z,xlim=xlims,ylim=ylims,\
            clim=self.clim,**kwargs)
            
    def to_latlon(self, x, y):
        """
        Convert x,y to lat/lon using the current grid projection

        Need to set these attributes manually:
            - projstr (leave as none for utm)
            - utmzone (integer)
            - isnorth (true/false)
        """
        if self.projstr is None and self.utmzone is None:
            if self.VERBOSE:
                print('Warning - no cartographic projection specified')
            lon = self._FillValue*np.ones_like(x)
            lat = self._FillValue*np.ones_like(y)
        else:
            P = MyProj(self.projstr, utmzone=self.utmzone, isnorth=self.isnorth)
            lon, lat = P.to_ll(x, y)


        return lon, lat

    def to_wgs84(self, utmzone, isnorth):
        """
        Convert the node coordinates to WGS84 map projection
        """
        # Convert the nodes
        #ll = utm2ll(np.column_stack([self.xp, self.yp]),\
        #        utmzone, north=isnorth)
        #self.xp = ll[:,0]
        #self.yp = ll[:,1]

        ## Convert the cell centers
        #ll = utm2ll(np.column_stack([self.xv, self.yv]),\
        #        utmzone, north=isnorth)
        #self.xv = ll[:,0]
        #self.yv = ll[:,1]
        P = MyProj(None, utmzone=utmzone, isnorth=isnorth)
        self.xp, self.yp, lat = P.to_ll(self.xp, self.yp)
        self.xv, self.yv, lat = P.to_ll(self.xv, self.yv)



        # Update the plot coordinates
        self.xy = self.cellxy(self.xp, self.yp)

    def cellxy(self, xpin, ypin):
        """ 
        Returns a list of Nx2 vectors containing the grid cell node coordinates
            
        Used by spatial ploting routines 
        """
        xp = np.zeros((self.Nc,self.maxfaces+1))
        yp = np.zeros((self.Nc,self.maxfaces+1))
        
        cells=self.cells.copy()
        #cells[self.cells.mask]=0
        cellmask = cells==int(self._FillValue)
        cells[cellmask]=0
        #print(self._FillValue)

        xp[:,:self.maxfaces]=xpin[cells]
        xp[list(range(self.Nc)),self.nfaces]=xpin[cells[:,0]]
        yp[:,:self.maxfaces]=ypin[cells]
        yp[list(range(self.Nc)),self.nfaces]=ypin[cells[:,0]]

        #xp[self.cells.mask]==0
        #yp[self.cells.mask]==0

        xy = np.zeros((self.maxfaces+1,2))
        def _closepoly(ii):
            nf=self.nfaces[ii]+1
            xy[:nf,0]=xp[ii,:nf]
            xy[:nf,1]=yp[ii,:nf]
            return xy[:nf,:].copy()

        return [_closepoly(ii) for ii in range(self.Nc)]


        # Old Method
        #return [closePoly(self.xp[self.cells[ii,0:self.nfaces[ii]]],\
        #    self.yp[self.cells[ii,0:self.nfaces[ii]]]) for ii in range(self.Nc)]
        
    def saveBathy(self,filename):
        """
            Saves the grid bathymetry to an xyz ascii file
        """
        f = open(filename,'w')
        
        for x,y,z in zip(self.xv,self.yv,self.dv):
            f.write('%10.6f %10.6f %10.6f\n'%(x,y,z))
            
        f.close()
    
    def loadBathy(self,filename):
        """
        Loads depths from a text file into the attribute 'dv'
        """
        depths = readTXT(filename)
        if len(depths) != self.Nc:
            print('Error - number of points in depth file (%d) does not match Nc (%d)'%(len(depths),self.Nc))
        else:
            dv=depths[:,2]

        self.dv = dv
        return dv
    
    def saveVertspace(self,filename):
        """
        Saves vertspace.dat
        """
        f = open(filename,'w')
        
        for dz in self.dz:
            f.write('%10.6f\n'%(dz))
            
        f.close()
        
    def saveCells(self,filename):
        """
        Save cells.dat into hybrid grid format
        """

        f = open(filename,'w')

        for ii in range(self.Nc):
            outstr = '%d %10.6f %10.6f '%(self.nfaces[ii],self.xv[ii],self.yv[ii])
            for nn in range(self.nfaces[ii]):
                outstr += '%d '%self.cells[ii,nn]

            for nn in range(self.nfaces[ii]):
                outstr += '%d '%self.neigh[ii,nn]

            outstr += '\n'
            f.write(outstr)
            

        f.close()

    def saveEdges(self,filename):
        """
        Saves the edges.dat data to a text file
        
        Used e.g. when the boundary markers have been updated
        """

        f = open(filename,'w')
#        
#        for e1,e2,m,g1,g2 in zip(self.edges[:,0],self.edges[:,1],self.mark,self.grad[:,0],self.grad[:,1]):
#            f.write('%d %d  %d  %d  %d\n'%(e1,e2,m,g1,g2))

        # Write an extra column that has the boundary edge segment flag    
        if 'edge_id' in self.__dict__:
            for e1,e2,m,g1,g2,ef1 in zip(self.edges[:,0],self.edges[:,1],self.mark,self.grad[:,0],self.grad[:,1],self.edge_id):
                f.write('%d %d  %d  %d  %d  %d\n'%(e1,e2,m,g1,g2,ef1))
        else:
            for e1,e2,m,g1,g2 in zip(self.edges[:,0],self.edges[:,1],self.mark,self.grad[:,0],self.grad[:,1]):
                f.write('%d %d  %d  %d  %d\n'%(e1,e2,m,g1,g2))
            
        f.close()

    def saveGrid(self,outpath):
        """
        Saves the suntnas grid ascii files to the specified path
        """
        self.saveCells(outpath+'/cells.dat')
        self.saveEdges(outpath+'/edges.dat')

        # Save to points.dat
        f = open(outpath+'/points.dat','w')
    
        for x,y in zip(self.xp,self.yp):
            f.write('%10.6f %10.6f  0\n'%(x,y))
    
        f.close()

        print('Complete - grid saved to folder: %s'%outpath)


        
    def find_nearest(self,xy,NNear=1):
        """
        Returns the grid indices of the closest points to the nx2 array xy
        
        Uses the scipy KDTree routine
        """
        
        if 'kd' not in self.__dict__:
            self._kd = spatial.cKDTree(np.vstack((self.xv,self.yv)).T)
    
        # Perform query on all of the points in the grid
        dist,ind=self._kd.query(xy,k=NNear)
        
        return dist, ind

    def find_nearest_boundary(self,markertype=1):
        """
        Returns the index 
        """
        if 'xe' not in self.__dict__:
            self.calc_edgecoord()

        bdyidx = self.mark==markertype
        kd = spatial.cKDTree(np.vstack((self.xe[bdyidx],self.ye[bdyidx])).T)

        xy = np.vstack((self.xv,self.yv)).T

        dist,edgeidx=kd.query(xy,k=1)

        return dist, edgeidx

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

        
    def calc_dg(self):
        """
        Manually calculate the distance between voronoi points, 'dg'
        """
        if self.VERBOSE:
            print('Calculating dg...')
            print(np.shape(self.grad))
        
        grad = self.grad
        Ne = len(grad)
        for ii in range(Ne):
            if grad[ii,0]==-1:
                grad[ii,0]=grad[ii,1]
            elif grad[ii,1]==-1:
                grad[ii,1]=grad[ii,0]
                
                
        x1 = self.xv[grad[:,0]]
        x2 = self.xv[grad[:,1]]
        y1 = self.yv[grad[:,0]]
        y2 = self.yv[grad[:,1]]
        
        dx=x1-x2
        dy=y1-y2
        
        self.dg = np.sqrt( dx*dx + dy*dy )

    def count_cells(self):
        """
        Count the total number of 3-D cells
        """
        return np.sum(self.Nk+1)
        
    def calc_tangent(self):
        """
        Calculate the tangential vector for the edges of each cell
        """
        if '_tx' not in self.__dict__:
            dx = np.zeros(self.cells.shape)    
            dy = np.zeros(self.cells.shape)  
    
            dx[:,0:-1] = self.xp[self.cells[:,1::]] - self.xp[self.cells[:,0:-1]]               
            dy[:,0:-1] = self.yp[self.cells[:,1::]] - self.yp[self.cells[:,0:-1]]               
    
            for ii in range(self.Nc):
                dx[ii,self.nfaces[ii]-1] = self.xp[self.cells[ii,0]] - self.xp[self.cells[ii,self.nfaces[ii]-1]]  
                dy[ii,self.nfaces[ii]-1] = self.yp[self.cells[ii,0]] - self.yp[self.cells[ii,self.nfaces[ii]-1]]  
   
            
            mag = np.sqrt(dx*dx + dy*dy)
            
            self._tx = dx/mag
            self._ty = dy/mag
            self._mag = mag
                     
        return self._tx, self._ty, self._mag

    def calc_edgecoord(self):
        """
        Manually calculate the coordinates of the edge points
        """
        self.xe = np.mean(self.xp[self.edges],axis=1)
        self.ye = np.mean(self.yp[self.edges],axis=1)    

    def get_facemark(self):
        """
        Finds the cells next to type-2 or type-3 boundaries
        """
        mask = self.face.mask
        face = self.face.copy()
        face[mask]=0
        facemark = self.mark[face]
        facemark[mask]=0
        return np.min(np.max(facemark,axis=-1),3)



    def setDepth(self,vertspace):
        """
        Calculates and sets the depth variables based on the vertspace vector
        """
        self.dz=vertspace
        self.Nkmax=np.size(self.dz)
        
        # Calculate the mid-point depth
        if not self.Nkmax == 1:
            #z_bot = np.cumsum(self.dz)
            #z_top = np.hstack((0.0,z_bot[:-1]))
            #self.z_r = 0.5*(z_bot+z_top)
            self.z_r = calc_z(self.dz)
        else:
            self.z_r=np.array([0.0])
            
    def calcVertSpace(self, Nkmax, r, depthmax):
        """
        Calculates the vertical spacing based on an exponential stretching function
        """
        
        vertspace = np.zeros((Nkmax,))
        
        if r < 1.0 or r > 1.1:
            print('r must be between 1.0 and 1.1')
            
        if Nkmax == 0:
            vertspace[0] = depthmax
        else:
            if r == 1.0:
                vertspace[0] = depthmax/Nkmax
            else:
                vertspace[0] = depthmax * (r-1.0) / (r**float(Nkmax) - 1.0)
            for k in range(1,Nkmax):
                vertspace[k] = r*vertspace[k-1]
                
        return vertspace
            
    def pnt2cells(self,pnt_i):
        """
        Returns the cell indices for a point, pnt_i
        
        (Stolen from Rusty's TriGrid class)
        """
        if '_pnt2cells' not in self.__dict__:
            # build hash table for point->cell lookup
            self._pnt2cells = {}
            for i in range(self.Nc):
                for j in range(3):
                    if self.cells[i,j] not in self._pnt2cells:
                        #self._pnt2cells[self.cells[i,j]] = set()
                        self._pnt2cells[self.cells[i,j]] = []
                    #self._pnt2cells[self.cells[i,j]].add(i)
                    self._pnt2cells[self.cells[i,j]].append(i)
        return self._pnt2cells[pnt_i]
        
    def cell2node(self,cell_scalar):
        """
        Map a cell-based scalar onto a node
        
        This is calculated via a mean of the cells connected to a node(point)
        """
        # Simple mean
        #node_scalar = [np.mean(cell_scalar[self.pnt2cells(ii)]) for ii in range(self.Np)]
        
        # Area weighted interpolation
        node_scalar = [np.sum(cell_scalar[self.pnt2cells(ii)]*self.Ac[self.pnt2cells(ii)])\
            / np.sum( self.Ac[self.pnt2cells(ii)]) for ii in range(self.Np)]
        return np.array(node_scalar)


    def interpLinear(self,cell_scalar,xpt,ypt,cellind,k=0):
        """
        Linearly interpolates onto the coordinates 'xpt' and 'ypt'
        located at cell index: 'cellind'.capitalize
        
        Use trisearch to get the cell index of xpt/ypt.
        """
        
        # Find the spatial shift
        dx = xpt - self.xv[cellind]
        dy = ypt - self.yv[cellind]
        
        # Find the gradient at the cell centres
        dphi_dx, dphi_dy = self.gradH(cell_scalar,cellind=cellind,k=k)
        
        return cell_scalar[cellind] + dphi_dx*dx + dphi_dy*dy
        

    def gradH(self,cell_scalar,k=0,cellind=None):
        """
        Compute the horizontal gradient of a cell-centred quantity

        """
        if self.maxfaces==3:
            dX,dY=self.gradHplane(cell_scalar,k=k,cellind=cellind)
        else:
            dX,dY=self.gradHdiv(cell_scalar,k=k)

        return dX,dY
 

    def gradHplane(self,cell_scalar,k=0,cellind=None):
        """
        Compute the horizontal gradient of a cell-centred quantity
        
        Finds the equation for the plane of the 3 points then takes the gradient
        of this equation
        
        Returns: d(phi)/dx, d(phi)/dy
        
        Derived from Phil's C code (diffusion.c):
          //PJW get corresponding points
          REAL xA = grid->xp[pA];
          REAL yA = grid->yp[pA];
          REAL xB = grid->xp[pB];
          REAL yB = grid->yp[pB];
          REAL xC = grid->xp[pC];
          REAL yC = grid->yp[pC];
        
          //PJW form the vectors
          //PJW first we need to get the points to build the vectors
          REAL ABx = xB - xA;
          REAL ABy = yB - yA;
          REAL ABz = z[pB][alevel] - z[pA][alevel];
        
          REAL ACx = xC - xA;
          REAL ACy = yC - yA;
          REAL ACz = z[pC][alevel] - z[pA][alevel];
        
          //PJW now we can take the cross product
          REAL mx = ABy*ACz - ABz*ACy;
          REAL my = ABz*ACx - ABx*ACz;
          REAL mz = ABx*ACy - ACx*ABy;
        
        //  printf("mx = %.3e my = %.3e mz = %.3e\n",mx,my,mz);
          //PJW now we can get the slopes and return them
          duv_over_dxj[0] = -mx/mz;
          duv_over_dxj[1] = -my/mz;
        """
        if cellind is None:
            cellind = np.arange(self.Nc,dtype=np.int)
            
        node_scalar = self.cell2nodekind(cell_scalar,cellind,k=k)
        #self.nc.variables[varname].dimensions
        xA = self.xp[self.cells[cellind,0]]
        yA = self.yp[self.cells[cellind,0]]
        xB = self.xp[self.cells[cellind,1]]
        yB = self.yp[self.cells[cellind,1]]
        xC = self.xp[self.cells[cellind,2]]
        yC = self.yp[self.cells[cellind,2]]
        
        zA = node_scalar[self.cells[cellind,0]]
        zB = node_scalar[self.cells[cellind,1]]
        zC = node_scalar[self.cells[cellind,2]]
        
        ABx = xB - xA
        ABy = yB - yA
        ABz = zB - zA
        
        ACx = xC - xA
        ACy = yC - yA
        ACz = zC - zA
        
        mx = ABy*ACz - ABz*ACy
        my = ABz*ACx - ABx*ACz
        mz = ABx*ACy - ACx*ABy
        
        return -mx/mz, -my/mz

        
    def cell2nodekind(self,cell_scalar,cellind,k=0):
        """
        Map a cell-based scalar onto a node
        
        Only does it for nodes connected to cells in: 'cellind'. This is faster and more generic.
        
        The node_scalar array still has size(Np) although nodes that aren't connected
        to cells in 'cellind' are simply zero.
        
        Uses sparse matrices to do the heavy lifting
        """
        

        Nc = cellind.shape[0]
        ii = np.arange(0,Nc)
        i = cellind[ii]

        # Find the row and column indices of the sparse matrix
        rowindex = self.cells[i,:] 
        colindex = np.repeat(i.reshape((Nc,1)), self.maxfaces,axis=1)

        # This is a bit of a hack for when cells is not a MaskedArray object
        try:
            mask = (k <= self.Nk[colindex]) & (self.cells.mask==False)
        except:
            mask = (k <= self.Nk[colindex])
        
        cell_scalar3d = np.repeat(cell_scalar[i].reshape((Nc,1)), self.maxfaces,axis=1)
        area = np.repeat(self.Ac[i].reshape((Nc,1)), self.maxfaces,axis=1)
        
        #Build the sparse matrices
        Asparse = sparse.coo_matrix((area[mask],(rowindex[mask],colindex[mask])),shape=(self.Np,self.Nc),dtype=np.double)
        datasparse = sparse.coo_matrix((cell_scalar3d[mask],(rowindex[mask],colindex[mask])),shape=(self.Np,self.Nc),dtype=np.double)

        # This step is necessary to avoid summing duplicate elements
        datasparse=sparse.csr_matrix(datasparse).tocoo()
        Asparse=sparse.csr_matrix(Asparse).tocoo()
        
        node_scalar = datasparse.multiply(Asparse).sum(axis=1) / Asparse.sum(axis=1)
         
        # Return a masked array
        
        #node_scalar=np.ma.array(node_scalar,mask=mask).squeeze()
        node_scalar=np.array(node_scalar).squeeze()
        mask = np.isnan(node_scalar)
        node_scalar[mask]=0.
        
        return node_scalar
 
    def cell2nodekindold(self,cell_scalar,cellind,k=0):
        """
        Map a cell-based scalar onto a node
        
        Only does it for nodes connected to cells in: 'cellind'. This is faster and more generic.
        
        The node_scalar array still has size(Np) although nodes that aren't connected
        to cells in 'cellind' are simply zero.
        
        Uses sparse matrices to do the heavy lifting
        """
        
        Nc = cellind.shape[0]
        
        if '_datasparse' not in self.__dict__:
            self._datasparse=[]
            self._Asparse=[]
            for kk in range(self.Nkmax):
                self._datasparse.append(sparse.dok_matrix((self.Np,self.Nc),dtype=np.double))
                self._Asparse.append(sparse.dok_matrix((self.Np,self.Nc),dtype=np.double))
                
                for ii in range(Nc):
                    i = cellind[ii]
                    for j in range(3):
                        if kk <= self.Nk[i]:
                            self._Asparse[kk][self.cells[i,j],i] = self.Ac[i]
                
                self._Asparse[kk]=self._Asparse[kk].tocoo() # COO sparse format is faster for multiplication

        for ii in range(Nc):
            i = cellind[ii]
            for j in range(3):
                if k <= self.Nk[i]:
                    self._datasparse[k][self.cells[i,j],i] = cell_scalar[i]
                
        node_scalar = self._datasparse[k].tocoo().multiply(self._Asparse[k]).sum(axis=1) / self._Asparse[k].sum(axis=1)
        
        return np.array(node_scalar).squeeze()
                                            
    
    def gradHdiv(self,phi,k=0):
        """
        Computes the horizontal gradient of variable, phi, along layer, k.
        
        Uses divergence (Green's) theorem to compute the gradient. This can be noisy 
        for triangular control volumes.
        
        Returns: d(phi)/dx, d(phi)/dy
        
        Based on MATLAB code sungradient.m
        """
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
            indk = operator.or_(k>=self.Nk[nc1], k>=self.Nk[nc2])
            ind3 = operator.and_(indk, self.Nk[nc2]>self.Nk[nc1])
            nc1[ind3]=nc2[ind3]
            ind4 = operator.and_(indk, self.Nk[nc1]>self.Nk[nc2])
            nc2[ind4]=nc1[ind4]
            
            # Calculate gradient across face            
            return (phi[nc1]-phi[nc2]) / self.dg[jj]
            
        ne = self.face #edge-indices
        mask = ne.mask.copy()
        ne[mask]=0
        
        Gn_phi = _GradientAtFace(phi,ne,k)
        Gn_phi[mask]=0
        
        Gx_phi = Gn_phi * self.n1[ne] * self.DEF * self.df[ne]
        Gy_phi = Gn_phi * self.n2[ne] * self.DEF * self.df[ne]
        dX = np.sum(Gx_phi,axis=1)/self.Ac;
        dY = np.sum(Gy_phi,axis=1)/self.Ac;

        return dX, dY

    def spatialfilter(self,phi,dx, ftype='low'):
        """
        Perform a gaussian spatial lowpass filter on the
        variable, phi. 
        dx is the filter length (gaussian simga parameter).
        """

        if '_spatialfilter' not in self.__dict__:
            print('Building the filter matrix...')
            xy = np.vstack((self.xv,self.yv)).T
            self._spatialfilter = ufilter(xy,dx)

        if ftype=='low':
            return self._spatialfilter(phi)
        elif ftype=='high':
            return phi - self._spatialfilter(phi)
        else:
            raise Exception('unknown filter type %s. Must be "low" or "high".'%ftype)

    
    def writeNC(self,outfile):
        """
        Export the grid variables to a netcdf file
        """
        
        nc = Dataset(outfile, 'w', format='NETCDF4_CLASSIC')
        nc.Description = 'SUNTANS History file'
        nc.Author = ''
        nc.Created = datetime.now().isoformat()

        nc.createDimension('Nc', self.Nc)
        nc.createDimension('Np', self.Np)
        try:
            nc.createDimension('Ne', self.Ne)
        except:
            print('No dimension: Ne')
        nc.createDimension('Nk', self.Nkmax)
        nc.createDimension('Nkw', self.Nkmax+1)
        nc.createDimension('numsides', self.maxfaces)
        nc.createDimension('two', 2)
        nc.createDimension('time', 0) # Unlimited
        
        # Write the grid variables
        def write_nc_var(var, name, dimensions, attdict, dtype='f8'):
            tmp=nc.createVariable(name, dtype, dimensions)
            for aa in list(attdict.keys()):
                tmp.setncattr(aa,attdict[aa])
            nc.variables[name][:] = var
    
        #gridvars = ['suntans_mesh','cells','face','nfaces','edges','neigh','grad','xp','yp','xv','yv','xe','ye',\
        #    'normal','n1','n2','df','dg','def','Ac','dv','dz','z_r','z_w','Nk','Nke','mark']
        
        
        self.Nk += 1 # Set to one-base in the file (reset to zero-base after)
        self.suntans_mesh=[0]  
        for vv in self.gridvars:
            if vv in self.__dict__ and vv != 'time':
                if self.VERBOSE:
                    print('Writing variables: %s'%vv)
                write_nc_var(self[vv],vv,ugrid[vv]['dimensions'],ugrid[vv]['attributes'],dtype=ugrid[vv]['dtype'])
            
            # Special treatment for "def"
            if vv == 'def' and 'DEF' in self.__dict__:
                if self.VERBOSE:
                    print('Writing variables: %s'%vv)
                write_nc_var(self['DEF'],vv,ugrid[vv]['dimensions'],ugrid[vv]['attributes'],dtype=ugrid[vv]['dtype'])

        nc.close()
        self.Nk -= 1 # set back to zero base

    def create_nc_var(self,outfile, name, dimensions, attdict, dtype='f8',zlib=False,complevel=0,fill_value=999999.0):
        
        nc = Dataset(outfile, 'a')
        tmp=nc.createVariable(name, dtype, dimensions,zlib=zlib,complevel=complevel,fill_value=fill_value)
        for aa in list(attdict.keys()):
            tmp.setncattr(aa,attdict[aa])
        #nc.variables[name][:] = var	
        nc.close()
        
    def __del__(self):
        if 'nc' in self.__dict__:
            self.nc.close()
        
    def __openNc(self):
        #nc = Dataset(self.ncfile, 'r', format='NETCDF4') 
        if self.VERBOSE:
            print('Loading: %s'%self.ncfile)
        try: 
            self.nc = MFDataset(self.ncfile,aggdim='time')
        except:
            if type(self.ncfile)==list:
                self.ncfile = self.ncfile[0]
            self.nc = Dataset(self.ncfile, 'r')
 
    def __getitem__(self,y):
        x = self.__dict__.__getitem__(y)
        return x
        

#################################################

class Spatial(Grid):
    
    """ Class for reading SUNTANS spatial netcdf output files """
    
    # Set some default parameters
    tstep=0
    klayer=[0] # -1 get seabed
    # Note that if j is an Nx2 array of floats the nearest cell will be found 
    j=None

    variable='eta'
    
    # Plotting parmaters
    clim=None
    
    def __init__(self,ncfile, **kwargs):
        
        self.ncfile = ncfile
        

        self.__dict__.update(kwargs)

        # Open the netcdf file
        #self.__openNc()
        
        # Load the grid (superclass)
        Grid.__init__(self, ncfile)  
        
        #self.xy = self.cellxy()
        
        # Load the time variable
        try:
            self.loadTime()
            # Update tstep 
            self.updateTstep()
        except:
            print('No time variable.')
        
        # Load the global attributes
        self.loadGlobals()


    def loadData(self, variable=None):
        """
        High-level wrapper to load different variables into the 'data' attribute
        
        """
        
        if variable is None:
            variable=self.variable
            
        if variable=='speed':
            return self.loadSpeed()
        elif variable=='vorticity':
            self.data = self.vorticity_circ(k=self.klayer[0])
            self.long_name = 'Vertical vorticity'
            self.units = 's-1'
            return
        elif variable=='div_H':
            self.data = self.calc_divergence()
            self.long_name = 'Horizontal divergence'
            self.units = 's-1'
            return
        elif variable=='streamfunction':
            self.data = self.calc_streamfunction()
            self.long_name = 'Horizontal streamfunction'
            self.units = 'm2 s-1'
            return
        elif variable=='PEanom':
            self.data = self.calc_PEanom()
        elif variable=='agemean':
        
            self.data = self.agemean()
            self.long_name = 'Mean age'
            self.units = 'days'
            return self.data
        elif variable=='dzz':
            if self.hasVar('dzz'):
                return self.loadDataRaw(variable='dzz')
            else:
                # Calculate dzz internally
                eta = self.loadDataRaw(variable='eta',setunits=False)
                return self.getdzz(eta)

        elif variable=='dzf':
            if self.hasVar('dzf'):
                return self.loadDataRaw(variable=variable)
            else:
                eta = self.loadDataRaw(variable='eta',setunits=False)
                return self.getdzf(eta)

        elif variable=='ctop':
            eta = self.loadDataRaw(variable='eta',setunits=False)
            return self.getctop(eta)

        elif variable=='etop':
            eta = self.loadDataRaw(variable='eta',setunits=False)
            etop,etaedge = self.getetop(eta)
            return etop

        elif variable=='buoyancy':
            self.data = self.calc_buoyancy()
            return self.data

        elif variable=='KE':
            self.data = self.calc_KE()
            return self.data

        elif variable=='PE':
            self.data = self.calc_PE()
            return self.data

        else:
            return self.loadDataRaw(variable=variable)
        
    
    def loadDataRaw(self,variable=None,setunits=True):
        """ 
            Load the specified suntans variable data as a vector
            
        """
        if variable is None:
            variable=self.variable
	
        # Get the indices of the horizontal dimension
        if self.hasDim(variable,self.griddims['Ne']) and self.j is None:
            j=list(range(self.Ne))
        elif self.hasDim(variable,self.griddims['Nc']) and self.j is None:
            j=list(range(self.Nc))
        else:
            if isinstance(self.j, np.ndarray):
                j = [self.j[0]]
            else:
                j = [self.j]


        nc = self.nc
        ndim = nc.variables[variable].ndim

        # Get the indices of the vertical dimension
        if ndim>2:
            if self.klayer[0] in [-1,-99,'seabed']:
                klayer = np.arange(0,self.Nkmax)
            elif self.klayer[0] == 'surface':
                #eta = self.loadDataRaw(variable='eta',setunits=False)
                eta = nc.variables['eta'][self.tstep,j]
                ctop=self.getctop(eta)
                klayer = list(range(0,ctop.max()+1))
            else:
                klayer = self.klayer
            
        # Set the long_name and units attribute
        if setunits:
            try:
                self.long_name = nc.variables[variable].long_name
            except:
                self.long_name = ''
            self.units= nc.variables[variable].units

        #        ndims = len(nc.variables[variable].dimensions)
        if ndim==1:
            self.data=nc.variables[variable][j]
        elif ndim==2:
            #print self.j
            self.data=nc.variables[variable][self.tstep,j]
        else: # 3D array
            data = nc.variables[variable][self.tstep,klayer,j]
            if self.klayer[0]==-99:
                self.data=data
            elif self.klayer[0]=='surface':
                self.data = self.get_surfacevar(data.squeeze(),ctop)
            elif self.klayer[0] in [-1,'seabed']:
                self.data = self.get_seabedvar(data.squeeze())
            else:
                self.data=data

            #if self.klayer[0]==-1: # grab the seabed values
            #    klayer = np.arange(0,self.Nkmax)

            #    #if type(self.tstep)==int:
            #    if isinstance(self.tstep,(int,long)):
            #        data=nc.variables[variable][self.tstep,klayer,j]
            #        self.data = data[self.Nk[j],j]
            #    else: # need to extract timestep by timestep for animations to save memory
            #        self.data=np.zeros((len(self.tstep),len(j)))
            #        i=-1
            #        for t in self.tstep:
            #            i+=1
            #            data=nc.variables[variable][t,klayer,j]
            #            self.data[i,:] = data[self.Nk[j],j]
            #elif self.klayer[0]==-99: # Grab all layers
            #    klayer = np.arange(0,self.Nkmax) 
            #    self.data=nc.variables[variable][self.tstep,klayer,j]
            #    
            #
            #else:
            #    klayer = self.klayer
            #    self.data=nc.variables[variable][self.tstep,klayer,j]
        
        # Mask the data
#        try:
#            fillval = nc.variables[variable]._FillValue
#        except:
        
        #self.mask = self.data==self._FillValue
        #if isinstance(self.data, np.ma.MaskedArray):
        #    self.data.mask[self.mask] = True
        #else:
        #    self.data[self.mask]=0.
        self.data = self.data.squeeze()
        
        return self.data
    
    def loadDataBar(self,variable=None):
        """
        Load a 3D variable and depth-average i.e. u
        """
        if variable is None:
            variable=self.variable
            
        ndim = self.nc.variables[variable].ndim
        if not ndim==3:
            print('Warning only 3D arrays can be used')
            return -1
        
        # Load time step by time step
        nt = len(self.tstep)
        databar = np.zeros((nt,self.Nc))
        
        self.klayer=[-99]

        tstep = self.tstep
        for ii,t in enumerate(tstep):
            self.tstep=[t]
            data = self.loadData(variable=variable)
            databar[ii,:] = self.depthave(data)            
            
        self.tstep=tstep
        
        return databar
        
    def loadSpeed(self):
        """
        Load the velocity magnitude into the data variable
        """
        u,v,w = self.getVector()
        
        self.long_name = 'Current speed'
        self.units = 'm s-1'
        
        speed = np.sqrt(u**2 + v**2 + w**2)
        
        self.data = speed
        
        return speed
        
    def loadTime(self):
         """
            Load the netcdf time as a vector datetime objects
         """
         #nc = Dataset(self.ncfile, 'r', format='NETCDF4') 
         nc = self.nc
         t = nc.variables[self.gridvars['time']]
         self.time = num2pydate(t[:],t.units)
         self.timeraw = t[:]
         
         self.Nt = self.time.shape[0]

    def loadGlobals(self):
        """
        Loads the global attributes into a dictionary
        """
        nc=self.nc
        self.globalatts={}
        for name in nc.ncattrs():
            self.globalatts.update({name:getattr(nc,name)})

    def listCoordVars(self):
        """
        List all of the variables that have the 'coordinate' attribute
        """
        
        vname=[]
        for vv in list(self.nc.variables.keys()):
            if hasattr(self.nc.variables[vv],'coordinates'):
                vname.append(vv)
                
        # Derived variables that can also be computed
        if 'vc' in vname and 'uc' in vname:
            vname.append('speed')
            vname.append('vorticity')
            vname.append('div_H') # Horizontal strain
            vname.append('streamfunction') # Horizontal streamfunction
            vname.append('KE')
        if 'rho' in vname:
            vname.append('PE')
            vname.append('buoyancy')
        if 'agec' in vname:
            vname.append('agemean')

            #vname.append('PEanom')
        return vname

    def build_tri(self):
        """
        Create a matplotlib triangulation object from the grid
        
        Used primarily to contour the data       
        """
        if '_tri' not in self.__dict__:
            if self.maxfaces==3:
                self._tri =tri.Triangulation(self.xp,self.yp,self.cells)
            else:

                # Need to compute a delaunay triangulation for mixed meshes
                if self.VERBOSE:
                    print('Computing delaunay triangulation and computing mask...')
                pts = np.vstack([self.xp,self.yp]).T
                D = spatial.Delaunay(pts)
                self._tri =tri.Triangulation(self.xp,self.yp,D.simplices)

                # Compute a mask by removing triangles outside of the polygon
                xy = D.points
                cells=D.simplices
                xv,yv = circumcenter(xy[cells[:,0],0],xy[cells[:,0],1],\
                    xy[cells[:,1],0],xy[cells[:,1],1],xy[cells[:,2],0],xy[cells[:,2],1])
                mask = self.find_cell(xv,yv)
                self._tri.set_mask(mask==-1)

    def build_tri_voronoi(self):
        """
        Create a matplotlib triangulation object from the grid voronoi points
        
        Used primarily to contour the data       
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



    def interpolate(self, z, x, y,  kind='linear'):
        """
        Interpolate data in 'z' on current grid onto points 'x' and 'y'

        kind = 'linear' or 'cubic'
        """
        self.build_tri_voronoi()

        if kind == 'linear':
            F = tri.LinearTriInterpolator(self._triv, z)
        elif kind == 'cubic':
            F = tri.CubicTriInterpolator(self._triv, z)

        return F(x, y)
        

 
    
    def plot(self,z=None, xlims=None, ylims=None, titlestr=None,
        vector_overlay=False, scale=1e-4, subsample=10, **kwargs):
        """
          Plot the unstructured grid data
        """
            
        if z is None:
            # Load the data if it's needed
            if 'data' not in self.__dict__:
                self.loadData() 
            z=self.data.ravel()
        
        # Find the colorbar limits if unspecified
        if self.clim is None:
            self.clim=[]
            self.clim.append(np.min(z))
            self.clim.append(np.max(z))
        # Set the xy limits
        if xlims is None or ylims is None:
            xlims=self.xlims 
            ylims=self.ylims
        
        self.fig,self.ax,self.patches,self.cb=unsurf(self.xy,z,xlim=xlims,ylim=ylims,\
            clim=self.clim,**kwargs)
            
        if titlestr is None:
            plt.title(self.genTitle())
        else:
            plt.title(titlestr)
            
        if vector_overlay:
             u,v,w = self.getVector()
             plt.quiver(self.xv[1::subsample],self.yv[1::subsample],u[0,1::subsample],v[0,1::subsample],scale=scale,scale_units='xy')
            #print 'Elapsed time: %f seconds'%(time.clock()-tic)
            
    def plotedgedata(self,z=None,xlims=None,ylims=None,titlestr=None,**kwargs):
        """
          Plot the unstructured grid edge data
        """
            
        if z is None:
            # Load the data if it's needed
            if 'data' not in self.__dict__:
                self.loadData() 
            z=self.data.ravel()

        assert(z.shape[0] == self.Ne),\
            ' length of z scalar vector not equal to number of edges, Ne.'
        
        # Find the colorbar limits if unspecified
        if self.clim is None:
            self.clim=[]
            self.clim.append(np.min(z))
            self.clim.append(np.max(z))
        # Set the xy limits
        if xlims is None or ylims is None:
            xlims=self.xlims 
            ylims=self.ylims
        
        xylines = [self.xp[self.edges],self.yp[self.edges]]
        self.fig,self.ax,self.collection,self.cb=edgeplot(xylines,z,xlim=xlims,ylim=ylims,\
            clim=self.clim,**kwargs)
            
        if titlestr is None:
            plt.title(self.genTitle())
        else:
            plt.title(titlestr)
 
    def contourf(self, clevs=20,z=None,xlims=None,ylims=None,filled=True,\
            k=0,cellind=None,vector_overlay=False,colorbar=True,\
            scale=1e-4,subsample=10,titlestr=None,\
            ax=None, fig=None, **kwargs):
        """
        Filled contour plot of  unstructured grid data
        """
        #if not self.__dict__.has_key('data'):
        #    self.loadData()
            
        if z is None:
            z=self.data
            
        # Find the colorbar limits if unspecified
        if self.clim is None:
            self.clim=[]
            self.clim.append(np.min(z))
            self.clim.append(np.max(z))
        # Set the xy limits
        if xlims is None or ylims is None:
            xlims=self.xlims 
            ylims=self.ylims
        
        # Build the triangulation object
        self.build_tri()
       
        if fig is None:
            fig = plt.gcf()
        if ax is None:
            ax = fig.gca()
        
        # Calculate the nodal data and mask if necessary
        #if k==0:
        #    zdata = self.cell2node(z)
        #else:
        if cellind is None:
            cellind = np.arange(self.Nc)
            
        zdata = self.cell2nodekind(z,cellind=cellind,k=k)
        
        # Amplitude plot (note that the data must be on nodes for tricontourf)
        if type(clevs)==type(1): # is integer
            V = np.linspace(self.clim[0],self.clim[1],clevs)
        else:
            V = clevs
            
        if filled:
            camp = ax.tricontourf(self._tri, zdata, V, **kwargs)
        else:
            camp = ax.tricontour(self._tri, zdata, V, **kwargs)
                
        ax.set_aspect('equal')
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)
        
        if filled and colorbar:
            self.cb = fig.colorbar(camp)
        
        if titlestr is None:
            plt.title(self.genTitle())
        else:
            plt.title(titlestr)
            
        if vector_overlay:
             u,v,w = self.getVector()
             plt.quiver(self.xv[1::subsample],self.yv[1::subsample],u[1::subsample],v[1::subsample],scale=scale,scale_units='xy')

        return camp
             
            
    def plotvtk(self,vector_overlay=False,scale=1e-3,subsample=1,**kwargs):
        """
          Plot the unstructured grid data using vtk libraries
        """
        # Mayavi libraries
        from mayavi import mlab
        
        # Load the data if it's needed
        if 'data' not in self.__dict__:
            self.loadData()
        
        points = np.column_stack((self.xp,self.yp,0.0*self.xp))
                
        self.fig,h,ug,d,title = unsurfm(points,self.cells,self.data,clim=self.clim,title=self.genTitle(),**kwargs)
        
        if vector_overlay:
             u,v,w = self.getVector()
             # Add vectorss to the unctructured grid object
             # This doesn't work ???       
             #vector = np.asarray((u,v,w)).T
             #ug.cell_data.vectors=vector
             #ug.cell_data.vectors.name='vectors'
             #ug.modified()
             #d.update()
             #d = mlab.pipeline.add_dataset(ug)
             #h2=mlab.pipeline.vectors(d,color=(0,0,0),mask_points=subsample,scale_factor=1./scale,scale_mode='vector')
             # This works             
             vec=mlab.pipeline.vector_scatter(self.xv, self.yv, self.yv*0, u, v, w)
             h2=mlab.pipeline.vectors(vec,color=(0,0,0),mask_points=subsample,scale_factor=1./scale,scale_mode='vector')
             
             # try out the streamline example
#             magnitude = mlab.pipeline.extract_vector_norm(vec)
#             pdb.set_trace()
#             h2 = mlab.pipeline.streamline(magnitude)
            
    def plotTS(self,j=None,**kwargs):
        """
         Plots a time series of the data at a given grid cell given by:
             self.j, self.klayer
        """

        # Load the time-series
        if j == None:
            j = self.j
            
        self.tstep = np.arange(0,len(self.time))
        self.j=j
        self.loadData()
        
        plt.ioff()
        self.fig = plt.gcf()
        ax = self.fig.gca()
        h = plt.plot(self.time,self.data,**kwargs) 
        ax.set_title(self.genTitle())
        
        return h
        
       
        
    def savefig(self,outfile,dpi=150):
        mod=self.fig.__class__.__module__
        name=self.fig.__class__.__name__
        
        if mod+'.'+name == 'matplotlib.figure.Figure':
            self.fig.savefig(outfile,dpi=dpi)
        else:
            self.fig.scene.save(outfile,size=dpi)
            
        if self.VERBOSE:
            print('SUNTANS image saved to file:%s'%outfile)
    
    def animate(self,xlims=None, ylims=None,\
            vector_overlay=False, scale=1e-4, subsample=10,\
            cbar=None, **kwargs):
        """
        Animates a spatial plot over all time steps
        
        Animation object is stored in the 'anim' property
        """
        #anim = unanimate(self.xy,self.data,self.tstep,xlim=self.xlims,ylim=self.ylims,clim=self.clim,**kwargs)
        #return anim
        
        # Load the vector data
        if vector_overlay:
            U,V,W = self.getVector()
            
        # Need to reload the data
        self.loadData()
    
        # Create the figure and axes handles
        #plt.ion()
        fig = plt.gcf()
        ax = fig.gca()
        #ax.set_animated('True')
        
        # Find the colorbar limits if unspecified
        if self.clim is None:
            self.clim=[]
            self.clim.append(np.min(self.data.ravel()))
            self.clim.append(np.max(self.data.ravel()))
           
        # Set the xy limits
        if xlims is None or ylims is None:
            xlims=self.xlims 
            ylims=self.ylims
        
            
        collection = PolyCollection(self.xy,**kwargs)
        collection.set_array(np.array(self.data[0,:]))
        collection.set_clim(vmin=self.clim[0],vmax=self.clim[1])
        ax.add_collection(collection)    
        ax.set_aspect('equal')
        ax.set_xlim(xlims)
        ax.set_ylim(ylims)

        ax.set_xlabel('Easting [m]')
        ax.set_ylabel('Northing [m]')

        title=ax.set_title("")
        if cbar is None:
            fig.colorbar(collection)
        
        qh=plt.quiver([],[],[],[])
        if vector_overlay:
            u=U[0,:]
            v=V[0,:]
            qh=plt.quiver(self.xv[1::subsample],self.yv[1::subsample],u[1::subsample],v[1::subsample]\
             ,scale=scale,scale_units='xy')
  
        def init():
            collection.set_array([])
            title.set_text("")
            qh.set_UVC([],[])
            return (title, collection)
               
        def updateScalar(i):
            if self.VERBOSE:
                print('\tStep %d of %d...'%(i,len(self.tstep)))

            collection.set_array(np.array(self.data[i,:]))
            collection.set_edgecolors(collection.to_rgba(np.array((self.data[i,:])))) 
            title.set_text(self.genTitle(i))
            if vector_overlay:
                qh.set_UVC(U[i,1::subsample],V[i,1::subsample])

            return (title,collection,qh)
  
        self.anim = animation.FuncAnimation(fig, updateScalar,\
                init_func=init, frames=len(self.tstep), interval=50, blit=True)

    def saveanim(self,outfile,fps=15):
        """
        Save the animation object to an mp4 movie
        """

        print('Building animation sequence...')

        ext = outfile[-4:]

        if ext=='.gif':
            self.anim.save(outfile,writer='imagemagick',fps=6)
        elif ext=='.mp4':
            print('Saving html5 video...')
            # Ensures html5 compatibility
            self.anim.save(outfile,writer='mencoder',fps=6,\
                bitrate=3600,extra_args=['-ovc','x264']) # mencoder options
                #bitrate=3600,extra_args=['-vcodec','libx264'])
        else:
            self.anim.save(outfile,writer='mencoder',fps=6,bitrate=3600)

        print('Complete - animation saved to: %s'%outfile)


        #try:
        #    print 'Building animation sequence...'
        #    self.anim.save(outfile, fps=fps,bitrate=3600)
        #    print 'Complete - animation saved to: %s'%outfile
        #except:
        #    print 'Error with animation generation - check if either ffmpeg or mencoder are installed.'
            
    def animateVTK(self,figsize=(640,480),vector_overlay=False,scale=1e-3,subsample=1):
        """
        Animate a scene in the vtk window
        
        """
        from mayavi import mlab
        
        # Load all the time steps into memory
        self.tstep=np.arange(0,len(self.time))
        nt = len(self.tstep)
        if vector_overlay:
             u,v,w = self.getVector()        
        
        self.loadData()
        
        
        # Find the colorbar limits if unspecified
        if self.clim is None:
            self.clim=[]
            self.clim.append(np.min(self.data))
            self.clim.append(np.max(self.data))
            
        # Generate the initial plot
        points = np.column_stack((self.xp,self.yp,0.0*self.xp))
        
        titlestr='%s [%s]\n Time: %s'%(self.long_name,self.units,\
                datetime.strftime(self.time[self.tstep[0]],'%d-%b-%Y %H:%M:%S'))
        
        mlab.figure(size=figsize)        
        self.fig,self.h,ug,d,title = unsurfm(points,self.cells,np.array(self.data[0,:]),clim=self.clim,title=titlestr)        
        
        if vector_overlay:
             # Add vectorss to the unctructured grid object
             #ug.cell_data.vectors=np.asarray((u[0,:],v[0,:],w[0,:])).T
             #ug.cell_data.vectors.name='vectors'
             #d.update()
             #h2=mlab.pipeline.vectors(d,color=(0,0,0),mask_points=1,scale_factor=1./scale)
             vec=mlab.pipeline.vector_scatter(self.xv, self.yv, self.yv*0, u[0,:], v[0,:], w[0,:])
             h2=mlab.pipeline.vectors(vec,color=(0,0,0),mask_points=subsample,scale_factor=1./scale,scale_mode='vector')

             
       # Animate the plot by updating the scalar data in the unstructured grid object      
#        for ii in range(nt):
#            print ii
#            # Refresh the unstructured grid object
#            ug.cell_data.scalars = self.data[ii,:]
#            ug.cell_data.scalars.name = 'suntans_scalar'
#
#            # Update the title
#            titlestr='%s [%s]\n Time: %s'%(self.long_name,self.units,\
#            datetime.strftime(self.time[self.tstep[ii]],'%d-%b-%Y %H:%M:%S'))
#            title.text=titlestr
#
#            self.fig.scene.render()
#            self.savefig('tmp_vtk_%00d.png'%ii)
#
#        mlab.show()
        
        @mlab.animate
        def anim():
            ii=-1
            while 1:
                if ii<nt-1:
                    ii+=1
                else:
                    ii=0
                ug.cell_data.scalars = self.data[ii,:]
                ug.cell_data.scalars.name = 'suntans_scalar'
                if vector_overlay:
                    #ug.cell_data.vectors=np.asarray((u[ii,:],v[ii,:],w[ii,:])).T
                    #ug.cell_data.vectors.name='vectors'
                    #d.update()
                    #vec=mlab.pipeline.vector_scatter(self.xv, self.yv, self.yv*0, u[ii,:], v[ii,:], w[ii,:])
                    #h2=mlab.pipeline.vectors(vec,color=(0,0,0),mask_points=subsample,scale_factor=1./scale,scale_mode='vector')
                    #h2.update_data()
                    h2.update_pipeline()
                    vectors=np.asarray((u[ii,:],v[ii,:],w[ii,:])).T
                    h2.mlab_source.set(vectors=vectors)
                    
                titlestr='%s [%s]\n Time: %s'%(self.long_name,self.units,\
                datetime.strftime(self.time[self.tstep[ii]],'%d-%b-%Y %H:%M:%S'))
                title.text=titlestr
                
                self.fig.scene.render()
                yield
        
        a = anim() # Starts the animation.

    def getVector(self):
        """
        Retrieve U and V vector components
        """
        tmpvar = self.variable
        
        u=self.loadDataRaw(variable='uc',setunits=False)
        
        v=self.loadDataRaw(variable='vc',setunits=False)
        
        try:
            w=self.loadDataRaw(variable='w',setunits=False)
        except:
            w=u*0.0
                               
        self.variable=tmpvar
        # Reload the original variable data
        #self.loadData()
        
        return u,v,w
        
    def getTstep(self,tstart,tend,timeformat='%Y%m%d.%H%M'):
        """
        Returns a vector of the time indices between tstart and tend
        
        tstart and tend can be string with format=timeformat ['%Y%m%d.%H%M' - default]
        
        Else tstart and tend can be datetime objects
        """
        
        try:
            t0 = datetime.strptime(tstart,timeformat)
            t1 = datetime.strptime(tend,timeformat)
        except:
            # Assume the time is already in datetime format
            t0 = tstart
            t1 = tend
                        
        n1 = othertime.findNearest(t0,self.time)
        n2 = othertime.findNearest(t1,self.time)
        
        if n1==n2:
            return [n1,n2]
        else:
            return list(range(n1,n2))
        
    def updateTstep(self):
        """
        Updates the tstep variable: -99 all steps, -1 last step
        """
        try:
            if self.tstep.any()==-99:
                self.tstep=np.arange(0,len(self.time))
            elif self.tstep.any()==-1:
                self.tstep=len(self.time)-1
        except:
            if self.tstep==-99:
                self.tstep=np.arange(0,len(self.time))
            elif self.tstep==-1:
                self.tstep=len(self.time)-1
       
            
    def checkIndex(self):
        """
        Ensure that the j property is a single or an array of integers
        """
        
        #shp = np.shape(self.j)
        #shp = self.j

        if isinstance(self.j[0], int):
            return self.j
        else:
            print('x/y coordinates input instead of cell index. Finding nearest neighbours.')
            dd, j = self.findNearest(self.j)
            print('Nearest cell: %d, xv[%d]: %6.10f, yv[%d]: %6.10f'%(j,j,self.xv[j],j,self.yv[j]))
            #self.j = j.copy()

        return j
     
    def hasVar(self,varname):
        """
        Tests if a variable exists in the file
        """
        try:
            self.nc.variables[varname][0]
            return True
        except:
            return False

    def hasDim(self,varname,dimname):
        """
        Tests if a variable contains a dimension
        """
        derivedvars  =[\
                      'speed',\
                      'vorticity',\
                      'div_H',\
                      'streamfunction',\
                      'KE',\
                      'PE',\
                      'buoyancy',\
                      'PEanom',\
                      'agemean',\
                      ]
        if varname in derivedvars:
            dimensions = ['Nt','Nc','Nk']
        else:
            dimensions = self.nc.variables[varname].dimensions
        
        return dimname in dimensions
        #try:
        #    self.nc.dimensions[dimname].__len__()
        #    return True
        #except:
        #    return False
     
    def calc_streamfunction(self):
        """
        Calculate the horizontal streamfunction
        """

        self.klayer = [-99] # Make sure this is set
        uc=self.loadDataRaw(variable='uc')
        vc=self.loadDataRaw(variable='vc')

        return self.depthint(uc) + self.depthint(vc)

    def calc_divergence(self):
        """
        Calculate the horizontal divergence
        """
        
        u,v,w = self.getVector()
            
        sz = u.shape
         
        if len(sz)==1: # One layer
            du_dx,du_dy = self.gradH(u,k=self.klayer[0])
            dv_dx,dv_dy = self.gradH(v,k=self.klayer[0])
            
            #data = du_dy + dv_dx
            data = du_dx + dv_dy
            
        else: # 3D
            data = np.zeros(sz)
            
            for k in self.klayer:
                du_dx,du_dy = self.gradH(u[:,k],k=k)
                dv_dx,dv_dy = self.gradH(v[:,k],k=k)
            
                data[:,k] = du_dy + dv_dx
                
        return data
    
    def vorticity_circ(self,k=0):
        """
        Calculate vertical vorticity component using the 
        circulation method
        """
        # Load the velocity
        u,v,w = self.getVector()
                             
        def _AverageAtFace(phi,jj,k):
            
            grad1 = self.grad[:,0]
            grad2 = self.grad[:,1]
            #Apply mask to jj
            jj[self.cellmask]=0
            nc1 = grad1[jj]
            nc2 = grad2[jj]
                    
            # check for edges (use logical indexing)
            ind1 = nc1==-1
            nc1[ind1]=nc2[ind1]
            ind2 = nc2==-1
            nc2[ind2]=nc1[ind2]
            
            # check depths (walls)
            indk = operator.or_(k>=self.Nk[nc1], k>=self.Nk[nc2])
            ind3 = operator.and_(indk, self.Nk[nc2]>self.Nk[nc1])
            nc1[ind3]=nc2[ind3]
            ind4 = operator.and_(indk, self.Nk[nc1]>self.Nk[nc2])
            nc2[ind4]=nc1[ind4]
            
            # Average the values at the face          
            return 0.5*(phi[nc1]+phi[nc2]) 
            
        # Calculate the edge u and v
        ne = self.face #edge-indices
        
        ue = _AverageAtFace(u,ne,k)
        ve = _AverageAtFace(v,ne,k)
        ue[self.cellmask]=0
        ve[self.cellmask]=0
        
        tx,ty,mag = self.calc_tangent()
        
        tx[self.cellmask]=0
        ty[self.cellmask]=0
        
        # Now calculate the vorticity
        return np.sum( (ue*tx + ve*ty )*mag,axis=-1)/self.Ac
        
    
        # Plot the result
#        fig = plt.gcf()
#        ax = fig.gca()
#        self.plotmesh(edgecolors='b')
#        scale = 1e-3
#        ind = 330
#        plt.quiver(self.xe[self.face[ind,:]],self.ye[self.face[ind,:]],tx[ind,:],ty[ind,:],scale=scale,scale_units='xy',color='r')
#        ax.set_aspect('equal')
#        plt.show()
#        
#        pdb.set_trace()

        
        
        
    def vorticity(self):
        """
        Calculate the vertical vorticity component
        
        Uses gradient method
        """
        
        u,v,w = self.getVector()
            
        sz = u.shape
         
        if len(sz)==1: # One layer
            du_dx,du_dy = self.gradH(u,k=self.klayer[0])
            dv_dx,dv_dy = self.gradH(v,k=self.klayer[0])
            
            data = dv_dx - du_dy
            
        else: # 3D
            data = np.zeros(sz)
            
            for k in self.klayer:
                du_dx,du_dy = self.gradH(u[:,k],k=k)
                dv_dx,dv_dy = self.gradH(v[:,k],k=k)
            
                data[:,k] = dv_dx - du_dy
                
        return data
        
    def agemean(self):
        """
        Calculates the mean age of a tracer

        This is already calculated in the time-averaged output files
        """
        sec2day = 1.0/86400.0

        if self.hasVar('agemean'):
            agemean = self.loadDataRaw(variable='agemean')
            return agemean*sec2day

        elif self.hasVar('agec') and self.hasVar('agealpha'):
            agec=self.loadDataRaw(variable='agec')
            agealpha=self.loadDataRaw(variable='agealpha')

        else:
            raise Exception("variables: 'agec' and 'agealpha' are not present in file.")


        mask = agec>=1e-10

        agemean = np.zeros_like(agec)

        agemean[mask] = agealpha[mask]/agec[mask]
            
        # Returns the mean age in days
        return agemean*sec2day

    def calc_buoyancy(self):
        """
        Returns the buoyancy of the fluid
        """

        # density is stored as rho' / rho_0
        rho=self.loadDataRaw(variable='rho')

        self.long_name='buoyancy'
        self.units='m s-2'
        return -GRAV*rho

    def calc_PE(self,b=None):
        """
        Calculate the potential energy of the fluid
        """
        if self.klayer[0]==-99:
            z = -self.z_r
        else:   
            z = -self.z_r[self.klayer]

        if b == None:
            b = self.calc_buoyancy()

        self.long_name = 'Potential energy'
        self.units = 'm2 s-2'

        return (b.swapaxes(0,1)*z).swapaxes(0,1)

    def calc_KE(self,u=None,v=None):
        """
        Calculate the kinetic energy
        """
        if u is None:
            u=self.loadDataRaw(variable='uc')
        if v is None:
            v=self.loadDataRaw(variable='vc')

        self.long_name = 'Kinetic energy'
        self.units = 'm2 s-2'
        return 0.5*(u*u + v*v)

    def calc_PEanom(self):
        """
        Calculate the potential energy anomaly as defined by Simpson et al., 1990
        """
        self.klayer=[-99]
        
        # Load the base data
        #self.loadDataRaw(variable='eta')
        #eta=self.data.copy()
        self.loadDataRaw(variable='rho')
        rho=self.data.copy()*1000.0+1000.0
        mask =rho>1500.0
      
        try:
            self.loadDataRaw(variable='dzz')
            dzz=self.data.copy()
            mask = dzz.mask.copy()
            dzz[mask]=0.
            h=np.sum(dzz,axis=0)
            #h = self.dv+eta
        except: # No dzz variable
            dz2=np.reshape(self.dz,(self.dz.shape[0],1))
            dzz = np.tile(dz2,(1,rho.shape[1]))
            h=self.dv
                
        rho[mask]=0.

        rho_bar = self.depthave(rho,dz=dzz,h=h)
        rho_bar = np.tile(rho_bar,(self.Nkmax,1))
        
        z = -1.0*np.cumsum(dzz,axis=0) 
        
        rho_pr = rho_bar - rho
        rho_pr[mask]=0.
        
        self.long_name='Potential Energy Anomaly'
        self.units = 'J m-2'
        
        return self.depthint(rho_pr*z*9.81,dz=dzz) 
#        return self.depthave(rho_pr*z*9.81,dz=dzz,h=h) 
    
    def getWind(self,nx=10,ny=10):
        """
        Interpolates the wind data onto a regular grid
        """
        from interpXYZ import interpXYZ
        
        x = np.linspace(self.xv.min(),self.xv.max(),nx)  
        y = np.linspace(self.yv.min(),self.yv.max(),ny) 
        
        X,Y = np.meshgrid(x,y)
        
        Uwind = self.loadDataRaw(variable='Uwind')
        Vwind = self.loadDataRaw(variable='Vwind')
        
        XYout = np.array((X.ravel(),Y.ravel())).T
        
        F = interpXYZ(np.array((self.xv,self.yv)).T, XYout, method='idw')
        
        return X,Y,F(Uwind).reshape((nx,ny)),F(Vwind).reshape((nx,ny))
        
        
    def depthave(self,data,dz=None, h=None):
        """ Calculate the depth average of a variable
        Variable should have dimension: [nz*nx] or [nt*nz*nx]
        
        """
        ndim = np.ndim(data)
        
        if h is None:
            h = self.dv
            
        if ndim == 2:
            return self.depthint(data,dz=dz) / h
        elif ndim == 3:
            nt = np.size(data,0)
            h = np.tile(h,(nt,1))
            return self.depthint(data,dz=dz) / h
                
    def depthint(self,data,dz=None,cumulative=False):
        """
        Calculates the depth integral of a variable: data
        Variable should have dimension: [nz*nx] or [nt*nz*nx]
        
        """
        ndim = np.ndim(data)
        
        nz = np.size(self.dz)
                
        if ndim == 1:
            return data*self.dv

        elif ndim == 2:
            nx = np.size(data,1)
            if dz is None:
                dz2=np.reshape(self.dz,(nz,1))
                dz2 = np.tile(dz2,(1,nx))
            else:
                dz2=dz
                
            if cumulative:
                return np.cumsum(data*dz2,axis=0)
            else:
                return np.sum(data*dz2,axis=0)
                
        elif ndim == 3:
            nt = np.size(data,0)
            nx = np.size(data,2)
            if dz is None:
                dz3 = np.reshape(self.dz,(1,nz,1))
                dz3 = np.tile(dz3,(nt,1,nx))
            else:
                dz3=dz
            if cumulative:
                return np.cumsum(data*dz3,axis=0)
            else:
                return np.sum(data*dz3,axis=0)

    def gradZ(self, data, zvar='z_r'):
        """
        Vertical gradient calculation on an unevenly spaced grid
        
        Variable data should have shape: [nz], [nz,nx] or [nt,nz,nx]
        """
                
        ndim = np.ndim(data)
        zin = getattr(self, zvar)

        if ndim == 1:
            axis=0
            z = -zin
        elif ndim == 2:
            axis=0
            z = -zin[:,None]
        elif ndim == 3:
            axis=1
            z = -zin[None, :, None]

        return grad_z(data, z, axis=axis)

        #rho = dsc['rho'][0:2,:,5:10].values
        #grad_z(rho, -z[None,:,None], axis=1) # 3D rho

        #rho = dsc['rho'][0,:,5].values
        #grad_z(rho, -z, axis=0) # 1D rho

        #rho = dsc['rho'][0,:,5:10].values
        #print(rho.shape)
        #grad_z(rho, -z[:,None], axis=0) # 2D rho

        ###
        # Old
        #if ndim == 1:
        #    phi = np.hstack((data[0],data,data[-1])) # nz+2
        #    phi_npm = (phi[1:]+phi[0:-1])*0.5 # nz+1
        #    return (phi_npm[0:-1] - phi_npm[1:])/self.dz
        #    
        #elif ndim == 2:
        #    nx = np.size(data,1)
        #    dz2=np.reshape(self.dz,(nz,1))
        #    dz2 = np.tile(dz2,(1,nx))
        #    phi = np.concatenate((data[[0],:],data,data[[-1],:]),axis=0) # nz+2
        #    phi_npm = (phi[1:,:]+phi[0:-1,:])*0.5 # nz+1
        #    dphi_dz = (phi_npm[0:-1,:] - phi_npm[1:,:])/dz2
        #    
        #    # Take care of the near bed gradient
        #    bedbot = np.ravel_multi_index((self.Nk,list(range(self.Nc))),(self.Nkmax,self.Nc))
        #    Nk1=self.Nk-1
        #    Nk1[Nk1<0]=0
        #    bedtop = np.ravel_multi_index((Nk1,list(range(self.Nc))),(self.Nkmax,self.Nc))
        #    bedbot = np.ravel_multi_index((self.Nk,list(range(self.Nc))),(self.Nkmax,self.Nc))
        #    dphi_dz.flat[bedbot] = (data.flat[bedtop]-data.flat[bedbot]) / dz2.flat[bedbot]
        #    
        #    return dphi_dz
        #
        #elif ndim == 3:
        #    nt = np.size(data,0)
        #    nx = np.size(data,2)
        #    dz3=np.reshape(self.dz,(1,nz,1))
        #    dz3 = np.tile(dz3,(nt,1,nx)) 
        #    phi = np.concatenate((data[:,[0],:],data[:,:,:],data[:,[-1],:]),axis=1) # nz+2
        #    phi_npm = (phi[:,1:,:]+phi[:,0:-1,:])*0.5 # nz+1
        #    
        #    return (phi_npm[:,0:-1,:] - phi_npm[:,1:,:])/dz3

    def areaint(self,phi,mask=None):
        """
        Calculate the area-integral of data at phi points 
        """
        if mask is None:#no mask
            mask = np.ones_like(phi)

        phiA = np.sum(phi*self.Ac*mask)
        area = np.sum(self.Ac*mask)
        
        return phiA, area
 
    def areaintold(self,phi,xypoly):
        """
        Calculate the area-integral of data at phi points 
        """
        try:
            mask = nxutils.points_inside_poly(np.array((self.xv,self.yv)).T, xypoly)
        except:
            import matplotlib.nxutils as nxutils #inpolygon equivalent lives here
            mask = nxutils.points_inside_poly(np.array((self.xv,self.yv)).T, xypoly)
  
        mask = nxutils.points_inside_poly(np.array((self.xv,self.yv)).T, xypoly)
        
        phiA = np.sum(phi[mask]*self.Ac[mask])
        area = np.sum(self.Ac[mask])
        
        return phiA, area
        
    def getdzz(self,eta):
        """
        Calculate the cell-centred vertical grid spacing based
        on the free surface height only
        """

        z = np.cumsum(self.dz)
        dzz = np.repeat(self.dz[:,np.newaxis],self.Nc,axis=1)

        # Find dzz of the top cell
        ctop = self.getctop(eta)
        #dztop = self.dz[ctop]+eta
        dztop = z[ctop]+eta
        dzz[ctop,list(range(self.Nc))]=dztop

        # Find dzz at the bottom
        dzbot = self.dv - z[self.Nk-1] 
        dzz[self.Nk,list(range(self.Nc))]=dzbot

        # Mask the cells
        Nk=self.Nk+1
        for ii in range(self.Nc):
            dzz[0:ctop[ii],ii]=0.0
            dzz[Nk[ii]::,ii]=0.0

        return dzz

    def getdzf(self,eta,U=None,method='max',j=None):
        """
        Calculate the edge-centred vertical grid spacing based
        on the free surface height only
        """

        etop,etaedge = self.getetop(eta,U=U,method=method,j=j)
        Ne = etop.shape[0]

        # Find dzz of the top cell
        dztop = self.dz[etop]+etaedge

        dzf = np.repeat(self.dz[:,np.newaxis],Ne,axis=1)

        dzf[etop,list(range(Ne))]=dztop

        # Mask the cells
        #mask = self.get_zmask(etop,self.Nke)
        for ii in range(Ne):
            dzf[0:etop[ii],ii]=0.0
            dzf[self.Nke[ii]::,ii]=0.0

        return dzf


    def getctop(self,eta):
        """
        Return the layer of the top cell
        """
        
        # Find the layer of the top cell
        ctop = np.searchsorted(self.z_w,-eta)
        ctop[ctop>0] -= 1
        return ctop

    def getetop(self,eta,method='max',U=None,j=None):
        """
        Return the layer of the top edge
        """
        if j is None:
            eta_edge = self.get_edgevar(eta,U=U,method=method)
        else:
            eta_edge = eta
        
        # Find the layer of the top cell
        etop = np.searchsorted(self.z_w,-eta_edge)
        etop[etop>0] -= 1
        return etop, eta_edge


    def get_zmask(self,ktop,nk):
        """
        Return a mask array (1/0) that is false (0) for points outside 
        of ktop and nk

        This is a good candidate for a cython function

        ##TBC##
        """
        sz = nk.shape
        mask = np.zeros(sz)

        nc = ktop.shape[0]
        cols = [list(range(ktop[ii],nk[ii])) for ii in range(nc)]


    def get_surfacevar(self,phi,ctop):
        """
        Retrieves the surface layer of a 3d array [Nk,Nc]
        """
        #assert phi.shape==(self.Nkmax,self.Nc),'size mismatch'
        #ctop = self.getctop(eta)
        if phi.ndim==1: # This happens if the surface is all contained in one
            #layer
            return phi
        else:
            j = list(range(self.Nc))
            return phi[ctop[j],j]

    def get_seabedvar(self,phi):
        """
        Retrieves the seabed layer of a 3d array [Nk,Nc]
        """
        assert phi.shape==(self.Nkmax,self.Nc),'size mismatch'

        j = list(range(self.Nc))
        return phi[self.Nk[j],j]

    def get_edgevar(self,phi,k=0,U=None,method='max'):
        """
        Return the edge value of a cell-based variable

        Method can be one of:
            'max' - maximum value either side
            'min' - minimum value either side
            'mean' - average value
            'linear' - linear interpolated value
            'upwind' - upwind value. Requires 'U'. 
        """
        nc1 = self.grad[:,0].copy()
        nc2 = self.grad[:,1].copy()
        Ne = self.Ne
                
        # check for edges (use logical indexing)
        ind1 = nc1==-1
        nc1[ind1]=nc2[ind1]
        ind2 = nc2==-1
        nc2[ind2]=nc1[ind2]

        
        # check depths (walls)
        #indk = operator.or_(k>=self.Nk[nc1], k>=self.Nk[nc2])
        indk = operator.or_(k>self.Nk[nc1], k>self.Nk[nc2]) # Nk is zero-based
        ind3 = operator.and_(indk, self.Nk[nc2]>self.Nk[nc1])
        nc1[ind3]=nc2[ind3]
        ind4 = operator.and_(indk, self.Nk[nc1]>self.Nk[nc2])
        nc2[ind4]=nc1[ind4]

        
        if method=='max':
            tmp = np.zeros((Ne,2),dtype=phi.dtype)
            tmp[:,0]=phi[nc1]
            tmp[:,1]=phi[nc2]
            return tmp.max(axis=-1)

        elif method=='min':
            tmp = np.zeros((Ne,2),dtype=phi.dtype)
            tmp[:,0]=phi[nc1]
            tmp[:,1]=phi[nc2]
            return tmp.min(axis=-1)

        elif method=='mean':
            # Average the values at the face          
            return 0.5*(phi[nc1]+phi[nc2]) 

        elif method=='upwind':
            if U is None:
                raise Exception('U must be set to use upwind method')

            ind = U>0
            nc1[ind]=nc2[ind]
            return phi[nc1]

        else:
            raise Exception('Method: %s not implemented.'%method)

    def genTitle(self,tt=None):
        
        if tt  is None:
            if type(self.tstep)==int:
                tt = self.tstep
            else:
                tt = self.tstep[0]
            
        if self.klayer[0] in [-1,'seabed']:
            zlayer = 'seabed'
        elif self.klayer[0] =='surface':
            zlayer = 'surface'
        elif self.klayer[0]==-99:
            zlayer = 'depth-averaged'
        elif self.klayer[0]>=0:
            zlayer = '%3.1f [m]'%self.z_r[self.klayer[0]]

        if 'time' in self.__dict__:
            tstr = datetime.strftime(self.time[tt],\
                '%Y-%m-%d %H:%M:%S')
            tstr = 'Time: %s'%tstr
        else:
            tstr = ''
        titlestr='%s [%s]\n z: %s, %s'%(self.long_name,self.units,zlayer,tstr )
                
        return titlestr

class TimeSeries(timeseries, Spatial):
    """
    Time series class for suntans output
    """    
    
    def __init__(self,ncfile,XY,Z=None,klayer=None,**kwargs):
        
        Spatial.__init__(self,ncfile,**kwargs)
        
        self.XY = XY
        self.klayer = klayer
        if not Z is None:
            self.Z = np.abs(Z)
        else:
            self.Z=Z
        self.tstep = list(range(0,len(self.time))) # Load all time steps
        
        self.update()

    def __call__(self, XY, Z=None, varname=None):
        """
        Returns an updated version of the object
        """

        if not varname is None:
            self.variable = varname


        if not Z is None:
            self.Z = Z

        self.XY = XY

        self.update()

        # Convert the suntans-TimeSeries object to an xray.DataArray
        coords = {'time':self.t.copy()}
        if self.y.ndim == 1:
            dims = ('time',)
        else:
            dims = ('depth','time')
            coords.update({'depth':self.Z})

        return xray.DataArray( self.y.copy(), \
                    dims=dims,\
                    name=self.variable,\
                    #attrs = nc[varname].attrs,\
                    coords = coords,\
               )
        
    def update(self):
        
        # Update the j position
        dist, j = self.find_nearest(self.XY)
        self.j = j
        
        # Find the klayer
        if isinstance(self.Z,type(np.array(()))):
            self.klayer=[]
            for zz in self.Z.tolist():
                k0=self.findNearestZ(zz)
                self.klayer.append(k0[0])

        else:
            if self.Z == None and self.klayer== None:
                self.klayer = [-99]
            elif not self.Z == None or self.klayer == None:
                self.klayer = self.findNearestZ(self.Z)
            else:
                self.klayer=self.klayer

        # Loads in a time series object                     
        data = self.loadDataRaw().reshape((self.Nt, len(self.klayer)))
        timeseries.__init__(self, self.time[self.tstep], data)
        
    def contourf(self,clevs=20,**kwargs):
        """
        Filled contour plot
        """
        if self.klayer[0]==-99:
            z = -self.z_r
        else:
            z = -self.z_r[self.klayer]

        h1=plt.contourf(self.time,z,self.y,clevs,**kwargs)
        plt.xticks(rotation=17)        
        return h1
        
        
    def findNearestZ(self,Z):
        
        dist = np.abs(Z-self.z_r)
        
        return np.where(dist == dist.min())[0].tolist()
    
    def __setitem__(self,key,value):
        
        if key == 'variable':
            self.variable=value
            self.update()
            
        elif key == 'XY':
            self.XY = value
            self.update()   
        elif key == 'Z':
            self.Z = np.abs(value)
            self.update()
        else:
            self.__dict__[key]=value
        
        
                  
class Profile(object):
    """
        Class for handling SUNTANS profile netcdf files
    """        
    
    def __init__(self,ncfile ,**kwargs):
        
       
       self.ncfile = ncfile
       
       self.__loadMeta()
       
       # Set defaults
       self.indices = np.arange(0,self.Np)
       self.tstep = 0 # -99 all steps, -1 last step
       self.klayer = np.arange(0,self.Nz)
       self.variable = 'u'
       self.xcoord = 'xp' # 'yp', 'time' or 'dist'
       self.ycoord = 'z'
       self.clim = None
       self.clevels = 12 # Number of contour levels
       
       # Linear EOS for files with no rho
       self.beta=1.0
       self.S0 = -1.0271
       
       # Distance calculation stuff
       self.smoothdist=True # "smooth" the transect by taking out jagged bits
       self.nsmooth = 10

       
       # Update user-defined properties
       self.__dict__.update(kwargs)
        
       # Update tstep 
       self.__updateTstep()

    
#    def __setattr__(self, name, value):
#        """
#        Call on other methods when certain attributes are set
#        """
#
#        self.__dict__[name] = value
#        
#        if name in ['xplot','yplot']:
#            self.__loadXY()
            
    def __loadMeta(self):
        """
        Loads the metadata from the profile netcdf file
        """
        #nc = Dataset(self.ncfile, 'r', format='NETCDF4') 
        try: 
            nc = MFDataset(self.ncfile, 'r')
        except:
            nc = Dataset(self.ncfile, 'r')
        
        # Note that some of these variable names may change
        try:
            self.xp = nc.variables['x'][:]
        except:
            self.xp = nc.variables['xv'][:]
        try:    
            self.yp = nc.variables['y'][:]
        except:
            self.yp = nc.variables['yv'][:]
        try:
            self.dv = nc.variables['h'][:]
        except:
            self.dv = nc.variables['dv'][:]
        try:
            self.dz = nc.variables['vertspace'][:]
        except:
            self.dz = nc.variables['dz'][:]
        try:
            self.z =  nc.variables['vertdepth'][:]
        except:
            self.z =  nc.variables['z_r'][:]
        try:
            self.Nk = nc.variables['klayers'][:]
        except:
            self.Nk = nc.variables['Nk'][:]
            
        self.Np = len(self.xp)
        self.Nz = len(self.dz)
        
        self.xlims = [self.xp.min(),self.xp.max()]
        self.ylims = [self.yp.min(),self.yp.max()]
        
        try:
            t = nc.variables['suntime']
        except:
            t = nc.variables['time']
        self.time = num2pydate(t[:],t.units)
        
        #print nc.variables.keys()
        nc.close()
        
    def loadData(self):
        """
        Loads the actual data for the given variable, indices, tsteps and zlayers
        """ 
        
        # Open the dataset    
        try: 
            nc = MFDataset(self.ncfile, 'r')
        except:
            nc = Dataset(self.ncfile, 'r')
                  
        # "Higher order" variable stuff
        tmpvar = self.variable
        if tmpvar == 'ubar':
            self.variable = 'u'
        if tmpvar == 'vbar':
            self.variable = 'v'
        if tmpvar in  ['rho','bvf2']:
            if 'rho' not in nc.variables:
                self.variable='S'
            else:
                self.variable='rho'
            
        
        # Load the data
        self.long_name = nc.variables[self.variable].long_name
        self.units= nc.variables[self.variable].units
        #        ndims = len(nc.variables[self.variable].dimensions)
        #print nc.variables[self.variable].dimensions
        ndim = nc.variables[self.variable].ndim
        if ndim==1:
            self.data=nc.variables[self.variable][self.indices]
        elif ndim==2:
            self.data=nc.variables[self.variable][self.tstep,self.indices]
        else:
            self.data=nc.variables[self.variable][self.tstep,self.klayer,self.indices]
    
        nc.close()
        
        # Calculate the higher order variables from the raw data
        # To Add:
        #   uprime, vprime (baroclinic velocity)
        #   time mean variables
        #   seabed values of variables
        #   hydrostatic pressure perturbation?
        #   buoyancy frequency        
        if tmpvar in ['ubar','vbar']:
            self.data=depthave(self.data,self.dz,np.abs(self.dv[self.indices]))
        if tmpvar in ['rho','bvf2'] and 'rho' not in nc.variables:
            self.data = linearEOS(self.data,S0=self.S0,beta=self.beta)
        if tmpvar == 'bvf2':
            self.data = calcN2(self.data,self.dz)
            
        self.variable = tmpvar
        
        
        # Update the colorbar limits if not set
        if self.clim is None:
            self.clim=[np.min(self.data),np.max(self.data)]
        
    
    def __loadXY(self):
        """
        Loads the X-Y coordinates used by 2-D plotting commands
        """
        if self.xcoord=='xp':
            self.xplot = self.xp[self.indices]
        elif self.xcoord=='yp':
            self.xplot = self.yp[self.indices]
        elif self.xcoord=='dist':
            # Calculate the distance along the transect
            self.xplot=self.calcDistAxis()            
        elif self.xcoord=='time':
            self.xplot = self.time[self.tstep]
            
        if self.ycoord=='z':
            self.yplot = self.z[self.klayer]
        elif self.ycoord=='time':
            self.yplot = self.time[self.tstep]
        elif self.ycoord=='dist':
            # Calculate the distance along the transect
            self.yplot=self.calcDistAxis()     
    
    def __updateTstep(self):
        """
        Updates the tstep variable: -99 all steps, -1 last step
        """
        try:
            if self.tstep.any()==-99:
                self.tstep=np.arange(0,len(self.time))
            elif self.tstep.any()==-1:
                self.tstep=len(self.time)-1
        except:
            if self.tstep==-99:
                self.tstep=np.arange(0,len(self.time))
            elif self.tstep==-1:
                self.tstep=len(self.time)-1
    def __checkDims(self):
        """
        Check that the dimensions sizes match for plotting
        
        If not transpose the data
        """        
        rc = np.shape(self.data)
        nx = len(self.xplot)
        ny = len(self.yplot)
        
        if ny!=rc[0] or nx !=rc[1]:
            self.data=np.transpose(self.data)
    
    
            
    def plotIndices(self):
        """
        Plots the locations of the points with the index numbers overlaid
        """
        offset=20
        plt.figure()
        plt.plot(self.xp,self.yp,'b.')
        plt.plot(self.xp[self.indices],self.yp[self.indices],'o',markeredgecolor='r',markerfacecolor=None)
        for s in range(self.Np):
            plt.text(self.xp[s]+offset,self.yp[s]+offset,'%d'%s)
            
        plt.axis('equal')
        plt.show()
    
    def closetTime(self,t):
        """
        Find the index of the closest time to the datetime object "t"
        """
        dtall = []
        for tt in self.time:
            dt = tt - t
            dtall.append(np.abs(dt.total_seconds()))
            
        dtall = np.asarray(dtall)

        return np.argwhere(dtall == dtall.min())
    
    def calcDistAxis(self):
        """
        Calculates distance along the transect
        """
        
        print('Setting x-axis to distance...')
        
        x = self.xp[self.indices]
        y = self.yp[self.indices]
        nx = len(x)
        
        if self.smoothdist:
            from scipy import interpolate
            F = interpolate.UnivariateSpline(x[1:-1:self.nsmooth],y[1:-1:self.nsmooth])
            xnew = np.linspace(x[0],x[-1],nx)
            ynew = F(xnew)
            x=xnew
            y=ynew
            
        dxdy = np.sqrt( (x[1:]-x[0:-1])**2 + (y[1:]-y[0:-1])**2 )
        dxdy = np.concatenate(([0.0],dxdy))
        return np.cumsum(dxdy)
        
        
        
    def pcolor(self,data=None,**kwargs):
        """
        Pcolor plot of the given variable (doesn't like time variable)
        """     
        if 'xplot' not in self.__dict__:
            self.__loadXY()
        if 'data' not in self.__dict__:
            self.loadData()
            self.__checkDims()  
        if data == None:
            data=self.data
          
        #plt.ion()
        self.fig=plt.gcf().s
        self.ax =plt.gca()
        self.h = plt.pcolor(self.xplot,self.yplot,data,**kwargs)
        self.cb = plt.colorbar()
        
        
    def contourf(self,data=None,V=None,**kwargs):
        """
        Filled contour plot of the given variable
        """
        if 'xplot' not in self.__dict__:
            self.__loadXY()
        if 'data' not in self.__dict__:
            self.loadData()
            self.__checkDims()  
        if data == None:
            data=self.data
        if V == None:
            V = np.linspace(self.clim[0],self.clim[1],num=self.clevels)    

        #plt.ion()
        self.fig=plt.gcf()
        self.ax =plt.gca()
        self.h = plt.contourf(self.xplot,self.yplot,data,V,**kwargs)
        #self.cb = plt.colorbar()  
        
        
    def contour(self,data=None,V=None,**kwargs):
        """
        Contour plot of the given variable
        """
        if 'xplot' not in self.__dict__:
            self.__loadXY()
        if 'data' not in self.__dict__:
            self.loadData()
            self.__checkDims()  
        if data == None:
            data=self.data
        if V == None:
            V = np.linspace(self.clim[0],self.clim[1],num=self.clevels)
            
            
        #plt.ion()
        self.fig=plt.gcf()
        self.ax =plt.gca()
        self.h = plt.contour(self.xplot,self.yplot,data,V,**kwargs)
        #self.cb = plt.colorbar()  


    def savefig(self,outfile,dpi=150):
        """ Saves a figure to file (matplotlib only)"""
        
        self.fig.savefig(outfile,dpi=dpi)
         
        print('SUNTANS image saved to file:%s'%outfile)
        
    def animate(self,fig=None,ax=None,h=None,cb=None,tsteps=None):
        """
        Method to animate the current method for multiple time steps
        """
         
        if tsteps == None:
            tsteps = np.arange(self.tstep,len(self.time))
            
        
        
        # Set the tsteps and grab the data
        tstepold = self.tstep
        self.tstep = tsteps
        self.loadData()
        
        if (ax is None) or (h  is None) or (fig is None):
            fig, ax, h, cb = self.pcolor(data=np.squeeze(self.data[0,:,:]))
         
        
        # Update the image
        ax.set_animated(True)
        h.set_clim(vmin=self.clim[0],vmax=self.clim[1])
        cb.set_clim(vmin=self.clim[0],vmax=self.clim[1])
        cb.draw_all()
        for ii in tsteps:
            zdata = np.squeeze(self.data[ii,:,:])
            h.set_array(np.ravel(zdata))
            ax.set_title('%d'%ii)
            fig.canvas.draw()
#        def updateanim(ii):
#            zdata = np.squeeze(self.data[ii,:,:])
#            h.set_array(np.ravel(zdata))
#            ax.set_title('%d'%ii)
#            
#        ani = animation.FuncAnimation(fig, updateanim, tsteps, interval=10)                
        
        # Set the data array to pre-animation array
        self.tstep = tstepold
        self.loadData()
        
    
    

####################################################################
#
# General functions to be used by all classes
#
####################################################################        
def closePoly(x,y):

    """ 
    Returns an Nx2 closed polygon for vectors x and y
        
    This output is required for plotting by unsurf. 
    """
    
    nx = len(x)
    ny = len(y)
    if nx != ny:
        print("Error: The lengths of vector x and y must be equal")
        return

    x = np.reshape(x,(nx,1))
    y = np.reshape(y,(ny,1))

    x = np.vstack((x,x[0]))
    y = np.vstack((y,y[0]))
    
    return np.hstack((x,y))


def linearEOS(S,S0=1.0271,beta=1.0,RHO0=1000.0):
    """
    Linear equation of state
    
    Returns density from salinity and/or temperature
    """    

    return RHO0 * ( beta * (S-S0) )
        
def calcN2(rho,dz):
    """
    Calculate the buoyancy frequency squared
    """
    g=9.81
    rho0=1024
    return   - g/rho0 * gradZ(rho,dz) 
    
def readTXT(fname,sep=None):
    """
    Reads a txt file into an array of floats
    """
    
    fp = open(fname,'rt')
    data = np.array([list(map(float,line.split(sep))) for line in fp])
    fp.close()
    
    return data

def unsurfm(points, cells, z,clim=None,title=None,**kwargs):
    """
    Plot cell-centred data using the mayavi/tvtk libraries
    
    """        
    if clim is None:
        clim=[]
        clim.append(np.min(z))
        clim.append(np.max(z))
    
    try:    
        tri_type = tvtk.Triangle().cell_type
    except:
        # Load tvtk libraries here as they slow the script down
        from tvtk.api import tvtk
        from mayavi import mlab
        tri_type = tvtk.Triangle().cell_type
        
    ug = tvtk.UnstructuredGrid(points=points)
    ug.set_cells(tri_type, cells)
    
    ug.cell_data.scalars = z
    ug.cell_data.scalars.name = 'suntans_scalar'
    
    d = mlab.pipeline.add_dataset(ug)
    h=mlab.pipeline.surface(d,vmin=clim[0],vmax=clim[1],**kwargs)
    f=mlab.gcf()
    f.scene.background = (0.,0.,0.)
    mlab.colorbar(object=h,orientation='vertical')
    mlab.view(0,0)
    title=mlab.title(title,height=0.95,size=0.15)
    
    return f, h, ug, d, title
    
def unsurf(xy,z,xlim=[0,1],ylim=[0,1],clim=None,colorbar=True,**kwargs):
    """
    Plot cell-centred data on an unstructured grid as a series of patches
        
    Similar functionality to the suntans matlab unsurf function.
        
    Inputs:
        xy -list[Nc] of N*2 arrays of closed polygons
        z - scalar vector[Nc]
            
    """     
    
    # Create the figure and axes handles
    plt.ioff()
    fig = plt.gcf()
    ax = fig.gca()
    

    # Find the colorbar limits if unspecified
    if clim is None:
        clim=[]
        clim.append(np.min(z))
        clim.append(np.max(z))
    
    collection = PolyCollection(xy,**kwargs)
    #collection.set_array(np.array(z))
    collection.set_array(z)
    collection.set_clim(vmin=clim[0],vmax=clim[1])
    #collection.set_linewidth(0)
    #collection.set_edgecolors(collection.to_rgba(np.array(z)))    
    collection.set_edgecolors(collection.to_rgba(z))    
    
    ax.add_collection(collection)

    ax.set_aspect('equal')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    if colorbar:
        axcb = fig.colorbar(collection)
    else:
        axcb = None
    
    return fig, ax, collection, axcb

def edgeplot(xylines,edata,xlim=[0,1],ylim=[0,1],clim=None,**kwargs):
    """
    Plot edge data as a series of colored lines

    Inputs:
        xylines - list with two Ne x 2 arrays containing x/y end points
        edata - Ne vector with edge values
    """
    # Create the figure and axes handles
    plt.ioff()
    fig = plt.gcf()
    ax = fig.gca()
    

    # Find the colorbar limits if unspecified
    if clim is None:
        clim=[]
        clim.append(np.min(z))
        clim.append(np.max(z))
 
    # Create the inputs needed by line collection
    Ne = xylines[0].shape[0]

    # Put this into the format needed by LineCollection
    linesc = [list(zip(xylines[0][ii,:],xylines[1][ii,:])) for ii in range(Ne)]

    collection = LineCollection(linesc,array=edata,**kwargs)

    collection.set_clim(vmin=clim[0],vmax=clim[1])
    
    ax.add_collection(collection)

    ax.set_aspect('equal')
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    axcb = fig.colorbar(collection)
    
    return fig, ax, collection, axcb

def unanimate(xy,z,tsteps,xlim=[0,1],ylim=[0,1],clim=None,**kwargs):
    """
    Plot cell-centred data on an unstructured grid as a series of patches
        
    Similar functionality to the suntans matlab unsurf function.
        
    Inputs:
        xy -list[Nc] of N*2 arrays of closed polygons
        z - scalar array [nt x Nc]
        t - vector of time steps
            
    """     
    from matplotlib import animation
    
    # Create the figure and axes handles
    #plt.ion()
    fig = plt.gcf()
    ax = fig.gca()
    #ax.set_animated('True')
    
    # Find the colorbar limits if unspecified
    if clim is None:
        clim=[]
        clim.append(np.min(z))
        clim.append(np.max(z))
    

        
    collection = PolyCollection(xy)
    collection.set_array(np.array(z[0,:]))
    collection.set_clim(vmin=clim[0],vmax=clim[1])
    ax.add_collection(collection)    

    
    ax.axis('equal')
    ax.set_aspect('equal',fixLimits=True)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    #ax.axis([xlim[0],xlim[1],ylim[0],ylim[1]])
    
    #plt.axes().set_aspect('equal')
    title=ax.set_title("")
    fig.colorbar(collection)
    
    def init():
        collection.set_array([])
        title.set_text("")
        return (collection,title)

        
    def updateScalar(i):
        ts = i
        collection.set_array(np.array(z[i,:]))
        collection.set_edgecolors(collection.to_rgba(np.array((z[i,:]))))    
        title.set_text('%d'%ts)
        return (title,collection)

    
    anim = animation.FuncAnimation(fig, updateScalar, frames=200, interval=50, blit=True)
    return anim
#    print 'Building animation sequence...'
#    anim.save('C:/Projects/GOMGalveston/CODE/PYTHON/SUNTANS/test_animation.mp4', fps=15)


def usage():
    """
    Command line usage output
    """
    print("--------------------------------------------------------------")
    print("sunpy.py   -h                 # show this help message      ")
    print("         -f suntans.nc        # SUNTANS output netcdf file  ")
    print("         -v varname           # variable name to plot [default: u]      ")
    print("         -k  N                # vertical layer to plot [default: 0]")
    print("         -t  N                # time step to plot. -1 = last step. [default: 0]")
    print("         -j  N                # grid cell index to plot (timeseries only) [default: 0]")
    print('         -c "N N"             # Color bar limits !! IN DOUBLE QUOTES !! [default: None]')
    print("         -s figure.png        # Save to a figure")
    print("         --animate            # Animate plot through all time steps")
    print("         --timeseries         # Plot time series of individual point")
    print("         --profile            # time-depth contour plot at cell: j")
    print("         --vtk                # Use the vtk plotting libraries")
    print("\n Example Usage:")
    print("--------")
    print(" python sunpy.py -f suntans.nc -v temp -k 5 -t 10 -c '10 29' ")
    print("")

    
if __name__ == '__main__':
    """
    Command line call to sunpy
    """        
    
    # Defaults
    k = [0]
    t = 0
    j = 0
    varname = 'eta'
    plottype = 0
    usevtk = False
    clim = None
    save = False
    
    try:
            opts,rest = getopt.getopt(sys.argv[1:],'hf:v:t:k:j:c:s:',
                                      ['animate',
                                       'timeseries',
                                       'profile',
                                       'vtk'])
    except getopt.GetoptError as e:
        print(e)
        print("-"*80)
        usage()
        exit(1)

    for opt,val in opts:
        if opt == '-h':
            usage()
            exit(1)
        elif opt == '-f':
            ncfile = str(val)
            print(ncfile)
        elif opt == '-v':
            varname = val
        elif opt == '-t':
            t = int(val)
        elif opt == '-k':
            k = [int(val)]
        elif opt == '-j':
            j = int(val)
        elif opt == '-c':
            clim = np.asarray(val.split(),dtype=float)
        elif opt == '-s':
            outfile = val
            save = True
        elif opt == '--animate':
            plottype = 1
        elif opt == '--timeseries':
            plottype = 2
        elif opt == '--profile':
            plottype = 3
        elif opt == '--vtk':
            usevtk = True
            
    # Load the class and plot
    if plottype in [0,1,2]:
        sun = Spatial(ncfile,klayer=k,tstep=t,variable=varname,clim=clim)
    
    if plottype == 0:
        # Spatial Plot
        if usevtk:
            #mlab.figure(size=(640,480))
            sun.plotvtk()
            if save:
                sun.savefig(outfile)
            #mlab.show()

        else:
            plt.figure(figsize=(10,7))
            sun.plot() 
            if save:
                sun.savefig(outfile)
            #plt.show()
        
    elif plottype == 1:
        # Animation
        sun.tstep=np.arange(0,len(sun.time))
        sun.loadData()
        plt.figure(figsize=(10,7))
        sun.animate()
        if save:
            sun.saveanim(outfile)
    
    elif plottype == 2:
        # Time series plot
        plt.figure(figsize=(10,7))
        sun.plotTS(j=j)
        if save:
            sun.savefig(outfile)
        #plt.show()
    
    elif plottype == 3:
        sun = Profile(ncfile,indices=j,tstep=-99,xcoord='time',variable=varname,clim=clim)

        plt.figure(figsize=(10,5))
        sun.contourf()
        titlestr='%s [%s]'%(sun.long_name,sun.units)
        plt.title(titlestr)
        sun.ax.set_xlabel('Time')
        sun.ax.set_ylabel('z [m]')

        if save:
            sun.savefig(outfile)


############################
# Testing stuff -spatial plots and animation
#ncfile = 'C:/Projects/GOMGalveston/MODELLING/GalvestonSquare/rundata/suntan_output.nc.0'
#sun = Spatial(ncfile,klayer=0,tstep=20)
#plottype = 1 # 0 - spatial, 1 - animation, 2 - time series
#sun.variable='Tair'
#
#if plottype == 0:
#    # Spatial Plot
##    plt.figure(figsize=(10,7))
##    sun.plot()
##    plt.show()
#    
#    mlab.figure(size=(640,480))
#    sun.plotvtk()
#    
#    sun.savefig('C:/Projects/GOMGalveston/MODELLING/GalvestonSquare/rundata/test.png')
#    
#elif plottype == 1:
#    # Animation
##    sun.tstep=np.arange(0,len(sun.time))
##    sun.loadData()
##    ani=unanimate(sun.xy,sun.data,sun.tstep,xlim=sun.xlims,ylim=sun.ylims)
#    sun.animateVTK()
#
#elif plottype == 2:
#    # Time series plot
#    sun.plotTS(j=80)
    

############################
# Profile testing
#ncfile =  'E:/Projects/ScottReef/MODELLING/SUNTANS/netcdf/ScottRfTVD3_prof.nc';
#
#sun = Profile(ncfile,indices=np.arange(26,528),tstep=100,variable='ubar')
#sun.loadData()
#u=sun.data
#fig, ax, h, cb 2= sun.contourf()
#sun.plotIndices()
