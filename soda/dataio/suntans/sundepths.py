# -*- coding: utf-8 -*-
"""
SUNTANS bathymetry interpolation tools

Created on Fri Oct 05 11:24:10 2012

@author: mrayson
"""
import numpy as np
from . import sunpy
import matplotlib.pyplot as plt

from .sunpy import Grid
from soda.dataio.ugrid.gridsearch import GridSearch
from soda.utils.interpXYZ import Inputs, interpXYZ
from soda.utils.ufilter import ufilter
from soda.dataio.conversion.dem import DEM

import time

import pdb

class DepthDriver(object):
    """
    Driver class for interpolating depth data onto a suntans grid
    """
    
    # Interpolation method
    interpmethod='idw' # 'nn', 'idw', 'kriging', 'griddata'
    
    # Type of plot
    plottype='mpl' # 'mpl', 'vtk2' or 'vtk3'
    
    # Interpolation options
    NNear=3

    p = 1.0 #  power for inverse distance weighting
    # kriging options
    varmodel = 'spherical'
    nugget = 0.1
    sill = 0.8
    vrange = 250.0
    
    # Projection conversion info for input data
    clip=False
    convert2utm=False
    CS='NAD83'
    utmzone=15
    isnorth=True
    projstr=None
    vdatum = 'MSL'
    shapefieldname='contour'
    
    # Smoothing options
    smooth=False
    smoothmethod='kriging' # USe kriging, idw or gaussian for smoothing
    smoothnear=5 # No. of points to use for smoothing
    smoothrange = 50.

    # Set if the input data is a DEM
    isDEM = False
    
    def __init__(self,depthfile,**kwargs):
        
        self.__dict__.update(kwargs)
        
        # Parse the depth data into an object
        if self.isDEM:
            self.indata = DEM(depthfile, convert2utm=self.convert2utm,\
                utmzone=self.utmzone, projstr=self.projstr, isnorth=self.isnorth,\
            )
        else:
            self.indata = Inputs(depthfile, convert2utm=self.convert2utm,\
                CS=self.CS, utmzone=self.utmzone, projstr=self.projstr,\
                isnorth=self.isnorth, vdatum=self.vdatum,\
                shapefieldname=self.shapefieldname)


    def __call__(self, suntanspath, depthmax=0.0, depthmin=1e4,\
        scalefac=-1.0, interpnodes=True):
        
        self.suntanspath=suntanspath
        self.interpnodes = interpnodes
        self.scalefac = scalefac
        
        # Initialise the interpolation points
        print('Loading suntans grid points...')
        self.grd = sunpy.Grid(self.suntanspath)


        if self.interpnodes:
            print('Interpolating depths onto nodes and taking min...')
            self.xy = np.column_stack((self.grd.xp,self.grd.yp))


        else:
            print('Interpolating depths straight to cell centres...')
            self.xy = np.column_stack((self.grd.xv,self.grd.yv))

        # Call the interpolation routine
        self.interp_depths()
           

        # Smooth
        if self.smooth:
            self.smoothDepths()
        
        # Cap the maximum and minimum depth
        ind = self.grd.dv<=depthmax
        self.grd.dv[ind]=depthmax

        ind = self.grd.dv>=depthmin
        self.grd.dv[ind]=depthmin

        
        # Write the depths to file
        print('Writing depths.dat...')
        self.grd.saveBathy(suntanspath+'/depths.dat-voro')
        print('Data saved to %s.'%suntanspath+'/depths.dat-voro')
        
        # Plot
        if self.plottype=='mpl':
            self.plot()
        elif self.plottype=='vtk2':
            self.plotvtk()
        elif self.plottype=='vtk3':
            self.plotvtk3D()
        elif self.plottype==None:
            print('Not plotting')
            
        print('Finished depth interpolation.')
        
    def interp_depths(self):

        # Initialise the Interpolation class
        if not self.isDEM:
            print('Building interpolant class...')
            self.F = interpXYZ(self.indata.XY,\
                self.xy,\
                clip=self.clip,\
                method=self.interpmethod,\
                NNear=self.NNear,\
                p=self.p,\
                varmodel=self.varmodel,\
                nugget=self.nugget,\
                sill=self.sill,\
                vrange=self.vrange)

            # Interpolate the data
            print('Interpolating data...')
            dv = self.F(self.indata.Zin)*self.scalefac

        else:
            print('Interpolating DEM data...')
            dv = self.indata.interp(self.xy[:,0],self.xy[:,1])*self.scalefac


        if self.interpnodes:
            self.grd.dv = np.zeros_like(self.grd.xv)
            for nn in range(self.grd.Nc):
                self.grd.dv[nn] = np.max(dv[self.grd.cells[nn,0:self.grd.nfaces[nn] ] ])
                #self.grd.dv[nn] = np.mean(dv[self.grd.cells[nn,0:self.grd.nfaces[nn] ] ])
                
        else:
            self.grd.dv = dv
 
    def smoothDepths(self):
        """ 
        Smooth the data by running an interpolant over the model grid points
        """
        print('Smoothing the data...')
        if self.smoothmethod=='gaussian':
            print('\tUsing Gaussian filter...')
            Fsmooth = ufilter(self.xy, self.smoothrange)

        else:
            Fsmooth =interpXYZ(self.xy,\
                self.xy,\
                method=self.smoothmethod,\
                NNear=self.smoothnear,\
                vrange=self.smoothrange)

        self.grd.dv = Fsmooth(self.grd.dv)

    def plot(self):
        """
        Plot using matplotlib
        """
        fig=plt.figure()
        self.grd.plot(cmap=plt.cm.gist_earth)
        outfile = self.suntanspath+'/depths.png'
        fig.savefig(outfile,dpi=150)
        print('Figure saved to %s.'%outfile)
        
    def plotvtk(self):
        """
        2D plot using the vtk libraries
        """
        self.grd.plotvtk()
        
        outfile = self.suntanspath+'/depths.png'
        self.grd.fig.scene.save(outfile)      
        print('Figure saved to %s.'%outfile)
        
    def plotvtk3D(self):
        """
        3D plot using the vtk libraries
        """
        from tvtk.api import tvtk
        from mayavi import mlab
        # Plot the data on a 3D grid
        xy = np.column_stack((self.grd.xp,self.grd.yp))
        dvp = self.F(xy)
        vertexag = 50.0
        
        points = np.column_stack((self.grd.xp,self.grd.yp,dvp*vertexag))
        tri_type = tvtk.Triangle().cell_type
        #tet_type = tvtk.Tetra().cell_type
        ug = tvtk.UnstructuredGrid(points=points)
        ug.set_cells(tri_type, self.grd.cells)
        
        ug.cell_data.scalars = self.grd.dv
        ug.cell_data.scalars.name = 'depths'
        
        f=mlab.gcf()
        f.scene.background = (0.,0.,0.)
        d = mlab.pipeline.add_dataset(ug)
        h=mlab.pipeline.surface(d,colormap='gist_earth')
        mlab.colorbar(object=h,orientation='vertical')
        mlab.view(0,0)
        
        outfile = self.suntanspath+'/depths.png'
        f.scene.save(outfile)      
        print('Figure saved to %s.'%outfile)
        
        #mlab.show()


class AverageDepth(Grid):
    """
    Returns the average of all depths inside each cells
    """
    # Projection conversion info for input data
    convert2utm=False
    CS='NAD83'
    utmzone=15
    isnorth=True
    vdatum = 'MSL'

    def __init__(self,suntanspath,**kwargs):
        
        self.__dict__.update(kwargs)
        Grid.__init__(self,suntanspath)
        
        # Initialise the trisearch object
        self.tsearch =  GridSearch(self.xp,self.yp,self.cells)
        
    def __call__(self,depthfile,**kwargs):
        
        self.__dict__.update(kwargs)
        
        # Parse the depth data into an object
        self.indata = Inputs(depthfile,convert2utm=False,CS=self.CS,utmzone=self.utmzone,\
            isnorth=self.isnorth,vdatum=self.vdatum)
            
        tic = time.clock()
        print('Performing triangle search...')
        cells = self.tsearch(self.indata.XY[:,0],self.indata.XY[:,1])
        
        toc = time.clock()
        print('Search time: %f seconds.'%(toc-tic))
            
def adjust_channel_depth(grd,shpfile,lcmax=500.):
    """
    Adjusts the depths of a suntans grid object using a line shapefile.

    The shapefile must have an attribute called "contour"
    """
    from shapely import geometry, speedups
    from maptools import readShpPointLine

    if speedups.available:
        speedups.enable()

    print('Adjusting depths in channel regions with a shapefile...')
        
    # Load the shapefile
    xyline,contour = readShpPointLine(shpfile,FIELDNAME='contour')
    
    # Load all of the points into shapely type geometry
    
    # Distance method won't work with numpy array
    #P = geometry.asPoint(xy)
    
    P = [geometry.Point(grd.xv[i],grd.yv[i]) for i in range(grd.Nc)]
    
    L=[]
    for ll in xyline:
        L.append(geometry.asLineString(ll))
     
    nlines = len(L)
    weight_all = np.zeros((grd.Nc,nlines))
    for n in range(nlines):
        print('Calculating distance from line %d...'%n)
        
        dist = [L[n].distance(P[i]) for i in range(grd.Nc)]
        dist = np.array(dist)

        # Calculate the weight from the distance
        weight = -dist/lcmax+1.
        weight[dist>=lcmax]=0.
        
        weight_all[:,n] = weight

    # Now go through and re-calculate the depths
    dv =  grd.dv*(1-weight_all.sum(axis=-1))
    for n in range(nlines):
        dv += weight_all[:,n]*contour[n]
        
    grd.dv=dv   
    return grd


