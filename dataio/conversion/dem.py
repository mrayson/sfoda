# -*- coding: utf-8 -*-
"""

General class for parsing different DEMs
Created on Mon Sep 10 14:43:26 2012

@author: mrayson
"""

import numpy as np
from netCDF4 import Dataset
from scipy import spatial, interpolate
import matplotlib.pyplot as plt
import time
import shutil
import gdal
from gdalconst import * 

from scipy.ndimage import gaussian_filter, generic_filter

from soda.utils.interpXYZ import tile_vector, nn
from soda.utils.myproj import MyProj

import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
import pdb

######
def weight_inverse_dist(dist, maxdist):
    # Compute the actual weight
    w = dist/maxdist
    #w = (dist-self.maxdist)/self.maxdist # Linear
    w[dist>maxdist]=1.0
    return w

def weight_tanh_dist(dist,maxdist):
    return 0.5+0.5*np.tanh( (dist-0.5*maxdist) / (0.17*maxdist) )

def despike(values):
    centre = int(values.size / 2)
    avg = np.mean([values[:centre], values[centre+1:]])
    std = np.std([values[:centre], values[centre+1:]])
    if (avg + 2 * std < values[centre]) | (avg - 2*std > values[centre]):
        #return np.nan
        return avg
    else:
        return values[centre]



########
# Main class

class DEM(object):
    """
        General DEM class
    """    
    
    W = 1.0 # Weight
    maxdist = 250.0

    convert2utm=False
    utmzone=None
    isnorth=True
    projstr=None

    meshgrid=True
 
    def __init__(self,infile,**kwargs):
        self.infile=infile
        self.__dict__.update(kwargs)
        if self.infile[-3:]=='.nc':
            xgrd,ygrd,self.Z = self.loadnc()
        elif self.infile[-3:] in ['dem','asc']:
            xgrd,ygrd,self.Z = self.readraster()

        # Define a projection:
        if self.utmzone is not None or self.projstr is not None:
            self.P = MyProj(self.projstr, utmzone=self.utmzone, isnorth=self.isnorth)
        else:
            self.P = None

        self.update_grid(xgrd, ygrd)

        if self.convert2utm:
            self.to_xy(self.P)

    def update_grid(self, xgrd, ygrd):
        # Generate the grid
        self.x0 = xgrd.min()
        self.y0 = ygrd.min()
        self.x1 = xgrd.max()
        self.y1 = ygrd.max()
        self.dx= xgrd[1]-xgrd[0]
        self.dy= ygrd[1]-ygrd[0]
#        

        self.nx = len(xgrd)
        self.ny = len(ygrd)
        self.npts = self.nx*self.ny
#        
        if self.meshgrid:
            self.X,self.Y = np.meshgrid(xgrd,ygrd)

        self.x, self.y = xgrd, ygrd

        #if self.convert2utm:                     
        #    print 'Transforming the DEM coordinates...'
        #    # Define a projection
        #    self.X, self.Y = self.P(self.X, self.Y)
        #    self.x = self.X[0,:]
        #    self.y = self.Y[:,0]
        
    def to_xy(self, P):
        """
        Projects the coordinates using the projection object P
        """
        X, Y = P.to_xy(self.X, self.Y)
        x = X[0,:]
        y = Y[:,0]

        xgrd, ygrd = np.meshgrid(x, y)
        longrd,latgrd = P.to_ll(xgrd, ygrd)

        # Need to interpolate the depths onto this new grid
        self.Z = self.interp(longrd, latgrd, method='linear')

        #xyin = np.array([X.ravel(), Y.ravel()]).T
        #    
        #ny,nx = xgrd.shape
        #xyout = np.array([xgrd.ravel(), ygrd.ravel()]).T

        #F = nn(xyin, xyout)

        #Z = F(self.Z.ravel())
        #self.Z = Z.reshape((ny,nx))

        self.update_grid(x, y)
 
    def to_ll(self):
        """
        Projects the coordinates using the projection object P
        """
        X, Y = self.P.to_ll(self.X, self.Y)
        x = X[0,:]
        y = Y[:,0]
        self.meshgrid=False
        self.update_grid(x,y)
        self.X = X
        self.Y = Y
        
    def loadnc(self):
        """ Load the DEM data from a netcdf file"""        
        nc = Dataset(self.infile, 'r')

        xopts = ['x','X','lon','longitude']
        yopts = ['y','Y','lat','latitude']
        zopts = ['z','Z','depth','topo','elevation']

        for xx in xopts:
            if xx in list(nc.variables.keys()):
                xvar = xx

        for yy in yopts:
            if yy in list(nc.variables.keys()):
                yvar = yy

        for zz in zopts:
            if zz in list(nc.variables.keys()):
                zvar = zz

        X = nc.variables[xvar][:]
        Y = nc.variables[yvar][:]
        Z = nc.variables[zvar][:]


        #print nc.variables.keys()
        #try:
        #    X = nc.variables['X'][:]
        #    Y = nc.variables['Y'][:]
        #    Z = nc.variables['topo'][:]
        #except:
        #    try:
        #        X = nc.variables['x'][:]
        #        Y = nc.variables['x'][:]
        #        Z = nc.variables['z'][:]
        #    except:
        #        X = nc.variables['lon'][:]
        #        Y = nc.variables['lat'][:]
        #        Z = nc.variables['topo'][:]
                
                
        nc.close()
        return X,Y,Z
        
    def interp(self, x, y, method='linear'):
        """
        Interpolate DEM data onto scattered data points using
        scipy.interpolate.RectBivariateSpline
        """
        #if not self.__dict__.has_key('_Finterp'):
        if method == 'spline':
            self._Finterp = interpolate.RectBivariateSpline(self.y, self.x, self.Z,
                    bounds_error=False,fill_value=0)
        elif method == 'linear':
            self._Finterp = interpolate.RegularGridInterpolator((self.y, self.x), self.Z,
                    bounds_error=False,fill_value=0)

        if method == 'spline':
            return self._Finterp(y, x, grid=False)
        elif method == 'linear':
            return self._Finterp((y, x))

    def despike(self, size=5):
        """
        Removes spikes with an area-averaged mean
        """
        self.Z = generic_filter(self.Z, despike, size=size)

    def regrid(self, x, y, meshgrid=False):
        """
        Regrid the data onto a different constant grid

        Uses interp1d (linear)
        """

        ## Interp1d - only works on 1d grids
        if meshgrid:
            kind='linear'
            Fy = interpolate.interp1d(self.y, self.Z, axis=0, \
                    bounds_error=False, fill_value=np.nan, kind=kind)
            Zytmp = Fy(y)

            Fx = interpolate.interp1d(self.x, Zytmp, axis=1,\
                    bounds_error=False, fill_value=np.nan, kind=kind)
            self.Z =  Fx(x)

            ## Inter2d - Usually crashes due to memory
            #F = interpolate.interp2d(self.X, self.Y, self.Z)
            #self.Z = F(x, y)

            self.update_grid(x, y)
            self.X,self.Y = x,y
        else: # 2D structured grid

            ## Use nearest neighbour
            self.X,self.Y = np.meshgrid(x,y)
            xyin = np.array([self.X.ravel(), self.Y.ravel()]).T
            
            ny,nx = x.shape
            xyout = np.array([x.ravel(), y.ravel()]).T

            F = nn(xyin, xyout)

            Z = F(self.Z.ravel())
            self.Z = Z.reshape((ny,nx))

            ##
            self.update_grid(x[0,:], y[:,0])
            self.X,self.Y = x,y

    def clip(self, x0, x1, y0 , y1):
        """
        Clips the domain
        """
        #
        xidx = np.argwhere( (self.x > x0) & (self.x < x1))
        i1 = xidx[0][0]
        i2 = xidx[-1][0]

        yidx = np.argwhere( (self.y > y0) & (self.y < y1))
        j1 = yidx[0][0]
        j2 = yidx[-1][0]

        self.Z = self.Z[j1:j2,i1:i2]
        self.update_grid(self.x[i1:i2], self.y[j1:j2])

    def readraster(self):
        """ Loads the data from a DEM raster file"""
        # register all of the drivers
        gdal.AllRegister()
        # open the image
        ds = gdal.Open(self.infile, GA_ReadOnly)
        
        # Read the x and y coordinates
        cols = ds.RasterXSize
        rows = ds.RasterYSize
        bands = ds.RasterCount
        
        geotransform = ds.GetGeoTransform()
        originX = geotransform[0]
        originY = geotransform[3]
        pixelWidth = geotransform[1]
        pixelHeight = geotransform[5]
        
        x = originX + np.linspace(0,cols-1,cols)*pixelWidth
        y = originY + np.linspace(0,rows-1,rows)*pixelHeight
        
        # Read the actual data
        data = ds.ReadAsArray(0,0,cols,rows)
        
        # Remove missing points
        data[data==-32767]=np.nan
        
        return x, y, data
    
        
    def ravel(self):
        """ Returns the grid coordinates as a vector"""
        return np.concatenate( (np.reshape(np.ravel(self.X),(self.npts,1)),\
            np.reshape(np.ravel(self.Y),(self.npts,1))),axis=1)
            
    def nanxy(self):
        """
            Returns the x,y locations of the nan points
        """
        ind = np.isnan(self.Z)
        nc = np.sum(ind)
        xy = np.zeros((nc,2)) 
        xy[:,0] = self.X[ind]
        xy[:,1] = self.Y[ind]
        #n = -1
        #for jj in range(0,self.ny):  
        #    for ii in range(0,self.nx):  
        #        if ind[jj,ii]:
        #            n+=1
        #            xy[n,0]=self.X[jj,ii]
        #            xy[n,1]=self.Y[jj,ii]
        
        return xy
        
#        ind = np.isnan(np.ravel(self.Z))
#        nc = np.sum(ind)
#        
#        x=np.ravel(self.X)
#        y=np.ravel(self.Y)
#        
#        return np.concatenate((np.reshape(x[ind],(nc,1)),np.reshape(y[ind],(nc,1))),axis=1)

        
    def nonnanxy(self):
        """
            Returns the x,y locations of the non-nan points
        """
        ind = np.isnan(self.Z)
        ind = ind==False
        nc = np.sum(ind)
        xy = np.zeros((nc,2)) 
        xy[:,0] = self.X[ind]
        xy[:,1] = self.Y[ind]
        
        #n = -1
        #for jj in range(0,self.ny):  
        #    for ii in range(0,self.nx):  
        #        if ind[jj,ii]:
        #            n+=1
        #            xy[n,0]=self.X[jj,ii]
        #            xy[n,1]=self.Y[jj,ii]
        
        return xy
        
#        ind = np.isnan(np.ravel(self.Z))
#        ind = ind==False
#        nc = np.sum(ind)
#        print nc
#        x=np.ravel(self.X)
#        y=np.ravel(self.Y)
#        
#        return np.concatenate((np.reshape(x[ind],(nc,1)),np.reshape(y[ind],(nc,1))),axis=1)

    def returnij(self,x,y):
        """def contourf(self,Z,vv=range(-10,0),**kwargs):
        fig= plt.figure(figsize=(9,8))
        plt.contourf(self.X,self.Y,Z,vv,**kwargs)
        plt.colorbar()
        plt.axis('equal')
        return fig
        Returns the grid cell indices that points x,y reside inside of.
        
        """
        I = np.ceil( (x-self.x0)/self.dx)
        J =np.ceil( (y-self.y0)/self.dy)
        
        J = np.array(J,dtype=int)
        I = np.array(I,dtype=int)
        
        # blank out bad cells
        J[J<0]=-1
        J[J>self.ny-1]=-1
        I[I<0]=-1
        I[I>self.nx-1]=-1
        
        return J,I
        
    def calc_weight_convolve(self):
        """

        """
        weight = self.W*np.ones((self.ny,self.nx))
        weight[np.isnan(self.Z)] = 0.


        sigma = self.maxdist/self.dx
        weightf = gaussian_filter(weight, sigma)
        weightf[np.isnan(self.Z)] = 0.


        #plt.imshow(weightf[::4,::4])
        #plt.show()
        return weightf

    def calc_weight(self, weightfunc=weight_tanh_dist, xynan=None):
        
        """ Calculate the weight at each point """
        MAXPOINTS=20e6
        weight = np.zeros((self.ny,self.nx))
        
        # Calculate the distance from each point to a nan point
        xy = self.nonnanxy()

        if xynan is None:
            xynan = self.nanxy()
        
        # If there are no nan's return W
        if xynan.shape[0] == 0:
            weight[:] = self.W
            return weight

        # Compute the spatial tree
        print('Building KD-tree...')
        kd = spatial.cKDTree(xynan)
        
        nxy = len(xy)
        
        if nxy <= MAXPOINTS:
            # Perform query on all of the points in the grid
            dist,ind=kd.query(xy, distance_upper_bound=1.1*self.maxdist, n_jobs=-1)
            
            w = weightfunc(dist, self.maxdist)
            w=self.W*w
            
            # Map onto the grid
            J,I=self.returnij(xy[:,0],xy[:,1])
            weight[J,I]=w
        else:
            print('Dataset too large - calculating weights for chunks...')
            nchunks = np.ceil(len(xy)/MAXPOINTS)
            pt1,pt2=tile_vector(len(xy),int(nchunks))
            for p1,p2 in zip(pt1,pt2):
                print('Calculating points %d to %d of %d...'%(p1,p2,nxy))
                dist,ind=kd.query(xy[p1:p2,:])
                w = weightfunc(dist, self.maxdist)
                w=self.W*w
                
                # Map onto the grid
                J,I=self.returnij(xy[p1:p2,0],xy[p1:p2,1])
                weight[J,I]=w
        
        return weight   
        
    def contourf(self,vv=list(range(-10,0)),**kwargs):
        #fig= plt.figure(figsize=(9,8))
        plt.contourf(self.X,self.Y,self.Z,vv,**kwargs)
        plt.colorbar()
        plt.hold(True)
        plt.contour(self.X,self.Y,self.Z,[0.0,0.0],colors='k',linewidths=0.02)
        plt.axis('equal')
        
    def contour(self,vv=list(range(-10,0)),**kwargs):
        #fig= plt.figure(figsize=(9,8))
        C = plt.contour(self.X,self.Y,self.Z,vv,colors='k',linestyles='-')
        plt.axis('equal')
        return C
        
    def plot(self, shading=True, ve=1, cmap=plt.cm.gist_earth, vmin=-5000, vmax=0, **kwargs):

        # Illuminate the scene from the northwest
        Zplot = 1*self.Z

        if self.y[0] < self.y[-1]:
            Zplot = Zplot[::-1,:]

        if shading:
            Zplot[np.isnan(Zplot)]=0.
            ls = LightSource(azdeg=315, altdeg=25)
            rgb = ls.shade(Zplot, cmap=cmap, vert_exag=ve, blend_mode='overlay',\
                    vmin=vmin, vmax=vmax)
            #h= plt.figure(figsize=(9,8))
            #h.imshow(np.flipud(self.Z),extent=[bbox[0],bbox[1],bbox[3],bbox[2]])
            h = plt.imshow(rgb, extent=[self.x0,self.x1,self.y0,self.y1],**kwargs)

        else:
            h = plt.imshow(Zplot, extent=[self.x0,self.x1,self.y0,self.y1],\
                vmin=vmin, vmax=vmax, cmap=cmap, **kwargs)
        plt.colorbar(h)
        
    def savenc(self,outfile='DEM.nc'):
        """ Saves the DEM to a netcdf file"""
        
        # Create the global attributes
        
        globalatts = {'title':'DEM model',\
        'history':'Created on '+time.ctime(),\
        'Input dataset':self.infile}
        
        
        nc = Dataset(outfile, 'w', format='NETCDF4')
        # Write the global attributes
        for gg in list(globalatts.keys()):
            nc.setncattr(gg,globalatts[gg])
            
        # Create the dimensions
        dimnamex = 'nx'
        dimlength = self.nx
        nc.createDimension(dimnamex,dimlength)
        dimnamey = 'ny'
        dimlength = self.ny
        nc.createDimension(dimnamey,dimlength)
        
        # Create the lat lon variables
        tmpvarx=nc.createVariable('X','f8',(dimnamex,))
        tmpvary=nc.createVariable('Y','f8',(dimnamey,))
        #tmpvarx[:] = self.X[0,:]
        #tmpvary[:] = self.Y[:,0]
        tmpvarx[:] = self.x
        tmpvary[:] = self.y
        # Create the attributes
        tmpvarx.setncattr('long_name','Easting')
        tmpvarx.setncattr('units','metres')
        tmpvary.setncattr('long_name','Northing')
        tmpvary.setncattr('units','metres')
        
        # Write the topo data
        tmpvarz=nc.createVariable('topo','f8',(dimnamey,dimnamex),zlib=True,least_significant_digit=1)
        tmpvarz[:] = self.Z
        tmpvarz.setncattr('long_name','Topographic elevation')
        tmpvarz.setncattr('units','metres')
        tmpvarz.setncattr('coordinates','X, Y')
        tmpvarz.setncattr('positive','up')
        
        nc.close()
        
        print('DEM save to %s.'%outfile)

    
def blendDEMs(ncfile,outfile,W,maxdist):
    ### Combine multiple files ###   
    
    #Calculate the weights for each file
    nfiles = len(ncfile)
    ii=-1
    for nc in ncfile:
        ii+=1
        d = DEM(infile=nc,W=W[ii],maxdist=maxdist[ii])

        print('Calculating weights for %s...'%nc)
        print('Weight = %6.3f, maxdist = %f'%(W[ii],maxdist[ii]))
        w=d.calc_weight()

        ny = d.ny
        nx = d.nx
        
#        if ii == 1:
#            f=d.contourf(w,vv=np.linspace(0,W[ii],10))
#            f.savefig('%s_Weights.pdf'%outfile[:-2])
#            del f
        del d
        if ii == 0:
            Wall = np.zeros((ny,nx,nfiles))
        Wall[:,:,ii]=w
        del w
        
    # Normalise the weights
    print('Normalising the weights...')
    Wsum = np.sum(Wall,axis=2)
    for ii in range(0,nfiles):
        Wall[:,:,ii] = np.squeeze(Wall[:,:,ii]) / Wsum
        
    # Re-load in the depths from each file and sum
    print('Writing to an output file...')
    Zout = np.zeros((ny,nx))
    filestr = ''
    ii=-1
    for infile in ncfile:
        ii+=1
        nc = Dataset(infile, 'r')
        Zin = nc.variables['topo'][:]
        nc.close()
        
        Zin[np.isnan(Zin)]=0.0
        Zout +=  np.squeeze(Wall[:,:,ii]) * Zin 
        filestr +='%s, '%infile
    
    # Check for NaNs
    idx = np.isnan(Zout)
    if np.any(idx):
        print('Warning NaNs found. Zeroing...')
        Zout[idx] = 0.

    # Copy the data to a new netcdf file
    shutil.copyfile(ncfile[-1],outfile)
    nc = Dataset(outfile, 'r+')
    nc.variables['topo'][:]=Zout
    
    globalatts = {'title':'DEM model',\
        'history':'Created on '+time.ctime(),\
        'Input datasets':filestr}
    # Write the global attributes
    for gg in list(globalatts.keys()):
        nc.setncattr(gg,globalatts[gg])
        
    nc.close()
    
    print('Completed write to %s.'%outfile)


#ncfile = [\
#'C:/Projects/GOMGalveston/DATA/Bathymetry/DEMs/USACELIDAR_dx25_blockavg.nc',\
#'C:/Projects/GOMGalveston/DATA/Bathymetry/DEMs/NOAADEM_dx25_IDW_dist100_NNEar3.nc',\
#'C:/Projects/GOMGalveston/DATA/Bathymetry/DEMs/NOAASoundingsDEM_dx25_KRIG_dist500_Nnear3_range200.nc',\
#'C:/Projects/GOMGalveston/DATA/Bathymetry/DEMs/TNRIS_dx25_GridData.nc'\
#]
#W = [50.0,1.0,10.0,0.1]
#maxdist = [100.0,100.,1500.,1000.0]
#outfile = 'C:/Projects/GOMGalveston/DATA/Bathymetry/DEMs/Blended/NOAA_Blended_All.nc'
#
#blendDEMs(ncfile,outfile,W,maxdist)

#d = DEM(infile=outfile)
#f=d.contourf(d.Z,vv=range(-15,3),vmin=-15,vmax=4,cmap=plt.cm.gist_earth)
#f.savefig(outfile[:-2]+'pdf')
#plt.show()

## Load in other formats
#infile = 'C:/Projects/GOMGalveston/DATA/Bathymetry/NGDCCoastalRelief/galveston_tx.asc'
#print 'Loading %s...'%infile
#d = DEM(infile=infile)
#
#print 'Saving to an image...'
##f=d.contourf(d.Z,vv=range(-15,3),vmin=-15,vmax=4,cmap=plt.cm.gist_earth)
#f=d.plot(d.Z,vmin=-20,vmax=5,cmap=plt.cm.gist_earth)
#f.savefig(infile[:-3]+'png',dpi=1200)
#d.savenc(outfile='C:/Projects/GOMGalveston/DATA/Bathymetry/DEMs/NOAA_10m_DEM.nc')

#infile='C:/Projects/GOMGalveston/DATA/Bathymetry/DEMs/NOAA_10m_DEM.nc'
#print 'Loading %s...'%infile
#d = DEM(infile=infile)
#
#print 'Saving to an image...'
#f=d.contourf(d.Z+0.14,vv=range(-20,3),vmin=-20,vmax=4,cmap=plt.cm.gist_earth)
#f.savefig(infile[:-3]+'_MSL.'+'pdf')
