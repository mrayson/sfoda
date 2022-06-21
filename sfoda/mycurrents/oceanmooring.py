"""
Module for storing and manipulating ocean mooring data: z, t

Wrapper around soda->timeseries object
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator, interp1d
import xarray as xray
import pandas as pd

from sfoda.utils.timeseries import timeseries, othertime
from sfoda.utils.isoslice import isoslice
import sfoda.utils.mynumpy as mynp

import pdb

class OceanMooring(timeseries):
    """
    Ocean mooring data class
    """

    verbose = False
    positive='up'

    def __init__(self, t, y, z, zvar=None, **kwargs):
        """
        Initialize the object using a time series

        Main difference to the time series object is that
        the z-coordinate must be specified
        """
        timeseries.__init__(self, t, y, **kwargs)

        # Use datetime64
        self.t = othertime.datetimetodatetime64(self.t)

        # Input the z-coordinate
        self.is2Dz, self.Z, self.zvar, self.Nz = self._check_z(z, zvar)
        

    #############
    # Construction methods
    #############
    def vstack(self, other, clip=True):
        """
        Stack the vertical coordinate of two OceanMooring objects

        Inputs:
                other - other OceanMooring object
                clip - True, clip to overlapping time region only
                       False, nearest neighbour interpolation onto time in self
        
        """
        # Only do this if is2Dz is False
        assert not self.is2Dz

        if not clip:
            use_self = True
        #    raise NotImplementedError, 'clip must be true (for now)'

        # Check the time dimensions
        #if clip:
        #    # time step must match if clipping
        #    assert self.dt == other.dt

        if clip:
            tstart = max(self.t[0], other.t[0])
            tend = min(self.t[-1], other.t[-1])

            # Check who has the greater time step
            # (always downsample)
            if other.dt > self.dt:
                print('Warning: time steps do not match\
                (dt1 = %3.1f, dt2 = %3.1f). Downsampling...'%\
                    (self.dt, other.dt))
                use_self = False
            else:
                use_self = True
        else:
            tstart = self.t[0]
            tend = self.t[-1]

        if self.verbose:
            print(tstart, tend)

        # Clip (or interpolate) as necessary
        if clip:
            if use_self:
                tself, yself = self.subset(tstart, tend)
            else:
                tother, yother = other.subset(tstart, tend)

        else:
            tself, yself = self.t, self.y

        # This is safer than simply subsetting
        if use_self:
                tother, yother = other.interp(tself, method='nearest')
        else:
            tself, yself = self.interp(tother, method='nearest')

        # Check that the two arrays match
        # (Most times array sizes are out by one)
        if tother.shape > tself.shape:
            tother = tother[1:]
            yother = yother[1:]
        elif tother.shape < tself.shape:
            tself = tself[1:]
            yself = yself[...,1:]
        
        assert tself.shape == tother.shape,\
	    'Failed: Vectors could not be matched.'
        
        # Stack the arrays y and vector z
        if self.ndim==1 and other.ndim ==1:
            self.y = np.ma.column_stack([\
                    yself,\
                    yother,\
            ]).T
        elif self.ndim==2 and other.ndim==1:
            #print yself.shape, yother.shape
            #print yself.mask, yother.mask
            self.y = np.ma.vstack([\
                    yself,\
                    yother[np.newaxis,...],\
            ])
        else:
            self.y = np.ma.vstack([\
                    yself,\
                    yother,\
            ])


        # Sort by z axis
        Zout = np.hstack([self.Z, other.Z])
        idx = np.argsort(Zout)

        self.Z = Zout[idx]
        self.y = self.y[idx, :]

        # Update the mask (check for nans)
        nanidx = np.isnan(self.y)
        self.y.mask[nanidx] = True

        # Update the time
        self.t = tself
        self.tsec = self._get_tsec(tself)

        self.ndim = 2

        self.Nt = self.t.shape[0]
        self.Nz = self.Z.shape[0]

        return self

    def concat(self, other):
        """
        Concatenate two OceanMooring objects along the time dimesion
        """
        # check that depths are the same
        if self.is2Dz:
            zvar = self.get_2d_z()
            assert self.get_2d_z() == other.get_2d_z()
        else:
            zvar = self.Z

        t = np.hstack([self.t, other.t])
        y = np.hstack([self.y, other.y])


        return self.copy_like(t, y, zvar=zvar)


    def clip(self,time1,time2):
        """
        Returns a subset of the array between time1 and time2

        Returns an ocean mooring object
        """
        t0, t1 = self.get_tslice(time1, time2)

        if self.is2Dz:
            zvar = self.get_2d_z()[...,t0:t1]
        else:
            zvar = None

        return self.copy_like(self.t[t0:t1], self.y[...,t0:t1], zvar=zvar)
        #zvar=None
        #if self.is2Dz:
        #    zvar = self.get_2d_z()[...,t0:t1]
        #   
        #return OceanMooring(self.t[t0:t1],\
        #    self.y[...,t0:t1],\
        #    self.Z,\
        #    zvar,\
        #    units=self.units,\
        #    long_name=self.long_name,\
        #    StationID = self.StationID,\
        #    StationName = self.StationName,\
        #    varname = self.varname,\
        #    X= self.X,\
        #    Y= self.Y,\
        #    Z= self.Z,\
        #)

    def copy_like(self, t, y, zvar=None):
        """
        Create a copy of the object
        """
        if self.is2Dz and zvar is None:
            zvar = self.get_2d_z()
 
        return OceanMooring(t,\
		    y,\
		    self.Z,\
		    zvar,\
		    units=self.units,\
		    long_name=self.long_name,\
		    StationID = self.StationID,\
		    StationName = self.StationName,\
		    varname = self.varname,\
		    X= self.X,\
		    Y= self.Y,\
		    Z= self.Z,\
                    positive=self.positive,\
        )

    ############
    # Calculation methods
    ############
    def depthint(self, ztop=None, zbed=None, cumulative=False):
        """
        Depth integrate the variable

        Inputs (optional):
            ztop : (scalar) set to top of water column if different from z[-1]
            zbed : (scalar) set to bottom of water column if different from z[0]
        """
        if self.is2Dz:
            #raise NotImplementedError, 'depth integration only works for 1D depth array'
            # Calculate the vertical grid spacing
            z = self.get_2d_z()

            zmid = np.zeros((self.Nz+1, self.Nt))
            zmid[1:-1,...] = 0.5*(z[1:,...]+z[0:-1,...])
            if ztop is None:
                zmid[-1,...] = z[-1,...]
            else:
                zmid[-1,...] = ztop

            if zbed is None:
                zmid[0,...] = z[0,...]
            else:
                zmid[0,...] = zbed

            dz = np.abs(zmid[1:,...] - zmid[0:-1,...])

            # Zero out masked cells
            mask = self.y.mask
            dz[mask] = 0.

            # Perform the integration
            if cumulative:
                return np.cumsum(self.y*dz, axis=0), dz
            else:
                return np.sum(self.y*dz, axis=0), dz

        else:
            # Calculate the vertical grid spacing
            z = self.Z.copy()
            zmid = np.zeros((self.Nz+1,))
            zmid[1:-1] = 0.5*(z[1:]+z[0:-1])
            if ztop is None:
                zmid[-1] = z[-1]
            else:
                zmid[-1] = ztop

            if zbed is None:
                zmid[0] = z[0]
            else:
                zmid[0] = zbed

            dz = np.abs(zmid[1:] - zmid[0:-1])

            # Perform the integration
            if cumulative:
                return np.cumsum(self.y*dz[:,np.newaxis], axis=0), dz
            else:
                return np.sum(self.y*dz[:,np.newaxis], axis=0), dz

    def depthavg(self, ztop=None, zbed=None):
        """
        Depth average the variable

        Inputs (optional):
            ztop : (scalar) set to top of water column if different from z[-1]
            zbed : (scalar) set to bottom of water column if different from z[0]
        """
        #if self.is2Dz:
        #    #TODO
        #    raise NotImplementedError, 'depth integration only works for 1D depth array'

        y_dz, dz = self.depthint(ztop=ztop, zbed=zbed)
        H = dz.sum(axis=0)

        return y_dz/H

    def grad_z(self):
        """
        Compute the vertical gradient
        """

        zhat = self.get_2d_z()

        return mynp.grad_z(self.y, zhat)

        #if self.is2Dz:
        #    #TODO
        #    raise NotImplementedError, 'depth integration only works for 1D depth array'

        #dy_dz = np.zeros_like(self.y)
        #
        ## Second-order accurate for mid-points
        #ymid = 0.5*(self.y[1:,...]+self.y[0:-1,...])

        #zmid = 0.5*(self.Z[1:]+self.Z[0:-1])

        #dzmid  = zmid[1:] - zmid[0:-1] 

        #dy_dz[1:-1, ...] = (ymid[1:,...] - ymid[0:-1,...])/\
        #        dzmid[:,np.newaxis]

        ## First-order accurate for top and bottom cells
        #dy_dz[0,...] = (self.y[0,...] - self.y[1,...])/dzmid[0]
        #dy_dz[-1,...] = (self.y[-1,...] - self.y[-2,...])/dzmid[-1]

        #return dy_dz

    def interp_z(self, zout, method='pchip', fill_value='extrapolate'):
        """
        Interpolate in the vertical direction

        If zvar is 2D will remove NaNs before interpolation (fills gaps)
        """
        if self.is2Dz:
            zhat = self.get_2d_z()
            # Isolslice only works on scalar zout (not vector)
            #raise NotImplementedError, 'depth interpolation only works for 1D depth array'
            nz = np.size(zout)
            if nz == 1:
                return isoslice(self.y, zhat, zout, masking=True)
            else:
                y = self.y
                z = zhat
                yout = np.zeros((nz, y.shape[1]))

                for ii in range(y.shape[1]):
                    idx = y.mask[:,ii]
                    if (~idx).sum()>2:
                        if method == 'pchip':
                            F = PchipInterpolator(z[~idx,ii], y.data[~idx,ii],\
                                extrapolate=True)
                        else:
                            F = interp1d(z[~idx,ii], y[~idx,ii], kind=method,\
                                bounds_error=False, fill_value=fill_value)

                        yout[:,ii] = F(zout)

                mask = np.isnan(yout)
                return np.ma.MaskedArray(yout, mask=mask)
                   
        if method == 'pchip':
            Fi = PchipInterpolator(self.Z, self.y.data, axis=0, extrapolate=True)
        else:
            Fi = interp1d(self.Z, self.y, axis=0, kind=method,\
                bounds_error=False, fill_value=fill_value) #self.y[-1,:])

        return Fi(zout)

    def isoslice(self, isoval):
        """
        Find the vertical height of an iso-value
        """
        zhat = self.get_2d_z()
        

        return isoslice(zhat, self.y, isoval, masking=False)
    
    def get_2d_z(self):
        """
        Return a 2D z-array

        Needed for plotting, isoslice, ...
        """
        # Time and depth need to be 2D arrays
        if self.is2Dz:
            zhat = self.zvar
        else:
            zhat = self.Z[:, np.newaxis].repeat(self.Nt, axis=1)

        return zhat

    def fill_gaps_z(self, kind='nearest'):
        """     
        Fills missing values in the vertical by interpolation
        Nearest neighbour by default

        """
        y = self.y
        z = self.Z
        yout = np.zeros(y.shape)

        for ii in range(y.shape[1]):
            idx = y.mask[:,ii]
            if (~idx).sum() > 2:
                F = interp1d(z[~idx], y[~idx,ii], kind=kind,\
                    bounds_error=False, fill_value='extrapolate')
                yout[idx,ii] = F(z[idx])
                yout[~idx,ii] = y[~idx,ii]

        return self.copy_like(self.t, yout)

    def fill_gaps_z_cheby(self, order=3):
        """     
        Fills missing values in the vertical using chebyshev
        polynomials

        """
        y = self.y
        z = self.Z
        yout = np.zeros(y.shape)

        for ii in range(y.shape[1]):
            idx = y.mask[:,ii]
            if idx.sum() > 2:
                coefs = np.polynomial.chebyshev.chebfit(z[~idx], y[~idx,ii], order)
                yout[idx,ii] = np.polynomial.chebyshev.chebval(z[idx],coefs)
                yout[~idx,ii] = y[~idx,ii]

        return self.copy_like(self.t, yout)

    ###########
    # Time series overloaded functions
    ###########

    ############
    # Plotting methods
    ############
    def contourf(self, clevs, cbar=True, filled=True,
    		t0=0, t1=-1, ax=None, fig=None,
                xangle=17., cmap='Spectral_r', 
                colors='k', linewidths=0.25,
                **kwargs):
        """
        Contourf plot of the data
        """
        zhat = self.get_2d_z()

        time = self.t[np.newaxis,t0:t1].repeat(self.Nz, axis=0)

        if ax is None:
            ax=plt.gca()

        if fig is None:
            fig=plt.gcf()

        if filled:
            C = ax.contourf(time, zhat[:,t0:t1], self.y[:,t0:t1],\
	    	clevs, cmap=cmap, **kwargs)

        else:
            C = ax.contour(time, zhat[:,t0:t1], self.y[:,t0:t1],\
	    	clevs, colors=colors, linewidths=linewidths)

        # Adds some label
        if cbar:
            cb = fig.colorbar(C, ax=ax)
        else:
            cb = None
        
        plt.xticks(rotation=xangle)

        ax.set_xlim(self.t[t0], self.t[t1])

        zlim, zunits = self._check_zunits(positive=self.positive)
        ax.set_ylim(zlim)
        ax.set_ylabel(zunits)

        return C, cb


    ############
    # Conversion methods
    ############
    def to_xray(self):
        """
        Convert to an xray object
        """
        if self.is2Dz:
            da = self.to_dataarray()
            daz = self.to_dataarray(self.zvar)
            return xray.Dataset({self.varname:da,
                'zhat':daz})
        else:
            return self.to_dataarray()
    
    def to_dataarray(self, ydata=None):
        """
        Convert to an xray dataArray object
        """
        #if self.is2Dz:
        #    print 'Warning: 2D vertical variable will be lost in conversion to xray'
        if ydata is None:
            ydata = self.y

        coords = {'time':self.t}
        dims = ('time',)
        coordstr = 'time'

        if self.ndim == 2:
            coords = {'depth':self.Z, 'time':self.t}
            dims = ('depth', 'time')
        
        # Coordinate attribute for cf-compliance
        if self.is2Dz and self.ndim == 2:
            coordstr = 'zhat, time'
        if not self.is2Dz and self.ndim == 2:
            coordstr = 'depth, time'

        attrs=self.get_meta()
        attrs.update({'coordinates':coordstr})
        
        return xray.DataArray(ydata,
            dims = dims,\
            coords = coords,\
            attrs = attrs,\
        )
    
    def to_Dataset(self, ds=None):
        """
        Convert to an xray dataset
        """
        dsdict = {self.varname:self.to_xray()}
        if ds is None:
            ds = xray.Dataset(dsdict)
        else:
            ds.update(dsdict)
        
        return ds
     
    def to_netcdf(self, outfile, ds=None, group=None,\
    		encoding=None, mode='w'):
        """
        Write to a netcdf file
        """
        if self.is2Dz:
           ds = self.to_xray() # Returns a Dataset object
        else:
            ds = self.to_Dataset(ds=ds)

        ds.to_netcdf(outfile, group=group, mode=mode, encoding=encoding)

    def to_pandas(self):
        """
        Converts to a pandas DataFrame
        """
        return pd.DataFrame(self.y.T, index=self.t, columns=self.Z)

    def _check_z(self, z, zvar):
        """
        Check the z-array if it matches the array size of data in "y"

        zvar is specified if one wants a 2D z-variable
            attribute:
                Z - z coordinate variable 1D
                Zvar - either same as Z or 2D array
        """
        is2Dz = False
        if isinstance(z, float):
            return is2Dz, z, z, 1

        Nz = z.shape[0]
        if zvar is None:
            zvar = z
            is2Dz == False
            return is2Dz, z, zvar, Nz

        # Make sure the 2D zvar dimensions are in the correct order
        if zvar.ndim == 2 and zvar.shape[0] == self.y.shape[0]: 
            is2Dz = True
            #zvar = zvar.T


        #if self.Z.ndim==1 and self.Z.shape[0] == self.y.shape[0]:
        #   is2Dz = False
        #elif self.Z.shape == self.y.shape:
        #   is2Dz = True
        #elif self.Z.shape[0] == self.y.shape[0]: 
        #   is2Dz = True
        #   self.Z = self.Z.T
        #else:
        #    raise Exception, "Z array does not match shape of y"

        return is2Dz, z, zvar, Nz

    def _check_zunits(self, positive='up'):

        zhat = self.get_2d_z()

        # Work out if depths are increasing or decreasing
        minz = np.min(zhat)
        maxz = np.max(zhat)

        if positive == 'up' and minz<0.:
            zlim = [minz,0]
            zunits = 'Depth [m below free-surface]'
        elif positive == 'up' and maxz > 0:
            zlim = [0, maxz]
            zunits = 'Depth [m above seabed]'
        elif positive == 'down':
             zlim = [maxz, 0]
             zunits = 'Depth [m above seabed]'

        return zlim, zunits



    def get_meta(self):
        """
        Returns the metadata dictionary
        """
        attrvars = [
           'X',\
           'Y',\
           'Z',\
           'long_name',\
           'units',\
           'StationID',\
           'StationName',\
           'positive',\
        ]
        meta = {}
        for vv in attrvars:
            try:
                #meta.update({vv: self.__dict__[vv]})
                meta.update({vv: '%s'%getattr(self,vv)})
            except:
                if self.verbose:
                    print('Warning: Missing attribute: %s'%vv)

        return meta

    def __repr__(self):
         
        return 'OceanMooring(Nz=%d, Nt=%d)\n\tvariable: %s [%s]'\
                %(self.Nz, self.Nt, self.varname, self.units)
    
################
# Interface functions
################
def from_xray(xray_da):
    """
    Converts an xray data array to an OceanMooring object
    """
    # Load the attributes
    attrs = xray_da.attrs

    attrs.update({'varname':xray_da.name})

    ## Get the depth
    depthvars = ['DepthHeight','depth','height','Z']
    Z = 0.
    for dd in depthvars:
        if hasattr(xray_da, dd):
            try:
                Z = getattr(xray_da, dd).values
            except:
                Z = getattr(xray_da, dd)

    timevars = ['time','TIME']
    for tt in timevars:
        if hasattr(xray_da, tt):
            try:
                time = getattr(xray_da, tt).values
            except:
                time = getattr(xray_da, tt)

    # Load
    return OceanMooring(time,\
        xray_da.values, Z, **attrs)

def from_sunprofile(sunTS):
    """
    Convert a suntans profile object (xray-like object)
    """
    return OceanMooring(sunTS.time.values,  \
        sunTS.values.squeeze()[:,:],\
        -sunTS.z_r.values[:])



def from_netcdf(ncfile, varname, group=None):
    """
    Loads an OceanMooring data object directly from a netcdf file

    Uses xray
    """
    da = xray.open_dataset(ncfile, group=group, decode_coords=False)
    data = da[varname]

    # Include all of the attributes
    attrs = da.attrs
    attrs.update(data.attrs)
    attrs.update({'varname':varname})

    # possible names for the vertical coordinate
    zvars = [
	'depth',
	'distance',
	'height',
        #'zhat',
    	]
    # Work out the depth coordinate
    zcoordname = None
    for vv in list(data.coords.keys()):
        #print vv
        for zz in zvars: # loop through possible names
            if zz in vv.lower():
                #print 'zcoord name = ', vv
                zcoordname = vv
                z = data[zcoordname].values

    # Check to see if 'zhat' is in the coordinates
    try:
        if 'zhat' in data.coordinates:
            zvar=da['zhat'].values
            if 'positive' in da['zhat'].attrs:
                attrs.update({'positive':da['zhat'].attrs['positive']})
        else:
            zvar = None
    except:
        zvar=None

    #z=data.Z
    
    if zcoordname is None:
        z = float(data.Z) # Read from the attributes
    #
    #if zvar is None:
    #    z = data.Z # Read from the attributes
    #else:
    #    z = da[zvar].values
    #if zvar is not None:
    #    zvar = da[zvar].values

    return OceanMooring(data.time.values, data.values, z, zvar=zvar, **attrs)


##############
# Testing
##############


######
# Test 1: Load a series of individual files and stack them
######
#from soda.dataio import netcdfio
#dbfile = '/home/suntans/Share/ScottReef/DATA/FIELD/ScottReef_AllObs.db'
#tablename = 'observations'
#station = 'SCR400'
#
## Load the data
#temp = netcdfio.load_sql_ncstation(dbfile, station, "temperature",\
#        otherquery='StationID LIKE "%SBE39%"')

#
## Test the stacking function
##data1 = from_xray(temp[0])
##data2 = from_xray(temp[1])
##data1.vstack(data2)
#
#for ii, tt in enumerate(temp):
#    if ii == 0:
#        data = from_xray(tt)
#    else:
#        data.vstack( from_xray(tt), clip=True )
#
#datax = data.to_xray()


#########
# Test 2: Load a gridded file from a netcdf file
#########
#Tfile = 'ProcessedData/FK150410_Gridded_Mooring_TP.nc'
#
#stationT = 'SCR400'
##stationUV = 'SCR200_RDI150'
#
#daT = xray.open_dataset(Tfile, group=stationT)   
#
## Variable vertical coordinate
##data = from_netcdf(Tfile ,'temperature', group=stationT, zvar='zhat')
#
## Invariant vertical coordinate
#data = from_netcdf(Tfile ,'temperature', group=stationT, zvar=None)
#
## Test the vertical integration etc
##Tz, dz = data.depthint()
##Tavg = data.depthavg()
#dT_dz = data.grad_z()
