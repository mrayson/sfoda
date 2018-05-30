"""
Tools for handling SUNTANS profile bindary data

Usage:
    # Save data to netcdf
    save_profile_nc(basedir, outfile, basetime, varnames=None, **kwargs):

    # Load object
    p = Profile(ncfile)
    eta = p(x, y, z, 'eta') # load the eta variable

"""

from soda.dataio.netcdfio import dict_toxray
from soda.dataio.suntans.sunpy import calc_z
from soda.utils.othertime import datetimetodatetime64

import xray
from netCDF4 import num2date
from datetime import datetime
import numpy as np
from scipy.interpolate import PchipInterpolator
import os

import pdb

# Defaults
vardict = {
    'salt':{
        'file':'s.dat.prof',
        'is3D':True,
        'nvars':1
        },
    'temp':{
        'file':'T.dat.prof',
        'is3D':True,
        'nvars':1
        },
    'kappa_t':{
        'file':'kappat.dat.prof',
        'is3D':True,
        'nvars':1
        },
    'uvw':{
        'file':'u.dat.prof',
        'is3D':True,
        'nvars':3
        },
    'eta':{
        'file':'fs.dat.prof',
        'is3D':False,
        'nvars':1
        },
}

FILL_VALUE = -999999.

class Profile(object):

    def __init__(self, ncfile):
        self.ds = xray.open_dataset(ncfile)

    def __call__(self, x ,y, z, varname, zinterp='nearest', nnearest=1):
        """
        Returns the profile data at the nearest point

        Inputs:
        ------
                x, y, z: coordinates to extract. Set z=None to return all layers
                varname: string of variable name
                zinterp: vertical interpolation 'nearest'(default) or 'pchip'
                nnearest: number of nearest neighbours (default=1)

        Returns:
        ------
                data - xray.DataArray of the nearest point
                if nnearest>1 returns an average of the nnearest points
        """
        # Find the nearest x index
        dist = np.sqrt( (x-self.ds.xv.values)**2 + (y-self.ds.yv.values)**2)
        jj = np.argsort(dist)[0:nnearest]


        # Return the data with x,y,z coordinates as attributes
        ndim = self.ds[varname].ndim
        if ndim == 2:
            data = self.ds[varname][:,jj]
        elif ndim == 3:
            #data = self.ds[varname][:,kk,jj]
            if not z is None:
                data = self._vertical_interpolation(self.ds[varname][:,:,jj], z, zinterp)
            else:
                data = self.ds[varname][:,:,jj]

        #data.attrs.update({'X':self.ds.xv[jj],\
        #                   'Y':self.ds.yv[jj],\
        #                   'Z':self.ds.z_r[kk]})
        if nnearest > 1:
            # Zero nans
            #data.values[np.isnan(data.values)] = 0.

            # Take the mean of the three points
            dataout = data[...,0]
            #dataout.values = data.values.mean(axis=-1)
            dataout.values = np.nanmean(data.values, axis=-1)
            dataout['xv'] = data.xv.mean()
            dataout['yv'] = data.yv.mean()
            #print dataout.values.max()

        else:
            dataout = data

        return dataout

    def _vertical_interpolation(self, data2d, z, zinterp):
        """

        """
        kk, znear = self.find_nearest_depth(z)

        if zinterp == 'nearest':

            return data2d[:,kk]

        elif zinterp == 'pchip':
            
            Fi = PchipInterpolator(self.ds.z_r, data2d, axis=1, extrapolate=True)

            data = data2d[:,kk]
            data.values = Fi(z)

            return data

        else:
             raise Exception('unknown value for zinterp: '%zinterp)

    def find_nearest_depth(self, z):
        """
        Find the nearest z index
        """
        nk = self.ds.z_r.shape[0]
        Zw = np.zeros((nk+1,))
        Zw[1:] = self.ds.dz.values.cumsum()
        kk = np.searchsorted(Zw, z) - 1 

        return kk, self.ds.z_r[kk]

    def close(self):
        """
        Close the opened data file
        """
        self.ds.close()

    def __del__(self):
        self.ds.close()


###
# Read the profdata
class ProfData(object):
    """
    Object for storing the profile data file header info
    """

    profdatafile = 'profdata.dat'

    fill_value = 999999.

    def __init__(self, basedir, **kwargs):
        
        self.__dict__.update(**kwargs)

        f = open('%s/%s'%(basedir, self.profdatafile))

        self.Np = np.fromfile(f, '<i4', count=1)[0] # Number of profile stations
        self.Ni = np.fromfile(f, '<i4', count=1)[0] # Number of interpolation points
        self.Nkmax = np.fromfile(f, '<i4', count=1)[0]
        self.Nt = np.fromfile(f, '<i4', count=1)[0]
        self.Ntout = np.fromfile(f, '<i4', count=1)[0]
        self.dt = np.fromfile(f, '<f8', count=1)[0]
        self.dz = np.fromfile(f, '<f8', count=self.Nkmax)

        self.indices = np.fromfile(f, '<i4', count=self.Np)
        #self.xyorig = np.fromfile(f, '<f8', count=2*self.Np*self.Ni)
        self.xyorig = np.fromfile(f, '<f8', count=2*self.Np)

        self.xv = np.fromfile(f, '<f8', count=self.Np*self.Ni)
        self.yv = np.fromfile(f, '<f8', count=self.Np*self.Ni)

        f.close()

        self.z_r = calc_z(self.dz)

class ReadProf(ProfData):
    
    def __init__(self, basedir, varname, basetime):


        ProfData.__init__(self,basedir)

        self.varname = varname

        # Get the file size to determine the number of steps
        self.datafile = '%s/%s'%(basedir, vardict[varname]['file'])

        self.nvars = vardict[varname]['nvars']

        if vardict[varname]['is3D']:
            self.nk = self.Nkmax
            self.is3D = True
        else:
            self.nk = 1 
            self.is3D = False

        f_size = os.stat(self.datafile).st_size

        # Size of variable for one time step
        self.varsize = self.Np*self.Ni*self.nvars*self.nk
        self.nsteps = f_size// (self.varsize*8)

        # Calculate the time
        self.calc_time(basetime)

    def calc_time(self, basetime):
        """
        calculate the time
        """
        dtout = self.dt*self.Ntout
        self.tsec = np.arange(0, self.nsteps*dtout, dtout)

        # Convert the basetime
        try:
            t0 = datetime.strptime(basetime, '%Y%m%d.%H%M%S')
        except:
            raise Exception('time format must be "YYYYmmdd.HHMMSS"')

        self.tunits = 'seconds since %s'%(datetime.strftime(t0, '%Y-%m-%d:%H%M%S'))

        time = num2date( self.tsec, self.tunits)
        self.time = datetimetodatetime64(time)

    def __call__(self, tstep=-1):

        # Read all of the data
        f = open(self.datafile,'r')

        data = np.fromfile(f, '<f8', count = self.nsteps*self.varsize)

        f.close()

        if self.is3D and self.nvars==1:
            #order = (self.Np*self.Ni, self.nk , self.nsteps)
            order = (self.nsteps , self.Np*self.Ni, self.nk )
            data = data.reshape(order)
            data = data.swapaxes(1,2)

            return {self.varname:np.ma.MaskedArray(data, mask = data==self.fill_value)}

        elif self.is3D and self.nvars==3: # varname=='uvw'
            data = data.reshape((self.nsteps, 3, self.Np*self.Ni, self.nk))
            data = data.swapaxes(2,3)

            u = data[:,0,...].squeeze()
            v = data[:,1,...].squeeze()
            w = data[:,2,...].squeeze()

            return {'uc':np.ma.MaskedArray(u, mask = u==self.fill_value),\
                'vc':np.ma.MaskedArray(v, mask = v==self.fill_value),\
                'w':np.ma.MaskedArray(w, mask = w==self.fill_value),\
                }

        else:
            data = data.reshape((self.nsteps,self.Np*self.Ni ))
            #data = data.swapaxes(0,1)

            return {self.varname:np.ma.MaskedArray(data, mask = data==self.fill_value)}


        #return data


def save_profile_nc(basedir, outfile, basetime, varnames=None, **kwargs):
    """
    Convert the profile binary file to netcdf
    """
    if varnames is None:
        varnames = list(vardict.keys())

    # Loop through and populate a dataset
    ds = {}
    for varname in varnames:
        print('Loading variable %s...'%varname)

        V = ReadProf(basedir, varname, basetime, **kwargs)
        data = V()

        if V.is3D:
            dims =('time','Nk','Nc') 
            coords = {'time':V.time,\
                'Nk':list(range(0,V.Nkmax)),\
                'Nc':list(range(0,V.Np*V.Ni))}
        else:
            dims =('time','Nc') 
            coords = {'time':V.time,\
                'Nc':list(range(0,V.Np*V.Ni))}

        encoding = {
                'complevel':5,
                'zlib':True,
            }

        ds = dict_toxray(data, ds=ds,\
            dims=dims,\
            coords = coords,\
            encoding = encoding)


    # Add in the coordinates last
    ds = dict_toxray({'xv':V.xv, 'yv':V.yv},ds=ds,\
        dims=('Nc'))

    ds = dict_toxray({'dz':V.dz, 'z_r':V.z_r},ds=ds,\
        dims=('Nk'))

    # Save to netcdf
    ds.to_netcdf('%s/%s'%(basedir,outfile), format='NETCDF4_CLASSIC')

    print('\tProfile data saved to %s.'%outfile)

###########
## Testing
###########

