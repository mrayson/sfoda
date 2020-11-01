"""
Lightweight SUNTANS class that takes advantage of other unstructured grid
classes and xray/dask.
"""

import xarray as xr
from netCDF4 import Dataset

import numpy as np

import dask.array as da
import dask
from dask import delayed
import glob

import pdb

# matplotlib
from sfoda.ugrid.uplot import Plot as UPlot

## VTK stuff
#from sfoda.dataio.ugrid.uplotvtk import PlotVTK as UPlot
#from mayavi import mlab

# Test manually setting the scheduler
#dask.config.set(scheduler='processes')

class Sunxray(UPlot):
    """
    SUNTANS-xray wrapper

    UPlot inherits everything from the HybridGrid class

    xray object is stored in the "_ds" attribute
    """
    
    # Default chunking
    #chunks={'Nk':50,'Nc':-1}
    chunks=None

    def __init__(self, ncfile, lazy=False, **kwargs):
        self.__dict__.update(kwargs)

        try:
            self._ds = xr.open_dataset(ncfile, \
                    chunks=self.chunks,\
                    mask_and_scale=True, decode_times=True)
        except:
            self._ds = xr.open_mfdataset(ncfile, \
                    concat_dim = 'time',
                    data_vars='minimal',
                    chunks=self.chunks,\
                    mask_and_scale=True, decode_times=True)
        

        if not lazy:
            self._init_grid(**kwargs)
            xp = self._ds.xp.values
            yp = self._ds.yp.values
            self.xlims = [xp.min(), xp.max()]
            self.ylims = [yp.min(), yp.max()]

            # Calculate the voronoi
            self.calc_centroids()


    def _init_grid(self, **kwargs):
        nfaces = np.sum(self._ds.cells.values >-1,axis=1)
        UPlot.__init__(self, self._ds.xp.values, self._ds.yp.values,\
            self._ds.cells.values,\
            nfaces=nfaces,\
            _FillValue=-999999,\
                **kwargs)

    def has_dim(self,varname, dimname):
        """
        Tests if a variable contains a dimension
        """

        dimensions = self._ds[varname].dims
        
        return dimname in dimensions
        #try:
        #    self.nc.dimensions[dimname].__len__()
        #    return True
        #except:
        #    return False
 
    def list_coord_vars(self):
        """
        List all of the variables that have the 'mesh' attribute
        """
        vname=[]
        for vv in list(self._ds.variables.keys()):
            # "mesh" attribute is standard in the ugrid convention
            if hasattr(self._ds[vv], 'mesh'):
                vname.append(vv)
                
        return vname


    def load_step(self, varname, tstep=slice(None,), klayer=slice(None,)):
        """
        Load the data from the xray.Dataset object
        """
        ndim = self._ds[varname].ndim
        if ndim==2:
            return self._ds[varname][tstep,:]
        elif ndim==3:
            return self._ds[varname][tstep,klayer,:]

    def loadfull(self, varname, axis=0):
        """
        Load data from a complete/full variable
        """
        return self._ds[varname][:]

    def __repr__(self):
        return self._ds.__repr__()


class Sundask(UPlot):
    """
    Parallel SUNTANS reader
    Goal is to have the same functionality to the user as Sunxray
    i.e. not worry about the parallel io, dask, etc
    """

    timedim = 'time'
    client = None
    _FillValue = 999999
    isparralel = True
    chunks = {}

    def __init__(self, ncfiles, **kwargs):
        self.__dict__.update(kwargs)

        # Load all of the files into list as xray objects
        #self._myfiles = [Sunxray(url, lazy=True) for url in filenames]
        self.filenames = sorted(glob.glob(ncfiles))
        print('Time dimension %s'%self.timedim)

        def openfile(url):
            try:
                #print(url)
                # Chunking is necessary to use dask
                #ds = xr.open_dataset(url, chunks={'Nk':-1,'Nc':-1,'time':-1})
                ds = xr.open_dataset(url, chunks=self.chunks)
                #ds = xr.open_dataset(url, chunks={})

            except:
                print('Failed to open file %s'%url)
            
            return ds
        
        self._myfiles = [openfile(url) \
            for url in self.filenames]

        # Keep this for compatibility with methods from the superclass
        self._ds = self._myfiles[0]

        try:
            self.Nt = self._ds[self.timedim].shape[0]
        except:
            print('No time dimesion')

        # Load the grid variables 

        ### Grid does not need re-sorting... but we need ghost points
        self.ghost = self.find_ghosts()
        
        # Each files stores all node coordinates (luckily)
        xp = self._myfiles[0]['xp']
        yp = self._myfiles[0]['yp']

        # Load the actual data
        cells = self.stack_var_2d('cells', axis=0)[self.ghost,...].compute()
        cells[cells == self._FillValue] = -1

        nfaces = self.stack_var_2d('nfaces', axis=0)[self.ghost].compute()
        
        # Finish initializing the class
        UPlot.__init__(self, xp, yp, cells, nfaces=nfaces,\
            _FillValue=self._FillValue,\
                **kwargs)

        # Optional variables (not necessary but nice)
        self.xv  = self.stack_var_2d('xv', axis=0)[self.ghost].compute() 
        self.yv  = self.stack_var_2d('yv', axis=0)[self.ghost].compute() 
        self.dv  = self.stack_var_2d('dv', axis=0)[self.ghost].compute() 
        self.Nk  = self.stack_var_2d('Nk', axis=0)[self.ghost].compute() 

        self.xlims = [xp.min(), xp.max()]
        self.ylims = [yp.min(), yp.max()]

    def list_coord_vars(self):
        """
        List all of the variables that have the 'mesh' attribute
        """
        vname=[]
        for vv in self._ds.variables.keys():
            # "mesh" attribute is standard in the ugrid convention
            if hasattr(self._ds[vv], 'mesh'):
                vname.append(vv)
                
        return vname


    def has_dim(self,varname, dimname):
        """
        Tests if a variable contains a dimension
        """

        dimensions = self._ds[varname].dims
        
        return dimname in dimensions

    def find_ghosts(self):

        ## Serial file type
        #if ~hasattr(self._ds, 'mnptr'):
        #    allghost = np.ones((self._ds['xv'].shape), dtype=np.bool)
        #    return allghost

        # find ghost cells in each dataset
        mnptr = self.stack_var_2d('mnptr',axis=0).compute()
        cells_unique = np.unique(mnptr)
        #Nc = cells_unique.shape[0]
        Nc = mnptr.shape[0]
        allghost = np.ones((Nc,),dtype=np.bool)
        allcells = set()

        for ii,cc in enumerate(mnptr):
            if cc in allcells:
                allghost[ii] = False
            else:
                allcells.add(cc)
                
        return allghost
    
    def find_ghosts_old(self):
        # find ghost cells in each dataset
        allcells = set()
        allghost = np.array((0,),dtype=np.bool)

        for myfile in self._myfiles:
            #positions = [myfile['mnptr'].values for myfile in myfiles]
            myghost = np.ones((myfile.dims['Nc'],), dtype=np.bool)
            for ii, cc in enumerate(myfile['mnptr'].values):
                if cc in allcells:
                    myghost[ii] = False
                else:
                    allcells.add(cc)

            #allghost.append(myghost)
            #allghost = np.hstack([np.array(myghost),allghost])
            allghost = np.hstack([allghost,np.array(myghost)])


        return allghost


    ### Loading functions
    def load_data(self, varname, tstep=slice(None,), klayer=slice(None,)):
        """
        Load the data from the xray.Dataset object
        """
        data = self.get_data(varname, tstep=tstep, klayer=klayer)
        return data.compute()[...,self.ghost]
    
    def get_data(self, varname, tstep=slice(None,), klayer=slice(None,)):
        """
        Get the data from the xray.Dataset object
        
        This does not load the actual data into memory see: `load_data`
        """
        myvar = self._myfiles[0][varname]
        ndim = myvar.ndim
        arrays=[]
        if ndim==1:
            for ii, data in enumerate(self._myfiles):
                mydata = data[varname].data[:]
                arrays.append(mydata)
        if ndim==2:
            for ii, data in enumerate(self._myfiles):
                mydata = data[varname].data[tstep,:]
                arrays.append(mydata)
        elif ndim==3:   
            for ii, data in enumerate(self._myfiles):
                mydata = data[varname].data[tstep,klayer,:]
                arrays.append(mydata)
        
        # Stack all small Dask arrays into one
        return dask.array.concatenate(arrays, axis=-1)#[...,ghost]
    
    def stack_var_3d(self, varname, klayer, axis=-1):
        arrays=[]
        for ii, data in enumerate(self._myfiles):
            mydata = data[varname].data[:,klayer,:]
            arrays.append(mydata)

        return dask.array.concatenate(arrays, axis=axis)

    def stack_var_2d(self, varname, axis=-1):
        arrays=[]
        for ii, data in enumerate(self._myfiles):
            mydata = data[varname].data[:]
            arrays.append(mydata)

        return dask.array.concatenate(arrays, axis=axis)
    
    #########
    # Apply functions
    def apply_func(self, func, varname, *args, **kwargs):
        """
        Apply the function to each block of data (doesn't use xarray)

        See here:
                http://dask.pydata.org/en/latest/delayed-best-practices.html
        """
        @delayed
        def load_single_nc(ncfile, varname):
            with Dataset(ncfile) as nc:
                # Load the data
                X = nc.variables[varname][:]
                X[np.isnan(X)] = 0.
            return X

        @delayed
        def lazy_func(func, X, *args, **kwargs):
            return func(X, *args, **kwargs)

        def f(func, ncfiles, varname, *args, **kwargs):
            output = []
            for ncfile in ncfiles:
                X = load_single_nc(ncfile, varname)
                output.append(lazy_func(func, X, *args, **kwargs))

            return output

        stack = dask.persist(f(func, self.filenames, varname, *args, **kwargs))
        return np.concatenate([ii.compute() for ii in stack[0]], axis=-1)[...,self.ghost]

    #########
    # Old
    def loadfull(self, varname, axis=0):
        """
        Load a full array (no indexing) into a dask array
        """

        def localload(ds, varname):
            return ds[varname][:].values

        loadfull_d = dask.delayed(localload, pure=True) 
        all_files = self._myfiles

        # Lazily evaluate loadfull on each file
        lazy_data = [loadfull_d(url._ds, varname) for url in all_files]

        # Construct a Dask array
        arrays = [da.from_delayed(lazy_value,           
                      dtype=sun._ds[varname].dtype,   # for every lazy value
                      shape=sun._ds[varname].shape,
                      ) for sun, lazy_value in zip(all_files, lazy_data)]

        # Stack all small Dask arrays into one
        return dask.array.concatenate(arrays, axis=axis)

    def load(self, varname, tstep, klayer):
        """
        Load the data from the xray.Dataset object
        """
        myvar = self.loadfull(varname, axis=-1)
        ndim = myvar.ndim
        if ndim==1:
            return myvar[:].compute()
        if ndim==2:
            return myvar[tstep,:].compute().ravel()
        elif ndim==3:   
            return myvar[tstep,klayer,:].compute().ravel()
        

#class Sundask(Sunxray):
#    """
#    Parallel SUNTANS reader
#
#    Goal is to have the same functionality to the user as Sunxray
#    i.e. not worry about the parallel io, dask, etc
#    """
#
#    def __init__(self, ncfiles, loadgrid=True, **kwargs):
#
#        print ncfiles
#
#        # Get the files from the first file
#        filenames = sorted(glob.glob(ncfiles))
#        
#        #for ff in filenames:
#        #    print ff
#        # Load all of the files into list as xray objects
#        self._myfiles = [Sunxray(url, lazy=True, **kwargs) \
#                for url in filenames]
#
#        # Keep this for compatibility with methods from the superclass
#        self._ds = self._myfiles[0]._ds
#
#        # Load the grid variables and re-sort
#        if loadgrid:
#            # Each files stores all node coordinates (luckily)
#            xp = self._myfiles[0].loadfull('xp')
#            yp = self._myfiles[0].loadfull('yp')
#
#            # Load the actual data
#            cells = self.loadfull('cells').compute()
#            nfaces = self.loadfull('nfaces').compute()
#
#            ### Grid does not need re-sorting...
#
#            # Finish initializing the class
#            UPlot.__init__(self, xp, yp, cells, nfaces=nfaces,\
#                _FillValue=-999999,\
#                    **kwargs)
#
#            self.xlims = [xp.min(), xp.max()]
#            self.ylims = [yp.min(), yp.max()]
#
#
#
#    def loadfull(self, varname, axis=0):
#        """
#        Load a full array (no indexing) into a dask array
#        """
#
#        def localload(ds, varname):
#            return ds[varname][:].values
#
#        loadfull_d = dask.delayed(localload, pure=True) 
#        all_files = self._myfiles
#
#        # Lazily evaluate loadfull on each file
#        lazy_data = [loadfull_d(url._ds, varname) for url in all_files]
#
#        # Construct a Dask array
#        arrays = [da.from_delayed(lazy_value,           
#                      dtype=sun._ds[varname].dtype,   # for every lazy value
#                      shape=sun._ds[varname].shape,
#                      ) for sun, lazy_value in zip(all_files, lazy_data)]
#
#        # Stack all small Dask arrays into one
#        return da.concatenate(arrays, axis=axis)
#
#    def load(self, varname, tstep, klayer):
#        """
#        Load the data from the xray.Dataset object
#        """
#        myvar = self.loadfull(varname, axis=-1)
#        ndim = myvar.ndim
#        if ndim==1:
#            return myvar[:].compute()
#        if ndim==2:
#            return myvar[tstep,:].compute().ravel()
#        elif ndim==3:
#            return myvar[tstep,klayer,:].compute().ravel()
#
