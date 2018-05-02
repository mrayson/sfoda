"""
Lightweight SUNTANS class that takes advantage of other unstructured grid
classes and xray/dask.
"""

import xray

import matplotlib.pyplot as plt
import numpy as np

import glob
import dask.array as da
import dask

import pdb

# matplotlib
from soda.dataio.ugrid.uplot import Plot as UPlot

## VTK stuff
#from soda.dataio.ugrid.uplotvtk import PlotVTK as UPlot
#from mayavi import mlab

class Sunxray(UPlot):
    """
    SUNTANS-xray wrapper

    UPlot inherits everything from the HybridGrid class

    xray object is stored in the "_ds" attribute
    """
    def __init__(self, ncfile, lazy=False, **kwargs):

        self._ds = xray.open_dataset(ncfile, \
                mask_and_scale=True, decode_times=True)

        if not lazy:
            self._init_grid(**kwargs)
            xp = self._ds.xp.values
            yp = self._ds.yp.values
            self.xlims = [xp.min(), xp.max()]
            self.ylims = [yp.min(), yp.max()]


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
        for vv in self._ds.variables.keys():
            # "mesh" attribute is standard in the ugrid convention
            if hasattr(self._ds[vv], 'mesh'):
                vname.append(vv)
                
        return vname


    def load(self, varname, tstep, klayer):
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
        return self._ds[varname][:].values

    def __repr__(self):
        return self._ds.__repr__()


class Sundask(Sunxray):
    """
    Parallel SUNTANS reader

    Goal is to have the same functionality to the user as Sunxray
    i.e. not worry about the parallel io, dask, etc
    """

    def __init__(self, ncfiles, **kwargs):

        print ncfiles

        # Get the files from the first file
        filenames = sorted(glob.glob(ncfiles))
        
        #for ff in filenames:
        #    print ff
        # Load all of the files into list as xray objects
        self._myfiles = [Sunxray(url, lazy=True) for url in filenames]

        # Keep this for compatibility with methods from the superclass
        self._ds = self._myfiles[0]._ds

        # Load the grid variables and re-sort

        # Each files stores all node coordinates (luckily)
        xp = self._myfiles[0].loadfull('xp')
        yp = self._myfiles[0].loadfull('yp')

        # Load the actual data
        cells = self.loadfull('cells').compute()
        nfaces = self.loadfull('nfaces').compute()

        ### Grid does not need re-sorting...

        # Finish initializing the class
        UPlot.__init__(self, xp, yp, cells, nfaces=nfaces,\
            _FillValue=-999999,\
                **kwargs)

        self.xlims = [xp.min(), xp.max()]
        self.ylims = [yp.min(), yp.max()]



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
        return da.concatenate(arrays, axis=axis)

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



