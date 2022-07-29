"""
SUNTANS xarray wrapper

Example usage:
```
ds = open_parallel(ncfiles)
ds.suntans.plot_celldata(ds['eta'][0,...])
```
"""
import xarray as xr
import numpy as np
from sfoda.ugrid.uplot import Plot as UPlot

@xr.register_dataset_accessor("suntans")
class SUNTANS(UPlot):
    """SUNTANS xarray subclass"""
    def __init__(self, xarray_obj):
        self._obj = xarray_obj
        UPlot.__init__(self, self._obj['xp'].values, 
                       self._obj['yp'].values, self._obj['cells'].values, 
                       nfaces=self._obj['nfaces'].values,
            _FillValue=-999999)
        
    def testing(self):
        print('this is a test')
        
def load_sun(ds):
    """
    Preprocessing function: currently drops variables with the dimension 'Ne'
    """
    #print(ds.encoding['source'])
    return ds.drop_dims('Ne')

def find_ghost_cells(mnptr):
    """
    Finds ghost cells by identifying non-unique cells from 'mnptr'
    
    'mnptr' is the original cell ordering index before parallel partitioning
    """
    cells_unique = np.unique(mnptr)
    Nc = mnptr.shape[0]
    allghost = np.ones((Nc,),dtype=bool)
    allcells = set()

    for ii,cc in enumerate(mnptr.values):
        if cc in allcells:
            allghost[ii] = False
        else:
            allcells.add(cc)
            
    return allghost

def open_parallel(ncfiles, parallel=True):
    """
    Main function for loading parallel suntans output netcdf files
    
    Uses xarray.open_mfdataset but removes edge-based variables and
    ghost points.
    """
    ds = xr.open_mfdataset(ncfiles, preprocess=load_sun, 
                           concat_dim='Nc', 
                           combine='nested',
                           data_vars='minimal',
                          parallel=parallel)

    ghost = find_ghost_cells(ds['mnptr'])

    return ds.isel(Nc=ghost)
    