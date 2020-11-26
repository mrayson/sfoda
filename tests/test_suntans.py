"""
Test the suntans loading functions

"""

from sfoda.suntans.sunpy import Spatial
from sfoda.suntans.sunslice import SliceEdge
from sfoda.suntans.sunxray import Sunxray
from netCDF4 import Dataset
import numpy as np

datafile = '/home/suntans/Projects/testdata/ScottReef3D_AVG_0000.nc'

print('Loading suntans file: ',datafile)

#with Dataset(datafile) as nc:
#    print(nc)
#    print(nc.variables['cells'][:])
#sun = Spatial(datafile, VERBOSE=True)

sunx = Sunxray(datafile)

print('Done')

xpt = np.array([364000., 388000])
ypt = np.array([8450000, 8460000.])

print('Loading slice object...')
sunslc = SliceEdge(datafile, xpt=xpt, ypt=ypt)

sunslc.tstep = [1,2,3,4,5]

temp = sunslc.loadData(variable='temp')

print(temp.shape)



