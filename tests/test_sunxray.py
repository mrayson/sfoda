"""
Test a standalone suntans xarray like object
"""
from sfoda.suntans.sunxray import Sunxray
datafile = '/home/suntans/code/iwatlas/DATA/NWS_2km_GLORYS_hex_2013_2014_InternalWave_Atlas.nc'
print('Loading file', datafile)
sun=Sunxray(datafile)
print('Done')
