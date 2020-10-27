"""
Load an example suntans netcdf file
"""

import matplotlib.pyplot as plt
from sfoda.suntans import sunxray

# get the sfoda path
import os
import sfoda
_path = os.path.abspath(sfoda.__file__)
_dirpath = os.path.dirname(_path)


ncfile = "{}/../tmpdata/NWS_2km_GLORYS_hex_2013_2014_SSHBC_Harmonics.nc".format(_dirpath)
print('Loading suntans data from:\n\t{}'.format(ncfile))

sun = sunxray.Sunxray(ncfile)
print(sun)

# Xarray object is stored in the ._ds attribute
# To interface with the data
data = sun._ds['SSH_BC_var']
# Plot some data
print(data)
plt.figure()
sun.plotcelldata(data)
plt.show()


