"""
Test the GDAL map packages on some data
"""

from sfoda.utils.maptools import readShpPoly

datafile = '/home/suntans/Projects/testdata/suntans_bc_markers_type2.shp'

print('Reading file: ',datafile)
data = readShpPoly(datafile)
print(data)
print('Done')

