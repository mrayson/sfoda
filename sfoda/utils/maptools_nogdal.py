"""
Map tools wrapper
"""

try:
    from .maptools_gdal import *
except:
    print('Warning: cannot import GDAL libraries')

    def readShpPoly():
        raise Exception

    def readShpBathy():
        raise Exception

    def readraster():
        raise Exception

    def readDEM():
        raise Exception

    def ll2lcc():
        raise Exception

    def ll2utm():
        raise Exception
