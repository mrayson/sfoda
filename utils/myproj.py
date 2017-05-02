"""
Wrapper for the pyproj library
"""

from pyproj import Proj
import pdb

class MyProj(object):
   
    def __init__(self, projstr, utmzone=51, isnorth=False):
        """
        Wrapper for Proj class

        Assists with creation of commonly used projections (UTM and Mercator)

        Provides convenient methods for converting between xy and ll
        """

        # Create a UTM projection string
        if projstr is None:
            projstr = "+proj=utm +zone=%d, +ellps=WGS84 +datum=WGS84 +units=m +no_defs"%utmzone
            if not isnorth:
                projstr += ' +south'

        elif projstr.lower() == 'merc':
            # Mercator string
            projstr = '+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +no_defs'
    
        # Projection object (super-classing doesn't work...)
        self.P = Proj(projstr)
        
        # Create the inverse projection here
        #self.inverseProj = self.P.to_latlong()

    def __call__(self, lon, lat):
        return self.P(lon,lat)
    ###
    def to_xy(self, lon, lat):
        return self.P(lon, lat)

    def to_ll(self, x, y):
        return self.P(x, y, inverse=True)

    #def __new__(self, projstr, **kwargs):
    #    if projstr is None:
    #        return
    #    else:
    #        return Proj.__new__(self, projstr)

################
#### Testing ###
#utmzone = 51
#isnorth = False
#projstr = 'merc'
#P = MyProj(projstr, utmzone=utmzone, isnorth=isnorth)
#
#print P([124.,124.1],[-12.,-12.1])
#
