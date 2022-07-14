"""
Wrapper for the pyproj library
"""
import pyproj 
from pyproj import Proj
from packaging import version

class MyProjOld(object):
   
    def __init__(self, projstr, utmzone=51, isnorth=False, init=None):
        """
        Wrapper for Proj class
        Assists with creation of commonly used projections (UTM and Mercator)
        Provides convenient methods for converting between xy and ll
        """
        self.isnorth = isnorth

        
        # Create a UTM projection string
        if projstr is None:
            projstr = "+proj=utm +zone=%d, +ellps=WGS84 +datum=WGS84 +units=m +no_defs"%utmzone
            if not isnorth:
                projstr += ' +south'

        elif projstr.lower() == 'merc':
            # Mercator string
            projstr = '+proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +no_defs'

        # Projection object (super-classing doesn't work...)
        #self.P = Proj(projstr, init=init)
        self.P = Proj(projstr )

        
        # Create the inverse projection here
        #self.inverseProj = self.P.to_latlong()

    def __call__(self, lon, lat):
        return self.P(lon,lat)
    ###
    def to_xy(self, lon, lat):
        return self.P(lon, lat)

    def to_ll(self, x, y):
        return self.P(x, y, inverse=True)
    
class MyProjNew(object):
    
    def __init__(self, projstr, utmzone=51, isnorth=False):
        """
        Wrapper for Proj class
        Assists with creation of commonly used projections (UTM and Mercator)
        Provides convenient methods for converting between xy and ll
        """
        self.isnorth = isnorth
        
        if projstr is None:
            self.P = Proj(proj='utm', zone=utmzone, ellps='WGS84', north=False)
        else:
            raise(Exception('Not implemented'))
            
    def __call__(self, **args):
        
        return self.to_xy(**args)
    ###
    def to_xy(self, lon, lat):
        x, y = self.P(lon, lat)
        if not self.isnorth:
            y += 1e7
        return x, y 

    def to_ll(self, x, y):
        if not self.isnorth:
            y -= 1e7
        return self.P(x, y, inverse=True)
    
if version.parse(pyproj.__version__) > version.parse("3.0.1"):
    MyProj = MyProjNew
else:
    MyProj = MyProjOld