"""
Parse a Vemco temperature logger csv file to an xray.Dataset
"""

import xray
import numpy as np
from datetime import datetime


def hours2sec(x):
    return datetime.strptime(x,'%H:%M:%S') - datetime(1900,1,1)

def date2num(x):
    return datetime.strptime(x,'%Y-%m-%d')

def read_vemco(infile):
    """
    Read the csv file and return an xray.DataArray
    """
    headerlines = 8

    rec = np.genfromtxt(infile,
	    delimiter=',', skip_header=headerlines,
	    converters={1:hours2sec, 0:date2num})

    time = np.array([d+dt for d,dt,T in rec])
    temp = np.array([T for d,dt,T in rec])

    attrs = {'long_name':'Water temperature',
	    'units':'degC',
	    #'filename':infile,
	    }

    da = xray.DataArray(temp,
	    dims=('time',),
	    coords={'time':time},
	    attrs=attrs)
	
    return xray.Dataset({'temperature':da})


