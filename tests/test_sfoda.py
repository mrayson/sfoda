"""
Test imports scripts
"""

# Try loading some of the troubling packages...

print("Running sfoda test imports...")
from sfoda.dbase import netcdfio
from sfoda.utils import timeseries 
from sfoda.utils import maptools_nogdal

from sfoda.ugrid import hybridgrid
from sfoda.suntans import sunxray #, sunplotpyqt
from sfoda.suntans import sunboundary, sundepths
from sfoda.suntans import suntides

from sfoda.dataio.conversion import demBuilder
from sfoda.tides import readotps
from sfoda.dataio.datadownload import get_metocean_dap

print("Done")
