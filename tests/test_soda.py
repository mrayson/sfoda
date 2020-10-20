"""
Test imports scripts
"""

# Try loading some of the troubling packages...

print("Running sfoda test imports...")
from sfoda.dataio import netcdfio
from sfoda.utils import timeseries, maptools

from sfoda.dataio.ugrid import hybridgrid
from sfoda.dataio.suntans import sunxray, sunplotpyqt
from sfoda.dataio.suntans import sunboundary, sundepths

from sfoda.dataio.conversion import demBuilder
from sfoda.dataio.conversion import readotps
from sfoda.dataio.datadownload import get_metocean_dap

print("Done")
