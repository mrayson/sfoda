"""
Test imports scripts
"""

# Try loading some of the troubling packages...

print("Running soda test imports...")
from soda.dataio import netcdfio
from soda.utils import timeseries, maptools

from soda.dataio.ugrid import hybridgrid
from soda.dataio.suntans import sunxray, sunplotpyqt
from soda.dataio.suntans import sunboundary, sundepths

from soda.dataio.conversion import demBuilder
from soda.dataio.conversion import readotps
from soda.dataio.datadownload import get_metocean_dap

print("Done")
