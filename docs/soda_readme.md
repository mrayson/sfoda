SODA - **S**erious **O**cean **D**ata **A**nalysis
=====

# Overview

A suite of python tools for pre- and post-processing ocean model data. Originally designed for use with the [SUNTANS model](https://github.com/ofringer/suntans).

Not everything in this toolbox is specific to the SUNTANS model though. Other uses include:
     
    - Processing [ROMS](www.myroms.org) model data
    - Download ocean and atmosphere model data via opendap e.g., CFSR, HYCOM.
    - Downloading observations from e.g., NOAA PORTS, NWS, USGS
    - Time-series analysis and signal processing (tidal harmonic fitting, power spectra, ...)
    - Interfacing with GIS data
    - Interfacing with SQL data
    - ...

*Created*:
   
   - Matt Rayson
   - Stanford University
   - 2012

# Installation

Set the *PYTHONPATH* environment variable to point to the path of this package.

# Usage

There is some incomplete documentation [here](http://suntanspy.readthedocs.org/en/latest/).

The python files are(were) organized into the following directories based on their general usage:

* **DataDownload** Scripts for downloading observations and model data from different web servers.

* **DataIO** General data input-output functions. Uses include netcdf, hdf and sql database files read/write/querying. 

* **GIS** GIS and other mapping tools. Reading and writing various GIS formats, creating digital elevation models and some plotting

* **SUNTANS** Python tools specific to the SUNTANS model. Includes classes for parsing model output and constructing model input data. 
	
* **Utils** Miscellaneous utilities for performing general tasks like signal processing, interpolation and other data manipulation processes.


