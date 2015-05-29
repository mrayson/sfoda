Sunpy
=====

A suite of python tools for pre- and post-processing numerical model data. Specifically designed for use with the [SUNTANS model](https://github.com/ofringer/suntans).

The python files are organized into the following directories based on their general usage:

* **DataDownload** Scripts for downloading observations and model data from different web servers.

* **DataIO** General data input-output functions. Uses include netcdf, hdf and sql database files read/write/querying. 

* **GIS** GIS and other mapping tools. Reading and writing various GIS formats, creating digital elevation models and some plotting

* **SUNTANS** Python tools specific to the SUNTANS model. Includes classes for parsing model output and constructing model input data. 
	
* **Utils** Miscellaneous utilities for performing general tasks like signal processing, interpolation and other data manipulation processes.

*Created* Matt Rason, Stanford University 2012
