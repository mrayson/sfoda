# Blueprint document for SFODA python package

*SFODA* is Stuff For Ocean Data Analysis

*SODA* is/was Serious Ocean Data Analysis

# Questions

## Why create a completely new package from SODA?

 - SODA was written between 2012 and 2020 so many of the packages I now use (e.g. xarray, pandas, cartopy) did not exist when I started the project
 - SODA was written before I knew about many python coding conventions like naming classes
 - SODA has many dependencies (like GDAL) that are only used by certain functions but are still a pain to install on some systems. Moving these isolated functions out into their own sub-library would help ease this burden.
 - I now have a few other packages that would be good to migrate back into the SFODA framework (e.g. mycurrents, cymetis, oi)

## What will be different in SFODA?

 - pip installable
 - Libraries will be broken up into sub-packages. Each sub-package has a distinct hierarchy so there are no more weird cross-dependencies
 - Sub-packages to vary
 - Use of consistent naming conventions and class names
 - Replace cython code with numba to avoid compiling
 - Testing (wishful)
 - Appropriate versioning (again, wishful)
 - Documentation (dreaming)

# Structure

- sfoda: base class with existing `utils` folder e.g. time-series class

  - sfoda-ugrid: unstructured grid tools
    - sfoda-suntans: suntans libraries
  - sfoda-roms: ROMS ocean model libraries
  - sfoda-dataio: (largely redundant IO routines written prior to xarray)
  - sfoda-dbase: database tools

  (NEW STUFF)

  - sfoda-tides: tidal analysis tools (migrate some iwatlas pacakge into here)
  - sfoda-oi: optimal interpolation library
  - sfoda-oceanobs: migration of the `mycurrents` library from bitbucket
