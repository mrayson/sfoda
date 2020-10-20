"""
Build with:
     python setup.py build_ext --inplace

See this site for building on windows-64:
	https://github.com/cython/cython/wiki/CythonExtensionsOnWindows

See this example for packaging cython modules:
        https://github.com/thearn/simple-cython-example/blob/master/setup.py
"""

#from distutils.core import setup
#from distutils.extension import Extension
try:
    from setuptools import setup
    from setuptools import Extension
    from setuptools.dist import Distribution
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension

from Cython.Build import cythonize
#from Cython.Distutils import build_ext
import os
import numpy

class BinaryDistribution(Distribution):
    def is_pure(self):
        return False


os.environ["CC"]='cc'

extensions =[
    Extension("sfoda.ugrid.ugridutils",["sfoda/ugrid/ugridutils.pyx"],
        include_dirs=[numpy.get_include()],
        extra_compile_args=
            ['-shared', '-pthread', '-fPIC', '-fwrapv', '-O2', '-Wall',
            '-fno-strict-aliasing'],),
    Extension("sfoda.ugrid.searchutils",["sfoda/ugrid/searchutils.pyx"],
        include_dirs=[numpy.get_include()],
        extra_compile_args=['-O3','-ffast-math','-march=native','-fopenmp'],
        extra_link_args=['-fopenmp'],),
]

setup(
    name = "sfoda",
    packages=[
        'sfoda',
        'sfoda.utils',
        'sfoda.ugrid',
        'sfoda.suntans',
        'sfoda.roms',
        'sfoda.tides',    
        'sfoda.dbase',
        'sfoda.dataio.conversion',
        'sfoda.dataio.datadownload',
        ],
    #package_dir={'sfoda':''},
    ext_modules = cythonize(extensions, language_level=3),
    #cmdclass={"build_ext": build_ext},
    #ext_modules = extensions,
    version="0.X.X",
    description='Stuff For Ocean Data Analysis',
    author='Matt Rayson',
    author_email='matt.rayson@uwa.edu.au',
    #packages=find_packages(exclude=["*.tests", "*.tests.*", "tests.*", "tests"]),
    install_requires=[
      'numpy',
      'scipy',
      'matplotlib',
      'netcdf4',
      'xarray',
      'pyyaml',
      'numba',
      'pyproj',
      'dask[complete]',
      #'gdal',
      #'shapely',
      ],
    license='LICENSE',
    include_package_data=True,
    distclass=BinaryDistribution,

)
