"""
Build with:
     python setup.py build_ext --inplace

See this site for building on windows-64:
	https://github.com/cython/cython/wiki/CythonExtensionsOnWindows
"""

#from distutils.core import setup
#from distutils.extension import Extension
try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension

from Cython.Build import cythonize
import os
import numpy

os.environ["CC"]='cc'

extensions =[
    Extension("ugridutils",["dataio/ugrid/ugridutils.pyx"],
        include_dirs=[numpy.get_include()],
        extra_compile_args=
            ['-shared', '-pthread', '-fPIC', '-fwrapv', '-O2', '-Wall',
            '-fno-strict-aliasing'],),
    Extension("searchutils",["dataio/ugrid/searchutils.pyx"],
        include_dirs=[numpy.get_include()],
        extra_compile_args=['-O3','-ffast-math','-march=native','-fopenmp'],
        extra_link_args=['-fopenmp'],),
]

setup(
    name = "SODA utilities",
    ext_modules = cythonize(extensions, language_level=3)
)
