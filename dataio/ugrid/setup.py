"""
Build with:
     python setup.py build_ext --inplace
"""

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import os
import numpy

os.environ["CC"]='cc'

extensions =[
    Extension("ugridutils",["ugridutils.pyx"],
        include_dirs=[numpy.get_include()],
        extra_compile_args=['-O3'],),
    Extension("searchutils",["searchutils.pyx"],
        include_dirs=[numpy.get_include()],
        extra_compile_args=['-O3','-ffast-math','-march=native','-fopenmp'],
        extra_link_args=['-fopenmp'],),
]

setup(
    name = "Shallow water utilities",
    ext_modules = cythonize(extensions)
)
