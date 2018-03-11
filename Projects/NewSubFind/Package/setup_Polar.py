from distutils.core import setup, Extension
import numpy
from Cython.Distutils import build_ext

setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=[Extension("SPack_NewSubfind_Polar",
                 sources=["_Main_Polar.pyx", "Main_Polar.c"],
                 include_dirs=[numpy.get_include()])],
  )
#python setup.py build_ext --inplace
#python setup_Polar.py install --user
