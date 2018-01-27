from distutils.core import setup, Extension
import numpy
from Cython.Distutils import build_ext

setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=[Extension("SPack_NewSubfind",
                 sources=["_Main.pyx", "Main.c"],
                 include_dirs=[numpy.get_include()])],
  )
#python setup.py build_ext --inplace
#python setup.py install --user
