from distutils.core import setup, Extension
import numpy
from Cython.Distutils import build_ext

setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=[Extension("SPack_NewSubfind_iso",
                 sources=["_Main_iso.pyx", "Main_iso.c"],
                 include_dirs=[numpy.get_include()])],
  )
#python setup.py build_ext --inplace
#python setup_iso.py install --user
