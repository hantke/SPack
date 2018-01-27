""" Example of wrapping a C function that takes C double arrays as input using
    the Numpy declarations from Cython """

# cimport the Cython declarations for numpy
cimport numpy as np
# if you want to use the Numpy-C-API from Cython
# (not strictly necessary for this example, but good practice)
np.import_array()

# cdefine the signature of our c function
cdef extern from "Main.h":
	double Tdf_NFW_C(double M_c, double M_s, double jmin,double epsilon,double alpha, double dt, double MAX, int verbose)
# create the wrapper code, with numpy type annotations

def Tdf_NFW(
    double M_c, 
    double M_s, 
    double jmin,
    double epsilon,
    double alpha, 
    double dt, 
    double MAX, 
    int verbose):
	Val = Tdf_NFW_C( M_c,  M_s,  jmin, epsilon, alpha,  dt,  MAX, verbose)
	return Val
