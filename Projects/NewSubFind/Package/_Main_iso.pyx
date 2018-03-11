""" Example of wrapping a C function that takes C double arrays as input using
    the Numpy declarations from Cython """

# cimport the Cython declarations for numpy
cimport numpy as np
# if you want to use the Numpy-C-API from Cython
# (not strictly necessary for this example, but good practice)
np.import_array()

# cdefine the signature of our c function
cdef extern from "Main_iso.h":
	double Tdf_ISO_C(double M_c, double M_s, double jmin,double epsilon,double alpha, double dt, double MAX, double dt_print, int verbose)
	double * jM_eV_ISO_C(double M_c, double M_s, double jmin,double epsilon,double alpha, double dt, double MAX, double dt_print, int verbose)
# create the wrapper code, with numpy type annotations

def Tdf_ISO(
    double M_c, 
    double M_s, 
    double jmin,
    double epsilon,
    double alpha, 
    double dt, 
    double MAX, 
    int verbose):
	dt_print = 0.1
	Val = Tdf_ISO_C( M_c,  M_s,  jmin, epsilon, alpha,  dt,  MAX, dt_print, verbose)
	return Val

def jM_eV_ISO(
    double M_c, 
    double M_s, 
    double jmin,
    double epsilon,
    double alpha, 
    double dt, 
    double MAX, 
    int verbose):
	dt_print = 0.1
	Val = jM_eV_ISO_C( M_c,  M_s,  jmin, epsilon, alpha,  dt,  MAX, dt_print, verbose)
	A = []
	Len = int(MAX/dt_print)
	for i in range(3*Len):
		if Val[i%Len] >= 0: A.append(Val[i])
	return A
