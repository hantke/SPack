""" Example of wrapping a C function that takes C double arrays as input using
    the Numpy declarations from Cython """

# cimport the Cython declarations for numpy
cimport numpy as np
# if you want to use the Numpy-C-API from Cython
# (not strictly necessary for this example, but good practice)
np.import_array()

# cdefine the signature of our c function
cdef extern from "Main_Polar.h":
	double Tdf_C(double M_c, double M_s,double epsilon,double alpha, double eta, double dt, double MAX, double dt_print, double sigma_r, double sigma_rp, double sigma_j, int verbose)
	double * jM_eV_C(double M_c, double M_s,double epsilon,double alpha, double eta, double dt, double MAX, double dt_print, double sigma_r, double sigma_rp, double sigma_j, int verbose)
# create the wrapper code, with numpy type annotations

def Tdf(
    double M_c, 
    double M_s, 
    double epsilon,
    double alpha, 
    double eta,
    double dt, 
    double MAX,
    double dt_print,
    double sigma_r, 
    double sigma_rp, 
    double sigma_j,
    int verbose):
	Val = Tdf_C( M_c,  M_s, epsilon, alpha, eta,  dt,  MAX, dt_print, sigma_r,  sigma_rp, sigma_j, verbose)
	return Val

def jM_eV(
    double M_c, 
    double M_s, 
    double epsilon,
    double alpha, 
    double eta,
    double dt, 
    double MAX,
    double dt_print,
    double sigma_r, 
    double sigma_rp, 
    double sigma_j,
    int verbose):
	Val = jM_eV_C( M_c,  M_s, epsilon, alpha, eta,  dt,  MAX, dt_print, sigma_r,  sigma_rp, sigma_j, verbose)
	A = []
	Len = int(MAX/dt_print)
	for i in range(6*Len):
		if Val[i%Len] >= 0: A.append(Val[i])
	return A
