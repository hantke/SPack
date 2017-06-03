""" Example of wrapping a C function that takes C double arrays as input using
    the Numpy declarations from Cython """

# cimport the Cython declarations for numpy
cimport numpy as np
# if you want to use the Numpy-C-API from Cython
# (not strictly necessary for this example, but good practice)
np.import_array()

# cdefine the signature of our c function
cdef extern from "Main.h":
	double Line_C(double x,double X1,double X2,double Y1,double Y2)
	
	double FuncSingle_C(double x0,double * X,double * Y,int N)
	
	double * FuncArr_C(double * X0,double * X,double * Y,int N,int N_Arr)
# create the wrapper code, with numpy type annotations

def Line(
	double x,
	double X1, 
	double X2, 
	double Y1, 
	double Y2):
	
	Val = Line_C(x,X1,X2,Y1,Y2)
	return Val

def FuncSingle(
	double x0,
	np.ndarray[double, ndim=1, mode="c"] X not None,
	np.ndarray[double, ndim=1, mode="c"] Y not None):
	
	Val = FuncSingle_C(
	x0,
	<double*> np.PyArray_DATA(X),
	<double*> np.PyArray_DATA(Y),
	len(X))
	return Val

def FuncArr(
	np.ndarray[double, ndim=1, mode="c"] X0 not None,
	np.ndarray[double, ndim=1, mode="c"] X not None,
	np.ndarray[double, ndim=1, mode="c"] Y not None):
	
	Val = FuncArr_C(
	<double*> np.PyArray_DATA(X0),
	<double*> np.PyArray_DATA(X),
	<double*> np.PyArray_DATA(Y),
	len(X),len(X0))
	A = []
	for i in range(len(X0)): A.append(Val[i])
	return A
	




