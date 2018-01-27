""" Example of wrapping a C function that takes C double arrays as input using
    the Numpy declarations from Cython """

# cimport the Cython declarations for numpy
cimport numpy as np
# if you want to use the Numpy-C-API from Cython
# (not strictly necessary for this example, but good practice)
np.import_array()

# cdefine the signature of our c function
cdef extern from "Main.h":
	int Index_C(double x, double Xmin, double Xmax, int NBin)	
	
	double Line_C(double x,double X1,double X2,double Y1,double Y2)
	
	double FuncSingle_C(double x0,double * X,double * Y,int N)
	
	double * FuncArr_C(double * X0,double * X,double * Y,int N,int N_Arr)
	
	double * AcumMassFunction_C(double * M, double Mmin, double Mmax, double Volume, int NBin, int NGal)
	
	double * MassFunction_C(double * M, double Mmin, double Mmax, double Volume, int NBin, int NGal)
	void * Histo2D_C(long * Arr, double * X, double * Y,double Xa_min, double Xa_max,double Xb_min, double Xb_max, int NBin_a, int NBin_b, long N)
	
	double Dist3D_C(double X1,double Y1,double Z1,double X2,double Y2,double Z2)
	
	double LogDist3D_C(double X1,double Y1,double Z1,double X2,double Y2,double Z2)
	
	double Dist2D_C(double X1,double Y1,double X2,double Y2)

	double LogDist2D_C(double X1,double Y1,double X2,double Y2)

# create the wrapper code, with numpy type annotations

def Index(double x, double Xmin, double Xmax, int NBin):
	return Index_C(x, Xmin, Xmax, NBin)

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
	

def AcumMassFunction(
	np.ndarray[double, ndim=1, mode="c"] M not None,
	double Mmin, double Mmax, double Volume, int NBin):
	
	Val = AcumMassFunction_C(
	<double*> np.PyArray_DATA(M),
	Mmin,Mmax,Volume,NBin,len(M))
	A = []
	for i in range(NBin): A.append(Val[i])
	return A

def MassFunction(
	np.ndarray[double, ndim=1, mode="c"] M not None,
	double Mmin, double Mmax, double Volume, int NBin):
	
	Val = MassFunction_C(
	<double*> np.PyArray_DATA(M),
	Mmin,Mmax,Volume,NBin,len(M))
	A = []
	for i in range(NBin): A.append(Val[i])
	return A

def Histo2D(
	np.ndarray[long, ndim=1, mode="c"] Arr not None,
	np.ndarray[double, ndim=1, mode="c"] X not None,
	np.ndarray[double, ndim=1, mode="c"] Y not None,
	double Xa_min, double Xa_max,double Xb_min, double Xb_max, int NBin_a, int NBin_b):
	
	Histo2D_C(
	<long*> np.PyArray_DATA(Arr),
	<double*> np.PyArray_DATA(X),
	<double*> np.PyArray_DATA(Y),
	Xa_min,Xa_max,Xb_min,Xb_max, NBin_a, NBin_b, len(X))

def Dist3D(double X1,double Y1,double Z1,double X2,double Y2,double Z2):	return Dist3D_C(X1,Y1,Z1,X2,Y2,Z2)

def LogDist3D(double X1,double Y1,double Z1,double X2,double Y2,double Z2):	return LogDist3D_C(X1,Y1,Z1,X2,Y2,Z2)

def Dist2D(double X1,double Y1,double X2,double Y2):	return Dist2D_C(X1,Y1,X2,Y2)

def LogDist2D(double X1,double Y1,double X2,double Y2):	return LogDist2D_C(X1,Y1,X2,Y2)
#Histo2D_C(double * X,double Xa_min, double Xa_max,double Xb_min, double Xb_max, int NBin_a, int NBin_b, long N)
