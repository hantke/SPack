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
	
	void Vel_1D_C(double * X, double * Y, double * Z, double * VZ, double * gv,double * gv_pc, double * gg,double * gg_pc, long * Prim_ID , long * Sec_ID, long * Main_Des_ID , int Prim_N, int Sec_N,int N, double NBin, int iNBin, double LIMIT, double Xmax, double Xmin, double Lbox, int NCPU, int CPU)
	
	void Halo_Population_C(long * Halo_NGal_All,long * Halo_NGal_Cen,long * Halo_NGal_Sat1,long * Halo_NGal_Sat2, long * Halo_ID,long * Gal_HaloID, double * Gal_Prop,long * Gal_type, double Cut, int NGal)
	
	void HOD_Group_C(double * HOD_Full, double * HOD_Cen, double * HOD_Sat1, double * HOD_Sat2,double * HOD_NHalo,long * Halo_NGal_Full,long * Halo_NGal_Cen,long * Halo_NGal_Sat1,long * Halo_NGal_Sat2, double * Halo_Mass, int NHalos, double Xmin, double Xmax, double NBin)
	
# create the wrapper code, with numpy type annotations

def Line(
	double x,
	double X1, 
	double X2, 
	double Y1, 
	double Y2):
	Val = Line_C(x,X1,X2,Y1,Y2)
	return Val

def Vel_1D(
	np.ndarray[double, ndim=1, mode="c"] X  not None,
	np.ndarray[double, ndim=1, mode="c"] Y  not None,
	np.ndarray[double, ndim=1, mode="c"] Z  not None,
	np.ndarray[double, ndim=1, mode="c"] VZ not None,
	np.ndarray[double, ndim=1, mode="c"] gv not None,
	np.ndarray[double, ndim=1, mode="c"] gv_pc not None,
	np.ndarray[double, ndim=1, mode="c"] gg  not None,
	np.ndarray[double, ndim=1, mode="c"] gg_pc  not None,
	
	np.ndarray[long, ndim=1, mode="c"] Prim_ID  not None,
	np.ndarray[long, ndim=1, mode="c"] Sec_ID  not None,
	np.ndarray[long, ndim=1, mode="c"] Main_Des_ID  not None,
	double NBin, int iNBin, double LIMIT,double Xmax, 
	double Xmin, double Lbox, int NCPU, int CPU
	):
	Vel_1D_C(
        <double*> np.PyArray_DATA(X),
        <double*> np.PyArray_DATA(Y),
        <double*> np.PyArray_DATA(Z),
        <double*> np.PyArray_DATA(VZ),
        <double*> np.PyArray_DATA(gv),
        <double*> np.PyArray_DATA(gv_pc),
        <double*> np.PyArray_DATA(gg),
        <double*> np.PyArray_DATA(gg_pc),
        <long*> np.PyArray_DATA(Prim_ID),
        <long*> np.PyArray_DATA(Sec_ID),
        <long*> np.PyArray_DATA(Main_Des_ID),
        len(Prim_ID), len(Sec_ID),len(X), NBin, iNBin, LIMIT, Xmax, Xmin, Lbox, NCPU, CPU
        )

def Halo_Population(
	np.ndarray[long, ndim=1, mode="c"] Halo_NGal_All  not None,
	np.ndarray[long, ndim=1, mode="c"] Halo_NGal_Cen  not None,
	np.ndarray[long, ndim=1, mode="c"] Halo_NGal_Sat1 not None,
	np.ndarray[long, ndim=1, mode="c"] Halo_NGal_Sat2 not None,
	np.ndarray[long, ndim=1, mode="c"] Halo_ID        not None,
	np.ndarray[long, ndim=1, mode="c"] Gal_HaloID     not None,
	np.ndarray[double, ndim=1, mode="c"] Gal_Prop     not None,
	np.ndarray[long, ndim=1, mode="c"] Gal_type       not None,
	double Cut,
	int NGal):
	Halo_Population_C(
		<long*>    np.PyArray_DATA(Halo_NGal_All),
		<long*>    np.PyArray_DATA(Halo_NGal_Cen),
		<long*>    np.PyArray_DATA(Halo_NGal_Sat1),
		<long*>    np.PyArray_DATA(Halo_NGal_Sat2),
		<long*>   np.PyArray_DATA(Halo_ID),
		<long*>   np.PyArray_DATA(Gal_HaloID),
		<double*> np.PyArray_DATA(Gal_Prop),
		<long*>    np.PyArray_DATA(Gal_type),
		Cut,NGal)

def HOD_Group(
	np.ndarray[double, ndim=1, mode="c"] HOD_Full       not None,
	np.ndarray[double, ndim=1, mode="c"] HOD_Cen        not None,
	np.ndarray[double, ndim=1, mode="c"] HOD_Sat1       not None,
	np.ndarray[double, ndim=1, mode="c"] HOD_Sat2       not None,
	np.ndarray[double, ndim=1, mode="c"] HOD_NHalo       not None,
	np.ndarray[long, ndim=1, mode="c"] Halo_NGal_Full not None,
	np.ndarray[long, ndim=1, mode="c"] Halo_NGal_Cen  not None,
	np.ndarray[long, ndim=1, mode="c"] Halo_NGal_Sat1 not None,
	np.ndarray[long, ndim=1, mode="c"] Halo_NGal_Sat2 not None,
	np.ndarray[double, ndim=1, mode="c"] Halo_Mass    not None,
	int NHalos, double Xmin, double Xmax, double NBin):
	HOD_Group_C(
		<double*>    np.PyArray_DATA(HOD_Full),
		<double*>    np.PyArray_DATA(HOD_Cen),
		<double*>    np.PyArray_DATA(HOD_Sat1),
		<double*>    np.PyArray_DATA(HOD_Sat2),
		<double*>    np.PyArray_DATA(HOD_NHalo),
		<long*>    np.PyArray_DATA(Halo_NGal_Full),
		<long*>    np.PyArray_DATA(Halo_NGal_Cen),
		<long*>    np.PyArray_DATA(Halo_NGal_Sat1),
		<long*>    np.PyArray_DATA(Halo_NGal_Sat2),
		<double*>  np.PyArray_DATA(Halo_Mass),
		NHalos, Xmin, Xmax, NBin)


