void Vel_1D_C(double * X, double * Y, double * Z, double * VZ, double * gv,double * gv_pc, double * gg,double * gg_pc, long * Prim_ID , long * Sec_ID, long * Main_Des_ID , int Prim_N, int Sec_N,int N, 

double NBin, int iNBin, double LIMIT, double Xmax, double Xmin, double Lbox, int NCPU, int CPU);

void Halo_Population_C(long * Halo_NGal_All,long * Halo_NGal_Cen,long * Halo_NGal_Sat1,long * Halo_NGal_Sat2, long * Halo_ID,long * Gal_HaloID, double * Gal_Prop,long * Gal_type, double Cut, int NGal);

void HOD_Group_C(double * HOD_Full, double * HOD_Cen, double * HOD_Sat1, double * HOD_Sat2,double * HOD_NHalo,long * Halo_NGal_Full,long * Halo_NGal_Cen,long * Halo_NGal_Sat1,long * Halo_NGal_Sat2, double * Halo_Mass, int NHalos, double Xmin, double Xmax, double NBin);

void Gaussian_Smooth_CIC_C(double * A,double * B, double sigma, int sigmaMax, int N1, int N2, int N3);

void ACF_DD_C(double * X, double * Y, double * Z, double * JN_Random, double * gg,double N, double NBin, int iNBin, double LIMIT, double Xmin, double Xmax, double Lbox,int JN, int NCPU, int CPU);
