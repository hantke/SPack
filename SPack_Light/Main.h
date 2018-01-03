int Index(double x, double Xmin, double Xmax, int NBin);

double Line_C(double x,double X1,double X2,double Y1,double Y2);

double FuncSingle_C(double x0,double * X,double * Y,int N);

double * FuncArr_C(double * X0,double * X,double * Y,int N,int N_Arr);

double * AcumMassFunction_C(double * M, double Mmin, double Mmax, double Volume, int NBin, int NGal);

double * MassFunction_C(double * M, double Mmin, double Mmax, double Volume, int NBin, int NGal);

void * Histo2D_C(long * Arr, double * X, double * Y,double Xa_min, double Xa_max,double Xb_min, double Xb_max, int NBin_a, int NBin_b, long N);

double Dist3D_C(double X1,double Y1,double Z1,double X2,double Y2,double Z2);

double LogDist3D_C(double X1,double Y1,double Z1,double X2,double Y2,double Z2);
