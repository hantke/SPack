#include <math.h>


double Line_C(double x,double X1,double X2,double Y1,double Y2){
	return (Y2-Y1)/(X2-X1)*(x-X1)+Y1;
}

double FuncSingle_C(double x0,double * X,double * Y,int N){
	int i;
	if (x0 < X[0])	return Line_C(x0,X[0],X[1],Y[0],Y[1]);
	for(i=1;i<N;i++)	if (x0 < X[i])	return Line_C(x0,X[i-1],X[i],Y[i-1],Y[i]);
	return Line_C(x0,X[N-2],X[N-1],Y[N-2],Y[N-1]);
}

double * FuncArr_C(double * X0,double * X,double * Y,int N,int N_Arr){
	int i,j;
	double * Arr;
	Arr = (double*) calloc (N_Arr,sizeof(double));
// 	double Arr[N_Arr];
	for(j=1;j<N_Arr;j++){
		if (X0[j] < X[0])	Arr[j] = Line_C(X0[j],X[0],X[1],Y[0],Y[1]);
		else if	(X0[j] > X[N-1]) Arr[j] = Line_C(X0[j],X[N-2],X[N-1],Y[N-2],Y[N-1]);
		else	for(i=1;i<N;i++)	if (X0[j] < X[i]){
			Arr[j] = Line_C(X0[j],X[i-1],X[i],Y[i-1],Y[i]);
			i = N;
		}
	}
	return Arr;
}
