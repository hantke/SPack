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



// TEST THESE ONES!

long * Histo_C(double * X,double Xmin, double Xmax, int NBin, long N){
	int i,index;
	long * Arr;
	Arr = (long*) calloc (NBin,sizeof(long));
	
	for (i=0;i<NBin;i++) Arr[i] = 0;
	
	for (i=0;i<N;i++){
		index = (int) ((X[i] - Xmin)/(Xmax-Xmin)*(double) NBin);
		if ((index > -1) && (index < NBin))	Arr[index]++;
	}
	
	return Arr;
}

long * Histo2D_C(double * X,double Xa_min, double Xa_max,double Xb_min, double Xb_max, int NBin_a, int NBin_b, long N){
	int i,index,index_a,index_b;
	int NBin = NBin_a*NBin_b
	long * Arr;
	Arr = (long*) calloc (NBin,sizeof(long));
	
	for (i=0;i<NBin;i++) Arr[i] = 0;
	
	for (i=0;i<N;i++){
		index_a = (int) ((X[i] - Xa_min)/(Xa_max-Xa_min)*(double) NBin_a);
		index_b = (int) ((X[i] - Xb_min)/(Xb_max-Xb_min)*(double) NBin_b);
		index = index_a*NBin_b+NBin_b;
		if ((index > -1) && (index < NBin))	Arr[index]++;
	}
	
	return Arr;
}
