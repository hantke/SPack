#include <math.h>
// Extra Tools
int Index_C(double x, double Xmin, double Xmax, int NBin){
    return (int) (NBin* (x-Xmin) / (Xmax-Xmin));
}
double Line_C(double x,double X1,double X2,double Y1,double Y2){
	return (Y2-Y1)/(X2-X1)*(x-X1)+Y1;
}
int imin(int a,int b){
	if (a<b) return a;
	return b;
}
int imax(int a,int b){
	if (a>b) return a;
	return b;
}
float min(float a,float b){
	if (a<b) return a;
	return b;
}
float max(float a,float b){
	if (a>b) return a;
	return b;
}
int iLimit(int x,int iNBin){
	if (x < 0)	return x+iNBin;
	if (x > iNBin-1)	return x-iNBin;
	return x;
	
}
//

/////// Find elment of a table A in a table B

long Find_ainA_C(int a, long * A, int a0, int a1, int LenA){
	if(a1-a0<3){
		int t=a0;
		while(t <= a1 && t < LenA){
			if (a == A[t]) return t;
			t++;
		}
		return -99;
	}
	int med = (a1+a0)/2;
		
		
	
	
}

// def Find_ainA(a,A,a0,a1):#Look for a in Y-Table, initial a0,a1 = 0, len(A-1)
// 	if a1 - a0 < 3:
// 		t = a0
// 		while t <= a1 and t < len(A):
// 			if a == A[t]:	return t
// 			t += 1
// 		return -99	
// 	med = (a1+a0)/2
// 	if a < A[med]:	return Find_ainA(a,A,a0,med)
// 	if a > A[med]:	return Find_ainA(a,A,med,a1)
// 	return med	
// 
// def Find_AinB(A, B):# Only if A is complete in B
// 	b=max(0,Find_ainA(A[0],B,0,len(B)-1)-1)
// 	a = 0
// 	#b=0
// 	index = []
// 	while a<len(A):
// 		if (A[a] < B[b]) or (len(B)<b):
// 			print 'warning', a, b, A[a], B[b]
// 			return -1
// 		#if b%1000000==0: print a, b, A[a], B[b]
// 		if A[a]== B[b]:
// 			index.append(b)
// 			a+=1
// 		b+=1
// 	return np.array(index)

//////

//Used on Contreras et al. 2018 Zehavi et al. in Prep.
void Halo_Population_simple_C(long * Halo_NGal_All,long * Halo_NGal_Cen,long * Halo_NGal_Sat, double * Gal_Prop,long * Gal_type, double Cut, int NGal){	
	int i=0;
	int j=-1;
	for (i=0;i < NGal;i++){
		if (Gal_type[i] == 0) j++;
		if (Gal_Prop[i] >= Cut){
			if (Gal_type[i] == 0){
				Halo_NGal_Cen[j]++;
			}
			else	Halo_NGal_Sat[j]++;
			Halo_NGal_All[j]++;
		}
	}
}

void Halo_Population_C(long * Halo_NGal_All,long * Halo_NGal_Cen,long * Halo_NGal_Sat1,long * Halo_NGal_Sat2, long * Halo_ID,long * Gal_HaloID, double * Gal_Prop,long * Gal_type, double Cut, int NGal){	
	int i=0;
	int j=0;
	while (i < NGal){
		while (Halo_ID[j] < Gal_HaloID[i])	j++;
		if (Gal_Prop[i] >= Cut){
			Halo_NGal_All[j]++;
			if (Gal_type[i] == 0)	    Halo_NGal_Cen[j]++;
			else if (Gal_type[i] == 1)	Halo_NGal_Sat1[j]++;
			else                      	Halo_NGal_Sat2[j]++;
		}
		i++;
	}
}

void HOD_Group_C(double * HOD_Full, double * HOD_Cen, double * HOD_Sat1, double * HOD_Sat2,double * HOD_NHalo, long * Halo_NGal_Full,long * Halo_NGal_Cen,long * Halo_NGal_Sat1,long * Halo_NGal_Sat2, double * Halo_Mass, int NHalos, double Xmin, double Xmax, double NBin){
	int i,index;
	for(i = 0; i < NHalos; i++){
		
		index = Index_C(Halo_Mass[i], Xmin, Xmax, NBin);
		if (index > -1 && index < NBin){
			HOD_NHalo[index] = HOD_NHalo[index]+1;
			HOD_Full[index] = HOD_Full[index] + Halo_NGal_Full[i];
			HOD_Cen[index] = HOD_Cen[index] + Halo_NGal_Cen[i];
			HOD_Sat1[index] = HOD_Sat1[index] + Halo_NGal_Sat1[i];
			HOD_Sat2[index] = HOD_Sat2[index] + Halo_NGal_Sat2[i];
		}
	}
	for(index = 0; index < NBin; index++){
		HOD_Full[index] = HOD_Full[index] / (HOD_NHalo[index] + 1e-10);
		HOD_Cen[index] = HOD_Cen[index]   / (HOD_NHalo[index] + 1e-10);
		HOD_Sat1[index] = HOD_Sat1[index] / (HOD_NHalo[index] + 1e-10);
		HOD_Sat2[index] = HOD_Sat2[index] / (HOD_NHalo[index] + 1e-10);
		
	}
	
}

//Used on Contreras et al. 2018 Zehavi et al. in Prep.
void HOD_Group_simple_C(double * HOD_Full, double * HOD_Cen, double * HOD_Sat,double * HOD_NHalo, long * Halo_NGal_Full,long * Halo_NGal_Cen,long * Halo_NGal_Sat, double * Halo_Mass, int NHalos, double Xmin, double Xmax, double NBin){
	int i,index;
	for(i = 0; i < NHalos; i++){
		
		index = Index_C(Halo_Mass[i], Xmin, Xmax, NBin);
		if (index > -1 && index < NBin){
			HOD_NHalo[index] = HOD_NHalo[index]+1;
			HOD_Full[index] = HOD_Full[index] + Halo_NGal_Full[i];
			HOD_Cen[index] = HOD_Cen[index] + Halo_NGal_Cen[i];
			HOD_Sat[index] = HOD_Sat[index] + Halo_NGal_Sat[i];
		}
	}
	for(index = 0; index < NBin; index++){
		HOD_Full[index] = HOD_Full[index] / (HOD_NHalo[index] + 1e-10);
		HOD_Cen[index] = HOD_Cen[index]   / (HOD_NHalo[index] + 1e-10);
		HOD_Sat[index] = HOD_Sat[index] / (HOD_NHalo[index] + 1e-10);
		
	}
	
}

int from3Dto1Dindex_C(int i1,int i2, int i3, int N1, int N2, int N3){
	return i1 + i2*N1 + i3*N1*N2;
}

void Gaussian_Smooth_CIC_C(double * A,double * B, double sigma, int sigmaMax, int N1, int N2, int N3){
	int i,j,k,ii,jj,kk,imin,imax,jmin,jmax,kmin,kmax,index1,index2;
	float D2;
	float norm = sqrt(2*3.14159265);
	
	for (i=0;i<N1;i++){
		for(j=0;j<N2;j++){
			for(k=0;k<N3;k++){
				
				index1 = from3Dto1Dindex_C(i,j,k,N1,N2,N3);
				
				imin = (int) max(0,i-sigma*sigmaMax);
				jmin = (int) max(0,j-sigma*sigmaMax);
				kmin = (int) max(0,k-sigma*sigmaMax);
				
				imax = (int) min(N1,i+sigma*sigmaMax+1);
				jmax = (int) min(N2,j+sigma*sigmaMax+1);
				kmax = (int) min(N3,k+sigma*sigmaMax+1);
				if (A[index1] > 0){
					for (ii=imin;ii<imax;ii++){
						for(jj=jmin;jj<jmax;jj++){
							for(kk=kmin;kk<kmax;kk++){
								index2 = from3Dto1Dindex_C(ii,jj,kk,N1,N2,N3);
								D2 = (i-ii)*(i-ii)+(j-jj)*(j-jj)+(k-kk)*(k-kk);
								B[index2]+= A[index1]/(sigma*norm)*exp(-0.5*D2/sigma/sigma);
							}
						}
					}	
				}
			}
		}
	}
}

void ACF_DD_C(double * X, double * Y, double * Z, long * JN_Random, double * gg,double N, double NBin, int iNBin, double LIMIT, double Xmin, double Xmax, double Lbox,int JN, int NCPU, int CPU){
	printf("Main Information: %f %f %f %f %f %d %d %d\n\n",N,NBin,Xmax,Xmin,Lbox,NCPU,CPU,JN);
	/////////////////////////////// Init var
	int i,j,k,index,x,y,z,tx,ty,tz,iLIMIT2,iLIMIT3,JN_Index,first = -99;
	double LgDis,dx,dy,dz;
	double binsize_i = NBin/( Xmax - Xmin);
	float Lbox2 = Lbox/2;
	int iLIMIT = (int) (LIMIT/Lbox*iNBin);
	printf("%d\n",iLIMIT);
	int iMin = (int)( N*(CPU-1.)/ (float) NCPU);
	int iMax = (int)( N*(CPU)/ (float) NCPU);
	/////////////////////////////// 	Linking List
	int * ll;
	int * ix;
	int * iy;
	int * iz;
	ll = (int*) calloc (N,sizeof(int));
	ix = (int*) calloc (N,sizeof(int));
	iy = (int*) calloc (N,sizeof(int));
	iz = (int*) calloc (N,sizeof(int));
	int lfirst [iNBin][iNBin][iNBin];
	for (i = 0;i < iNBin; i++) for (j = 0;j < iNBin; j++) for (k = 0;k < iNBin; k++)	lfirst[i][j][k] = -99;
		for (i = 0;i < N; i++){
		ix[i] = (int) (X[i]/Lbox*iNBin);
		iy[i] = (int) (Y[i]/Lbox*iNBin);
		iz[i] = (int) (Z[i]/Lbox*iNBin);
		lfirst[ix[i]][iy[i]][iz[i]]=i;
	}
	for (i = 0;i < N; i++){
		ll[i] = lfirst[ix[i]][iy[i]][iz[i]];
		lfirst[ix[i]][iy[i]][iz[i]]=i;
	}
	/////////////////////////////// 	Linking List
	
	for(i = iMin; i < iMax; i++){
		JN_Index = JN_Random[i];
// 		if (i% (int) (N/100) == 0)	printf("%f\n",100.0* (float) i/(float) N );
		for(tx = ix[i]-iLIMIT; tx < ix[i]+iLIMIT+1;tx++){
			iLIMIT2 = (int) sqrt( (float) (iLIMIT*iLIMIT) - (ix[i] - tx)*(ix[i] - tx)  );
			for(ty = iy[i]-iLIMIT2; ty < iy[i]+iLIMIT2+1;ty++){
// 				iLIMIT3 = (int) sqrt( (float) (iLIMIT*iLIMIT) - (ix[i] - tx)*(ix[i] - tx) - (iy[i] - ty)*(iy[i] - ty));
				for(tz = iz[i]-iLIMIT2; tz < iz[i]+iLIMIT2+1;tz++){
					x = iLimit(tx,iNBin);
					y = iLimit(ty,iNBin);
					z = iLimit(tz,iNBin);
					j=lfirst[x][y][z];
					while (j != -99 && first != j){
						if (i>j){	
							dx = fabs( X[j] -  X[i]);
							dy = fabs( Y[j] -  Y[i]);
							dz = fabs( Z[j] -  Z[i]);
							if (dx > Lbox2)	dx = Lbox - dx;
							if (dy > Lbox2)	dy = Lbox - dy;
							if (dz > Lbox2)	dz = Lbox - dz;
							LgDis = log10(dx*dx+dy*dy+dz*dz+1e-10)/2.;
							index = (int) ((LgDis-Xmin)*binsize_i);
							if (index > -1 && index < NBin)	gg[index +JN_Index* (int) NBin] += 1;
						}
						j = ll[j];
						first = lfirst[x][y][z];
						
					}
				}
			}
		}
			
	}
}


// TEST!!!!!

void Vel_1D_C(double * X, double * Y, double * Z, double * VZ, double * gv,double * gv_pc, double * gg,double * gg_pc, long * Prim_ID , long * Sec_ID, long * Main_Des_ID , int Prim_N, int Sec_N,int N, double NBin, int iNBin, double LIMIT, double Xmax, double Xmin, double Lbox, int NCPU, int CPU){
/*	printf("Main Information: %d %f %f %f %f %d %d %d %d \n\n",N,NBin,Xmax,Xmin,Lbox,NCPU,CPU,iNBin);
	int i,j,k,i0,index,x,y,z,tx,ty,tz,dvx,dvy,dvz,iLIMIT2,JN_Index,index_dis,first = -99;
	/////////////////%TODO ADD INT LIST OF PRIM AND SEC
	int * ll;
	int * ix;
	int * iy;
	int * iz;
	ll = (int*) calloc (N,sizeof(int));
	ix = (int*) calloc (N,sizeof(int));
	iy = (int*) calloc (N,sizeof(int));
	iz = (int*) calloc (N,sizeof(int));
	int lfirst [iNBin][iNBin][iNBin];
	for (i = 0;i < iNBin; i++) for (j = 0;j < iNBin; j++) for (k = 0;k < iNBin; k++)	lfirst[i][j][k] = -99;
	for (i0 = 0;i0 < Sec_N; i0++){
		i = (int) Sec_ID[i0]; 
		ix[i] = (int) (X[i]/Lbox*iNBin);
		iy[i] = (int) (Y[i]/Lbox*iNBin);
		iz[i] = (int) (Z[i]/Lbox*iNBin);
		lfirst[ix[i]][iy[i]][iz[i]]=i;
	}
	for (i0 = 0;i0 < Sec_N; i0++){
		i = (int) Sec_ID[i0];
		ll[i] = lfirst[ix[i]][iy[i]][iz[i]];
		lfirst[ix[i]][iy[i]][iz[i]]=i;
	}
	double LgDis,dx,dy,dz,DX,DY,DZ,dl2;
	double binsize_i = NBin/( Xmax - Xmin);
	float Lbox2 = Lbox/2;
	int iLIMIT = (int) (LIMIT/Lbox*iNBin);
	for (i = 0; i < NBin*Prim_N; i++){
		gg[i] = 0;
		gv[i] = 0;
		gg_pc[i] = 0;
		gv_pc[i] = 0;
	}
	int iMin = (int)( Prim_N*(CPU-1.)/ (float) NCPU);
	int iMax = (int)( Prim_N*(CPU)/ (float) NCPU);
	for(i0 = iMin; i0 < iMax; i0++){
		i = (int) Prim_ID[i0];
// 		JN_Index = (int) (X[i] / Lbox * JN);
// 		if (i% (int) (N/100) == 0)	printf("%d,%d,%d,%d\n",i,i0,Prim_N,iMin );
		if (i0% (int) (Prim_N/100) == 0)	printf("%f\n",100.0* (float) i0/(float) Prim_N );
		for(tx = ix[i]-iLIMIT; tx < ix[i]+iLIMIT+1;tx++){
			iLIMIT2 = (int) sqrt( (float) (iLIMIT*iLIMIT) - (ix[i] - tx)*(ix[i] - tx)  );
			for(ty = iy[i]-iLIMIT2; ty < iy[i]+iLIMIT2+1;ty++){
				for(tz = iz[i]-iLIMIT2; tz < iz[i]+iLIMIT2+1;tz++){
					x = iLimit(tx,iNBin);
					y = iLimit(ty,iNBin);
					z = iLimit(tz,iNBin);
					j=lfirst[x][y][z];
					
					
					while (j != -99 && first != j){
						
						DX =  X[j] -  X[i];
						DY =  Y[j] -  Y[i];
						DZ =  Z[j] -  Z[i];
						
                        dvz = VZ[j] - VZ[i];
                        
						dx = fabs( DX);
						dy = fabs( DY);
						dz = fabs( DZ);
						
						if (dx > Lbox2)	dx = Lbox - dx;
						if (dy > Lbox2)	dy = Lbox - dy;
						if (dz > Lbox2)	dz = Lbox - dz;
						
						
						dl2 = dx*dx+dy*dy+dz*dz;

                        if (( DZ < 0) ^ (fabs( DZ)  > Lbox2))	dvz = -dvz;

						LgDis = log10(dl2+1e-10)/2.;
						index = (int) ((LgDis-Xmin)*binsize_i);
                        
						
						if (index > -1 && index < NBin){
							gg[index +(int) NBin*i0] += 1;
							gv[index +(int) NBin*i0] += dvz;
                            if (Main_Des_ID[i] == Main_Des_ID[j] && Main_Des_ID[i] != -99){
                                gg_pc[index +(int) NBin*i0] += 1;
                                gv_pc[index +(int) NBin*i0] += dvz;
                            }
						}
						j = ll[j];
						first = lfirst[x][y][z];
					}
				}
			}
		}
			
	}*/
}
/*
long * Histo2D_C(double * X,double Xa_min, double Xa_max,double Xb_min, double Xb_max, int NBin_a, int NBin_b, long N){
	int i,index,index_a,index_b;
	int NBin = NBin_a*NBin_b;
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
*/
