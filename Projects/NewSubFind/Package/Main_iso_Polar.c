#include <stdio.h>
#include <math.h>

// Global values:
double h,rho_c,G,pi;
// End of Global Values

// Basic Package
int imax(int a, int b){
	if (a < b) return b;
	return a;
}
double max(double a, double b){
	if (a < b) return b;
	return a;
}
double min(double a, double b){
	if (a > b) return b;
	return a;
}
double Dist3D(double * X, double * Y){
	return sqrt(pow((X[0]-Y[0]),2) + pow((X[1]-Y[1]),2) + pow((X[2]-Y[2]),2));
}
double Module3D(double * X){
	return sqrt(X[0]*X[0]+X[1]*X[1]+X[2]*X[2]);
}
void Cross_Product(double * X, double * Y, double * Z){
	Z[0] = X[1]*Y[2]-X[2]*Y[1];
	Z[1] = X[2]*Y[0]-X[0]*Y[2];
	Z[2] = X[0]*Y[1]-X[1]*Y[0];
}
// End Basic Package
void init_gobal_param(){
	h = 0.7; //0.76 Millennium Cosmology
	rho_c = 1.8788e-26*h*h * (h/1.989e30) / pow((3.24078e-23*h),3); // h factor will disapear, I just let it clear this are in Msun/h & Mpc/h
	G = 6.67408e-11 * pow((3.24078e-23*h),3) / (h/1.989e30) / pow((1e-9*3.17098e-8),2);
	pi = 3.1415926535897932384626433832795028841971693993751;
}



double c(double M){
	return 4.67*pow(M/1e14,-0.11); // Eq 11 Gan+10
}
double g(double x){
	return log(1.+x)-x/(1.+x);
}
double Rvir(double M){
    return pow(M*3/(800*pi*rho_c),(1./3.));
}
double Vvir(double M){
	return sqrt(G*M/Rvir(M)); //CAREFULL! in h^-1 MPC/s
}
double rs(double M){ //Virial values
	return Rvir(M)/c(M);
}
double M_in(double r,double M){
	return pow(Vvir(M),2)*r/G;
}
double Rm(double M, double Ms){
	return M/(Ms+1e-10);
}
double lnLamnda(double M,double Ms){
	return log(1.+Rm(M,Ms));
}
double rho(double r,double M){
	return pow(Vvir(M),2)/(4*pi*G*r*r);
}
double sigma(double r,double M){
	return Vvir(M)/sqrt(2); //Zentner 2003, eq. 6 
}
double Xi(double v,double r,double M){
	return v/(sqrt(2.)*sigma(r,M));
}


void init_PosVel(double M_c, double M_s, double * Sat_PosVel,double epsilon){
	Sat_PosVel[0] = Rvir(M_c); //Radial Distance
	Sat_PosVel[1] = 0; //Angle
	Sat_PosVel[2] = -sqrt(1-epsilon*epsilon)*Vvir(M_c); //Radial Velocity
	Sat_PosVel[3] = epsilon*Vvir(M_c)/Sat_PosVel[0]; //Angular velocity
}

double Velocity(double * Sat_PosVel){
	return sqrt(Sat_PosVel[2]*Sat_PosVel[2]+pow((Sat_PosVel[3]*Sat_PosVel[0]),2));
}

double rt_func(double M_c, double M_s, double r_t, double * Sat_PosVel){
	double r = Sat_PosVel[0];
	return G*M_in(r_t,M_s)/(pow(Sat_PosVel[3],2)+G*(2*M_in(r,M_c)*pow(r,-3)-4*pi*rho(r,M_c)))-pow(r_t,3);
}

double rt(double M_c, double M_s, double a,double b,double * Sat_PosVel,double err){
	double c = (a+b)/2.;
	while ((b-a)/2. > err){
		if (rt_func(M_c,M_s,a,Sat_PosVel)*rt_func(M_c,M_s,c,Sat_PosVel) < 0) b = c;
		else a = c;
		c = (a+b)/2.;
	}
	return c;
}

// double dm_ds(double M_c, double M_s,double * Sat_Pos,double * Sat_Vel, double alpha){
// 	double r_t = rt(M_c,M_s,1e-5,4*Rvir(M_s),Sat_Pos, Sat_Vel, 1e-5);
// 	double Torb = 2*pi/omega(Sat_Pos,Sat_Vel);
// 	return -alpha*(M_s-M_in(r_t,M_s))/Torb;
// }
/* Taylor & Babul (2001) assume no change in the internal density profile of the subhalo*/
double dm_ds(double M_c, double M_s, double Ms0,double * Sat_PosVel, double alpha){
	double r_t = rt(M_c,Ms0,1e-5,4*Rvir(M_s),Sat_PosVel, 1e-5);
	double Torb = 2*pi/Sat_PosVel[3];
	return -max(alpha*(M_s-M_in(r_t,Ms0))/Torb,0);
}
// double dm_ds(double M_c, double M_s, double Ms0,double * Sat_PosVel, double alpha){
// 	double r_t = rt(M_c,Ms0,1e-5,4*Rvir(M_s),Sat_PosVel, 1e-5);
// 	double Torb = 2*pi/Sat_PosVel[3];
// 	return -max(alpha*(M_s-Min_Sat(r_t,M_s,Ms0))/Torb,0);
// }
double A_g(double M_c, double * Sat_PosVel){
	double r = Sat_PosVel[0];
	return G*M_in(r,M_c)/r/r;
}
double A_d(double M_c,double M_s, double * Sat_PosVel){//Use M_s or Ms0?? CAREFULL!!!!!!
	double r = Sat_PosVel[0];
	double vt = Velocity(Sat_PosVel);
	if (vt == 0) return 0;
	double x = Xi(vt,r,M_c);
	return 4*pi*G*G*M_s*lnLamnda(M_c,M_s)*rho(r,M_c)*(erf(x) - 2*x/sqrt(pi)*exp(-x*x))/vt/vt;
}  

double RK4_integrator(double * Arr, double M_c, double M_s, double * Sat_PosVel, double jmin, double dt, double alpha, double MAX,  double dt_print, int verbose){
	int dim , k, rk, i=0, l = 0, Narr = (int) (MAX/dt_print);
	double r,j,Lim_rt,rp,Vtot,Msat,Ag,Ad,theta,t = 0;
		
	double r0  = Sat_PosVel[0];
	double j0  = Sat_PosVel[3]*Sat_PosVel[0]*Sat_PosVel[0];
	double Ms0 = M_s;
	
	j = j0;
	r = r0;
	theta = 0;
	
	double rk_mod[4]  = {0,0.5,0.5,1};
	double k_d[4][4]  = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
	double k_dm[4]    = {0,0,0,0};
	double SatPV[4]   = {0,0,0,0};
	
	while(j / j0 > jmin && M_s > 0 && r / r0 > 1e-5){
		
		if ( MAX <= t){
			printf("Warning, %lf>%lf, integration did not converge on MAX range\n",t,MAX);
			return MAX;
		}
		
		if (i % (int) (dt_print/dt) == 0 && l < Narr){
			if (verbose) printf("%lf\t%lf\t\t%G\t%G\t%G\t%G\t%G \n",t,j/j0,Sat_PosVel[0],Sat_PosVel[1],Sat_PosVel[2],Sat_PosVel[3],M_s/Ms0);
			Arr[l] = t;
			Arr[l+Narr] = M_s;
			Arr[l+2*Narr] = j;
			l++;
		}
		
		// RK4
		for (rk = 0; rk < 4; rk++){
				
// 			SatPV[4]    = (0,0,0,0);
		
			for (dim=0; dim < 4; dim++)	SatPV[dim]    = Sat_PosVel[dim]+k_d[dim][imax(0,rk-1)]*dt*rk_mod[rk];
			Msat = M_s+k_dm[imax(0,rk-1)]*dt*rk_mod[rk];
			
			Ag         = A_g(M_c,Sat_PosVel);
			Ad         = A_d(M_c,Msat,Sat_PosVel); 
// 			Ad = 0;
			
			theta = acos(fabs(SatPV[2]/Velocity(SatPV)));
			k_d[0][rk]  = SatPV[2];
			k_d[1][rk]  = SatPV[3];
			if (SatPV[2] == 0) k_d[2][rk]  = pow(SatPV[3],2)*SatPV[0]-Ag;
			else k_d[2][rk]  = pow(SatPV[3],2)*SatPV[0]-Ag-Ad*cos(theta)*SatPV[2]/fabs(SatPV[2]);
			if (SatPV[3] == 0) k_d[3][rk]  = -2*SatPV[2]*SatPV[3]/SatPV[0];
			else k_d[3][rk]  = -Ad*sin(theta)*SatPV[3]/fabs(SatPV[3])/SatPV[0]- 2*SatPV[2]*SatPV[3]/SatPV[0];
			if (alpha>0) k_dm[rk] = dm_ds(M_c, M_s, Ms0, Sat_PosVel, alpha);
			
// 			k_d[2][rk]  = SatPV[3]**2*SatPV[0]-Ag-Ad*np.cos(theta)
//             if SatPV[3] == 0: k_d[3][rk]  = -2*SatPV[2]*SatPV[3]/SatPV[0]
//             else: k_d[3][rk]  = -Ad*np.sin(theta)*SatPV[3]/abs(SatPV[3])/SatPV[0]- 2*SatPV[2]*SatPV[3]/SatPV[0]
			
		}
		if (alpha>0) M_s += dt/6.*(k_dm[0]+2*k_dm[1]+2*k_dm[2]+k_dm[3]);
		
		for (dim=0; dim < 4; dim++)	Sat_PosVel[dim] += dt/6.*(k_d[dim][0] +2*k_d[dim][1] +2*k_d[dim][2] +k_d[dim][3]);
// 		printf("%lf\t%lf\t\t%G\t%G\t%G\t\t%G\t%G \n",t,j/j0,Sat_PosVel[0],Sat_PosVel[1],Sat_PosVel[2],Sat_PosVel[3],M_s/Ms0);		
		t += dt;
		i += 1;
		r = Sat_PosVel[0];
		j = Sat_PosVel[3]*Sat_PosVel[0]*Sat_PosVel[0];
	}
	
	return t;
}

double Tdf_C(double M_c, double M_s, double jmin,double epsilon,double alpha, double dt, double MAX, double dt_print, int verbose){
	int i = 0;
	double * Arr;
	Arr = (double*)calloc((int) 3*(MAX/dt_print), sizeof(double));
	for (i=0; i < 3*(MAX/dt_print); i++) Arr[i] = -99;
	init_gobal_param();
    double Sat_PosVel[4] = {0,0,0,0};
	init_PosVel(M_c, M_s, Sat_PosVel, epsilon);
	double Tdf = RK4_integrator(Arr, M_c, M_s, Sat_PosVel, jmin, dt, alpha, MAX, dt_print, verbose);
	
	if (verbose) printf("%G\n",Tdf);
	return Tdf;
}


double * jM_eV_C(double M_c, double M_s, double jmin,double epsilon,double alpha, double dt, double MAX, double dt_print, int verbose){
	int i = 0;
	double * Arr;
	Arr = (double*)calloc((int) 3*(MAX/dt_print), sizeof(double));
	for (i=0; i < 3*(MAX/dt_print); i++) Arr[i] = -99;
	init_gobal_param();
	double Sat_PosVel[4] = {0,0,0,0};
	init_PosVel(M_c, M_s, Sat_PosVel, epsilon);
	double Tdf = RK4_integrator(Arr, M_c, M_s, Sat_PosVel, jmin, dt, alpha, MAX, dt_print, verbose);
	
	if (verbose) printf("%G\n",Tdf);
	return Arr;
}
