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
double rs(double M){ //Virial values
	return Rvir(M)/c(M);
}
double M_in(double r,double M){
	return M*g(r/rs(M))/g(c(M));// Eq 10 Gan+10
}
double Rm(double M, double Ms){
	return M/(Ms+1e-10);
}
double xm(double Ms, double Ms0){
	return log10(Ms/Ms0);
}
double r_te(double x, double Ms0){
	return pow(10,1.02+1.38*x+0.37*x*x)*rs(Ms0);
}
double ft(double x){
	return pow(10,-0.007+0.35*x+0.39*x*x+0.23*x*x*x);
}
double lnLamnda(double M,double Ms){
	return log(1.+Rm(M,Ms));
}
double rho_0(double M){
	return 200.*pow(c(M),3)/(3*g(c(M)));
}
double rho_nfw(double r,double M){
	return rho_0(M)*rho_c/(r/rs(M))/pow((1+r/rs(M)),2);
}
double rho_sat(double r, double Ms, double Ms0){
	double x = xm(Ms,Ms0);
	return rho_nfw(r,Ms0)*ft(x)/(1+pow(r/r_te(x,Ms0),3));
}
double Min_Sat(double r, double Ms, double Ms0){
	int i=0, Nint = 50;
	double dr,dv,rmin,rmed,rmax, Mass = 0;
	dr = r/Nint;
	for(i=0;i<Nint;i++){
		
		rmin = i*dr;
		rmax = (i+1)*dr;
		rmed = (i+0.5)*dr;
		Mass += 4./3.*pi*(pow(rmax,3)-pow(rmin,3))*rho_sat(rmed,Ms,Ms0);
	}
	return Mass;
}
double Vvir(double M){
	return sqrt(G*M/Rvir(M)); //CAREFULL! in h^-1 MPC/s
}
double Vmax(double M){
	return Vvir(M)*sqrt(0.216*c(M)/(g(c(M))));// #Prada 2011, eq. 9
}
double sigma(double r,double M){
	double x = r/rs(M);
	return Vmax(M)*(1.4393*pow(x,0.354))/(1+1.1756*pow(x,0.725)); //Zentner 2003, eq. 6 
}
double Xi(double v,double r,double M){
	return v/(sqrt(2.)*sigma(r,M));
}
//

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
	return G*M_in(r_t,M_s)/(pow(Sat_PosVel[3],2)+G*(2*M_in(r,M_c)*pow(r,-3)-4*pi*rho_nfw(r,M_c)))-pow(r_t,3);
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
	double r_t = rt(M_c,Ms0,1e-5,4*Rvir(Ms0),Sat_PosVel, 1e-5);
	double Torb = 2*pi/Sat_PosVel[3];
	return -max(alpha*(M_s-M_in(r_t,Ms0))/Torb,0);
}
// double dm_ds(double M_c, double M_s, double Ms0,double * Sat_PosVel, double alpha){
// 	double r_t = rt(M_c,Ms0,1e-5,4*Rvir(Ms0),Sat_PosVel, 1e-5);
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
	return 4*pi*G*G*M_s*lnLamnda(M_c,M_s)*rho_nfw(r,M_c)*(erf(x) - 2*x/sqrt(pi)*exp(-x*x))/vt/vt;
}  

double RK4_integrator(double * Arr, double M_c, double M_s, double * Sat_PosVel, double dt, double alpha, double MAX,  double dt_print, double sigma_r, double sigma_rp, double sigma_j, int verbose){
	int dim , k, rk, i=0, l = 0, Narr = (int) (MAX/dt_print);
	double r,j,Lim_rt,rp,Vtot,Msat,Ag,Ad,theta,t = 0;
	
	double v0  = Vvir(M_c);
	double r0  = Rvir(M_c);
	double j0  = v0*r0;
	double Ms0 = M_s;
	
	j = Sat_PosVel[3]*Sat_PosVel[0]*Sat_PosVel[0];
	theta = 0;
	
	double rk_mod[4]  = {0,0.5,0.5,1};
	double k_d[4][4]  = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
	double k_dm[4]    = {0,0,0,0};
	double SatPV[4]   = {0,0,0,0};
	
	while(sqrt( pow((j/j0/sigma_j),2) + pow((Sat_PosVel[0]/r0/sigma_r),2) + pow((Sat_PosVel[2]/v0/sigma_rp),2) ) > 1 && M_s > 0 ){
		
		if ( MAX <= t){
			printf("Warning, %lf>%lf, integration did not converge on MAX range\n",t,MAX);
			return MAX;
		}
		
		if (i % (int) (dt_print/dt) == 0 && l < Narr){
			if (verbose) printf("%lf\t%lf\t\t%G\t%G\t%G\t%G\t%G \n",t,j/j0,Sat_PosVel[0],Sat_PosVel[1],Sat_PosVel[2],Sat_PosVel[3],M_s/Ms0);
			Arr[l] = t;
			Arr[l+Narr] = Sat_PosVel[0];
			Arr[l+2*Narr] = Sat_PosVel[1];
			Arr[l+3*Narr] = Sat_PosVel[2];
			Arr[l+4*Narr] = Sat_PosVel[3];
			Arr[l+5*Narr] = M_s;
			l++;
		}
		
		// RK4
		for (rk = 0; rk < 4; rk++){
				
// 			SatPV[4]    = (0,0,0,0);
		
			for (dim=0; dim < 4; dim++)	SatPV[dim]    = Sat_PosVel[dim]+k_d[dim][imax(0,rk-1)]*dt*rk_mod[rk];
			Msat = M_s+k_dm[imax(0,rk-1)]*dt*rk_mod[rk];
			
			Ag         = A_g(M_c,Sat_PosVel);
			Ad         = A_d(M_c,Msat,Sat_PosVel); 
// // 			Ad = 0;
			
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
		j = Sat_PosVel[3]*Sat_PosVel[0]*Sat_PosVel[0];
	}
	
	return t;
}

double Tdf_C(double M_c, double M_s,double epsilon,double alpha, double eta, double dt, double MAX, double dt_print, double sigma_r, double sigma_rp, double sigma_j, int verbose){
	int i = 0;
	double * Arr;
	Arr = (double*)calloc((int) 6*(MAX/dt_print), sizeof(double));
	for (i=0; i < 6*(MAX/dt_print); i++) Arr[i] = -99;
	init_gobal_param();
    double Sat_PosVel[4] = {0,0,0,0};
	init_PosVel(M_c, M_s, Sat_PosVel, epsilon);
	double Tdf = RK4_integrator(Arr, M_c, M_s, Sat_PosVel, dt, alpha, MAX, dt_print, sigma_r, sigma_rp, sigma_j, verbose);
	
	if (verbose) printf("%G\n",Tdf);
	return Tdf;
}


double * jM_eV_C(double M_c, double M_s,double epsilon,double alpha, double eta, double dt, double MAX, double dt_print, double sigma_r, double sigma_rp, double sigma_j, int verbose){
	int i = 0;
	double * Arr;
	Arr = (double*)calloc((int) 6*(MAX/dt_print), sizeof(double));
	for (i=0; i < 6*(MAX/dt_print); i++) Arr[i] = -99;
	init_gobal_param();
	double Sat_PosVel[4] = {0,0,0,0};
	init_PosVel(M_c, M_s, Sat_PosVel, epsilon);
	double Tdf = RK4_integrator(Arr, M_c, M_s, Sat_PosVel, dt, alpha, MAX, dt_print, sigma_r,  sigma_rp, sigma_j, verbose);
	
	if (verbose) printf("%G\n",Tdf);
	return Arr;
}

