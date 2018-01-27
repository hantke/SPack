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
double lnLamnda(double M,double Ms){
	return log(1.+Rm(M,Ms));
}
double rho_0(double M){
	return 200.*pow(c(M),3)/(3*g(c(M)));
}
double rho(double r,double M){
	return rho_0(M)*rho_c/(r/rs(M))/pow((1+r/rs(M)),2);
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

void init_PosVel(double M_c, double M_s, double * Sat_Pos,double * Sat_Vel,double epsilon){
	Sat_Pos[0] = Rvir(M_c);
	Sat_Pos[1] = 0;
	Sat_Pos[2] = 0;
	Sat_Vel[0] = -sqrt(1-epsilon*epsilon)*Vvir(M_c);
	Sat_Vel[1] = epsilon*Vvir(M_c);
	Sat_Vel[2] = 0;
}

double omega(double * Sat_Pos,double * Sat_Vel){
	double r_mod = Module3D(Sat_Pos);
	double Pos_vec[3] = {Sat_Pos[0]/r_mod,Sat_Pos[1]/r_mod,Sat_Pos[2]/r_mod};
	
	double cross[3] = {0,0,0};
	Cross_Product(Pos_vec,Sat_Vel,cross);
	
    return Module3D(cross)/r_mod;
}

double rt_func(double M_c, double M_s, double r_t, double * Sat_Pos, double * Sat_Vel, double r){
	return G*M_in(r_t,M_s)/(pow(omega(Sat_Pos,Sat_Vel),2)+G*(2*M_in(r,M_c)*pow(r,-3)-4*pi*rho(r,M_c)))-pow(r_t,3);
}

double rt(double M_c, double M_s, double a,double b,double * Sat_Pos,double * Sat_Vel,double err){
	double r_mod = Module3D(Sat_Pos);
	double c = (a+b)/2.;
	while ((b-a)/2. > err){
		if (rt_func(M_c,M_s,a,Sat_Pos,Sat_Vel,r_mod)*rt_func(M_c,M_s,c,Sat_Pos,Sat_Vel,r_mod) < 0) b = c;
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
double dm_ds(double M_c, double M_s, double Ms0,double * Sat_Pos,double * Sat_Vel, double alpha){
	double r_t = rt(M_c,Ms0,1e-5,4*Rvir(M_s),Sat_Pos, Sat_Vel, 1e-5);
	double Torb = 2*pi/omega(Sat_Pos,Sat_Vel);
	return -max(alpha*(M_s-M_in(r_t,Ms0))/Torb,0);
}

double A_g(double M_c, double r){
	return G*M_in(r,M_c)/r/r;
}
double A_d(double M_c,double M_s,double r,double vt){
	if (vt == 0) return 0;
	double x = Xi(vt,r,M_c);
	return 4*pi*G*G*M_s*lnLamnda(M_c,M_s)*rho(r,M_c)*(erf(x) - 2*x/sqrt(pi)*exp(-x*x))/vt/vt;
}  

double RK4_integrator(double M_c, double M_s, double * Sat_Pos, double * Sat_Vel, double jmin, double dt, double alpha, double MAX, int verbose){
	int dim , k, rk, i=0;
	double r,j,Lim_rt,rp,Vtot,Msat,Ag,Ad,t = 0;
	
	double j_vec[3] = {0,0,0};
	Cross_Product(Sat_Pos,Sat_Vel,j_vec);
	
	double r0  = Module3D(Sat_Pos);
	double j0  = Module3D(j_vec);
	double Ms0 = M_s;
	
	j = j0;
	r = r0;
	Lim_rt = Rvir(M_c);
	
	double rk_mod[4]  = {0,0.5,0.5,1};
	double k_d[3][4]  = {{0,0,0,0},{0,0,0,0},{0,0,0,0}};
	double k_dd[3][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0}};
	double k_dm[4]    = {0,0,0,0};
	double SatV[3]    = {0,0,0};
	double SatP[3]    = {0,0,0};
	
	while(j / j0 > jmin && M_s > 0){
		
		if (i == MAX-1){
			printf("Warning, i=MAX, integration did not converge on MAX range\n");
			return 0;
		}
		
		// RK4
		for (rk = 0; rk < 4; rk++){
				
			SatV[3]    = (0,0,0);
			SatP[3]    = (0,0,0);
		
			for (dim=0; dim < 3; dim++){
				SatV[dim]    = Sat_Vel[dim]+k_dd[dim][imax(0,rk-1)]*dt*rk_mod[rk];
				SatP[dim]    = Sat_Pos[dim]+k_d[dim][imax(0,rk-1)] *dt*rk_mod[rk];
			}
			rp   = Module3D(SatP);
			Vtot = Module3D(SatV);
			
			Msat       = M_s+k_dm[imax(0,rk-1)]*dt*rk_mod[rk];
			Ag         = A_g(M_c,rp);
			Ad         = A_d(M_c,Msat,rp,Vtot); 
			for (dim=0; dim < 3; dim++){
				k_dd[dim][rk] = -Ag*SatP[dim]/rp-Ad*SatV[dim]/Vtot;
				k_d[dim][rk]  = SatV[dim];
			}
			if (alpha>0) k_dm[rk] = dm_ds(M_c,Msat,Ms0,SatP,SatV,alpha);
		}
		if (alpha>0) M_s += dt/6.*(k_dm[0]+2*k_dm[1]+2*k_dm[2]+k_dm[3]);
		
		for (dim=0; dim < 3; dim++){
			Sat_Vel[dim] += dt/6.*(k_dd[dim][0]+2*k_dd[dim][1]+2*k_dd[dim][2]+k_dd[dim][3]);
			Sat_Pos[dim] += dt/6.*(k_d[dim][0] +2*k_d[dim][1] +2*k_d[dim][2] +k_d[dim][3]);
		}
		
		if (verbose && i % (int) (0.1/dt) == 0) printf("%lf\t%lf\t%lf\t\t%G\t%G\t%G\t\t%G\t%G\t%G\t%G \n",t,r,j/j0,Sat_Pos[0],Sat_Pos[1],Sat_Pos[2],Sat_Vel[0],Sat_Vel[1],Sat_Vel[2],M_s/Ms0);
		
		Cross_Product(Sat_Pos,Sat_Vel,j_vec);
		t += dt;
		i += 1;
		j = Module3D(j_vec);
		r = Module3D(Sat_Pos);
	}
	
	return t;
}

double Tdf_NFW_C(double M_c, double M_s, double jmin,double epsilon,double alpha, double dt, double MAX, int verbose){
	init_gobal_param();
    double Sat_Pos[3] = {0,0,0};
	double Sat_Vel[3] = {0,0,0};
	init_PosVel(M_c, M_s, Sat_Pos, Sat_Vel, epsilon);
	double r = Module3D(Sat_Pos);
	double v = Module3D(Sat_Vel);
	double Tdf = RK4_integrator(M_c, M_s, Sat_Pos, Sat_Vel, jmin, dt, alpha, MAX, verbose);
	
	if (verbose) printf("%G\n",Tdf);
	return Tdf;
}

