#include "C_xi.h"

#include "gsl/gsl_spline.h"
#include <math.h>

#define rhomconst 2.775808e+11
//1e4*3.*Mpcperkm*Mpcperkm/(8.*PI*G); units are SM h^2/Mpc^3

double xi_nfw_at_R(double R, double Mass, double conc, int delta, double om){
  double rhom = om*rhomconst;//SM h^2/Mpc^3
  double Rdelta = pow(Mass/(1.33333333333*M_PI*rhom*delta), 0.33333333333);
  double Rscale = Rdelta/conc;
  double fc = log(1.+conc)-conc/(1.+conc); 
  return Mass/(4.*M_PI*Rscale*Rscale*Rscale*fc)/(R/Rscale*(1+R/Rscale)*(1+R/Rscale))/rhom - 1.0;
}

int calc_xi_nfw(double*R, int NR, double Mass, double conc, int delta, double om, double*xi_nfw){
  int i;
  for(i = 0; i < NR; i++)
    xi_nfw[i] = xi_nfw_at_R(R[i], Mass, conc, delta, om);
  return 0;
}

int calc_xi_2halo(int NR, double bias, double*xi_mm, double*xi_2halo){
  int i;
  for(i = 0; i < NR; i++){
    xi_2halo[i] = bias * xi_mm[i];
  }
}

int calc_xi_hm(int NR, double*xi_1h, double*xi_2h, double*xi_hm){
  int i;
  for(i = 0; i < NR; i++){
    if(xi_1h[i] >= xi_2h[i]) xi_hm[i] = xi_1h[i];
    else xi_hm[i] = xi_2h[i];
  }
  return 0;
}

double get_P(double,double,double*,double*,int,gsl_spline*,gsl_interp_accel*);

int calc_xi_mm(double*R, int NR, double*k, double*P, int Nk, double*xi, int N, double h){
  int i;
  for(i=0;i<NR;i++){
    xi[i] = xi_mm_at_R(R[i],k,P,Nk,N,h);
  }
  return 0;
}

///////Functions for calc_xi_mm
//The power spectrum
double get_P(double x,double R,double*k,double*P,int Nk,gsl_spline*Pspl,gsl_interp_accel*acc){
  double ki = x/R;
  double kmin = k[0];
  double kmax = k[Nk-1];
  double alpha,A;
  if (ki < kmin){
    alpha = log(P[1]/P[0])/log(k[1]/k[0]);
    A = P[0]/pow(k[0],alpha);
    return A*pow(ki,alpha);
  }else if (ki > kmax){
    alpha = log(P[Nk-1]/P[Nk-2])/log(k[Nk-1]/k[Nk-2]);
    A = P[Nk-1]/pow(k[Nk-1],alpha);
    return A*pow(ki,alpha);
  }// Assume power laws at ends
  return gsl_spline_eval(Pspl,ki,acc);
}

double xi_mm_at_R(double R, double*k, double*P, int Nk, int N, double h){
  double zero,psi,x,t,dpsi,f,PIsinht;
  double PI_h = M_PI/h;
  double PI_2 = M_PI*0.5;
  gsl_spline*Pspl = gsl_spline_alloc(gsl_interp_cspline,Nk);
  gsl_spline_init(Pspl,k,P,Nk);
  gsl_interp_accel*acc= gsl_interp_accel_alloc();

  double sum = 0;
  int i;
  for(i=0;i<N;i++){
    zero = i+1;
    psi = h*zero*tanh(sinh(h*zero)*PI_2);
    x = psi*PI_h;
    t = h*zero;
    PIsinht = M_PI*sinh(t);
    dpsi = (M_PI*t*cosh(t)+sinh(PIsinht))/(1+cosh(PIsinht));
    if (dpsi!=dpsi) dpsi=1.0;
    f = x*get_P(x,R,k,P,Nk,Pspl,acc);
    sum += f*sin(x)*dpsi;
  }

  gsl_spline_free(Pspl),gsl_interp_accel_free(acc);
  return sum/(R*R*R*M_PI*2);
}
