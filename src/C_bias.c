#include "C_bias.h"

#include "gsl/gsl_integration.h"
#include "gsl/gsl_spline.h"

#include <math.h>
#include <stdio.h>

#define BIAS_TOL 1e-6 //Used for the tinker bias
#define delta_c 1.686 //Critical collapse density
#define rhomconst 2.775808e+11
//1e4*3.*Mpcperkm*Mpcperkm/(8.*PI*G); units are Msun h^2/Mpc^3
#define workspace_size 8000

typedef struct integrand_params{
  gsl_spline *spline;
  gsl_interp_accel *acc;
  double r;
}integrand_params;

double M_to_R(double M, double Omega_m){
  //Lagrangian radius Mpc/h
  return pow(M/(1.33333333333*M_PI*rhomconst*Omega_m),0.3333333333);
}

double R_to_M(double R, double om){
  //Lagrangian mass Msun/h
  return R*R*R*1.33333333333*M_PI*rhomconst*om;
}

double integrand(double lk, void*params){
  integrand_params pars = *(integrand_params*)params;
  gsl_spline*spline = pars.spline;
  gsl_interp_accel*acc = pars.acc;
  double R = pars.r;
  double k = exp(lk);
  double x  = k*R;
  double P = gsl_spline_eval(spline, k, acc);
  double w = 3.0/x/x/x*(sin(x)-x*cos(x)); //Window function
  return k*k*k*P*w*w;
}

///////////// linar matter variance functions /////////////

double sigma2_at_R(double R, double*k, double*P, int Nk){
  gsl_spline*spline = gsl_spline_alloc(gsl_interp_cspline, Nk);
  gsl_spline_init(spline, k, P, Nk);
  gsl_interp_accel*acc = gsl_interp_accel_alloc();
  gsl_integration_workspace*workspace = gsl_integration_workspace_alloc(workspace_size);
  integrand_params *params = malloc(sizeof(integrand_params));
  params->spline = spline;
  params->acc = acc;
  params->r = R;
  gsl_function F;
  F.function = &integrand;
  F.params = params;
  double result, err;
  gsl_integration_qag(&F, log(k[0]), log(k[Nk-1]), BIAS_TOL ,BIAS_TOL/10., workspace_size, 6, workspace, &result, &err);
  gsl_spline_free(spline),gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  free(params);
  return result/(2*M_PI*M_PI);
}

double sigma2_at_M(double M, double*k, double*P, int Nk, double om){
  double R=M_to_R(M, om);
  return sigma2_at_R(R, k, P, Nk);
}

int sigma2_at_R_arr(double*R, int NR,  double*k, double*P, int Nk, double*s2){
  gsl_spline*spline = gsl_spline_alloc(gsl_interp_cspline,Nk);
  gsl_spline_init(spline,k,P,Nk);
  gsl_interp_accel*acc = gsl_interp_accel_alloc();
  gsl_integration_workspace*workspace = gsl_integration_workspace_alloc(workspace_size);
  integrand_params *params = malloc(sizeof(integrand_params));
  params->spline = spline;
  params->acc = acc;
  gsl_function F;
  F.function = &integrand;
  F.params = params;
  double lo = log(k[0]);
  double hi = log(k[Nk-1]);
  double result,abserr;
  int i;
  for(i = 0; i < NR; i++){
    params->r = R[i];
    gsl_integration_qag(&F, lo, hi, BIAS_TOL ,BIAS_TOL/10., workspace_size, 6, workspace, &result, &abserr);
    s2[i] = result/(2*M_PI*M_PI);
  }
  gsl_spline_free(spline),gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  free(params);
  return 0;
}

int sigma2_at_M_arr(double*M, int NM,  double*k, double*P, int Nk, double om, double*s2){
  int i;
  double*R = (double*)malloc(sizeof(double)*NM);
  for(i = 0; i < NM; i++)
    R[i] = M_to_R(M[i], om);
  sigma2_at_R_arr(R, NM, k, P, Nk, s2);
  return 0;
}

///////////PEAK HEIGHT FUNCTIONS///////////

double nu_at_R(double R, double*k, double*P, int Nk){
  return delta_c/sqrt(sigma2_at_R(R, k, P, Nk));
}

double nu_at_M(double M, double*k, double*P, int Nk, double om){
  double R=M_to_R(M, om);
  return nu_at_R(R, k, P, Nk);
}

int nu_at_R_arr(double*R, int NR, double*k, double*P, int Nk, double*nu){
  int i;
  double*s2 = (double*)malloc(sizeof(double)*NR);
  sigma2_at_R_arr(R, NR, k, P, Nk, s2);
  for(i = 0; i < NR; i++)
    nu[i] = delta_c/sqrt(s2[i]);
  return 0;
}

int nu_at_M_arr(double*M, int NM, double*k, double*P, int Nk, double om, double*nu){
  int i;
  double*R = (double*)malloc(sizeof(double)*NM);
  for(i = 0; i < NM; i++)
    R[i] = M_to_R(M[i], om);
  nu_at_R_arr(R, NM, k, P, Nk, nu);
  return 0;
}

///////////BIAS FUNCTIONS///////////

double bias_at_nu(double nu, int delta){
  double y = log10(delta);
  double xp = exp(-1.0*pow(4./y,4.));
  double A = 1.+0.24*y*xp, a = 0.44*y-0.88;
  double B = 0.183, b = 1.5;
  double C = 0.019+0.107*y+0.19*xp, c = 2.4;
  return 1 - A*pow(nu,a)/(pow(nu,a)+pow(delta_c,a)) + B*pow(nu,b) + C*pow(nu,c);
}

double bias_at_R(double R, int delta, double*k, double*P, int Nk){
  return bias_at_nu(nu_at_R(R, k, P, Nk), delta);
}

double bias_at_M(double M, int delta, double*k, double*P, int Nk, double om){
  return bias_at_nu(nu_at_M(M, k, P, Nk, om), delta);
}

int bias_at_nu_arr(double*nu, int Nnu, int delta, double*bias){
  double y = log10(delta);
  double xp = exp(-1.0*pow(4./y,4.));
  double A = 1.+0.24*y*xp, a = 0.44*y-0.88;
  double B = 0.183, b = 1.5;
  double C = 0.019+0.107*y+0.19*xp, c = 2.4;
  int i;
  for(i = 0; i < Nnu; i++)
    bias[i] = 1 - A*pow(nu[i],a)/(pow(nu[i],a)+pow(delta_c,a)) + B*pow(nu[i],b) + C*pow(nu[i],c);
  return 0;
}

int bias_at_R_arr(double*R, int NR, int delta, double*k, double*P, int Nk, double*bias){
  double*nu = (double*)malloc(sizeof(double)*NR);
  nu_at_R_arr(R, NR, k, P, Nk, nu);
  bias_at_nu_arr(nu, NR, delta, bias);
  free(nu);
  return 0;
}

int bias_at_M_arr(double*M, int NM, int delta, double*k, double*P, int Nk, double om, double*bias){
  double*nu = (double*)malloc(sizeof(double)*NM);
  nu_at_M_arr(M, NM, k, P, Nk, om, nu);
  bias_at_nu_arr(nu, NM, delta, bias);
  free(nu);
  return 0;
}
