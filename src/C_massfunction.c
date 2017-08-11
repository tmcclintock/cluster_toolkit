#include "C_massfunction.h"
#include "C_bias.h"

#include "gsl/gsl_integration.h"
#include "gsl/gsl_sf.h"
#include "gsl/gsl_spline.h"
#include <math.h>

#define TOL 1e-6 //Used for the tinker bias
#define rhomconst 2.775808e+11
//1e4*3.*Mpcperkm*Mpcperkm/(8.*PI*G); units are SM h^2/Mpc^3
#define workspace_size 8000


///////////////// dndM functions below ///////////////////

double dndM_at_M(double M, double*k, double*P, int Nk, double om, double d, double e, double f, double g){
  //Compute the prefactor B
  double d2 = 0.5*d;
  double gamma_d2 = gsl_sf_gamma(d2);
  double f2 = 0.5*f;
  double gamma_f2 = gsl_sf_gamma(f2);
  double B = 2./(pow(e, d)*pow(g, -d2)*gamma_d2 + pow(g, -f2)*gamma_f2);
  //Compute sigma(M)
  double sigma = sqrt(sigma2_at_M(M, k, P, Nk, om));
  //Compute g(sigma)
  double gsigma = B*exp(-g/(sigma*sigma))*(pow(sigma/e, -d)+pow(sigma, -f));
  double dM = 1e-6*M;
  double dlnsiginvdm = log(sqrt(sigma2_at_M(M-dM*0.5, k, P, Nk, om)/sigma2_at_M(M+dM*0.5, k, P, Nk, om)))/dM;
  return gsigma*om*rhomconst/M*dlnsiginvdm;
}

int dndM_at_M_arr(double*M, int NM, double*k, double*P, int Nk, double om, double d, double e, double f, double g, double*dndM){
  double rhom = om*rhomconst;
  //Compute the prefactor B
  double d2 = 0.5*d;
  double gamma_d2 = gsl_sf_gamma(d2);
  double f2 = 0.5*f;
  double gamma_f2 = gsl_sf_gamma(f2);
  double B = 2./(pow(e, d)*pow(g, -d2)*gamma_d2 + pow(g, -f2)*gamma_f2);
  //Compute all sigma(M) arrays
  double*M_top = (double*)malloc(sizeof(double)*NM);
  double*M_bot = (double*)malloc(sizeof(double)*NM);
  double*sigma2 = (double*)malloc(sizeof(double)*NM);
  double*sigma2_top = (double*)malloc(sizeof(double)*NM);
  double*sigma2_bot = (double*)malloc(sizeof(double)*NM);
  int i;
  double dM = 0;
  for(i = 0; i < NM; i++){
    dM = 1e-6*M[i];
    M_top[i] = M[i]-dM*0.5;
    M_bot[i] = M[i]+dM*0.5;
  }
  sigma2_at_M_arr(M, NM, k, P, Nk, om, sigma2);
  sigma2_at_M_arr(M_top, NM, k, P, Nk, om, sigma2_top);
  sigma2_at_M_arr(M_bot, NM, k, P, Nk, om, sigma2_bot);
  //Compute g(sigma)
  double*gsigma = (double*)malloc(sizeof(double)*NM);
  double*dlnsiginvdm = (double*)malloc(sizeof(double)*NM);
  double sigma;
  for(i = 0; i < NM; i++){
    dM = 1e-6*M[i];
    sigma = sqrt(sigma2[i]);
    gsigma[i] = B*exp(-g/sigma2[i])*(pow(sigma/e, -d)+pow(sigma, -f));
    dlnsiginvdm[i] = log(sqrt(sigma2_top[i]/sigma2_bot[i]))/dM;
    dndM[i] = gsigma[i]*rhom/M[i]*dlnsiginvdm[i];
  }
  free(M_top); free(M_bot);
  free(sigma2); free(sigma2_top); free(sigma2_bot);
  free(gsigma);
  free(dlnsiginvdm);
  return 0;
}

///////////////// N in bin functions below ///////////////////

typedef struct integrand_params{
  gsl_spline *spline;
  gsl_interp_accel *acc;
}integrand_params;

double dndm_integrand(double lM, void*params){
  double M = exp(lM);
  integrand_params pars = *(integrand_params*)params;
  gsl_spline*spline = pars.spline;
  gsl_interp_accel*acc = pars.acc;
  return M*gsl_spline_eval(spline, M, acc);
}

double N_in_bin(double*M, double*dndM, int NM, double volume, double Mlo, double Mhi){
  gsl_spline*spline = gsl_spline_alloc(gsl_interp_cspline, NM);
  gsl_spline_init(spline, M, dndM, NM);
  gsl_interp_accel*acc = gsl_interp_accel_alloc();
  gsl_integration_workspace*workspace = gsl_integration_workspace_alloc(workspace_size);
  integrand_params *params = malloc(sizeof(integrand_params));
  params->spline = spline;
  params->acc = acc;
  gsl_function F;
  F.function = &dndm_integrand;
  F.params = params;
  double result, err;
  gsl_integration_qag(&F, log(Mlo), log(Mhi), TOL ,TOL/10., workspace_size, 6, workspace, &result, &err);
  gsl_spline_free(spline),gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  free(params);
  //result contains number density
  return result*volume;
}

int N_in_bins(double*M, double*dndM, int NM, double volume, double*bins, int Nedges, double*N){
  //Note: N is one element less than bins
  gsl_spline*spline = gsl_spline_alloc(gsl_interp_cspline, NM);
  gsl_spline_init(spline, M, dndM, NM);
  gsl_interp_accel*acc = gsl_interp_accel_alloc();
  gsl_integration_workspace*workspace = gsl_integration_workspace_alloc(workspace_size);
  integrand_params *params = malloc(sizeof(integrand_params));
  params->spline = spline;
  params->acc = acc;
  gsl_function F;
  F.function = &dndm_integrand;
  F.params = params;
  double result, err;
  int i;
  for(i = 0; i < Nedges-1; i++){
    gsl_integration_qag(&F, log(bins[i]), log(bins[i+1]), TOL ,TOL/10., workspace_size, 6, workspace, &result, &err);
    N[i] = result*volume;
  }
  gsl_spline_free(spline),gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  free(params);
  return 0;
}
