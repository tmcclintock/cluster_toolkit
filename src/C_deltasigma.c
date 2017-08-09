#include "C_deltasigma.h"
#include "C_xi.h"

#include "gsl/gsl_integration.h"
#include "gsl/gsl_spline.h"
#include <math.h>

#define TOL 1e-5
#define workspace_size 8000
#define rhomconst 2.775808e+11
//1e4*3.*Mpcperkm*Mpcperkm/(8.*PI*G); units are SM h^2/Mpc^3

typedef struct integrand_params{
  gsl_spline*spline;
  gsl_interp_accel*acc;
  double Rp;
  double M;
  double conc;
  int delta;
  double om;
}integrand_params;

double integrand_small_scales(double lRz, void*params){
  double Rz = exp(lRz);
  integrand_params pars=*(integrand_params *)params;
  double Rp = pars.Rp;
  double om = pars.om;
  double M = pars.M;
  double conc = pars.conc;
  int delta = pars.delta;
  return Rz * xi_nfw_at_R(sqrt(Rz*Rz + Rp*Rp), M, conc, delta, om);
}

double integrand_medium_scales(double lRz, void*params){
  double Rz = exp(lRz);
  integrand_params pars=*(integrand_params*)params;
  gsl_spline*spline = pars.spline;
  gsl_interp_accel*acc = pars.acc;
  double Rp = pars.Rp;
  return Rz*gsl_spline_eval(spline, sqrt(Rz*Rz + Rp*Rp), acc);
}

double Sigma_at_R(double R, double*Rxi, double*xi, int Nxi, double M, double conc, int delta, double om){
  double rhom = om*rhomconst*1e-12; //SM h^2/pc^2/Mpc; integral is over Mpc/h
  double loR = Rxi[0];
  double hiR = Rxi[Nxi-1];
  double lnmax = log(sqrt(hiR*hiR - R*R));

  gsl_spline*spline = gsl_spline_alloc(gsl_interp_cspline, Nxi);
  gsl_spline_init(spline, Rxi, xi, Nxi);
  gsl_interp_accel*acc= gsl_interp_accel_alloc();
  gsl_integration_workspace*workspace = gsl_integration_workspace_alloc(workspace_size);

  integrand_params*params=malloc(sizeof(integrand_params));
  params->acc = acc;
  params->spline = spline;
  params->Rp = R;
  params->M = M;
  params->conc= conc;
  params->delta = delta;
  params->om = om;
  gsl_function F;
  F.params=params;
  double result1, err1;
  double result2, err2;
  F.function = &integrand_small_scales;
  gsl_integration_qag(&F, log(loR)-10, log(loR), TOL, TOL/10., workspace_size, 6, workspace, &result1, &err1);
  F.function = &integrand_medium_scales;
  gsl_integration_qag(&F, log(loR), lnmax, TOL, TOL/10., workspace_size, 6, workspace, &result2, &err2);
  gsl_spline_free(spline),gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  free(params);
  return (result1+result2)*rhom*2;
}

int Sigma_at_R_arr(double*R, int NR, double*Rxi, double*xi, int Nxi, double M, double conc, int delta, double om, double*Sigma){
  double rhom = om*rhomconst*1e-12; //SM h^2/pc^2/Mpc; integral is over Mpc/h
  double loR = Rxi[0];
  double hiR = Rxi[Nxi-1];
  double lnmax;

  gsl_spline*spline = gsl_spline_alloc(gsl_interp_cspline, Nxi);
  gsl_spline_init(spline, Rxi, xi, Nxi);
  gsl_interp_accel*acc= gsl_interp_accel_alloc();
  gsl_integration_workspace*workspace = gsl_integration_workspace_alloc(workspace_size);

  integrand_params*params=malloc(sizeof(integrand_params));
  params->acc = acc;
  params->spline = spline;
  params->M = M;
  params->conc= conc;
  params->delta = delta;
  params->om = om;
  gsl_function F;
  F.params=params;
  double result1, err1;
  double result2, err2;
  int i;
  for(i = 0; i < NR; i++){
    lnmax = log(sqrt(hiR*hiR - R[i]*R[i]));
    params->Rp = R[i];
    F.function = &integrand_small_scales;
    gsl_integration_qag(&F, log(loR)-10, log(loR), TOL, TOL/10., workspace_size, 6, workspace, &result1, &err1);
    F.function = &integrand_medium_scales;
    gsl_integration_qag(&F, log(loR), lnmax, TOL, TOL/10., workspace_size, 6, workspace, &result2, &err2);
    Sigma[i] = (result1+result2)*rhom*2;
  }
  gsl_spline_free(spline),gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  free(params);

  return 0;
}
