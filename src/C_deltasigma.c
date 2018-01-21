#include "C_deltasigma.h"
#include "C_xi.h"

#include "gsl/gsl_integration.h"
#include "gsl/gsl_spline.h"
#include <math.h>
#include <stdio.h>

#define ABSERR 0.0
#define RELERR 1e-4
#define workspace_size 8000
#define ulim 5.0
#define rhocrit 2.77533742639e+11
//1e4*3.*Mpcperkm*Mpcperkm/(8.*PI*G); units are SM h^2/Mpc^3

////////////// SIGMA(R) FUNCTIONS BELOW////////////////

double Sigma_nfw_at_R(double R, double M, double c, int delta, double om){
  double*Rs = (double*)malloc(sizeof(double));
  double*Sigma = (double*)malloc(sizeof(double));
  double result;
  Rs[0] = R;
  Sigma_nfw_at_R_arr(Rs, 1, M, c, delta, om, Sigma);
  result = Sigma[0];
  free(Rs);
  free(Sigma);
  return result;
}

int Sigma_nfw_at_R_arr(double*R, int NR, double M, double c, int delta, double om, double*Sigma){
  double rhom = om*rhocrit;//SM h^2/Mpc^3
  double deltac = delta*0.3333333333*c*c*c/(log(1.+c)-c/(1.+c));
  double Rdelta = pow(M/(1.333333333*M_PI*rhom*delta),0.333333333);//Mpc/h
  double Rscale = Rdelta/c;
  double gx = 0;
  double x;
  int i;
  for(i = 0; i < NR; i++){
    x = R[i]/Rscale;
    if(x<1){
      gx = (1 - 2./sqrt(1-x*x)*atanh(sqrt((1-x)/(1+x))))/(x*x-1);
    }else{// if(x>=1){
      gx = (1 - 2./sqrt(x*x-1)* atan(sqrt((x-1)/(1+x))))/(x*x-1);
    }
    Sigma[i] = 2*Rscale*deltac*rhom*gx*1.e-12; //SM h/pc^2
  }
  return 0;
}

typedef struct integrand_params{
  gsl_spline*spline;
  gsl_interp_accel*acc;
  double Rp;
  double M;
  double conc;
  int delta;
  double om;
  double slope;
  double intercept;
}integrand_params;

double integrand_small_scales(double lRz, void*params){
  double Rz = exp(lRz);
  integrand_params pars=*(integrand_params*)params;
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
  return Rz * gsl_spline_eval(spline, sqrt(Rz*Rz + Rp*Rp), acc);
}

double integrand_large_scales(double lRz, void*params){
  //Use a power law approximation. This is fine as long as xi_hm
  //doesn't have wiggles in it (i.e. it doesn't stop at BAO)
  double Rz = exp(lRz);
  integrand_params pars=*(integrand_params*)params;
  double Rp = pars.Rp;
  double slope = pars.slope;
  double intercept = pars.intercept;
  return Rz * intercept*pow(sqrt(Rz*Rz + Rp*Rp), slope);
}

double Sigma_at_R(double R, double*Rxi, double*xi, int Nxi, double M, double conc, int delta, double om){
  double*Rs = (double*)malloc(sizeof(double));
  double*Sigma = (double*)malloc(sizeof(double));
  double result;
  Rs[0] = R;
  Sigma_at_R_arr(Rs, 1, Rxi, xi, Nxi, M, conc, delta, om, Sigma);
  result = Sigma[0];
  free(Rs);
  free(Sigma);
  return result;
}

int Sigma_at_R_arr(double*R, int NR, double*Rxi, double*xi, int Nxi, double M, double conc, int delta, double om, double*Sigma){
  double rhom = om*rhocrit*1e-12; //SM h^2/pc^2/Mpc; integral is over Mpc/h
  double Rxi0 = Rxi[0];
  double Rxi_max = Rxi[Nxi-1];
  double ln_z_max;

  gsl_spline*spline = gsl_spline_alloc(gsl_interp_cspline, Nxi);
  gsl_spline_init(spline, Rxi, xi, Nxi);
  gsl_interp_accel*acc= gsl_interp_accel_alloc();
  gsl_integration_workspace*workspace = gsl_integration_workspace_alloc(workspace_size);

  integrand_params*params = malloc(sizeof(integrand_params));
  params->acc = acc;
  params->spline = spline;
  params->M = M;
  params->conc= conc;
  params->delta = delta;
  params->om = om;
  gsl_function F;
  F.params = params;
  double result1, err1, result2, err2;
  int i;
  
  for(i = 0; i < NR; i++){
    ln_z_max = log(sqrt(Rxi_max*Rxi_max - R[i]*R[i])); //Max distance to integrate to
    params->Rp = R[i];
    if(R[i] < Rxi0){
      F.function = &integrand_small_scales;
      gsl_integration_qag(&F, log(Rxi0)-10, log(sqrt(Rxi0*Rxi0-R[i]*R[i])), ABSERR, RELERR, workspace_size, 6, workspace, &result1, &err1);
      F.function = &integrand_medium_scales;
      gsl_integration_qag(&F, log(sqrt(Rxi0*Rxi0-R[i]*R[i])), ln_z_max, ABSERR, RELERR, workspace_size, 6, workspace, &result2, &err2);
    }else{ //R[i] > Rxi0
      result1 = 0;
      F.function = &integrand_medium_scales;
      gsl_integration_qag(&F, -10, ln_z_max, ABSERR, RELERR, workspace_size, 6, workspace, &result2, &err2);
    }
    Sigma[i] = (result1+result2)*rhom*2;
  }
  
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  free(params);

  return 0;
}

double Sigma_at_R_full(double R, double*Rxi, double*xi, int Nxi, double M, double conc, int delta, double om){
  double*Ra = (double*)malloc(sizeof(double));
  double*Sigma = (double*)malloc(sizeof(double));
  double result;
  Ra[0] = R;
  Sigma_at_R_full_arr(Ra, 1, Rxi, xi, Nxi, M, conc, delta, om, Sigma);
  result = Sigma[0];
  free(Ra);
  free(Sigma);
  return result;
}

int Sigma_at_R_full_arr(double*R, int NR, double*Rxi, double*xi, int Nxi, double M, double conc, int delta, double om, double*Sigma){
  Sigma_at_R_arr(R, NR, Rxi, xi, Nxi, M, conc, delta, om, Sigma);
  double rhom = om*rhocrit*1e-12; //Msun h^2/pc^2/Mpc; integral is over Mpc/h
  double Rxi_max = Rxi[Nxi-1];
  double ln_z_max;
  double result3, err3;
  gsl_integration_workspace*workspace = gsl_integration_workspace_alloc(workspace_size);
  integrand_params*params=malloc(sizeof(integrand_params));
  params->slope = log(xi[Nxi-1]/xi[Nxi-2])/log(Rxi[Nxi-1]/Rxi[Nxi-2]);
  params->intercept = xi[Nxi-1]/pow(Rxi[Nxi-1], params->slope);
  gsl_function F;
  F.params=params;
  F.function = &integrand_large_scales;
  int i;
  for(i = 0; i < NR; i++){
    ln_z_max = log(sqrt(Rxi_max*Rxi_max - R[i]*R[i]));
    params->Rp = R[i];
    gsl_integration_qag(&F, ln_z_max, ln_z_max+ulim, ABSERR, RELERR, workspace_size, 6, workspace, &result3, &err3);
    Sigma[i] += (result3*rhom*2);
  }
  gsl_integration_workspace_free(workspace);
  free(params);
  return 0;
}

////////////// DELTASIGMA FUNCTIONS BELOW////////////////

double DS_integrand_small_scales(double lR, void*params){
  double R = exp(lR);
  integrand_params pars = *(integrand_params*)params;
  double M = pars.M;
  double conc = pars.conc;
  int delta = pars.delta;
  double om = pars.om;
  return R * R * Sigma_nfw_at_R(R, M, conc, delta, om);
}

double DS_integrand_medium_scales(double lR, void*params){
  double R = exp(lR);
  integrand_params pars = *(integrand_params*)params;
  gsl_spline*spline = pars.spline;//Sigma(R) spline
  gsl_interp_accel*acc = pars.acc;
  return R * R * gsl_spline_eval(spline, R, acc);
}

double DeltaSigma_at_R(double R, double*Rs, double*Sigma, int Ns, double M, double conc, int delta, double om){
  double*Ra = (double*)malloc(sizeof(double));
  double*DS = (double*)malloc(sizeof(double));
  double result;
  Ra[0] = R;
  DeltaSigma_at_R_arr(Ra, 1, Rs, Sigma, Ns, M, conc, delta, om, DS);
  result = DS[0];
  free(Ra);
  free(DS);
  return result;
}

int DeltaSigma_at_R_arr(double*R, int NR, double*Rs, double*Sigma, int Ns, double M, double conc, int delta, double om, double*DeltaSigma){
  double lrmin = log(Rs[0]);

  gsl_spline*spline = gsl_spline_alloc(gsl_interp_cspline, Ns);
  gsl_spline_init(spline, Rs, Sigma, Ns);
  gsl_interp_accel*acc = gsl_interp_accel_alloc();
  gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(workspace_size);
  integrand_params*params=malloc(sizeof(integrand_params));
  params->spline = spline;
  params->acc = acc;
  params->M = M;
  params->conc = conc;
  params->delta = delta;
  params->om = om;
  double result1, result2, err1, err2;
  gsl_function F;
  F.params = params;
  F.function = &DS_integrand_small_scales;
  gsl_integration_qag(&F, lrmin-10, lrmin, ABSERR, RELERR, workspace_size, 6, workspace, &result1, &err1);
  F.function = &DS_integrand_medium_scales;
  int i;
  for(i = 0; i < NR; i++){
    gsl_integration_qag(&F, lrmin, log(R[i]), ABSERR, RELERR, workspace_size, 6, workspace, &result2, &err2);
    DeltaSigma[i] = (result1+result2)*2/(R[i]*R[i]) - gsl_spline_eval(spline, R[i], acc);
  }
  gsl_spline_free(spline),gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  free(params);
  return 0;
}
