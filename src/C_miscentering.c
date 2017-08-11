#include "C_miscentering.h"
#include "C_deltasigma.h"

#include "gsl/gsl_integration.h"
#include "gsl/gsl_spline.h"
#include <math.h>
#include <stdio.h>

#define TOL 1e-2 // Used for miscentering
#define TOL2 TOL*0.1
#define workspace_size 8000
#define rhomconst 2.775808e+11
//1e4*3.*Mpcperkm*Mpcperkm/(8.*PI*G); units are SM h^2/Mpc^3

////////////// SIGMA(R) FUNCTIONS BELOW////////////////

typedef struct integrand_params{
  gsl_spline*spline;
  gsl_interp_accel*acc;
  gsl_integration_workspace*workspace;
  gsl_integration_workspace*workspace2;
  double Rp;  //R_perp
  double Rp2; //R_perp^2
  double rmin;
  double rmax;
  double lrmin;
  double lrmax;
  double M;
  double conc;
  int delta;
  double om;
  double Rmis; //Miscentering length
  double Rmis2; //Rmis^2
  double Rp_cos_theta_2;
  double slope;
  double intercept;
}integrand_params;

double single_angular_integrand(double theta, void*params){
  integrand_params*pars = (integrand_params*)params;
  double Rp = pars->Rp;
  double Rmis = pars->Rmis;
  double arg = sqrt(Rp*Rp + Rmis*Rmis - 2*Rp*Rmis*cos(theta));
  double rmin = pars->rmin,rmax = pars->rmax;
  if (arg < rmin){
    return Sigma_nfw_at_R(arg, pars->M, pars->conc, pars->delta, pars->om);
  }else if(arg < rmax){
    gsl_spline*spline = pars->spline;
    gsl_interp_accel*acc = pars->acc;
    return gsl_spline_eval(spline, arg ,acc);
  }
  return 0;
}

double Sigma_mis_single_at_R(double R, double*Rs, double*Sigma, int Ns, double M, double conc, int delta, double om, double Rmis){
  gsl_spline*spline = gsl_spline_alloc(gsl_interp_cspline, Ns);
  gsl_spline_init(spline, Rs, Sigma, Ns);
  gsl_interp_accel*acc = gsl_interp_accel_alloc();
  gsl_integration_workspace*workspace = gsl_integration_workspace_alloc(workspace_size);
  gsl_integration_workspace*workspace2 = gsl_integration_workspace_alloc(workspace_size);
  integrand_params*params = malloc(sizeof(integrand_params));
  params->acc = acc;
  params->spline = spline;
  params->workspace = workspace;
  params->workspace2 = workspace2;
  params->Rp  = R;
  params->Rp2 = R * R;
  params->M = M;
  params->conc = conc;
  params->delta = delta;
  params->om = om;
  params->Rmis = Rmis;
  params->Rmis2 = Rmis*Rmis;
  params->rmin = Rs[0];
  params->rmax = Rs[Ns-1];
  params->lrmin = log(Rs[0]);
  params->lrmax = log(Rs[Ns-1]);
  gsl_function F;
  F.function=&single_angular_integrand;
  F.params=params;
  double result, err;
  gsl_integration_qag(&F, 0, M_PI, TOL, TOL2, workspace_size, 6, workspace, &result, &err);
  gsl_spline_free(spline),gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  gsl_integration_workspace_free(workspace2);
  free(params);
  return result/=M_PI; //Factor of PI from the angular integral
}

int Sigma_mis_single_at_R_arr(double*R, int NR, double*Rs, double*Sigma, int Ns, double M, double conc, int delta, double om, double Rmis, double*Sigma_mis){
  gsl_spline*spline = gsl_spline_alloc(gsl_interp_cspline, Ns);
  gsl_spline_init(spline, Rs, Sigma, Ns);
  gsl_interp_accel*acc = gsl_interp_accel_alloc();
  gsl_integration_workspace*workspace = gsl_integration_workspace_alloc(workspace_size);
  gsl_integration_workspace*workspace2 = gsl_integration_workspace_alloc(workspace_size);
  integrand_params*params = malloc(sizeof(integrand_params));
  params->acc = acc;
  params->spline = spline;
  params->workspace = workspace;
  params->workspace2 = workspace2;
  params->M = M;
  params->conc = conc;
  params->delta = delta;
  params->om = om;
  params->Rmis = Rmis;
  params->Rmis2 = Rmis*Rmis;
  params->rmin = Rs[0];
  params->rmax = Rs[Ns-1];
  params->lrmin = log(Rs[0]);
  params->lrmax = log(Rs[Ns-1]);
  gsl_function F;
  F.function=&single_angular_integrand;
  F.params=params;
  double result, err;
  int i;
  for(i = 0; i < NR; i++){
    params->Rp  = R[i];
    params->Rp2 = R[i] * R[i];
    gsl_integration_qag(&F, 0, M_PI, TOL, TOL2, workspace_size, 6, workspace, &result, &err);
    Sigma_mis[i] = result/M_PI;
  }
  gsl_spline_free(spline),gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  gsl_integration_workspace_free(workspace2);
  free(params);
  return 0;
}

/////////////////// SIGMA(R) CONVOLUTION BELOW /////////////

double g2d_radial_integrand(double lRc, void*params){
  integrand_params*pars = (integrand_params*)params;
  double Rc = exp(lRc);
  double Rc2 = Rc*Rc;
  double rmin = pars->rmin,rmax = pars->rmax;
  double Rp2 = pars->Rp2;
  double Rmis2 = pars->Rmis2;
  double Rp_cos_theta_2 = pars->Rp_cos_theta_2;
  double arg = sqrt(Rp2 + Rc2 - Rc*Rp_cos_theta_2);
  double answer = 0;
  /*Note: P_miscentering is a 2D gaussian.
    This is the Rc/Rmis^2 * exp() term.*/
  if(arg < rmin){
    answer = Sigma_nfw_at_R(arg, pars->M, pars->conc, pars->delta, pars->om);
  }else if(arg < rmax){
    gsl_spline*spline = pars->spline;
    gsl_interp_accel*acc = pars->acc;
    answer = gsl_spline_eval(spline, arg, acc);
  }
  return Rc2 * exp(-0.5 * Rc2/Rmis2) * answer;
}

double g2d_angular_integrand(double theta, void*params){
  integrand_params*pars = (integrand_params*)params;
  double Rmis2 = pars->Rmis2;
  double cos_theta = cos(theta);
  pars->Rp_cos_theta_2 = pars->Rp*cos_theta*2;
  double lrmin = pars->lrmin, lrmax = pars->lrmax;
  gsl_integration_workspace*workspace = pars->workspace2;
  gsl_function F;
  F.function = &g2d_radial_integrand;
  F.params = pars;
  double result, err;
  gsl_integration_qag(&F, lrmin-10, lrmax, TOL, TOL2, workspace_size, 6, workspace, &result, &err);
  return result/Rmis2;
}

//2D Gaussian kernel
double Sigma_mis_g2d_at_R(double R, double*Rs, double*Sigma, int Ns, double M, double conc, int delta, double om, double Rmis){
  gsl_spline*spline = gsl_spline_alloc(gsl_interp_cspline, Ns);
  gsl_spline_init(spline, Rs, Sigma, Ns);
  gsl_interp_accel*acc = gsl_interp_accel_alloc();
  gsl_integration_workspace*workspace = gsl_integration_workspace_alloc(workspace_size);
  gsl_integration_workspace*workspace2 = gsl_integration_workspace_alloc(workspace_size);
  integrand_params*params = malloc(sizeof(integrand_params));
  params->acc = acc;
  params->spline = spline;
  params->workspace = workspace;
  params->workspace2 = workspace2;
  params->Rp  = R;
  params->Rp2 = R * R;
  params->M = M;
  params->conc = conc;
  params->delta = delta;
  params->om = om;
  params->Rmis = Rmis;
  params->Rmis2 = Rmis*Rmis;
  params->rmin = Rs[0];
  params->rmax = Rs[Ns-1];
  params->lrmin = log(Rs[0]);
  params->lrmax = log(Rs[Ns-1]);
  gsl_function F;
  F.function=&g2d_angular_integrand;
  F.params=params;
  double result, err;
  gsl_integration_qag(&F, 0, M_PI, TOL, TOL2, workspace_size, 6, workspace, &result, &err);
  gsl_spline_free(spline),gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  gsl_integration_workspace_free(workspace2);
  free(params);
  return result/=M_PI; //Factor of PI from the angular integral
}

//2D Gaussian kernel
int Sigma_mis_g2d_at_R_arr(double*R, int NR, double*Rs, double*Sigma, int Ns, double M, double conc, int delta, double om, double Rmis, double*Sigma_mis){
  gsl_spline*spline = gsl_spline_alloc(gsl_interp_cspline, Ns);
  gsl_spline_init(spline, Rs, Sigma, Ns);
  gsl_interp_accel*acc = gsl_interp_accel_alloc();
  gsl_integration_workspace*workspace = gsl_integration_workspace_alloc(workspace_size);
  gsl_integration_workspace*workspace2 = gsl_integration_workspace_alloc(workspace_size);
  integrand_params*params = malloc(sizeof(integrand_params));
  params->acc = acc;
  params->spline = spline;
  params->workspace = workspace;
  params->workspace2 = workspace2;
  params->M = M;
  params->conc = conc;
  params->delta = delta;
  params->om = om;
  params->Rmis = Rmis;
  params->Rmis2 = Rmis*Rmis;
  params->rmin = Rs[0];
  params->rmax = Rs[Ns-1];
  params->lrmin = log(Rs[0]);
  params->lrmax = log(Rs[Ns-1]);
  gsl_function F;
  F.function=&g2d_angular_integrand;
  F.params=params;
  double result, err;
  int i;
  for(i = 0; i < NR; i++){
    params->Rp  = R[i];
    params->Rp2 = R[i] * R[i];
    gsl_integration_qag(&F, 0, M_PI, TOL, TOL2, workspace_size, 6, workspace, &result, &err);
    Sigma_mis[i] = result/M_PI;
  }
  gsl_spline_free(spline),gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  gsl_integration_workspace_free(workspace2);
  free(params);
  return 0;
}

//////////// DELTASIGMA(R) BELOW //////////////////

double DS_mis_integrand(double lR, void*params){
  double R = exp(lR);
  integrand_params pars = *(integrand_params*)params;
  gsl_spline*spline = pars.spline;//Sigma(R) spline
  gsl_interp_accel*acc = pars.acc;
  return R * R * gsl_spline_eval(spline, R, acc);
}

double DeltaSigma_mis_at_R(double R, double*Rs, double*Sigma, int Ns){
  double lrmin = log(Rs[0]);

  gsl_spline*spline = gsl_spline_alloc(gsl_interp_cspline, Ns);
  gsl_spline_init(spline, Rs, Sigma, Ns);
  gsl_interp_accel*acc = gsl_interp_accel_alloc();
  gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(workspace_size);
  integrand_params*params=malloc(sizeof(integrand_params));
  params->spline = spline;
  params->acc = acc;
  double slope = log(Sigma[0]/Sigma[1])/log(Rs[0]/Rs[1]);
  double intercept = Sigma[0]*pow(Rs[0], -slope);
  double low_part = intercept*pow(Rs[0], slope+2)/(slope+2);
  double result, err;
  gsl_function F;
  F.params = params;
  F.function = &DS_mis_integrand;
  gsl_integration_qag(&F, lrmin, log(R), TOL, TOL/10., workspace_size, 6, workspace, &result, &err);
  
  gsl_spline_free(spline),gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  free(params);
  return (low_part+result)*2/(R*R) - gsl_spline_eval(spline, R, acc);
}

int DeltaSigma_mis_at_R_arr(double*R, int NR, double*Rs, double*Sigma, int Ns, double*DeltaSigma_mis){
  double lrmin = log(Rs[0]);

  gsl_spline*spline = gsl_spline_alloc(gsl_interp_cspline, Ns);
  gsl_spline_init(spline, Rs, Sigma, Ns);
  gsl_interp_accel*acc = gsl_interp_accel_alloc();
  gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(workspace_size);
  integrand_params*params=malloc(sizeof(integrand_params));
  params->spline = spline;
  params->acc = acc;
  double slope = log(Sigma[0]/Sigma[1])/log(Rs[0]/Rs[1]);
  double intercept = Sigma[0]*pow(Rs[0], -slope);
  double low_part = intercept*pow(Rs[0], slope+2)/(slope+2);
  double result,  err;
  gsl_function F;
  F.params = params;
  F.function = &DS_mis_integrand;
  int i;
  for(i = 0; i < NR; i++){
    gsl_integration_qag(&F, lrmin, log(R[i]), TOL, TOL/10., workspace_size, 6, workspace, &result, &err);
    DeltaSigma_mis[i] = (low_part+result)*2/(R[i]*R[i]) - gsl_spline_eval(spline, R[i], acc);
  }
  gsl_spline_free(spline),gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  free(params);
  return 0;
}
