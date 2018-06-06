/** @file C_miscentering
 *  @brief Miscentering effects on projected cluster profiles.
 * 
 *  This file implements the functions that modify projected
 *  galaxy cluster weak lensing profiles. For stacked clusters,
 *  this includes different distributions of miscentering lengths,
 *  or the distribution of the amount that clusters are miscentered.
 *  
 *  @author Tom McClintock (tmcclintock)
 *  @bug No known bugs.
 */

#include "C_miscentering.h"
#include "C_deltasigma.h"

#include "gsl/gsl_integration.h"
#include "gsl/gsl_spline.h"
#include <math.h>
#include <stdio.h>

#define ABSERR 0.0
#define RELERR 1e-4 // Used for miscentering
#define workspace_size 8000
#define rhomconst 2.77533742639e+11
//1e4*3.*Mpcperkm*Mpcperkm/(8.*PI*G); units are SM h^2/Mpc^3
#define KEY 1 //Used for GSL QAG function

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
  gsl_function F_radial; // function for the radial part of the miscentering 
}integrand_params;

double single_angular_integrand(double theta, void*params){
  integrand_params*pars = (integrand_params*)params;
  gsl_spline*spline = pars->spline;
  gsl_interp_accel*acc = pars->acc;
  double Rp = pars->Rp;
  double Rmis = pars->Rmis;
  double arg = sqrt(Rp*Rp + Rmis*Rmis - 2*Rp*Rmis*cos(theta));
  double rmin = pars->rmin,rmax = pars->rmax;
  if (arg < rmin){
    return Sigma_nfw_at_R(arg, pars->M, pars->conc, pars->delta, pars->om);
  }else if(arg < rmax){
    return gsl_spline_eval(spline, arg ,acc);
  }
  return 0;
}

double Sigma_mis_single_at_R(double R, double*Rs, double*Sigma, int Ns, double M, double conc, int delta, double om, double Rmis){
  double*Ra = (double*)malloc(sizeof(double));
  double*Smis = (double*)malloc(sizeof(double));
  double result;
  Ra[0] = R;
  Sigma_mis_single_at_R_arr(Ra, 1, Rs, Sigma, Ns, M, conc, delta, om, Rmis, Smis);
  result = Smis[0];
  free(Ra);
  free(Smis);
  return result;
}

int Sigma_mis_single_at_R_arr(double*R, int NR, double*Rs, double*Sigma, int Ns, double M, double conc, int delta, double om, double Rmis, double*Sigma_mis){
  gsl_spline*spline = gsl_spline_alloc(gsl_interp_cspline, Ns);
  gsl_interp_accel*acc = gsl_interp_accel_alloc();
  gsl_integration_workspace*workspace = gsl_integration_workspace_alloc(workspace_size);
  gsl_integration_workspace*workspace2 = gsl_integration_workspace_alloc(workspace_size);
  integrand_params*params = malloc(sizeof(integrand_params));
  gsl_function F;
  double result, err;
  int i;
  gsl_spline_init(spline, Rs, Sigma, Ns);
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
  F.function=&single_angular_integrand;
  F.params=params;
  for(i = 0; i < NR; i++){
    params->Rp  = R[i];
    params->Rp2 = R[i] * R[i];
    gsl_integration_qag(&F, 0, M_PI, ABSERR, RELERR, workspace_size, KEY, workspace, &result, &err);
    Sigma_mis[i] = result/M_PI;
  }
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  gsl_integration_workspace_free(workspace2);
  free(params);
  return 0;
}

/////////////////// SIGMA(R) INTEGRANDS BELOW //////////////////////

double get_Sigma(double Rc, double Rc2, void*params){
  integrand_params*pars = (integrand_params*)params;
  gsl_spline*spline = pars->spline;
  gsl_interp_accel*acc = pars->acc;
  double rmin = pars->rmin,rmax = pars->rmax;
  double Rp2 = pars->Rp2;
  double Rp_cos_theta_2 = pars->Rp_cos_theta_2;
  double arg = sqrt(Rp2 + Rc2 - Rc*Rp_cos_theta_2);
  double Sigma = 0;
  if(arg > rmin && arg < rmax){
    Sigma = gsl_spline_eval(spline, arg, acc);
  }else if(arg < rmin){
    Sigma = Sigma_nfw_at_R(arg, pars->M, pars->conc, pars->delta, pars->om);
  }
  return Sigma;
}

double exp_radial_integrand(double lRc, void*params){
  integrand_params*pars = (integrand_params*)params;
  double Rc = exp(lRc);
  double Rc2 = Rc*Rc;
  double Rmis = pars->Rmis;
  double Sigma = get_Sigma(Rc, Rc2, params);
  return Rc2 * exp(-Rc/Rmis) * Sigma; //normalized outside
}

double g2d_radial_integrand(double lRc, void*params){
  integrand_params*pars = (integrand_params*)params;
  double Rc = exp(lRc);
  double Rc2 = Rc*Rc;
  double Rmis2 = pars->Rmis2;
  double Sigma = get_Sigma(Rc, Rc2, pars);
  return Rc2 * exp(-0.5 * Rc2/Rmis2) * Sigma; //normalized outside
}

double angular_integrand(double theta, void*params){
  integrand_params*pars = (integrand_params*)params;
  double cos_theta = cos(theta);
  double lrmin = pars->lrmin, lrmax = pars->lrmax;
  gsl_integration_workspace*workspace = pars->workspace2;
  gsl_function F = pars->F_radial;
  double result, err;
  pars->Rp_cos_theta_2 = pars->Rp*cos_theta*2;
  gsl_integration_qag(&F, lrmin-10, lrmax, ABSERR, RELERR, workspace_size, KEY, workspace, &result, &err);
  return result;
}

double Sigma_mis_at_R(double R, double*Rs, double*Sigma, int Ns, double M, double conc, int delta, double om, double Rmis, int integrand_switch){
  double*Ra = (double*)malloc(sizeof(double));
  double*Sigma_mis = (double*)malloc(sizeof(double));
  double result;
  Ra[0] = R;
  Sigma_mis_at_R_arr(Ra, 1, Rs, Sigma, Ns, M, conc, delta, om, Rmis, integrand_switch, Sigma_mis);
  result = Sigma_mis[0];
  free(Ra);
  free(Sigma_mis);
  return result;
}

int Sigma_mis_at_R_arr(double*R, int NR, double*Rs, double*Sigma, int Ns, double M, double conc, int delta, double om, double Rmis, int integrand_switch, double*Sigma_mis){
  gsl_spline*spline = gsl_spline_alloc(gsl_interp_cspline, Ns);
  gsl_interp_accel*acc = gsl_interp_accel_alloc();
  gsl_integration_workspace*workspace = gsl_integration_workspace_alloc(workspace_size);
  gsl_integration_workspace*workspace2 = gsl_integration_workspace_alloc(workspace_size);
  integrand_params*params = malloc(sizeof(integrand_params));
  gsl_function F;
  gsl_function F_radial;
  double result, err;
  int i;
  gsl_spline_init(spline, Rs, Sigma, Ns);
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
  F.function = &angular_integrand;
  F_radial.function = &g2d_radial_integrand; //integrand_switch == 0
  if(integrand_switch == 1){
    F_radial.function = &exp_radial_integrand;
  }
  F_radial.params = params;
  params->F_radial = F_radial;
  F.params = params;
  for(i = 0; i < NR; i++){
    params->Rp  = R[i];
    params->Rp2 = R[i] * R[i];
    gsl_integration_qag(&F, 0, M_PI, ABSERR, RELERR, workspace_size, KEY, workspace, &result, &err);
    Sigma_mis[i] = result/(M_PI*Rmis*Rmis);
  }
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
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
  double*Ra = (double*)malloc(sizeof(double));
  double*DSm = (double*)malloc(sizeof(double));
  double result;
  Ra[0] = R;
  DeltaSigma_mis_at_R_arr(Ra, 1, Rs, Sigma, Ns, DSm);
  result = DSm[0];
  free(Ra);
  free(DSm);
  return result;
}

int DeltaSigma_mis_at_R_arr(double*R, int NR, double*Rs, double*Sigma, int Ns, double*DeltaSigma_mis){
  double lrmin = log(Rs[0]);
  gsl_spline*spline = gsl_spline_alloc(gsl_interp_cspline, Ns);
  gsl_interp_accel*acc = gsl_interp_accel_alloc();
  gsl_integration_workspace* workspace = gsl_integration_workspace_alloc(workspace_size);
  integrand_params*params=malloc(sizeof(integrand_params));
  double slope = log(Sigma[0]/Sigma[1])/log(Rs[0]/Rs[1]);
  double intercept = Sigma[0]*pow(Rs[0], -slope);
  double low_part = intercept*pow(Rs[0], slope+2)/(slope+2);
  double result,  err;
  gsl_function F;
  int i;
  gsl_spline_init(spline, Rs, Sigma, Ns);
  params->spline = spline;
  params->acc = acc;
  F.params = params;
  F.function = &DS_mis_integrand;
  for(i = 0; i < NR; i++){
    gsl_integration_qag(&F, lrmin, log(R[i]), ABSERR, RELERR, workspace_size, KEY, workspace, &result, &err);
    DeltaSigma_mis[i] = (low_part+result)*2/(R[i]*R[i]) - gsl_spline_eval(spline, R[i], acc);
  }
  gsl_spline_free(spline),gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  free(params);
  return 0;
}
