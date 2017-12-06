#include "C_concentration.h"
#include "C_bias.h"

#include "gsl/gsl_roots.h"
#include "gsl/gsl_spline.h"

#include <math.h>
#include <stdio.h>

#define rhocrit 2.77533742639e+11
//1e4*3.*Mpcperkm*Mpcperkm/(8.*PI*G); units are Msun h^2/Mpc^3

//Structure to hold the parameters for the M-c root finder
typedef struct mc_params{
  double Mm;
  double Rm;
  double c;//will be the output
  double*k;
  double*P;
  int Nk;
  int delta;
  double Omega_m;
}mc_params;

double mass_comparison_function(double M, void*params){
  //M is M200c
  mc_params*pars = (mc_params*)params;
  double Mm = pars->Mm;
  //double Rm = pars->Rm;
  double*k = pars->k;
  double*P = pars->P;
  int Nk = pars->Nk;
  int delta = pars->delta;
  double Omega_m = pars->Omega_m;
  pars->c = DK15_concentration_at_Mcrit(M, k, P, Nk, delta, Omega_m);
  //double M_integral = stuff;
  return Mm - M;//WRONG for now
}

//This is for M200m(b)
double DK15_concentration_at_Mmean(double Mass, double*k, double*P, int Nk, int delta, double Omega_m){
  mc_params*pars = (mc_params*)malloc(sizeof(mc_params));
  double R = pow(Mass/(4./3.*M_PI*rhocrit*Omega_m*delta), 0.33333333); //R200m
  pars->Mm = Mass;
  pars->Rm = R;
  pars->k = k;
  pars->P = P;
  pars->Nk = Nk;
  pars->delta = delta;
  pars->Omega_m = Omega_m;
  /*
    Need to do root finding to solve the equation
    Mass - int_0^R200m dr r^2*rho_nfw(r|M200c, c200c) = 0.
    The solution to the right term is analytic. Check wolframalpha.
   */
  
  return R; //this is wrong right now
}

//We need to implement the M-c equation just for M200crit first
//This is the median M-c relation from DK15
double DK15_concentration_at_Mcrit(double Mass, double*k, double*P, int Nk, int delta, double Omega_m){
  double nu = nu_at_M(Mass, k, P, Nk, Omega_m);
  double R = M_to_R(Mass, Omega_m);
  double lnk_R = log(0.69 * 2*M_PI/R); //0.69 is DK15 kappa
  double*lnk = (double*)malloc(Nk*sizeof(double));
  double*lnP = (double*)malloc(Nk*sizeof(double));
  int i;
  for(i = 0; i < Nk; i++){
    lnk[i] = log(k[i]);
    lnP[i] = log(P[i]);
  }
  gsl_spline*lnPspl = gsl_spline_alloc(gsl_interp_cspline,Nk);
  gsl_spline_init(lnPspl, lnk, lnP, Nk);
  gsl_interp_accel*acc= gsl_interp_accel_alloc();
  double n = gsl_spline_eval_deriv(lnPspl, lnk_R, acc);
  //Free things we don't need anymore
  gsl_spline_free(lnPspl);
  gsl_interp_accel_free(acc);
  free(lnk);
  free(lnP);
  double phi0  = 6.58;
  double phi1  = 1.37;
  double eta0  = 6.82;
  double eta1  = 1.42;
  double alpha = 1.12;
  double beta  = 1.69;
  double c0 = phi0 + n * phi1;
  double nu0   = eta0 + n * eta1;
  return 0.5 * c0 * (pow(nu0/nu, -alpha) + pow(nu/nu0, beta));
}
