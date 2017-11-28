#include "C_concentration.h"
#include "C_bias.h"

#include "gsl/gsl_spline.h"

#include <math.h>
#include <stdio.h>

double DK15_concentration_at_M(double Mass, double*k, double*P, int Nk, double Omega_m){
  double nu = nu_at_M(Mass, k, P, Nk, Omega_m);
  //By default we use the median c-M relation presented in DK15
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
