/** @file C_profile_derivatives.c
 *  @brief Halo density profile derivatives
 *
 *  These functions are derivatives of different models of density profiles for halos.
 *
 *  @author Tom McClintock (tmcclintock)
 *  @bug No known bugs.
 */

#include "C_density.h"
#include "C_power.h"
#include "C_profile_derivatives.h"
#include "C_xi.h"

#include "gsl/gsl_integration.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_sf_gamma.h"
#include <math.h>
#include <stdio.h>

#define rhocrit 2.77533742639e+11
//1e4*3.*Mpcperkm*Mpcperkm/(8.*PI*G); units are Msun h^2/Mpc^3

double drho_nfw_dr_at_R(double R, double Mass, double conc, int delta, double Omega_m){
  double*Rarr = (double*)malloc(sizeof(double));
  double*drhodr  = (double*)malloc(sizeof(double));
  double result;
  Rarr[0] = R;
  drho_nfw_dr_at_R_arr(Rarr, 1, Mass, conc, delta, Omega_m, drhodr);
  result = drhodr[0];
  free(Rarr);
  free(drhodr);
  return result;
}

int drho_nfw_dr_at_R_arr(double*R, int NR, double Mass, double conc,
			 int delta, double Omega_m, double*drhodr){
  int i;
  double rhom = Omega_m*rhocrit;//Msun h^2/Mpc^3
  double Rdelta = pow(Mass/(1.33333333333*M_PI*rhom*delta), 0.33333333333);
  double Rscale = Rdelta/conc;
  double fc = log(1.+conc)-conc/(1.+conc);
  double rho0 = Mass/(4.*M_PI*Rscale*Rscale*Rscale*fc); //characteristic density
  double R_Rs;
  for(i = 0; i < NR; i++){
    R_Rs = R[i]/Rscale;
    drhodr[i] = - rho0*(1+3*R_Rs) / (Rscale * R_Rs*R_Rs * (1+R_Rs)*(1+R_Rs));
  }
  return 0;
}

typedef struct integrand_params_profile_derivs{
  gsl_spline*spline;
  gsl_interp_accel*acc;
  gsl_integration_workspace * workspace;
  double r; //3d r; Mpc/h, or inverse units of k
  double*kp; //pointer to wavenumbers
  double*Pp; //pointer to P(k) array
  int Nk; //length of k and P arrays
}integrand_params_profile_derivs;


double integrand_dxi_mm_dr_COSINE(double k, void*params){
  integrand_params_profile_derivs pars
    = *(integrand_params_profile_derivs*)params;
  gsl_spline*spline = pars.spline;
  gsl_interp_accel*acc = pars.acc;
  double*kp = pars.kp;
  double*Pp = pars.Pp;
  int Nk = pars.Nk;
  double R = pars.r;
  double x  = k*R;
  double P = get_P(x, R, kp, Pp, Nk, spline, acc);
  return P*k*k/R; //Note - cos(kR) is taken care of in the qawo table
}

double integrand_dxi_mm_dr_SINE(double k, void*params){
  integrand_params_profile_derivs pars
    = *(integrand_params_profile_derivs*)params;
  gsl_spline*spline = pars.spline;
  gsl_interp_accel*acc = pars.acc;
  double*kp = pars.kp;
  double*Pp = pars.Pp;
  int Nk = pars.Nk;
  double R = pars.r;
  double x  = k*R;
  double P = get_P(x, R, kp, Pp, Nk, spline, acc);
  return P*k/(R*R); //Note - sin(kR) is taken care of in the qawo table
}
