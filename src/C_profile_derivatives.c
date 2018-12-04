/** @file C_profile_derivatives.c
 *  @brief Halo density profile derivatives
 *
 *  These functions are derivatives of different models of density profiles for halos.
 *
 *  @author Tom McClintock (tmcclintock)
 *  @bug No known bugs.
 */

#include "C_density.h"
#include "C_profile_derivatives.h"
#include "C_xi.h"

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
