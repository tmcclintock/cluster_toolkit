/** @file C_triaxiality.c
 *  @brief Triaxiality and orientation angle affects on surface mass density profiles.
 *  
 *  These functions compute projected surface mass density profiles
 *  for galaxy clusters (halos) that can be ellipsoidal.
 *  
 *  @author Tom McClintock (tmcclintock)
 *  @bug No known bugs.
 */

#include "C_triaxiality.h"
#include "C_xi.h"

#include "gsl/gsl_integration.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_errno.h"

#include <math.h>
#include <stdio.h>

#define ABSERR 0.0
#define RELERR 1e-4
#define workspace_size 8000
#define ulim 5.0
#define rhocrit 2.77533742639e+11
//1e4*3.*Mpcperkm*Mpcperkm/(8.*PI*G); units are SM h^2/Mpc^3
#define KEY 5 //Used for GSL QAG function

typedef struct integrand_params{
  gsl_integration_workspace*workspace;
  gsl_integration_workspace*workspace2;
  double R;
  double M;
  double conc;
  double sin_i;
  double cos_i;
  double q;
  double s;
  double delta;
  double Omega_m;
  double z;
  gsl_function F_azimuthal; //Function for the azimuthal integrand
}integrand_params;



/**
 * \brief Radial vector norm in the frame of the triaxial halo,
 * mapped to the isodensity contour of a spherical halo
 * at r = \sqrt(z^2 + R^2).
 * Units are Msun/h comoving.
 */
double mapped_radii(double z, double R, double sin_i, double cos_i,
		    double sin_phi, double cos_phi, double q, double s){
  return sqrt(R * R * sin_phi * sin_phi / (q * q) +
	      pow(R * cos_phi * cos_i - z * sin_i, 2) / (s * s) +
	      pow(R * cos_phi * sin_i + z * cos_i, 2));
}

double azimuthal_integrand_Sigma_nfw(double phi, void*params){
  integrand_params pars=*(integrand_params*)params;
  double sin_phi = sin(phi);
  double cos_phi = cos(phi);
  double r = mapped_radii(pars.z, pars.R, pars.sin_i, pars.cos_i,
			  sin_phi, cos_phi, pars.q, pars.s);
  double xi = xi_nfw_at_r(r, pars.M, pars.conc, pars.delta, pars.Omega_m);
  //printf("(%.2e  %.2e  %.2e %.2e %.2e %.2e  r=%.2e  xi=%.2e)\n", pars.z, pars.R, pars.sin_i, pars.cos_i, pars.q, pars.s, r, xi);
  //exit(1);

  return xi;
}

double LOS_integrand_Sigma_nfw(double ln_z, void*params){
  double result, err;
  double z = exp(ln_z);
  integrand_params pars=*(integrand_params*)params;
  pars.z = z;
  gsl_integration_qag(&pars.F_azimuthal, 0, M_PI,
		      ABSERR, RELERR, workspace_size,
		      KEY, pars.workspace2, &result, &err);
  //printf("%e\n",result / );
  return z * result;
}

/**
 * \brief Surface mass density of an ellipsoidal NFW halo.
 * Units are hMsun/pc^2 comoving.
 */
double Ellipsoidal_Sigma_nfw_single_halo_at_R(double R, double M,
					      double c, double i,
					      double q, double s,
					      double delta, double Omega_m){
  double result = 0;
  Ellipsoidal_Sigma_nfw_single_halo_at_R_arr(&R, 1, M, c, i, q,
					     s, delta, Omega_m, &result);
  return result;
}

int Ellipsoidal_Sigma_nfw_single_halo_at_R_arr(double*R, int NR,
					       double M, double conc,
					       double i, double q, double s,
					       double delta, double Omega_m,
					       double*Sigma){

  //h^2Msun/pc^2/Mpc; integral is over Mpc/h
  double rhom = Omega_m*rhocrit*1e-12; 

  //Useful variables
  int j;
  double result, err;
  gsl_function F;
  gsl_function F_azimuthal;

  //Allocate things
  static int init_flag = 0;
  static gsl_integration_workspace*workspace = NULL;
  static gsl_integration_workspace*workspace2 = NULL;
  static integrand_params params;
  if (init_flag == 0){
    init_flag = 1;
    workspace = gsl_integration_workspace_alloc(workspace_size);
    workspace2 = gsl_integration_workspace_alloc(workspace_size);
  }

  //Set integrand parameters
  printf("%e %e %e  %e  %e   %e\n",M, conc, i, q, s, Omega_m);
  params.workspace = workspace;
  params.workspace2 = workspace2;
  params.M = M;
  params.sin_i = sin(i);
  params.cos_i = cos(i);
  params.q = q;
  params.s = s;
  params.conc = conc;
  params.delta = delta;
  params.Omega_m = Omega_m;
  F.function = &LOS_integrand_Sigma_nfw;
  F.params = &params;
  F_azimuthal.function = &azimuthal_integrand_Sigma_nfw;
  F_azimuthal.params = &params;
  params.F_azimuthal = F_azimuthal;

  for(j = 0; j < NR; j++){
    params.R = R[j];
    gsl_integration_qag(&F, -7, 7,
			ABSERR, RELERR, workspace_size,
			KEY, workspace, &result, &err);
    Sigma[j] = result * 2 * rhom / M_PI;
  }

  return 0;
}
