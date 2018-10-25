#include "C_exclusion.h"
#include "C_xi.h"
#include "C_peak_height.h"
#include "C_power.h"

#include "gsl/gsl_integration.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_sf.h" //for erfc
#include "gsl/gsl_errno.h"
#include <math.h>
#include <stdio.h>

#define rhocrit 2.77533742639e+11
//1e4*3.*Mpcperkm*Mpcperkm/(8.*PI*G); units are SM h^2/Mpc^3
#define sqrt2    1.4142135623730951 //sqrt(2), for computational efficiency
#define invsqrt2 0.7071067811865475 //1/sqrt(2), for computational efficiency

int xihm_exclusion_at_r(double*r, int Nr, double M, double c,
			double rt, double beta,
			double Ma, double ca, double Mb, double cb,
			double bias, double*ximm, int delta, double Omega_m,
			double*xihm){
  int i;
  double rhom = Omega_m * rhocrit; //SM h^2/Mpc^3 comoving
  double*xi_1h  = malloc(sizeof(double)*Nr);
  double*xi_2h  = malloc(sizeof(double)*Nr);
  double*xi_c   = malloc(sizeof(double)*Nr);
  xi_1h_at_r_arr(r, Nr, M, c, rt, beta, delta, Omega_m, xi_1h);
  //1halo and 2halo terms
  for(i = 0; i < Nr; i++){
    xi_2h[i] = bias*ximm[i];
    //xi_c[i]  = bias*otherstuff;
  }
  rhom = 0; //suppresses warning for now
  free(xi_1h);
  free(xi_2h);
  free(xi_c);
  return 0;
}

int xi_1h_at_r_arr(double*r, int Nr, double M, double c,
		   double rt, double beta, int delta, double Omega_m,
		   double*xi_1h){
  int i;
  double*thetas = malloc(sizeof(double)*Nr);
  calc_xi_nfw(r, Nr, M, c, delta, Omega_m, xi_1h);
  theta_erfc_at_r_arr(r, Nr, rt, beta, thetas);
  for(i = 0; i < Nr; i++){
    xi_1h[i] = xi_1h[i] * thetas[i];
  }
  return 0; //success
}

double xi_1h_at_r(double r, double M, double c,
		  double rt, double beta, int delta, double Omega_m){
  double*rs    = malloc(sizeof(double));
  double*xi_1h = malloc(sizeof(double));
  double result;
  rs[0] = r;
  xi_1h_at_r_arr(rs, 1, M, c, rt, beta, delta, Omega_m, xi_1h);
  result = xi_1h[0];
  free(rs);
  free(xi_1h);
  return result;
}

int utct_at_r_arr(double*r, int Nr, double rt, double M1, double M2, double conc, int delta, double Omega_m, double*utct){
  double rhom = Omega_m * rhocrit;//SM h^2/Mpc^3
  double r1 = pow(M1/(1.33333333333*M_PI*rhom*delta), 0.33333333333);
  double r2 = pow(M2/(1.33333333333*M_PI*rhom*delta), 0.33333333333);
  double rt1 = rt; //this line not really necessary...
  double con = rt1/r1;
  double rt2 = con*r2;
  rt2=0;
  //more stuff
  return 0;
}

int theta_erfc_at_r_arr(double*r, int Nr, double rt, double beta, double*theta){
  int i;
  double invbeta = 1./beta;
  for(i = 0; i < Nr; i++){
    theta[i] = 0.5*gsl_sf_erfc((r[i]-rt) * invbeta * invsqrt2);
  }
  return 0;
}

double theta_erfc_at_r(double r, double rt, double beta){
  double*rs     = malloc(sizeof(double));
  double*thetas = malloc(sizeof(double));
  double result;
  rs[0] = r;
  theta_erfc_at_r_arr(rs, 1, rt, beta, thetas);
  result = thetas[0];
  free(rs);
  free(thetas);
  return result;
}
