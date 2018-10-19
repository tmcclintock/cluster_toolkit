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
  double*xi_1h  = (double*)malloc(sizeof(double)*Nr);
  double*thetas = (double*)malloc(sizeof(double)*Nr);
  double*xi_2h  = (double*)malloc(sizeof(double)*Nr);
  double*xi_c   = (double*)malloc(sizeof(double)*Nr);

  //1halo and 2halo terms
  calc_xi_nfw(r, Nr, M, c, delta, Omega_m, xi_1h);
  theta_erfc_at_r_arr(r, Nr, rt, beta, thetas);
  for(i = 0; i < Nr; i++){
    xi_1h[i] = xi_1h[i] * thetas[i];
    xi_2h[i] = bias*ximm[i];
    //xi_c[i]  = bias*otherstuff;
  }
  rhom = 0; //suppresses warning for now
  free(xi_1h);
  free(thetas);
  free(xi_2h);
  free(xi_c);
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
  double*rs     = (double*)malloc(sizeof(double));
  double*thetas = (double*)malloc(sizeof(double));
  double result;
  rs[0] = r;
  theta_erfc_at_r_arr(rs, 1, rt, beta, thetas);
  result = thetas[0];
  free(rs);
  free(thetas);
  return result;
}
