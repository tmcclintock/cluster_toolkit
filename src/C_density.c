#include "C_density.h"
#include "C_xi.h"

#include "gsl/gsl_spline.h"
#include "gsl/gsl_sf_gamma.h"
#include <math.h>
#include <stdio.h>

#define rhomconst 2.77533742639e+11
//#define rhomconst 2.775808e+11 //old value used. Disregard
//1e4*3.*Mpcperkm*Mpcperkm/(8.*PI*G); units are Msun h^2/Mpc^3

double rho_nfw_at_R(double R, double Mass, double conc, int delta, double om){
  double rhom = om*rhomconst;//Msun h^2/Mpc^3
  double xi = xi_nfw_at_R(R, Mass, conc, delta, om);
  return rhom*(1+xi);
}

int calc_rho_nfw(double*R, int NR, double Mass, double conc, int delta, double om, double*rho_nfw){
  int i;
  double rhom = om*rhomconst;//Msun h^2/Mpc^3
  calc_xi_nfw(R, NR, Mass, conc, delta, om, rho_nfw); //rho_nfw actually holds xi_nfw here
  for(i = 0; i < NR; i++){
    rho_nfw[i] = rhom*(1+rho_nfw[i]);
  }
  return 0;
}

double rhos_einasto_at_M(double Mass, double rs, double alpha, int delta, double om){
  double rhom = om*rhomconst;//Msun h^2/Mpc^3
  // Rdelta in Mpc/h comoving
  double Rdelta = pow(Mass/(1.3333333333333*M_PI*rhom*delta), 0.333333333333);
  double x = 2./alpha * pow(Rdelta/rs, alpha); 
  double a = 3./alpha;
  double gam = gsl_sf_gamma(a) - gsl_sf_gamma_inc(a, x);
  double num = delta * rhom * Rdelta*Rdelta*Rdelta * alpha * pow(2./alpha, a);
  double den = 3. * rs*rs*rs * gam;
  return num/den;
}

double rho_einasto_at_R(double R, double Mass, double rhos, double rs, double alpha, int delta, double om){
  double*Rarr = (double*)malloc(sizeof(double));
  double*rhoe = (double*)malloc(sizeof(double));
  Rarr[0] = R;
  calc_rho_einasto(Rarr, 1, Mass, rhos, rs, alpha, delta, om, rhoe);
  double result = rhoe[0];
  free(Rarr);
  free(rhoe);
  return result;
}

int calc_rho_einasto(double*R, int NR, double Mass, double rhos, double rs, double alpha, int delta, double om, double*rho_einasto){
  if (rhos < 0)
    rhos = rhos_einasto_at_M(Mass, rs, alpha, delta, om);
  int i;
  double x;
  for(i = 0; i < NR; i++){
    x = 2./alpha * pow(R[i]/rs, alpha);
    rho_einasto[i] = rhos * exp(-x);
  }
  
  return 0;
}
