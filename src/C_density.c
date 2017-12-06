#include "C_density.h"

#include "gsl/gsl_spline.h"
#include "gsl/gsl_sf_gamma.h"
#include <math.h>
#include <stdio.h>

#define rhomconst 2.775808e+11
//1e4*3.*Mpcperkm*Mpcperkm/(8.*PI*G); units are Msun h^2/Mpc^3

double rho_nfw_at_R(double R, double Mass, double conc, int delta, double om){
  double rhom = om*rhomconst;//Msun h^2/Mpc^3
  double Rdelta = pow(Mass/(1.33333333333*M_PI*rhom*delta), 0.33333333333);
  double Rscale = Rdelta/conc;
  double fc = log(1.+conc)-conc/(1.+conc); 
  return Mass/(4.*M_PI*Rscale*Rscale*Rscale*fc)/(R/Rscale*(1+R/Rscale)*(1+R/Rscale));
}

int calc_rho_nfw(double*R, int NR, double Mass, double conc, int delta, double om, double*rho_nfw){
  int i;
  for(i = 0; i < NR; i++)
    rho_nfw[i] = rho_nfw_at_R(R[i], Mass, conc, delta, om);
  return 0;
}

double rho_einasto_at_R(double R, double Mass, double rhos, double rs, double alpha, int delta, double om){
  double x;
  double rhom = om*rhomconst;//Msun h^2/Mpc^3;
  if (rhos < 0){
    //Calculate rhos from Mass
    // Rdelta in Mpc/h comoving
    double Rdelta = pow(Mass/(1.3333333333333*M_PI*rhom*delta), 0.333333333333);
    x = 2./alpha * pow(Rdelta/rs, alpha); // Mpc/h comoving
    double a = 3./alpha;
    double gam = gsl_sf_gamma(a) - gsl_sf_gamma_inc(a, x);
    double num = delta * rhom * Rdelta*Rdelta*Rdelta * alpha * pow(2./alpha, a);
    double den = 3. * rs*rs*rs * gam;
    rhos = num/den;
  }
  x = 2./alpha * pow(R/rs, alpha);
  return rhos * exp(-x);  
}

int calc_rho_einasto(double*R, int NR, double Mass, double rhos, double rs, double alpha, int delta, double om, double*rho_einasto){
  if (rhos < 0){
    //Calculate rhos from Mass
    double rhom = om*rhomconst;//Msun h^2/Mpc^3
    // Rdelta in Mpc/h comoving
    double Rdelta = pow(Mass/(1.3333333333333*M_PI*rhom*delta), 0.333333333333);
    double x = 2./alpha * pow(Rdelta/rs, alpha); 
    double a = 3./alpha;
    double gam = gsl_sf_gamma(a) - gsl_sf_gamma_inc(a, x);
    double num = delta * rhom * Rdelta*Rdelta*Rdelta * alpha * pow(2./alpha, a);
    double den = 3. * rs*rs*rs * gam;
    rhos = num/den;
  }
  int i;
  for(i = 0; i < NR; i++)
    rho_einasto[i] = rho_einasto_at_R(R[i], Mass, rhos, rs, alpha, delta, om);

  return 0;
}
