#include "C_concentration.h"
#include "C_bias.h"

#include "gsl/gsl_math.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_roots.h"
#include "gsl/gsl_spline.h"

#include <math.h>
#include <stdio.h>

#define rhocrit 2.77533742639e+11
//1e4*3.*Mpcperkm*Mpcperkm/(8.*PI*G); units are Msun h^2/Mpc^3

/*
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

double Mm_from_Mc(double Mc, void*params){
  //M is M200c
  mc_params*pars = (mc_params*)params;
  double Mm = pars->Mm;
  double Rm = pars->Rm;
  double*k = pars->k;
  double*P = pars->P;
  int Nk = pars->Nk;
  int delta = pars->delta;
  double Omega_m = pars->Omega_m;
  double rhom = rhocrit*Omega_m;
  double cc = DK15_concentration_at_Mcrit(Mc, k, P, Nk, delta, Omega_m);
  pars->c = cc;
  //Figure out the total mass, but first figure out rho0 for Mcrit
  double rho0c = delta*rhom*cc*cc*cc/((cc+2)/(cc+1)+log(1+cc));
  double Rc = pow(Mc/(4./3.*M_PI*rhocrit*delta), 0.33333333); //R200c
  double Rs_c = Rc/cc; //Scale radius of Mcrit
  double Rm_Rsc = Rm/Rs_c;
  double Mout = 4*M_PI*rho0c*Rs_c*Rs_c*Rs_c*((Rm_Rsc+2)/(Rm_Rsc+1)+log(1+Rm_Rsc));
  return Mm - Mout;
}

//This is for M200m(b)
double DK15_concentration_at_Mmean(double Mass, double*k, double*P, int Nk, int delta, double n_s, double Omega_b, double Omega_m, double h, double T_CMB){
  int status;
  mc_params*pars = (mc_params*)malloc(sizeof(mc_params));
  double R = pow(Mass/(4./3.*M_PI*rhocrit*Omega_m*delta), 0.33333333); //R200m
  double M_lo = Mass/10;
  double M_hi = Mass*10;
  double Mm;
  double cm = -1; //will have result
  const gsl_root_fsolver_type*T=gsl_root_fsolver_bisection;
  gsl_root_fsolver*s = gsl_root_fsolver_alloc(T);
  gsl_function F;
  int iter = 0, max_iter = 100;
  pars->Mm = Mass;
  pars->Rm = R;
  pars->k = k;
  pars->P = P;
  pars->Nk = Nk;
  pars->delta = delta;
  pars->Omega_m = Omega_m;
  
  F.function = &Mm_from_Mc;
  F.params = pars;

  status = gsl_root_fsolver_set(s, &F, M_lo, M_hi);

  do{
    iter++;
    status = gsl_root_fsolver_iterate(s);
    Mm = gsl_root_fsolver_root(s);
    M_lo = gsl_root_fsolver_x_lower(s);
    M_hi = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(M_lo, M_hi, 0, 0.001);
  }while(status == GSL_CONTINUE && iter < max_iter);
  
  cm = pars->c;

  free(pars);
  gsl_root_fsolver_free(s);
  return cm;
}
*/

//Just the signature
double dlnP_dlnk(double kin, double n_s, double Omega_b, double Omega_m, double h, double T_CMB);


//We need to implement the M-c equation just for M200crit first
//This is the median M-c relation from DK15
double DK15_concentration_at_Mcrit(double Mass, double*k, double*P, int Nk, int delta, double n_s, double Omega_b, double Omega_m, double h, double T_CMB){
  double nu = nu_at_M(Mass, k, P, Nk, Omega_m);
  double R = M_to_R(Mass, Omega_m); //Lagrangian Radius
  double k_R = 0.69 * 2*M_PI/R;
  double n = dlnP_dlnk(k_R, n_s, Omega_b, Omega_m, h, T_CMB);  
  double phi0  = 6.58;
  double phi1  = 1.37;
  double eta0  = 6.82;
  double eta1  = 1.42;
  double alpha = 1.12;
  double beta  = 1.69;
  double c0 = phi0 + n * phi1;
  double nu0   = eta0 + n * eta1;
  return 0.5 * c0 * (pow(nu0/nu, alpha) + pow(nu/nu0, beta));
}

/*
This is an equation for P(k) of E&H with no BAO.
this is what the DK15 M-c relation uses in its derivative.
It may not be feasible to incorporate it into here, given the amount of
front end overhaul it would require.
*/
double transferFunc_EH98_zeroBaryon(double kin, double Omega_b, double Omega_m, double h, double T_CMB){


  double k, Tk;
  double omb, om0, omc;
  double obh2, omh2, och2;
  double ob_om;
  double theta2p7, s, q;
  double Gamma, alphaGamma, L0, C0;
  
  omb = Omega_b;
  om0 = Omega_m;
  omc = Omega_m - Omega_b;
  obh2 = Omega_b*h*h;
  omh2 = Omega_m*h*h;
  och2 = omc*h*h;
  ob_om = omb / om0;
  theta2p7 = T_CMB / 2.7;
  
  //convert k from hMpc^-1 to Mpc^-1
  k = kin * h;

  //eqn 26
  s = 44.5*log(9.83/om0/h/h)/sqrt(1.0 + 10.0*pow(omb*h*h,0.75));

  //eqn 31
  alphaGamma = 1.0 - 0.328*log(431.0*om0*h*h)*omb/om0 + 0.38*log(22.3*om0*h*h)*(omb/om0)*(omb/om0);

  //eqn 30
  Gamma = om0*h*(alphaGamma + (1.0 - alphaGamma)/(1.0 + pow(0.43*k*s,4.0)));

  //eqn 28
  q = kin * theta2p7 * theta2p7 / Gamma;

  //eqns 29
  C0 = 14.2 + 731.0 / (1.0 + 62.5 * q);
  L0 = log(2.0 * exp(1.0) + 1.8 * q);
  Tk = L0 / (L0 + C0 * q * q);

  return Tk;
}

double dlnP_dlnk(double kin, double n_s, double Omega_b, double Omega_m, double h, double T_CMB){
  //kin needs to have units of h/Mpc
  double dlnk = 1e-6;
  double dk = dlnk*kin;
  double T1 = transferFunc_EH98_zeroBaryon(kin+dk*0.5, Omega_b,  Omega_m,  h,  T_CMB);
  double T2 = transferFunc_EH98_zeroBaryon(kin-dk*0.5, Omega_b,  Omega_m,  h,  T_CMB);
  double P1 = pow(kin+dk*0.5, n_s)*T1*T1;
  double P2 = pow(kin-dk*0.5, n_s)*T2*T2;
  return log(P1/P2)/dlnk;
}
