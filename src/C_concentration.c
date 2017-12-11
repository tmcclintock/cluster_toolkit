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
double DK15_concentration_at_Mmean(double Mass, double*k, double*P, int Nk, int delta, double Omega_m){
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

    /*
    if (status == GSL_SUCCESS)
      printf("Converged:\n");
    printf ("%3d [%.2e, %.2e] %.2e %+.2e %.2e\n",
    iter, M_lo, M_hi,
    Mm, Mm - Mass, 
    M_hi - M_lo);
    */
  }while(status == GSL_CONTINUE && iter < max_iter);
  
  cm = pars->c;

  free(pars);
  gsl_root_fsolver_free(s);
  return cm;
}

//We need to implement the M-c equation just for M200crit first
//This is the median M-c relation from DK15
double DK15_concentration_at_Mcrit(double Mass, double*k, double*P, int Nk, int delta, double Omega_m){
  double nu = nu_at_M(Mass, k, P, Nk, Omega_m);
  double R = M_to_R(Mass, Omega_m); //Lagrangian Radius
  double lnk_R = log10(0.69 * 2*M_PI/R); //0.69 is DK15 kappa
  double*lnk = (double*)malloc(Nk*sizeof(double));
  double*lnP = (double*)malloc(Nk*sizeof(double));
  int i;
  for(i = 0; i < Nk; i++){
    lnk[i] = log10(k[i]);
    lnP[i] = log10(P[i]);
  }
  gsl_spline*lnPspl = gsl_spline_alloc(gsl_interp_cspline, Nk);
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
  return 0.5 * c0 * (pow(nu0/nu, alpha) + pow(nu/nu0, beta));
}
