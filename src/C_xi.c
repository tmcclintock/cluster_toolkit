#include "C_xi.h"
#include "C_power.h"

#include "gsl/gsl_integration.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_sf_gamma.h"
#include <math.h>
#include <stdio.h>

#define rhomconst 2.77533742639e+11
//1e4*3.*Mpcperkm*Mpcperkm/(8.*PI*G); units are SM h^2/Mpc^3

double xi_nfw_at_R(double R, double Mass, double conc, int delta, double om){
  double*Rarr = (double*)malloc(sizeof(double));
  double*xi  = (double*)malloc(sizeof(double));
  double result;
  Rarr[0] = R;
  calc_xi_nfw(Rarr, 1, Mass, conc, delta, om, xi);
  result = xi[0];
  free(Rarr);
  free(xi);
  return result;
}

int calc_xi_nfw(double*R, int NR, double Mass, double conc, int delta, double om, double*xi_nfw){
  int i;
  double rhom = om*rhomconst;//SM h^2/Mpc^3
  double Rdelta = pow(Mass/(1.33333333333*M_PI*rhom*delta), 0.33333333333);
  double Rscale = Rdelta/conc;
  double fc = log(1.+conc)-conc/(1.+conc);
  double R_Rs;
  for(i = 0; i < NR; i++){
    R_Rs = R[i]/Rscale;
    xi_nfw[i] = Mass/(4.*M_PI*Rscale*Rscale*Rscale*fc)/(R_Rs*(1+R_Rs)*(1+R_Rs))/rhom - 1.0;
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

double xi_einasto_at_R(double R, double Mass, double rhos, double rs, double alpha, int delta, double om){
  double*Rarr = (double*)malloc(sizeof(double));
  double*xi  = (double*)malloc(sizeof(double));
  double result;
  Rarr[0] = R;
  calc_xi_einasto(Rarr, 1, Mass, rhos, rs, alpha, delta, om, xi);
  result = xi[0];
  free(Rarr);
  free(xi);
  return result;
}

int calc_xi_einasto(double*R, int NR, double Mass, double rhos, double rs, double alpha, int delta, double om, double*xi_einasto){
  double rhom = om*rhomconst;//SM h^2/Mpc^3
  double x;
  int i;
  if (rhos < 0)
    rhos = rhos_einasto_at_M(Mass, rs, alpha, delta, om);
  for(i = 0; i < NR; i++){
    x = 2./alpha * pow(R[i]/rs, alpha);
    xi_einasto[i] = rhos/rhom * exp(-x) - 1;
  }
  return 0;
}

int calc_xi_2halo(int NR, double bias, double*xi_mm, double*xi_2halo){
  int i;
  for(i = 0; i < NR; i++){
    xi_2halo[i] = bias * xi_mm[i];
  }
  return 0;
}

int calc_xi_hm(int NR, double*xi_1h, double*xi_2h, double*xi_hm, int flag){
  //Flag specifies how to combine the two terms
  int i;
  if (flag == 0) { //Take the max
    for(i = 0; i < NR; i++){
      if(xi_1h[i] >= xi_2h[i]) xi_hm[i] = xi_1h[i];
      else xi_hm[i] = xi_2h[i];
    }
  } else if (flag == 1) { //Take the sum
    for(i = 0; i < NR; i++){
      xi_hm[i] = 1 + xi_1h[i] + xi_2h[i];
    }
  }
  return 0;
}

int calc_xi_mm(double*R, int NR, double*k, double*P, int Nk, double*xi, int N, double h){
  double zero,psi,x,t,dpsi,f,PIsinht;
  double PI_h = M_PI/h;
  double PI_2 = M_PI*0.5;
  gsl_spline*Pspl = gsl_spline_alloc(gsl_interp_cspline, Nk);
  gsl_interp_accel*acc= gsl_interp_accel_alloc();
  double sum;  
  int i,j;
  gsl_spline_init(Pspl, k, P, Nk);
  for(j = 0; j < NR; j++){
    sum = 0;
    for(i = 0; i < N; i++){
      zero = i+1;
      psi = h*zero*tanh(sinh(h*zero)*PI_2);
      x = psi*PI_h;
      t = h*zero;
      PIsinht = M_PI*sinh(t);
      dpsi = (M_PI*t*cosh(t)+sinh(PIsinht))/(1+cosh(PIsinht));
      if (dpsi != dpsi) dpsi=1.0;
      f = x*get_P(x,R[j],k,P,Nk,Pspl,acc);
      sum += f*sin(x)*dpsi;
    }
    xi[j] = sum/(R[j]*R[j]*R[j]*M_PI*2);
  }

  gsl_spline_free(Pspl);
  gsl_interp_accel_free(acc);
  return 0; //Note: factor of pi picked up in the quadrature rule
  //See Ogata 2005 for details, especially eq. 5.2
}

///////Functions for calc_xi_mm/////////

double xi_mm_at_R(double R, double*k, double*P, int Nk, int N, double h){
  double*Ra = (double*)malloc(sizeof(double));
  double*xi = (double*)malloc(sizeof(double));
  double result;
  Ra[0] = R;
  calc_xi_mm(Ra, 1, k, P, Nk, xi, N, h);
  result = xi[0];
  free(Ra);
  free(xi);
  return result;
}

//////////////////////////////////////////
//////////////Xi(R) exact below //////////
//////////////////////////////////////////

#define workspace_size 8000
#define workspace_num 100
#define ABSERR 0.0
#define RELERR 1.8e-2

typedef struct integrand_params_xi_mm_exact{
  gsl_spline*spline;
  gsl_interp_accel*acc;
  gsl_integration_workspace * workspace;
  double r; //3d r; Mpc/h, or inverse units of k
  double*kp; //pointer to wavenumbers
  double*Pp; //pointer to P(k) array
  int Nk; //length of k and P arrays
}integrand_params_xi_mm_exact;

double integrand_xi_mm_exact(double k, void*params){
  integrand_params_xi_mm_exact pars = *(integrand_params_xi_mm_exact*)params;
  gsl_spline*spline = pars.spline;
  gsl_interp_accel*acc = pars.acc;
  double*kp = pars.kp;
  double*Pp = pars.Pp;
  int Nk = pars.Nk;
  double R = pars.r;
  double x  = k*R;
  double P = get_P(x, R, kp, Pp, Nk, spline, acc);
  return P*k/R; //Note - sin(kR) is taken care of in the qawo table
}

int calc_xi_mm_exact(double*R, int NR, double*k, double*P, int Nk, double*xi){
  gsl_spline*Pspl = gsl_spline_alloc(gsl_interp_cspline, Nk);
  gsl_interp_accel*acc= gsl_interp_accel_alloc();
  gsl_integration_workspace*workspace = gsl_integration_workspace_alloc(workspace_size);
  gsl_integration_qawo_table*wf;
  integrand_params_xi_mm_exact*params=malloc(sizeof(integrand_params_xi_mm_exact));
  gsl_function F;
  double kmax = 4e3;
  double kmin = 5e-8;
  double result, err;
  int i;
  gsl_spline_init(Pspl, k, P, Nk);
  params->acc = acc;
  params->spline = Pspl;
  params->kp = k;
  params->Pp = P;
  params->Nk = Nk;

  F.function=&integrand_xi_mm_exact;
  F.params=params;

  wf = gsl_integration_qawo_table_alloc(R[0], kmax-kmin, GSL_INTEG_SINE, (size_t)workspace_num);
  for(i = 0; i < NR; i++){
    gsl_integration_qawo_table_set(wf, R[i], kmax-kmin, GSL_INTEG_SINE);
    params->r=R[i];
    gsl_integration_qawo(&F, kmin, ABSERR, RELERR, (size_t)workspace_num, workspace, wf, &result, &err);
    xi[i] = result/(M_PI*M_PI*2);
  }

  free(params);
  gsl_spline_free(Pspl);
  gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  gsl_integration_qawo_table_free(wf);

  return 0;
}

double xi_mm_at_R_exact(double R, double*k, double*P, int Nk){
  double*Ra = (double*)malloc(sizeof(double));
  double*xi = (double*)malloc(sizeof(double));
  double result;
  Ra[0] = R;
  calc_xi_mm_exact(Ra, 1, k, P, Nk, xi);
  result = xi[0];
  free(Ra);
  free(xi);
  return result;
}
