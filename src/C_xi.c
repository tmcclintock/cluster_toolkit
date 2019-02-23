#include "C_xi.h"
#include "C_peak_height.h"
#include "C_power.h"

#include "gsl/gsl_integration.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_errno.h"
#include <math.h>
#include <stdio.h>

#define rhomconst 2.77533742639e+11
//1e4*3.*Mpcperkm*Mpcperkm/(8.*PI*G); units are SM h^2/Mpc^3


/** @brief The NFW correlation function profile.
 * 
 *  The NFW correlation function profile of a halo a distance r from the center,
 *  assuming the halo has a given mass and concentration. It works 
 *  with any overdensity parameter and arbitrary matter fraction.
 *  This function calls xi_nfw_at_r_arr().
 * 
 *  @param r Distance from the center of the halo in Mpc/h comoving.
 *  @param M Halo mass in Msun/h.
 *  @param c Concentration.
 *  @delta Halo overdensity.
 *  @Omega_m Matter fraction.
 *  @return NFW halo correlation function.
 */
double xi_nfw_at_r(double r, double Mass, double conc, int delta, double om){
  double*rarr = (double*)malloc(sizeof(double));
  double*xi  = (double*)malloc(sizeof(double));
  double result;
  rarr[0] = r;
  calc_xi_nfw(rarr, 1, Mass, conc, delta, om, xi);
  result = xi[0];
  free(rarr);
  free(xi);
  return result;
}

int calc_xi_nfw(double*r, int Nr, double Mass, double conc, int delta, double om, double*xi_nfw){
  int i;
  double rhom = om*rhomconst;//SM h^2/Mpc^3
  double rho0_rhom = delta/(3.*(log(1.+conc)-conc/(1.+conc)));
  double rdelta = pow(Mass/(1.33333333333*M_PI*rhom*delta), 0.33333333333);
  double rscale = rdelta/conc;
  double r_rs;
  for(i = 0; i < Nr; i++){
    r_rs = r[i]/rscale;
    xi_nfw[i] = rho0_rhom/(r_rs*(1+r_rs)*(1+r_rs)) - 1.;
  }
  return 0;
}

double rhos_einasto_at_M(double Mass, double conc, double alpha, int delta, double om){
  double rhom = om*rhomconst;//Msun h^2/Mpc^3
  // rdelta in Mpc/h comoving
  double rdelta = pow(Mass/(1.3333333333333*M_PI*rhom*delta), 0.333333333333);
  double rs = rdelta/conc; //Scale radius in Mpc/h
  double x = 2./alpha * pow(conc, alpha); //pow(rdelta/rs, alpha); 
  double a = 3./alpha;
  double gam = gsl_sf_gamma(a) - gsl_sf_gamma_inc(a, x);
  double num = delta * rhom * rdelta*rdelta*rdelta * alpha * pow(2./alpha, a);
  double den = 3. * rs*rs*rs * gam;
  return num/den;
}

int calc_xi_einasto(double*r, int Nr, double Mass, double rhos, double conc, double alpha, int delta, double om, double*xi_einasto){
  double rhom = om*rhomconst;//SM h^2/Mpc^3
  double rdelta = pow(Mass/(1.3333333333333*M_PI*rhom*delta), 0.333333333333);
  double rs = rdelta/conc; //Scale radius in Mpc/h
  double x;
  int i;
  if (rhos < 0)
    rhos = rhos_einasto_at_M(Mass, conc, alpha, delta, om);
  for(i = 0; i < Nr; i++){
    x = 2./alpha * pow(r[i]/rs, alpha);
    xi_einasto[i] = rhos/rhom * exp(-x) - 1;
  }
  return 0;
}

int calc_xi_2halo(int Nr, double bias, double*xi_mm, double*xi_2halo){
  int i;
  for(i = 0; i < Nr; i++){
    xi_2halo[i] = bias * xi_mm[i];
  }
  return 0;
}

int calc_xi_hm(int Nr, double*xi_1h, double*xi_2h, double*xi_hm, int flag){
  //Flag specifies how to combine the two terms
  int i;
  if (flag == 0) { //Take the max
    for(i = 0; i < Nr; i++){
      if(xi_1h[i] >= xi_2h[i]) xi_hm[i] = xi_1h[i];
      else xi_hm[i] = xi_2h[i];
    }
  } else if (flag == 1) { //Take the sum
    for(i = 0; i < Nr; i++){
      xi_hm[i] = 1 + xi_1h[i] + xi_2h[i];
    }
  }
  return 0;
}

int calc_xi_mm(double*r, int Nr, double*k, double*P, int Nk, double*xi, int N, double h){
  int i,j;
  double PI_h = M_PI/h;
  double PI_2 = M_PI*0.5;
  gsl_spline*Pspl = gsl_spline_alloc(gsl_interp_cspline, Nk);
  gsl_interp_accel*acc= gsl_interp_accel_alloc();
  gsl_spline_init(Pspl, k, P, Nk);
  double sum;
  
  static int init_flag = 0;
  static double h_old = -1;
  static int N_old = -1;
  static double*x      = NULL;
  static double*xsinx  = NULL;
  static double*dpsi   = NULL;
  static double*xsdpsi = NULL;
  double t, psi, PIsinht;
  if ((init_flag == 0) || (h_old != h) || (N_old < N)){
    h_old = h;
    N_old = N;
    init_flag = 1; //been initiated
    
    if (x!=NULL)      free(x);
    if (xsinx!=NULL)  free(xsinx);
    if (dpsi!=NULL)   free(dpsi);
    if (xsdpsi!=NULL) free(xsdpsi);
    
    x      = malloc(N*sizeof(double));
    xsinx  = malloc(N*sizeof(double));
    dpsi   = malloc(N*sizeof(double));
    xsdpsi = malloc(N*sizeof(double));
    for(i = 0; i < N; i++){
      t = h*(i+1);
      psi = t*tanh(sinh(t)*PI_2);
      x[i] = psi*PI_h;
      xsinx[i] = x[i]*sin(x[i]);
      PIsinht = M_PI*sinh(t);
      dpsi[i] = (M_PI*t*cosh(t) + sinh(PIsinht))/(1+cosh(PIsinht));
      if (dpsi[i]!=dpsi[i]) dpsi[i]=1.0;
      xsdpsi[i] = xsinx[i]*dpsi[i];
    }
  }
  for(j = 0; j < Nr; j++){
    sum = 0;
    for(i = 0; i < N; i++){
      sum += xsdpsi[i] * get_P(x[i], r[j], k, P, Nk, Pspl, acc);
    }
    xi[j] = sum/(r[j]*r[j]*r[j]*M_PI*2);
  }
  
  /* //original code below, left in for testing for now
  double zero,psi,x,t,dpsi,f,PIsinht;
  for(j = 0; j < Nr; j++){
    sum = 0;
    for(i = 0; i < N; i++){
      zero = i+1;
      psi = h*zero*tanh(sinh(h*zero)*PI_2);
      x = psi*PI_h;
      t = h*zero;
      PIsinht = M_PI*sinh(t);
      dpsi = (M_PI*t*cosh(t)+sinh(PIsinht))/(1+cosh(PIsinht));
      if (dpsi != dpsi) dpsi=1.0;
      f = x*get_P(x,r[j],k,P,Nk,Pspl,acc);
      sum += f*sin(x)*dpsi;
    }
    xi[j] = sum/(r[j]*r[j]*r[j]*M_PI*2);
  }
  */

  gsl_spline_free(Pspl);
  gsl_interp_accel_free(acc);
  return 0; //Note: factor of pi picked up in the quadrature rule
  //See Ogata 2005 for details, especially eq. 5.2
}

///////Functions for calc_xi_mm/////////

//////////////////////////////////////////
//////////////Xi(r) exact below //////////
//////////////////////////////////////////

#define workspace_size 8000
#define workspace_num 100
#define ABSERR 0.0
#define RELERR 1.8e-4

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
  double r = pars.r;
  double x  = k*r;
  double P = get_P(x, r, kp, Pp, Nk, spline, acc);
  return P*k/r; //Note - sin(kr) is taken care of in the qawo table
}

int calc_xi_mm_exact(double*r, int Nr, double*k, double*P, int Nk, double*xi){
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
  int status;
  gsl_spline_init(Pspl, k, P, Nk);
  params->acc = acc;
  params->spline = Pspl;
  params->kp = k;
  params->Pp = P;
  params->Nk = Nk;

  F.function=&integrand_xi_mm_exact;
  F.params=params;

  wf = gsl_integration_qawo_table_alloc(r[0], kmax-kmin, GSL_INTEG_SINE, (size_t)workspace_num);
  for(i = 0; i < Nr; i++){
    status = gsl_integration_qawo_table_set(wf, r[i], kmax-kmin, GSL_INTEG_SINE);
    if (status){
      printf("Error in calc_xi_mm_exact, first integral.\n");
      exit(-1);
    }
    params->r=r[i];
    status = gsl_integration_qawo(&F, kmin, ABSERR, RELERR, (size_t)workspace_num, workspace, wf, &result, &err);
    if (status){
      printf("Error in calc_xi_mm_exact, second integral.\n");
      exit(-1);
    }

    xi[i] = result/(M_PI*M_PI*2);
  }

  free(params);
  gsl_spline_free(Pspl);
  gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(workspace);
  gsl_integration_qawo_table_free(wf);

  return 0;
}

double xi_mm_at_r_exact(double r, double*k, double*P, int Nk){
  double*ra = malloc(sizeof(double));
  double*xi = malloc(sizeof(double));
  double result;
  ra[0] = r;
  calc_xi_mm_exact(ra, 1, k, P, Nk, xi);
  result = xi[0];
  free(ra);
  free(xi);
  return result;
}

/*
 * Diemer-Kravtsov 2014 profiles below.
 */

int calc_xi_DK(double*r, int Nr, double M, double rhos, double conc, double be, double se, double alpha, double beta, double gamma, int delta, double*k, double*P, int Nk, double om, double*xi){
  double rhom = rhomconst*om; //SM h^2/Mpc^3
  //Compute rDeltam
  double rdelta = pow(M/(1.33333333333*M_PI*rhom*delta), 0.33333333333);
  xi[0] = 0;
  double*rho_ein = malloc(Nr*sizeof(double));
  double*f_trans = malloc(Nr*sizeof(double));
  double*rho_outer = malloc(Nr*sizeof(double));
  int i;
  double nu = nu_at_M(M, k, P, Nk, om);
  if (alpha < 0){ //means it wasn't passed in
    alpha = 0.155 + 0.0095*nu*nu;
  }
  if (beta < 0){ //means it wasn't passed in
    beta = 4;
  }
  if (gamma < 0){ //means it wasn't passed in
    gamma = 8;
  }
  if (rhos < 0){ //means it wasn't passed in
    rhos = rhos_einasto_at_M(M, conc, alpha, delta, om);
  }
  double g_b = gamma/beta;
  double r_t = (1.9-0.18*nu)*rdelta;
  calc_xi_einasto(r, Nr, M, rhos, conc, alpha, delta, om, rho_ein);
  //here convert xi_ein to rho_ein
  for(i = 0; i < Nr; i++){
    rho_ein[i] = rhom*(1+rho_ein[i]); //rho_ein had xi_ein in it
    f_trans[i] = pow(1+pow(r[i]/r_t,beta), -g_b);
    rho_outer[i] = rhom*(be*pow(r[i]/(5*rdelta), -se) + 1);
    xi[i] = (rho_ein[i]*f_trans[i] + rho_outer[i])/rhom - 1;
  }
  free(rho_ein);
  free(f_trans);
  free(rho_outer);
  return 0;
}

//////////////////////////////
//////Appendix version 1//////
//////////////////////////////
int calc_xi_DK_app1(double*r, int Nr, double M, double rhos, double conc, double be, double se, double alpha, double beta, double gamma, int delta, double*k, double*P, int Nk, double om, double bias, double*xi_mm, double*xi){
  double rhom = rhomconst*om; //SM h^2/Mpc^3
  //Compute r200m
  double rdelta = pow(M/(1.33333333333*M_PI*rhom*delta), 0.33333333333);
  //double rs = rdelta / conc; //compute scale radius from concentration
  xi[0] = 0;
  double*rho_ein = malloc(Nr*sizeof(double));
  double*f_trans = malloc(Nr*sizeof(double));
  double*rho_outer = malloc(Nr*sizeof(double));
  int i;
  double nu = nu_at_M(M, k, P, Nk, om);
  if (alpha < 0){ //means it wasn't passed in
    alpha = 0.155 + 0.0095*nu*nu;
  }
  if (beta < 0){ //means it wasn't passed in
    beta = 4;
  }
  if (gamma < 0){ //means it wasn't passed in
    gamma = 8;
  }
  if (rhos < 0){ //means it wasn't passed in
    rhos = rhos_einasto_at_M(M, conc, alpha, delta, om);
  }
  double g_b = gamma/beta;
  double r_t = (1.9-0.18*nu)*rdelta; //NEED nu for this
  calc_xi_einasto(r, Nr, M, rhos, conc, alpha, delta, om, rho_ein);
  //here convert xi_ein to rho_ein
  for(i = 0; i < Nr; i++){
    rho_ein[i] = rhom*(1+rho_ein[i]); //rho_ein had xi_ein in it
    f_trans[i] = pow(1+pow(r[i]/r_t,beta), -g_b);
    rho_outer[i] = rhom*(be*pow(r[i]/(5*rdelta), -se) * bias * xi_mm[i] + 1);
    xi[i] = (rho_ein[i]*f_trans[i] + rho_outer[i])/rhom - 1;
  }
  free(rho_ein);
  free(f_trans);
  free(rho_outer);
  return 0;
}

//////////////////////////////
//////Appendix version 2//////
//////////////////////////////
int calc_xi_DK_app2(double*r, int Nr, double M, double rhos, double conc, double be, double se, double alpha, double beta, double gamma, int delta, double*k, double*P, int Nk, double om, double bias, double*xi_mm, double*xi){
  double rhom = rhomconst*om; //SM h^2/Mpc^3
  //Compute r200m
  double rdelta = pow(M/(1.33333333333*M_PI*rhom*delta), 0.33333333333);
  //double rs = rdelta / conc; //compute scale radius from concentration
  xi[0] = 0;
  double*rho_ein = malloc(Nr*sizeof(double));
  double*f_trans = malloc(Nr*sizeof(double));
  double*rho_outer = malloc(Nr*sizeof(double));
  int i;
  double nu = nu_at_M(M, k, P, Nk, om);
  if (alpha < 0){ //means it wasn't passed in
    alpha = 0.155 + 0.0095*nu*nu;
  }
  if (beta < 0){ //means it wasn't passed in
    beta = 4;
  }
  if (gamma < 0){ //means it wasn't passed in
    gamma = 8;
  }
  if (rhos < 0){ //means it wasn't passed in
    rhos = rhos_einasto_at_M(M, conc, alpha, delta, om);
  }
  double g_b = gamma/beta;
  double r_t = (1.9-0.18*nu)*rdelta; //NEED nu for this
  calc_xi_einasto(r, Nr, M, rhos, conc, alpha, delta, om, rho_ein);
  //here convert xi_ein to rho_ein
  for(i = 0; i < Nr; i++){
    rho_ein[i] = rhom*(1+rho_ein[i]); //rho_ein had xi_ein in it
    f_trans[i] = pow(1+pow(r[i]/r_t,beta), -g_b);
    rho_outer[i] = rhom*((1+be*pow(r[i]/(5*rdelta), -se))*bias*xi_mm[i] + 1);
    xi[i] = (rho_ein[i]*f_trans[i] + rho_outer[i])/rhom - 1;
  }
  free(rho_ein);
  free(f_trans);
  free(rho_outer);
  return 0;
}
