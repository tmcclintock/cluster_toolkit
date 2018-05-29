/** @file C_bias.c
 *  @brief Halo bias functions.
 * 
 *  These functions are the halo bias for different
 *  types of inputs including halo mass, halo radius,
 *  and halo peak height.
 *  
 *  @author Tom McClintock (tmcclintock)
 *  @bug No known bugs.
 */

#include "C_bias.h"
#include "C_peak_height.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define delta_c 1.686 //Critical collapse density
#define rhocrit 2.77533742639e+11
//1e4*3.*Mpcperkm*Mpcperkm/(8.*PI*G); units are Msun h^2/Mpc^3

///////////BIAS FUNCTIONS///////////

double bias_at_nu(double nu, int delta){
  double*nus = (double*)malloc(sizeof(double));
  double*bias = (double*)malloc(sizeof(double));
  nus[0] = nu;
  bias_at_nu_arr(nus, 1, delta, bias);
  double result = bias[0];
  free(nus);
  free(bias);
  return result;
}

double bias_at_R(double R, int delta, double*k, double*P, int Nk){
  double*Rs = (double*)malloc(sizeof(double));
  double*bias = (double*)malloc(sizeof(double));
  Rs[0] = R;
  bias_at_R_arr(Rs, 1, delta, k, P, Nk, bias);
  double result = bias[0];
  free(Rs);
  free(bias);
  return result;
}

double bias_at_M(double M, int delta, double*k, double*P, int Nk, double om){
  return bias_at_nu(nu_at_M(M, k, P, Nk, om), delta);
}

int bias_at_nu_arr(double*nu, int Nnu, int delta, double*bias){
  double y = log10(delta);
  double xp = exp(-1.0*pow(4./y,4.));
  double A = 1.+0.24*y*xp, a = 0.44*y-0.88;
  double B = 0.183, b = 1.5;
  double C = 0.019+0.107*y+0.19*xp, c = 2.4;
  int i;
  for(i = 0; i < Nnu; i++)
    bias[i] = 1 - A*pow(nu[i],a)/(pow(nu[i],a)+pow(delta_c,a)) + B*pow(nu[i],b) + C*pow(nu[i],c);
  return 0;
}

int bias_at_R_arr(double*R, int NR, int delta, double*k, double*P, int Nk, double*bias){
  double*nu = (double*)malloc(sizeof(double)*NR);
  nu_at_R_arr(R, NR, k, P, Nk, nu);
  bias_at_nu_arr(nu, NR, delta, bias);
  free(nu);
  return 0;
}

int bias_at_M_arr(double*M, int NM, int delta, double*k, double*P, int Nk, double om, double*bias){
  double*nu = (double*)malloc(sizeof(double)*NM);
  nu_at_M_arr(M, NM, k, P, Nk, om, nu);
  bias_at_nu_arr(nu, NM, delta, bias);
  free(nu);
  return 0;
}

////EXTRA FUNCTION FOR EMULATOR DEVELOPMENT////
int bias_at_nu_arr_FREEPARAMS(double*nu, int Nnu, int delta, double A, double a, double B, double b, double C, double c, double*bias){
  int i;
  for(i = 0; i < Nnu; i++)
    bias[i] = 1 - A*pow(nu[i],a)/(pow(nu[i],a)+pow(delta_c,a)) + B*pow(nu[i],b) + C*pow(nu[i],c);
  return 0;
}
