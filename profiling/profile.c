#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "C_averaging.h"
#include "C_bias.h"
#include "C_boostfactors.h"
#include "C_concentration.h"
#include "C_deltasigma.h"
#include "C_density.h"
#include "C_massfunction.h"
#include "C_miscentering.h"
#include "C_xi.h"
#include "constants.h"


int main(){
  printf("Starting the profiling routines.\n");

  double M = 1e14; //Msun/h
  double c = 5;
  double om = 0.3;
  double B0 = 1.0;
  double Rs = 1.0;
  double alpha = -1.0;
  
  int i;
  double Rmin = 0.01; //Mpc/h
  double Rmax = 150; //Mpc/h
  int num = 100;
  double dlR = (log10(Rmax)-log10(Rmin))/num;
  double*R = (double*)malloc(num*sizeof(double));
  double*xi_nfw = (double*)malloc(num*sizeof(double));
  double*den_nfw = (double*)malloc(num*sizeof(double));
  double*boost_nfw = (double*)malloc(num*sizeof(double));
  double*boost_pl = (double*)malloc(num*sizeof(double));
  
  for(i = 0; i <num; i++){
    R[i] = pow(10, log10(Rmin)+i*dlR);
    xi_nfw[i] = xi_nfw_at_R(R[i], M, c, 200, om);
    den_nfw[i] = rho_nfw_at_R(R[i], M, c, 200, om);
    boost_nfw[i] = boost_nfw_at_R(R[i], B0, Rs);
    boost_pl[i] = boost_powerlaw_at_R(R[i], B0, Rs, alpha);
  }

  free(R);
  free(xi_nfw);
  free(den_nfw);
  free(boost_nfw);
  free(boost_pl);
  return 0;
}
