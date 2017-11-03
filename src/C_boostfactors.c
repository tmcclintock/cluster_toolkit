#include "C_boostfactors.h"

#include <math.h>
#include <stdio.h>

double boost_nfw_at_R(double R, double B0, double Rs){
  double x2m1 = R*R/(Rs*Rs)-1; //x^2 -1
  double sqx2m1;
  if (x2m1 > 0){ //x > 1
    sqx2m1 = sqrt(x2m1); //sqrt(x*x-1)
    return 1. + B0/x2m1 * (1-atan(sqx2m1)/sqx2m1);
  } else if (x2m1 < 0){ // x < 1
    sqx2m1 = sqrt(-x2m1); //sqrt(1-x*x)
    return 1. + B0/x2m1 * (1-atanh(sqx2m1)/sqx2m1);
  }  else{ // x = 1
    return 1;
  }
}

int boost_factor_nfw_at_R_arr(double*R, int NR, double B0, double Rs, double*boost){
  int i;
  for(i = 0; i < NR; i++){
    boost[i] = boost_nfw_at_R(R[i], B0, Rs);
  }
  return 0;
}
