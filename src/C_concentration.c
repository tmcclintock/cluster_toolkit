#include "C_concentration.h"
#include "C_bias.h"

#include <math.h>
#include <stdio.h>

double DK15_concentration_at_M(double Mass, double*k, double*P, int Nk, double Omega_m){
  double nu = nu_at_M(Mass, k, P, Nk, Omega_m);
  return 0;
}
