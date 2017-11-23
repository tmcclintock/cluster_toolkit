#include <stdio.h>
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

  printf("%e\n",boost_nfw_at_R(1.0, 1.0, 1.0));
  return 0;
}
