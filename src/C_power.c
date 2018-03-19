/** @file C_power.c
 *  @brief Function for integrands that require P(x/R).
 *
 * This file contains a function to access the power
 * spectrum P(k) in integrands that actually evaluate
 * P(x/R) where x=kR.
 * 
 * @author Tom McClintock
 */

#include "gsl/gsl_spline.h"
#include <math.h>

/**
 * \brief Evaluate the power spectrum P(x/R).
 *
 * Certain integrals require evaluating the power spectrum at P(k=x/R).
 * This function takes in a spline for a pre-computed power spectrum P(k)
 * and returns the evaluation of P(k) either within the region that
 * the spline is valid or with power law approximations at smaller or
 * larger scales.
 */
double get_P(double x, double R, double*k, double*P, int Nk, void*void_Pspl, void*void_acc){
  gsl_spline*Pspl = (gsl_spline*)void_Pspl;
  gsl_interp_accel*acc = (gsl_interp_accel*)void_acc;
  double ki = x/R;
  double kmin = k[0];
  double kmax = k[Nk-1];
  double alpha,A;
  if (ki < kmin){
    alpha = log(P[1]/P[0])/log(k[1]/k[0]);
    A = P[0]/pow(k[0],alpha);
    return A*pow(ki,alpha);
  }else if (ki > kmax){
    alpha = log(P[Nk-1]/P[Nk-2])/log(k[Nk-1]/k[Nk-2]);
    A = P[Nk-1]/pow(k[Nk-1],alpha);
    return A*pow(ki,alpha);
  }// Assume power laws at ends
  return gsl_spline_eval(Pspl,ki,acc);
}
