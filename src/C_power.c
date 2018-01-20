//////////////////////////////////////////////////////
//Some routines used to evaluate power spectra////////
//////////////////////////////////////////////////////

#include "gsl/gsl_spline.h"
#include <math.h>

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
