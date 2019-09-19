/** @file C_triaxiality.c
 *  @brief Triaxiality and orientation angle affects on surface mass density profiles.
 *  
 *  These functions compute projected surface mass density profiles
 *  for galaxy clusters (halos) that can be ellipsoidal.
 *  
 *  @author Tom McClintock (tmcclintock)
 *  @bug No known bugs.
 */

#include "C_xi.h"

#include "gsl/gsl_integration.h"
#include "gsl/gsl_spline.h"
#include "gsl/gsl_errno.h"

#include <math.h>
#include <stdio.h>

#define ABSERR 0.0
#define RELERR 1e-4
#define workspace_size 8000
#define ulim 5.0
#define rhocrit 2.77533742639e+11
//1e4*3.*Mpcperkm*Mpcperkm/(8.*PI*G); units are SM h^2/Mpc^3
#define KEY 3 //Used for GSL QAG function


/**
 * \brief Radial vector norm in the frame of the triaxial halo,
 * mapped to the isodensity contour of a spherical halo
 * at r = \sqrt(z^2 + R^2).
 * Units are Msun/h comoving.
 */
double mapped_radii(double z, double R, double sin_i, double cos_i,
		    double sin_phi, double cos_phi, double q, double s){
  return sqrt(R * R * sin_phi * sin_phi/(q * q) +
	      pow(R * cos_phi * cos_i - z * sin_i, 2)/(s * s) +
	      pow(R * cos_phi * sin_i + z * cos_i, 2));
}
