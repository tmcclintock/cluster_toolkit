int xihm_exclusion_at_r(double*r, int Nr, double M, double c, double rt, double beta, double Ma, double ca, double Mb, double cb, double bias, double*ximm, int delta, double Omega_m, double*xihm);

int xi_1h_at_r_arr(double*r, int Nr, double M, double c,
		   double rt, double beta, int delta, double Omega_m,
		   double*xi_1h);
double xi_1h_at_r(double r, double M, double c,
		  double rt, double beta, int delta, double Omega_m);

int xi_2h_at_r_arr(double*r, int Nr, double bias, double*ximm, double*xi2h);
double xi_2h_at_r(double r, double bias, double ximm);

int theta_erfc_at_r_arr(double*r, int Nr, double rt, double beta, double*theta);
double theta_erfc_at_r(double r, double rt, double beta);

