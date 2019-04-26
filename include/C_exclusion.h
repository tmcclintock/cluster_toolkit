int xi_hm_exclusion_at_r_arr(double*r, int Nr,
			    double M, double c, double alpha,
			    double rt, double beta,
			    double r_eff, double beta_eff,
			    double r_A, double r_B, double beta_ex,
			    double bias, double*ximm, int delta,
			     double Omega_m, double*xihm);

int xi_1h_at_r_arr(double*r, int Nr, double M, double c, double alpha,
		   double rt, double beta, int delta, double Omega_m,
		   double*xi_1h);

int xi_2h_at_r_arr(double*r, int Nr, double r_eff, double beta_eff,
		   double bias, double*ximm, double*xi2h);

int xi_C_at_r_arr(double*r, int Nr, double r_A, double r_B, double beta_2,
		  double*xi_2h, double*xi_C);

int theta_erfc_at_r_arr(double*r, int Nr, double rt, double beta,
			double*theta);

double r_exclusion(double r1, double r2, int scheme);
