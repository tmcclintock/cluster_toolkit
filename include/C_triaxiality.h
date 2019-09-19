double mapped_radii(double z, double R, double sin_i, double cos_i,
		    double sin_phi, double cos_phi, double q, double s);

double Ellipsoidal_Sigma_nfw_at_R(double R, double M, double c, double i,
				  double q, double s,
				  double delta, double Omega_m);

int Ellipsoidal_Sigma_nfw_at_R_arr(double*R, int NR, double M, double c,
				   double i, double q, double s,
				   double delta, double Omega_m, double*Sigma);
