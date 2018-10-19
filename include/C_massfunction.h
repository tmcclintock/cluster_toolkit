double G_at_M(double M, double*k, double*P, int Nk, double om, double d, double e, double f, double g);
int G_at_M_arr(double*M, int NM, double*k, double*P, int Nk, double om, double d, double e, double f, double g, double*G);

double G_at_sigma(double sigma, double d, double e, double f, double g);
int G_at_sigma_arr(double*sigma, int Ns, double d, double e, double f, double g, double*G);

double dndM_at_M(double M, double*k, double*P, int Nk, double om, double d, double e, double f, double g);
int dndM_at_M_arr(double*M, int NM, double*k, double*P, int Nk, double om, double d, double e, double f, double g, double*dndM);

double n_in_bin(double Mlo, double Mhi, double*M, double*dndM, int NM);
int n_in_bins(double*edges, int Nedges, double*M, double*dndM, int NM, double*N);

//int dndM_sigma2_precomputed(double*M, double*sigma2, double*sigma2_top, double*sigma2_bot, int NM, double Omega_m, double d, double e, double f, double g, double*dndM);
int dndM_sigma2_precomputed(double*M, double*sigma2, double*dsigma2dM, int NM, double Omega_m, double d, double e, double f, double g, double*dndM);
