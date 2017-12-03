double G_sigma(double sigma, double d, double e, double f, double g);
int G_sigma_arr(double*sigma, int Ns, double d, double e, double f, double g, double*G);

double dndM_at_M(double M, double*k, double*P, int Nk, double om, double d, double e, double f, double g);
int dndM_at_M_arr(double*M, int NM, double*k, double*P, int Nk, double om, double d, double e, double f, double g, double*dndM);

double n_in_bin(double Mlo, double Mhi, double*M, double*dndM, int NM);
int n_in_bins(double*edges, int Nedges, double*M, double*dndM, int NM, double*N);
