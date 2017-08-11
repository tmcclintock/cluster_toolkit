double dndM_at_M(double M, double*k, double*P, int Nk, double om, double d, double e, double f, double g);
int dndM_at_M_arr(double*M, int NM, double*k, double*P, int Nk, double om, double d, double e, double f, double g, double*dndM);

double N_in_bin(double*M, double*dndM, int NM, double volume, double Mlo, double Mhi);
int N_in_bins(double*M, double*dndM, int NM, double volume, double*bins, int Nedges, double*N);
