double Sigma_mis_single_at_R(double R, double*Rs, double*Sigma, int Ns, double M, double conc, int delta, double om, double Rmis);
int Sigma_mis_single_at_R_arr(double*R, int NR, double*Rs, double*Sigma, int Ns, double M, double conc, int delta, double om, double Rmis, double*Sigma_mis);

double Sigma_mis_at_R(double R, double*Rs, double*Sigma, int Ns, double M, double conc, int delta, double om, double Rmis, int integrand_switch);
int Sigma_mis_at_R_arr(double*R, int NR, double*Rs, double*Sigma, int Ns, double M, double conc, int delta, double om, double Rmis, int integrand_switch, double*Sigma_mis);

double DeltaSigma_mis_at_R(double R, double*Rs, double*Sigma, int Ns);
int DeltaSigma_mis_at_R_arr(double*R, int NR, double*Rs, double*Sigma, int Ns, double*DeltaSigma_mis);
