double Sigma_mis_single_at_R(double R, double*Rs, double*Sigma, int Ns, double M, double conc, int delta, double om, double Rmis);
int Sigma_mis_single_at_R_arr(double*R, int NR, double*Rs, double*Sigma, int Ns, double M, double conc, int delta, double om, double Rmis, double*Sigma_mis);

double Sigma_mis_g2d_at_R(double R, double*Rs, double*Sigma, int Ns, double M, double conc, int delta, double om, double Rmis); //2D Gaussian
int Sigma_mis_g2d_at_R_arr(double*R, int NR, double*Rs, double*Sigma, int Ns, double M, double conc, int delta, double om, double Rmis, double*Sigma_mis); //2D Gaussian

double DeltaSigma_mis_at_R(double R, double*Rs, double*Sigma, int Ns);
int DeltaSigma_mis_at_R_arr(double*R, int NR, double*Rs, double*Sigma, int Ns, double*DeltaSigma_mis);
