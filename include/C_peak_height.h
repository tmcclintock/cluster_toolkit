double M_to_R(double M, double Omega_m);
double R_to_M(double R, double Omega_m);

double sigma2_at_R(double R, double*k, double*P, int Nk);
double sigma2_at_M(double M, double*k, double*P, int Nk, double om);
int sigma2_at_R_arr(double*R, int NR,  double*k, double*P, int Nk, double*s2);
int sigma2_at_M_arr(double*M, int NM,  double*k, double*P, int Nk, double om, double*s2);

int dsigma2dR_at_R_arr(double*R, int NR, double*k, double*P, int Nk, double*ds2dR);
double dsigma2dR_at_R(double R, double*k, double*P, int Nk);
int dsigma2dM_at_M_arr(double*M, int NM, double*k, double*P, int Nk, double Omega_m, double*ds2dM);
double dsigma2dM_at_M(double M, double*k, double*P, int Nk, double Omega_m);

double nu_at_R(double R, double*k, double*P, int Nk);
double nu_at_M(double M, double*k, double*P, int Nk, double om);
int nu_at_R_arr(double*R, int NR, double*k, double*P, int Nk, double*nu);
int nu_at_M_arr(double*M, int NM, double*k, double*P, int Nk, double om, double*nu);
