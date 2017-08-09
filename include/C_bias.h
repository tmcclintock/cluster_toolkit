double nu_at_R(double R, double*k, double*P, int Nk);
double nu_at_M(double M, double*k, double*P, int Nk, double om);
double bias_at_nu(double nu, int delta);
double bias_at_R(double R, int delta, double*k, double*P, int Nk);
double bias_at_M(double M, int delta, double*k, double*P, int Nk, double om);
double bias_at_nu_arr(double*nu, int Nnu, int delta, double*bias);
double bias_at_R_arr(double*R, int NR, int delta, double*k, double*P, int Nk, double*bias);
double bias_at_M_arr(double*M, int NM, int delta, double*k, double*P, int Nk, double*bias, double om);

