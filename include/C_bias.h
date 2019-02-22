int bias_at_nu_arr(double*nu, int Nnu, int delta, double*bias);
int bias_at_R_arr(double*R, int NR, int delta, double*k, double*P, int Nk, double*bias);
int bias_at_M_arr(double*M, int NM, int delta, double*k, double*P, int Nk, double om, double*bias);

int bias_at_nu_arr_FREEPARAMS(double*nu, int Nnu, int delta, double A, double a, double B, double b, double C, double c, double*bias);
