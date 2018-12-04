double drho_nfw_dr_at_R(double R, double M, double c, int delta, double Omega_m);
int drho_nfw_dr_at_R_arr(double*R, int NR, double Mass, double conc, int delta, double Omega_m, double*drhodr);

double dxi_mm_dr_at_R(double R, double*k, double*P, int Nk);
int dxi_mm_dr_at_R_arr(double*R, int NR, double*k, double*P, int Nk, double*dxidr);
