double xi_nfw_at_R(double, double, double, int, double);
double xi_einasto_at_R(double, double, double, double, double, int, double);
double xi_mm_at_R(double,double*,double*,int,int,double);
int calc_xi_nfw(double*, int, double, double, int, double, double*);
int calc_xi_einasto(double*, int, double, double, double, double, int, double, double*);

double rhos_einasto_at_M(double Mass, double rs, double alpha, int delta, double om);

int calc_xi_2halo(int, double, double*, double*);
int calc_xi_hm(int, double*, double*, double*, int);
int calc_xi_mm(double*, int, double*, double*, int, double*, int, double);

double xi_mm_at_R_exact(double R, double*k, double*P, int Nk);
int calc_xi_mm_exact(double*R, int NR, double*k, double*P, int Nk, double*xi);

int calc_xi_DK(double*R, int NR, double M, double rhos, double rs, double alpha, double be, double se, double beta, double gamma, int delta, double*k, double*P, int Nk, double om, double*xi);
double xi_DK(double R, double M, double rhos, double rs, double alpha, double be, double se, double beta, double gamma, int delta, double*k, double*P, int Nk, double om);
