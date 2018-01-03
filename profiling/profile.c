#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "C_averaging.h"
#include "C_bias.h"
#include "C_boostfactors.h"
#include "C_concentration.h"
#include "C_deltasigma.h"
#include "C_density.h"
#include "C_massfunction.h"
#include "C_miscentering.h"
#include "C_xi.h"
#include "constants.h"

int read_file(FILE *fp, int N, double *data){
  int i;
  for(i = 0; i < N; i++)
    fscanf(fp,"%lf",&data[i]);
  return 0;
}

int main(){
  printf("Starting the profiling routines.\n");

  double M = 1e14; //Msun/h
  double c = 5;
  double om = 0.3;
  double B0 = 1.0;
  double Rs = 1.0;
  double alpha = -1.0;
  
  int i;
  double rmin = 0.01; //Mpc/h
  double rmax = 200; //Mpc/h
  double Rmin = 0.1; //Mpc/ projectd
  double Rmax = 80; //Mpc/h projected
  int num = 1000;
  double dlr = (log10(rmax)-log10(rmin))/num;
  double dlR = (log10(Rmax)-log10(Rmin))/num;
  double*r = (double*)malloc(num*sizeof(double));
  double*R = (double*)malloc(num*sizeof(double));
  double*xi_nfw = (double*)malloc(num*sizeof(double));
  double*den_nfw = (double*)malloc(num*sizeof(double));
  double*boost_nfw = (double*)malloc(num*sizeof(double));
  double*boost_pl = (double*)malloc(num*sizeof(double));
  double*Sigma_nfw = (double*)malloc(num*sizeof(double));
  
  for(i = 0; i <num; i++){
    r[i] = pow(10, log10(rmin)+i*dlr);
    R[i] = pow(10, log10(Rmin)+i*dlR);
    xi_nfw[i] = xi_nfw_at_R(r[i], M, c, 200, om);
    den_nfw[i] = rho_nfw_at_R(r[i], M, c, 200, om);
    boost_nfw[i] = boost_nfw_at_R(R[i], B0, Rs);
    boost_pl[i] = boost_powerlaw_at_R(R[i], B0, Rs, alpha);
    Sigma_nfw[i] = Sigma_nfw_at_R(R[i], M, c, 200, om);
  }
  printf("First arrays calculated.\n");

  printf("Reading in power spectra.\n");
  FILE*k_fp, *P_fp;
  char*line = NULL;
  size_t len = 0;
  int Nk = -1, Nklin = -1;
  int read;
  //Get the nl data
  k_fp = fopen("input_files/knl.txt","r");
  while ((read = getline(&line, &len, k_fp)) != -1){Nk++;}
  rewind(k_fp);
  getline(&line,&len,k_fp); //header line read off
  double*knl = (double*)malloc(Nk*sizeof(double));
  read_file(k_fp, Nk, knl);
  fclose(k_fp);
  P_fp = fopen("input_files/pnl.txt","r");
  getline(&line,&len,P_fp); //header line read off
  double*Pnl = (double*)malloc(Nk*sizeof(double));
  read_file(P_fp, Nk, Pnl);
  fclose(P_fp);
  //Get the lin data
  k_fp = fopen("input_files/klin.txt","r");
  while ((read = getline(&line, &len, k_fp)) != -1){Nklin++;}
  rewind(k_fp);
  getline(&line,&len,k_fp); //header line read off
  double*klin = (double*)malloc(Nklin*sizeof(double));
  read_file(k_fp, Nklin, klin);
  fclose(k_fp);
  P_fp = fopen("input_files/plin.txt","r");
  getline(&line, &len, P_fp); //header line read off
  double*Plin = (double*)malloc(Nklin*sizeof(double));
  read_file(P_fp, Nklin, Plin);
  fclose(P_fp);
  printf("Power spectra read in. Nlin=%d Nnl=%d\n",Nk, Nklin);

  printf("Freeing arrays.\n");
  free(r);
  free(R);
  free(xi_nfw);
  free(den_nfw);
  free(boost_nfw);
  free(boost_pl);
  free(Sigma_nfw);
  free(knl);
  free(klin);
  free(Pnl);
  free(Plin);
  return 0;
}
