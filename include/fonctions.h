#ifndef FONCTIONS_H
#define FONCTIONS_H

double flux(double E);
double flux_U_235(double E);
double flux_U_238(double E);
double flux_PU_239(double E);
double flux_PU_241(double E);

double probability(double L, char A, bool f);
double sigma (double E);
double visible_energy(double E);
double energy_positron(double E);
double integrale_spectre(double Emin, double Emax,double n);
double total_reactor_flux(double flux, double power);
double calcul_spectre(double tab_flux,double energy);
double standart_dev(double E);
double gauss_pdf(double E, double mu);
float* convolve(float E1[], float E2[], int lenE1, int lenE2, int* lenE3);
double* convert_Evis(double E[], int size);

#endif