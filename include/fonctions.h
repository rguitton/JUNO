#ifndef FONCTIONS_H
#define FONCTIONS_H

double flux(double E);
double flux_U_235(double E);
double flux_U_238(double E);
double flux_PU_239(double E);
double flux_PU_241(double E);

double probability(double L,  bool f,double sin2_teta_13, double sin2_teta_12, double delta2_m21, double delta2_m32,double delta2_m31);
double sigma (double E);
double visible_energy(double E);
double energy_positron(double E);
double total_reactor_flux(double flux, double power);
double calcul_spectre(double tab_flux,double energy,double sin2_teta_13, double sin2_teta_12, double delta2_m21, double delta2_m32, double delta2_m31);
double standart_dev(double E);
double gauss_pdf(double E, double Ei);
double probability_lenght(double E,bool f, double sin2_teta_13, double sin2_teta_12, double delta2_m21, double delta2_m32, double delta2_m31, double lenght );
//float* convolve(float E1[], float E2[], int lenE1, int lenE2, int* lenE3);
//double* convert_Evis(double E[], int size);
double calcul_spectre_lenght(double tab_flux,double energy, double sin2_teta_13, double sin2_teta_12, double delta2_m21, double delta2_m32,double delta2_m31, double lenght);

#endif