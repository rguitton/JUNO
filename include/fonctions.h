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


#endif