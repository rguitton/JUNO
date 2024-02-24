#ifndef PLOTTER_H
#define PLOTTER_H

void plotter_spectra(int size, const double* xaxis, const double** all_yaxis, const char* dir);
void plotter_flux(int size, const double* xaxis, const double** all_yaxis, const char* dir);

#endif