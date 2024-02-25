#ifndef PLOTTER_H
#define PLOTTER_H

#include <TGraphErrors.h>

// Fill ignoring nan values
void fillGraphIgnoringNaN(TGraphErrors &graph, const double *x, const double *y, const double *ex, const double *ey, int nPoints);

//Figure 2.4, article 1
void plotter_spectra(int size, const double* xaxis, const double** all_yaxis, const char* dir);

//Figure 2.6, article 1
void plotter_flux(int size, const double* xaxis, const double** all_yaxis, const char* dir);

//Figure 4, article 2: 
void plotter_energy_spectrum(int size, const double* xaxis, const double* yaxis,const char* dir);

void plotter_visible_energy_spectrum(int size, const double* xaxis, const double* yaxis,const char* dir);

void plotter_visible_energy_spectrum_six_years(int size, const double* xaxis, const double* yaxis,const char* dir);

void plotter_chi2(int nPoint, int nx, int ny, const double** all_xaxis, const double** all_yaxis, const char* dir);

void plotter_chi2_MH(int nPoint, const double** all_xaxis, const double** all_yaxis, const char* dir);

void plotter_ch2_Distance(int size, const double* xaxis, const double* yaxis,const char* dir);

#endif