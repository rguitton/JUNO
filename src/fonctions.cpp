#include "fonctions.h"
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
using namespace std;

double mass_neutron=939.565;//MeV
double mass_proton=938.272;//MeV
double mass_positron=0.510999;//MeV

double fission_energy_U_235=202.36;//MeV
double fission_energy_U_238=205.99;//MeV
double fission_energy_PU_239=211.12;//MeV
double fission_energy_PU_241=214.26;//MeV

double fission_fraction_U_235=0.58;//MeV
double fission_fraction_U_238=0.07;//MeV
double fission_fraction_PU_239=0.30;//MeV
double fission_fraction_PU_241=0.05;//MeV



double flux_U_235(double E){
    double phi;
    phi=exp(0.870-0.160*E-0.091*pow(E,2));
    return phi;
}

double flux_U_238(double E){
    double phi;
    phi=exp(0.896-0.239*E-0.0981*pow(E,2));
    return phi;
}

double flux_PU_239(double E){
    double phi;
    phi=exp(0.976-0.162*E-0.0790*pow(E,2));
    return phi;
}

double flux_PU_241(double E){//renommer en energy 
    double phi;
    phi=exp(0.793-0.080*E-0.1085*pow(E,2));
    return phi;
}

double flux(double E){
    double phi;
    phi=fission_fraction_U_235*flux_U_235(E)+fission_fraction_U_238*flux_U_238(E)+fission_fraction_PU_239*flux_PU_239(E)+fission_fraction_PU_241*flux_PU_241(E);
    return phi*1e-4;//conversion de m^-2 en cm^-2
}

double sigma (double E){
    double energy_positron=E-(mass_neutron-mass_proton);//Mev
    double moment_positron=sqrt(pow(energy_positron,2)-pow(mass_positron,2));
    return 0.0952e-42*energy_positron*moment_positron;//e-42
}

double probability(double E, bool f, double sin2_teta_13, double sin2_teta_12, double delta2_m21, double delta2_m32, double delta2_m31 ){
    double L=53e3;
    double teta_12_bis=asin(sqrt(sin2_teta_12));
    double teta_13_bis=asin(sqrt(sin2_teta_13));

    double P21=pow(cos(teta_13_bis),4)*pow(sin(2*teta_12_bis),2)*pow(sin(1.27*delta2_m21*L/E),2);
    if(f==0){ return 1-P21;}
    else {
            double P31=pow(cos(teta_12_bis),2)*0.1*pow(sin(1.27*delta2_m31*L/E),2);
            double P32=pow(sin(teta_12_bis),2)*0.1*pow(sin(1.27*delta2_m32*L/E),2);
            return 1-P21-P31-P32;
    }
    return 0;
}

double energy_positron(double E){
    double energy_positron=E-(mass_neutron-mass_proton);
    return energy_positron;
}


//use power in MeV
double total_reactor_flux(double tab_flux, double reactor_power){
    double energy_fraction=fission_energy_U_235*fission_fraction_U_235
                            +fission_energy_U_238*fission_fraction_U_238
                            +fission_energy_PU_239*fission_fraction_PU_239
                            +fission_energy_PU_241*fission_fraction_PU_241;//denominator
    double total_flux=reactor_power/(energy_fraction*1.602e-13)*tab_flux;//conversion Mev en j
    return total_flux;
}

double calcul_spectre(double tab_flux,double energy, double sin2_teta_13, double sin2_teta_12, double delta2_m21, double delta2_m32,double delta2_m31){
    double spectre=tab_flux*sigma(energy)*probability(energy, 1, sin2_teta_13, sin2_teta_12, delta2_m21, delta2_m32, delta2_m31)/(4*M_PI*pow(53e3,2));//on divise par 4*pi*L^2
    return spectre;
}

double standart_dev(double E){
    double dev=0.03*sqrt(E);
    return dev;
}

double gauss_pdf(double E, double Ei)
{   
    double mu=standart_dev(E);
	return 1.0 / (standart_dev(E) * sqrt(2.0 * M_PI)) * exp(-(pow((E-Ei)/standart_dev(E), 2)/2.0));
}



//partie longueur

double probability_lenght(double E, bool f, double sin2_teta_13, double sin2_teta_12, double delta2_m21, double delta2_m32, double delta2_m31, double lenght ){
    //double L=53e3;
    double teta_12_bis=asin(sqrt(sin2_teta_12));
    double teta_13_bis=asin(sqrt(sin2_teta_13));

    double P21=pow(cos(teta_13_bis),4)*pow(sin(2*teta_12_bis),2)*pow(sin(1.27*delta2_m21*lenght/E),2);
    if(f==0){ return 1-P21;}
    else {

            double P31=pow(cos(teta_12_bis),2)*0.1*pow(sin(1.27*delta2_m31*lenght/E),2);
            double P32=pow(sin(teta_12_bis),2)*0.1*pow(sin(1.27*delta2_m32*lenght/E),2);
            return 1-P21-P31-P32;
        }
    return 0;
}

    double calcul_spectre_lenght(double tab_flux,double energy, double sin2_teta_13, double sin2_teta_12, double delta2_m21, double delta2_m32,double delta2_m31, double lenght){
    double spectre=tab_flux*sigma(energy)*probability_lenght(energy, 1, sin2_teta_13, sin2_teta_12, delta2_m21, delta2_m32, delta2_m31, lenght)/(4*M_PI*pow(lenght,2));//on divise par 4*pi*L^2
    return spectre;
}