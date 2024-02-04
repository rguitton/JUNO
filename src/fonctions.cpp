#include "fonctions.h"
#include <math.h>
using namespace std;

double mass_neutron=939.565;//Mev
double mass_proton=938.272;//Mev
double mass_positron=0.510999;//Mev

//pensez à remettre les coeffs 
double flux_U_235(double E){
    double phi;
    phi=0.58*exp(0.870-0.160*E-0.091*pow(E,2));
    return phi;
}

double flux_U_238(double E){
    double phi;
    phi=0.3*exp(0.896-0.239*E-0.0981*pow(E,2));
    return phi;
}

double flux_PU_239(double E){
    double phi;
    phi=0.07*exp(0.976-0.162*E-0.0790*pow(E,2));
    return phi;
}

double flux_PU_241(double E){
    double phi;
    phi=0.05*exp(0.793-0.080*E-0.1085*pow(E,2));
    return phi;
}

double flux(double E){
    double phi;
    phi=flux_U_235(E)+flux_U_238(E)+flux_PU_239(E)+flux_PU_241(E);
    return phi;
}

double sigma (double E){
    double energy_positron=E-(mass_neutron-mass_proton);
    double moment_positron=sqrt(pow(energy_positron,2)-pow(mass_positron,2));
    return 0.0952e-42*energy_positron*moment_positron;//e-42
}

double sin2_teta_12=0.32;
double sin2_2teta_13=0.1;   

double teta_12=asin(sqrt(sin2_teta_12));
double teta_13=asin(sqrt(sin2_2teta_13))/2;


double delta2_m21=7.6e-5;
double delta2_m32=2.4e-3;

double delta2_m31_NH=delta2_m32+delta2_m21;
double delta2_m31_IH=delta2_m32-delta2_m21;

double probability(double E, char A, bool f){
    double L=53e3;
    double P21=pow(cos(teta_13),4)*pow(sin(2*teta_12),2)*pow(sin(1.27*delta2_m21*L/E),2);
    if(f==0){ return 1-P21;}
    else {
        
        if(A=='N'){
            double P31=pow(cos(teta_12),2)*0.1*pow(sin(1.27*delta2_m31_NH*L/E),2);
            double P32=pow(sin(teta_12),2)*0.1*pow(sin(1.27*delta2_m32*L/E),2);
            return 1-P21-P31-P32;
        }
        else if(A=='I'){
            double P31=pow(cos(teta_12),2)*0.1*pow(sin(1.27*delta2_m31_IH*L/E),2);
            double P32=pow(sin(teta_12),2)*0.1*pow(sin(1.27*(-delta2_m32)*L/E),2);
            return 1-P21-P31-P32;
        } 
    }
    return 0;
}