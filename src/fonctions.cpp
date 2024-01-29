#include "fonctions.h"
#include <math.h>
using namespace std;

double mass_neutron=939.565;//Mev
double mass_proton=938.272;//Mev
double mass_positron=0.510999;

double flux(double E){//energy en Mev 
    double phi;
    phi=0.58*exp(0.870-0.160*E-0.091*pow(E,2))
    +0.30*exp(0.896-0.239*E-0.0981*pow(E,2))
    +0.07*exp(0.976-0.162*E-0.0790*pow(E,2))
    +0.05*exp(0.793-0.080*E-0.1085*pow(E,2));
    return phi;
}


double sigma (double E){
    double energy_positron=E-(mass_neutron-mass_proton);
    double moment_positron=sqrt(pow(energy_positron,2)-pow(mass_positron,2));
    return 0.0952*pow(10,-0)*energy_positron*moment_positron;
}
/*
double sin2_teta_12=0.32;
double sin2_teta_13=
double sin2_teta_32=


double delta_m21=
double delta_m31=
double delta_m32=
*/

double probability(){

    return 1;
}