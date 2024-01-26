#include "fonctions.h"
#include <math.h>
using namespace std;


double flux(double E){
    double phi;
    phi=0.58*exp(0.870-0.160*E-0.091*pow(E,2))
    +0.30*exp(0.896-0.239*E-0.0981*pow(E,2))
    +0.07*exp(0.976-0.162*E-0.0790*pow(E,2))
    +0.05*exp(0.793-0.080*E-0.1085*pow(E,2));
    return phi;
}


double sigma (){
    return 0.0952*pow(10,-42);
}
/*
double sin2_teta_12=0.32;
double sin2_teta_13=
double sin2_teta_32=

double delta_21=
double delta_31=
double delta_32=
*/

double probability(){
    return 1;
}