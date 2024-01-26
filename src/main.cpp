#include "fonctions.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

int main()
{
    int size=100;
    double tableau_energy[size];
    double tableau_spectra[size];
    
    for (int i=0;i<size;i++){
        tableau_energy[i]=double(i)*35/size;
        printf("Voici l'élément du tableau %f \n",tableau_energy[i]);
        tableau_spectra[i]=flux(tableau_energy[i])*sigma()*probability();
        printf("Le spectre correspondant vaut  %f \n",tableau_spectra[i]);
    }

    return 0;
}