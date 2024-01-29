#include "fonctions.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main()
{
    int size=100;
    double tableau_lenght[size];
    double tableau_spectra[size];
    double L=60;
    double energy;


    for (int i=0;i<size;i++){
        tableau_lenght[i]=double(i)*35/size;
        energy=L/tableau_lenght[i];
        printf("Voici l'élément du tableau %f \n",tableau_lenght[i]);
        printf("L'énergie vaut donc %f Mev \n",energy);
        
        printf("Le flux correspondant vaut  %f \n",flux(energy));
        
        tableau_spectra[i]=flux(energy)*sigma(energy);
        printf("Le spectre correspondant vaut  %f \n",tableau_spectra[i]);
    }

        FILE *gnuplotPipe = popen("gnuplot -persist", "w");
    if (gnuplotPipe == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture de Gnuplot.\n");
        return 1;
    }

    fprintf(gnuplotPipe, "plot '-' with linespoints title 'Graphique'\n");
    for (int i = 0; i < size; i++) {
        fprintf(gnuplotPipe, "%lf %lf\n", tableau_lenght[i], tableau_spectra[i]);
    }
    fprintf(gnuplotPipe, "e\n");
    fflush(gnuplotPipe);

    // Attente de l'utilisateur avant de fermer la fenêtre Gnuplot
    printf("Appuyez sur Entrée pour fermer le graphique...\n");
    getchar();

    // Fermeture du pipe Gnuplot
    fclose(gnuplotPipe);
    
    return 0;
}