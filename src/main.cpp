#include "fonctions.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main()
{
    int size=1000;
    double tableau_lenght[size];
    double spectra_initial[size];
    double spectra_p21[size];
    double spectra_NH[size];
    double spectra_IH[size];

    double L=60;
    double energy;


    for (int i=0;i<size;i++){
        tableau_lenght[i]=double(i)*33/size;
        energy=L/tableau_lenght[i];
        //printf("Voici l'élément du tableau %f \n",tableau_lenght[i]);
        //printf("L'énergie vaut donc %f Mev \n",energy);
        
        //printf("Le flux correspondant vaut  %f \n",flux(energy));
        
        spectra_initial[i]=flux(energy)*sigma(energy);
        spectra_p21[i]=flux(energy)*sigma(energy)*probability(energy, 'I', 0);
        spectra_NH[i]=flux(energy)*sigma(energy)*probability(energy, 'N', 1);
        spectra_IH[i]=flux(energy)*sigma(energy)*probability(energy, 'I', 1);

        //printf("Le spectre correspondant vaut  %f \n",tableau_spectra[i]);
    }

   FILE *gnuplotPipe = popen("gnuplot -persist", "w");
    if (gnuplotPipe == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture de Gnuplot.\n");
        return 1;
    }

    fprintf(gnuplotPipe, "plot '-' with linespoints title 'Initial', '-' with linespoints title 'P21', '-' with linespoints title 'NH', '-' with linespoints title 'IH'\n");

    // Envoi des données de chaque tableau
    for (int i = 0; i < size; i++) {
        fprintf(gnuplotPipe, "%lf %lf\n", tableau_lenght[i], spectra_initial[i]);
    }
    fprintf(gnuplotPipe, "e\n");

    for (int i = 0; i < size; i++) {
        fprintf(gnuplotPipe, "%lf %lf\n", tableau_lenght[i], spectra_p21[i]);
    }
    fprintf(gnuplotPipe, "e\n");

    for (int i = 0; i < size; i++) {
        fprintf(gnuplotPipe, "%lf %lf\n", tableau_lenght[i], spectra_NH[i]);
    }
    fprintf(gnuplotPipe, "e\n");

    for (int i = 0; i < size; i++) {
        fprintf(gnuplotPipe, "%lf %lf\n", tableau_lenght[i], spectra_IH[i]);
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