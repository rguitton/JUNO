#include "fonctions.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main()
{
    int size=1000;
    double tableau_lenght[size];
    double tableau_energy[size];
    
    double tab_cross[size];
    
    double tab_flux_cross_U_235[size];
    double tab_flux_cross_U_238[size];
    double tab_flux_cross_PU_239[size];
    double tab_flux_cross_PU_241[size];

    double tab_flux_U_235[size];
    double tab_flux_U_238[size];
    double tab_flux_PU_239[size];
    double tab_flux_PU_241[size];


    double spectra_initial[size];
    double spectra_p21[size];
    double spectra_NH[size];
    double spectra_IH[size];

    double L=53;//à changer 
    double energy;


    for (int i=0;i<size;i++){
        tableau_lenght[i]=5+double(i)*28/size;
        tableau_energy[i]=1+double(i)*10/size;
    }

     for (int i=0;i<size;i++){
        double energy=L/tableau_lenght[i];
        spectra_initial[i]=flux(energy)*sigma(energy);
        spectra_p21[i]=flux(energy)*sigma(energy)*probability(energy, 'I', 0);
        spectra_NH[i]=flux(energy)*sigma(energy)*probability(energy, 'N', 1);
        spectra_IH[i]=flux(energy)*sigma(energy)*probability(energy, 'I', 1);

        tab_flux_U_235[i]=flux_U_235(tableau_energy[i]);
        tab_flux_U_238[i]=flux_U_238(tableau_energy[i]);
        tab_flux_PU_239[i]=flux_PU_239(tableau_energy[i]);
        tab_flux_PU_241[i]=flux_PU_241(tableau_energy[i]);

        tab_cross[i]=sigma(tableau_energy[i]);

        tab_flux_cross_U_235[i]=tab_cross[i]*tab_flux_U_235[i]*pow(10,2);
        tab_flux_cross_U_238[i]=tab_cross[i]*tab_flux_U_238[i]*pow(10,2);
        tab_flux_cross_PU_239[i]=tab_cross[i]*tab_flux_PU_239[i]*pow(10,2);
        tab_flux_cross_PU_241[i]=tab_cross[i]*tab_flux_PU_241[i]*pow(10,2);

    }

   FILE *gnuplotPipe1 = popen("gnuplot -persist", "w");
    if (gnuplotPipe1 == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture de Gnuplot.\n");
        return 1;
    }

    fprintf(gnuplotPipe1, "plot '-' with linespoints title 'Initial', '-' with linespoints title 'P21', '-' with linespoints title 'NH', '-' with linespoints title 'IH'\n");

    for (int i = 0; i < size; i++) {
        fprintf(gnuplotPipe1, "%g %g\n", tableau_lenght[i], spectra_initial[i]);
    }
    fprintf(gnuplotPipe1, "e\n");

    for (int i = 0; i < size; i++) {
        fprintf(gnuplotPipe1, "%g %g\n", tableau_lenght[i], spectra_p21[i]);
    }
    fprintf(gnuplotPipe1, "e\n");

    for (int i = 0; i < size; i++) {
        fprintf(gnuplotPipe1, "%g %g\n", tableau_lenght[i], spectra_NH[i]);
    }
    fprintf(gnuplotPipe1, "e\n");

    for (int i = 0; i < size; i++) {
        fprintf(gnuplotPipe1, "%g %g\n", tableau_lenght[i], spectra_IH[i]);
    }
    fprintf(gnuplotPipe1, "e\n");
        
    fflush(gnuplotPipe1);

     FILE *gnuplotPipe2 = popen("gnuplot -persist", "w");
    if (gnuplotPipe2 == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture de Gnuplot.\n");
        return 1;
    }

    fprintf(gnuplotPipe2, "plot '-' with linespoints lt 3 title 'flux U 235', '-' with linespoints lt 4 title 'flux U 238', '-' with linespoints lt 1 title 'flux PU 239', '-' with linespoints lt 2 title 'flux PU 241', '-' with linespoints title 'cross section','-' with linespoints lt 3 title 'U235 product','-' with linespoints lt 4 title 'U238 product','-' with linespoints lt 1 title 'PU239 product','-' with linespoints lt 2 title 'PU241 product'\n");
    for (int i = 0; i < size; i++) {
        fprintf(gnuplotPipe2, "%g %g\n",tableau_energy[i], tab_flux_U_235[i]);
    }
    fprintf(gnuplotPipe2, "e\n");

    for (int i = 0; i < size; i++) {
        fprintf(gnuplotPipe2, "%g %g\n",tableau_energy[i], tab_flux_U_238[i]);
    }
    fprintf(gnuplotPipe2, "e\n");

    for (int i = 0; i < size; i++) {
        fprintf(gnuplotPipe2, "%g %g\n",tableau_energy[i], tab_flux_PU_239[i]);
    }
    fprintf(gnuplotPipe2, "e\n");

    for (int i = 0; i < size; i++) {
        fprintf(gnuplotPipe2, "%g %g\n",tableau_energy[i], tab_flux_PU_241[i]);
    }
    fprintf(gnuplotPipe2, "e\n");
    
    for (int i = 0; i < size; i++) {
        fprintf(gnuplotPipe2, "%g %g\n",tableau_energy[i], tab_cross[i]);
    }
    fprintf(gnuplotPipe2, "e\n");
    
    for (int i = 0; i < size; i++) {
        fprintf(gnuplotPipe2, "%g %g\n",tableau_energy[i], tab_flux_cross_U_235[i]);
    }
    fprintf(gnuplotPipe2, "e\n");

    for (int i = 0; i < size; i++) {
        fprintf(gnuplotPipe2, "%g %g\n",tableau_energy[i], tab_flux_cross_U_238[i]);
    }
    fprintf(gnuplotPipe2, "e\n");

    for (int i = 0; i < size; i++) {
        fprintf(gnuplotPipe2, "%g %g\n",tableau_energy[i], tab_flux_cross_PU_239[i]);
    }
    fprintf(gnuplotPipe2, "e\n");

    for (int i = 0; i < size; i++) {
        fprintf(gnuplotPipe2, "%g %g\n",tableau_energy[i], tab_flux_cross_PU_241[i]);
    }
    fprintf(gnuplotPipe2, "e\n");

    fflush(gnuplotPipe2);

    // Attente de l'utilisateur avant de fermer la fenêtre Gnuplot
    printf("Appuyez sur Entrée pour fermer le graphique...\n");
    getchar();

    // Fermeture du pipe Gnuplot
    fclose(gnuplotPipe1);
    fclose(gnuplotPipe2);

    return 0;
}