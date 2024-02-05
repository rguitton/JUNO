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
    fprintf(gnuplotPipe1,"set title 'reactor antineutrino flux for different neutrino MHs'\n");
    fprintf(gnuplotPipe1, "set xlabel 'L/E(km/MeV)'\n");
    fprintf(gnuplotPipe1, "set ylabel 'Arbitrary unit'\n");
    fprintf(gnuplotPipe1, "plot '-' with linespoints title 'Non oscillation', '-' with linespoints title 'theta_21', '-' with linespoints title 'NH', '-' with linespoints title 'IH'\n");

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
    
    fprintf(gnuplotPipe2, "set multiplot layout 2, 1\n");
    fprintf(gnuplotPipe2, "set title 'Flux et spectre d émission'\n");
    fprintf(gnuplotPipe2, "set xlabel 'antineutrino energy'\n");
    fprintf(gnuplotPipe2, "set xrange [0:9]\n");

    fprintf(gnuplotPipe2, "set ylabel 'Number of antineutrino (MeV^-1.fission^-1)'\n");
    fprintf(gnuplotPipe2, "set ytics nomirror\n");
    fprintf(gnuplotPipe2, "set yrange [0:0.9]\n");

    fprintf(gnuplotPipe2, "set y2label 'Cross section(cm^-2)'\n");
    fprintf(gnuplotPipe2, "set y2tics nomirror\n");
    fprintf(gnuplotPipe2, "set y2range [0:5e-42]\n");

    //fprintf(gnuplotPipe2, "plot '-' with linespoints lt 3 title 'flux U 235', '-' with linespoints lt 4 title 'flux U 238', '-' with linespoints lt 1 title 'flux PU 239', '-' with linespoints lt 2 title 'flux PU 241', '-' with linespoints axes x1y2 title 'cross section','-' with linespoints axes x1y2 lt 3 title 'U235 product','-' with linespoints axes x1y2 lt 4 title 'U238 product','-' with linespoints axes x1y2 lt 1 title 'PU239 product','-' with linespoints axes x1y2 lt 2 title 'PU241 product'\n");
    fprintf(gnuplotPipe2, "plot '-' with linespoints lt 3 title 'flux U 235', '-' with linespoints lt 4 title 'flux U 238', '-' with linespoints lt 1 title 'flux PU 239', '-' with linespoints lt 2 title 'flux PU 241', '-' with linespoints axes x1y2 title 'cross section','-' with linespoints axes x1y2 lt 4 title 'U238 product'\n");

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
        fprintf(gnuplotPipe2, "%g %g\n",tableau_energy[i], tab_flux_cross_U_238[i]);
    }
    fprintf(gnuplotPipe2, "e\n");

    /*
    for (int i = 0; i < size; i++) {
        fprintf(gnuplotPipe2, "%g %g\n",tableau_energy[i], tab_flux_cross_U_235[i]);
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

    */    
    fprintf(gnuplotPipe2, "set title 'Donnée de l article'\n");
    fprintf(gnuplotPipe2, "set datafile separator ','\n plot '../data/fig_2_6.txt' using 1:2 with linespoints axes x1y1 title 'U235', '' using 3:4 with linespoints axes x1y1 title 'Pu239', '' using 5:6 with linespoints axes x1y1 title 'U238', '' using 7:8 with linespoints axes x1y1 title 'Pu241'\n");

    fprintf(gnuplotPipe2, "unset multiplot\n");
    fflush(gnuplotPipe2);

//Figure 5 article 2
//The kinetic energy of the positron, together with the typically two 0.511 MeV annihilation photons, is assumed to be fully deposited in the detector article 2 p 1515 and is defined as Edep.
double spectra_NH_dep[size];
double spectra_IH_dep[size];
double tableau_energy_deposited[size];
double Np=1.44e33;
 for (int i=0;i<size;i++){
        tableau_energy_deposited[i]=energy_positron(tableau_energy[i])+2*0.511;
        spectra_NH_dep[i]=flux(tableau_energy_deposited[i])*sigma(tableau_energy_deposited[i])*probability(tableau_energy_deposited[i], 'N', 1);
        spectra_IH_dep[i]=flux(tableau_energy_deposited[i])*sigma(tableau_energy_deposited[i])*probability(tableau_energy_deposited[i], 'I', 1);
 }
// sur 6 ans
 for (int i=0;i<size;i++){
        spectra_NH_dep[i]*=3600*24*365*6;
        spectra_IH_dep[i]*=3600*24*365*6;
 }
    FILE *gnuplotPipe3 = popen("gnuplot -persist", "w");
    if (gnuplotPipe3 == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture de Gnuplot.\n");
        return 1;
    }
    fprintf(gnuplotPipe3, "set title 'Juno 6 years data datking'\n");
    fprintf(gnuplotPipe2, "set xlabel 'Deposited energy [MeV]'\n");
    fprintf(gnuplotPipe2, "set ylabel 'Events per 20 keV'\n");

    fprintf(gnuplotPipe3, "plot '-' with linespoints title 'NH','-' with linespoints title 'IH'\n");

    for (int i = 0; i < size; i++) {    
        fprintf(gnuplotPipe3, "%g %g\n", tableau_energy_deposited[i], Np*spectra_NH_dep[i]/(20e3));
    }
    fprintf(gnuplotPipe3, "e\n");

    for (int i = 0; i < size; i++) {
        fprintf(gnuplotPipe3, "%g %g\n", tableau_energy_deposited[i], Np*spectra_IH_dep[i]/(20e3));
    }
    fprintf(gnuplotPipe3, "e\n");
        
    
    fflush(gnuplotPipe3);

    //Figure 4 article 2
//Evis ≃ E − 0.8MeV. p36 article 1
double spectra_NH_vis[size];
double spectra_IH_vis[size];
double tableau_energy_visible[size];
 for (int i=0;i<size;i++){
        tableau_energy_visible[i]=tableau_energy[i]-0.8;
        spectra_NH_vis[i]=flux(tableau_energy_visible[i])*sigma(tableau_energy_visible[i])*probability(tableau_energy_visible[i], 'N', 1);
        spectra_IH_vis[i]=flux(tableau_energy_visible[i])*sigma(tableau_energy_visible[i])*probability(tableau_energy_visible[i], 'I', 1);
 }
// sur 1 jour
 for (int i=0;i<size;i++){
        spectra_NH_vis[i]*=3600*24;
        spectra_IH_vis[i]*=3600*24;
 }
    FILE *gnuplotPipe4 = popen("gnuplot -persist", "w");
    if (gnuplotPipe4 == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture de Gnuplot.\n");
        return 1;
    }
    fprintf(gnuplotPipe4, "set title 'Juno 1 day data datking'\n");
    fprintf(gnuplotPipe4, "set xlabel 'Deposited energy [MeV]'\n");
    fprintf(gnuplotPipe4, "set ylabel 'Events/0.02 [Mev^-1.Day^-1]'\n");

    fprintf(gnuplotPipe4, "plot '-' with linespoints title 'NH','-' with linespoints title 'IH'\n");

    for (int i = 0; i < size; i++) {    
        fprintf(gnuplotPipe4, "%g %g\n", tableau_energy_visible[i], Np*spectra_NH_dep[i]/(0.02));
    }
    fprintf(gnuplotPipe4, "e\n");

    for (int i = 0; i < size; i++) {
        fprintf(gnuplotPipe4, "%g %g\n", tableau_energy_visible[i], Np*spectra_IH_dep[i]/(0.02));
    }
    fprintf(gnuplotPipe4, "e\n");
        
    
    fflush(gnuplotPipe4);


    //recherche du 83
    printf("il y a %g event par jour \n", Np*integrale_spectre(5.0,25.0,1000)*3600*24);
    // Attente de l'utilisateur avant de fermer la fenêtre Gnuplot
    printf("Appuyez sur Entrée pour fermer le graphique...\n");
    getchar();

    // Fermeture du pipe Gnuplot
    fclose(gnuplotPipe1);
    fclose(gnuplotPipe2);
    fclose(gnuplotPipe3);
    fclose(gnuplotPipe4);

    return 0;
}