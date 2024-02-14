#include "fonctions.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstdlib>
#include <fstream>

using namespace std;

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

    double fission_fraction_U_235=0.58;//MeV
    double fission_fraction_U_238=0.07;//MeV
    double fission_fraction_PU_239=0.30;//MeV
    double fission_fraction_PU_241=0.05;//MeV
    
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

    // Définition/ouverture de fichiers utiles aux plots ROOT, format %g %g %g %g
    // spectra
    fstream file_spectra_initial; file_spectra_initial.open("../results/spectra_initial.txt",ios::out);
    fstream file_spectra_p21; file_spectra_p21.open("../results/spectra_p21.txt",ios::out);
    fstream file_spectra_NH; file_spectra_NH.open("../results/spectra_NH.txt",ios::out);
    fstream file_spectra_IH; file_spectra_IH.open("../results/spectra_IH.txt",ios::out);
    // flux
    fstream file_flux_U_235; file_flux_U_235.open("../results/flux_U_235.txt",ios::out);
    fstream file_flux_U_238; file_flux_U_238.open("../results/flux_U_238.txt",ios::out);
    fstream file_flux_PU_239; file_flux_PU_239.open("../results/flux_PU_239.txt",ios::out);
    fstream file_flux_PU_241; file_flux_PU_241.open("../results/flux_PU_241.txt",ios::out);
    fstream file_flux_cross_U_238; file_flux_cross_U_238.open("../results/flux_cross_U_235.txt",ios::out);
    fstream file_cross; file_cross.open("../results/cross.txt",ios::out);


    for (int i=0;i<size;i++){
        tableau_lenght[i]=5+double(i)*28/size;
        tableau_energy[i]=1.81+double(i)*7.19/size;
    }

     for (int i=0;i<size;i++){
        double energy=L/tableau_lenght[i];
        spectra_initial[i]=flux(energy)*sigma(energy);
        spectra_p21[i]=flux(energy)*sigma(energy)*probability(energy, 'I', 0);
        spectra_NH[i]=flux(energy)*sigma(energy)*probability(energy, 'N', 1);
        spectra_IH[i]=flux(energy)*sigma(energy)*probability(energy, 'I', 1);

        tab_flux_U_235[i]=fission_fraction_U_235*flux_U_235(tableau_energy[i]);
        tab_flux_U_238[i]=fission_fraction_U_238*flux_U_238(tableau_energy[i]);
        tab_flux_PU_239[i]=fission_fraction_PU_239*flux_PU_239(tableau_energy[i]);
        tab_flux_PU_241[i]=fission_fraction_PU_241*flux_PU_241(tableau_energy[i]);

        tab_cross[i]=sigma(tableau_energy[i]);

        tab_flux_cross_U_235[i]=tab_cross[i]*tab_flux_U_235[i]*pow(10,2);
        tab_flux_cross_U_238[i]=tab_cross[i]*tab_flux_U_238[i]*pow(10,2);
        tab_flux_cross_PU_239[i]=tab_cross[i]*tab_flux_PU_239[i]*pow(10,2);
        tab_flux_cross_PU_241[i]=tab_cross[i]*tab_flux_PU_241[i]*pow(10,2);

        // Sauvegarde des données dans les fichiers .txt

        // spectra fig2.4, art[1]
        file_spectra_initial << tableau_lenght[i] << " " << spectra_initial[i] << " " << 0 << " " << 0 << endl;
        file_spectra_p21 << tableau_lenght[i] << " " << spectra_p21[i] << " " << 0 << " " << 0 << endl;
        file_spectra_NH << tableau_lenght[i] << " " << spectra_NH[i] << " " << 0 << " " << 0 << endl;
        file_spectra_IH << tableau_lenght[i] << " " << spectra_IH[i] << " " << 0 << " " << 0 << endl;
        // flux from reactor fig2.6, art[1]
        file_flux_U_235 << tableau_energy[i] << " " << tab_flux_U_235[i] << " " << 0 << " "<< 0 << endl;
        file_flux_U_238 << tableau_energy[i] << " " << tab_flux_U_238[i] << " " << 0 << " " << 0 << endl;
        file_flux_PU_239 << tableau_energy[i] << " " << tab_flux_PU_239[i] << " " << 0 << " " << 0 << endl;
        file_flux_PU_241 << tableau_energy[i] << " " << tab_flux_PU_241[i] << " " << 0 << " " << 0 << endl;
        file_flux_cross_U_238 << tableau_energy[i] << " " << tab_flux_cross_U_238[i] << " " << 0 << " " << 0 << endl;
        file_cross << tableau_energy[i] << " " << tab_cross[i] << " " << 0 << " " << 0 << endl;

    }

    // Fermeture des fichiers
    //spectra
    file_spectra_initial.close();
    file_spectra_p21.close();
    file_spectra_NH.close();
    file_spectra_IH.close();
    //flux
    file_flux_U_235.close();
    file_flux_U_238.close();
    file_flux_PU_239.close();
    file_flux_PU_241.close();
    file_flux_cross_U_238.close();
    file_cross.close();

    

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
/*

 
        spectra_NH_dep[i]=flux(tableau_energy_deposited[i])*sigma(tableau_energy_deposited[i])*probability(tableau_energy_deposited[i], 'N', 1);
        spectra_IH_dep[i]=flux(tableau_energy_deposited[i])*sigma(tableau_energy_deposited[i])*probability(tableau_energy_deposited[i], 'I', 1);
 }
// sur 6 ans
 for (int i=0;i<size;i++){
        spectra_NH_dep[i]*=3600*24*365*6;
        spectra_IH_dep[i]*=3600*24*365*6;
 }
 */

    double tableau_energy_deposited[size];
    double Np=1.44e33;
    for (int i=0;i<size;i++){
        tableau_energy_deposited[i]=energy_positron(tableau_energy[i])+2*0.511;}
        //printf("le spectre vaut %g \n",tableau_energy_deposited[i]);}
    
    double flux_per_fission_dep[size];
    for(int i=0;i<size;i++){
        flux_per_fission_dep[i]=flux(tableau_energy_deposited[i]);
    }

    double flux_total_dep[size];
    for(int i=0;i<size;i++){
        flux_total_dep[i]=total_reactor_flux(flux_per_fission_dep[i],36e9);
    }
    double spectre_final_dep[size];
    for(int i=0;i<size;i++){
        //total_spectre_final[i]=total_flux[i]*sigma(tableau_energy[i])*probability(tableau_energy[i], 'I', 1);     
        spectre_final_dep[i]=calcul_spectre(flux_total_dep[i],tableau_energy_deposited[i]);
        }
    FILE *gnuplotPipe3 = popen("gnuplot -persist", "w");
    if (gnuplotPipe3 == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture de Gnuplot.\n");
        return 1;
    }
    fprintf(gnuplotPipe3, "set title 'Juno 6 years data datking'\n");
    fprintf(gnuplotPipe3, "set xlabel 'Deposited energy [MeV]'\n");
    fprintf(gnuplotPipe3, "set ylabel 'Events per 20 keV'\n");

    fprintf(gnuplotPipe3, "plot '-' with linespoints title 'IH'\n");

    for (int i = 0; i < size; i++) {    
        fprintf(gnuplotPipe3, "%g %g\n", tableau_energy_deposited[i], Np*spectre_final_dep[i]*3600*24*365*6);
    }
    fprintf(gnuplotPipe3, "e\n");

    fflush(gnuplotPipe3);

    //Figure 4 article 2
//Evis ≃ E − 0.8MeV. p36 article 1

    double tableau_energy_vis[size];
    for (int i=0;i<size;i++){
        tableau_energy_vis[i]=tableau_energy[i]-0.8;}
            
    double flux_per_fission_vis[size];
    for(int i=0;i<size;i++){
        flux_per_fission_vis[i]=flux(tableau_energy_vis[i]);
    }

    double flux_total_vis[size];
    for(int i=0;i<size;i++){
        flux_total_vis[i]=total_reactor_flux(flux_per_fission_vis[i],36e9);
    }
    double spectre_final_vis[size];
    for(int i=0;i<size;i++){
        //total_spectre_final[i]=total_flux[i]*sigma(tableau_energy[i])*probability(tableau_energy[i], 'I', 1);     
        spectre_final_vis[i]=calcul_spectre(flux_total_vis[i],tableau_energy_vis[i]);
        }
    /*
    double event_per_10kev[size];
    double energy_per_10kev[size];

    double Emin=tableau_energy_vis[0];
    for(int i=0;i<size;i++){
        int u=0;
        double Emax=tableau_energy_vis[i];
        while(Emax-Emin<10e-3 && u<size){
            double Emax=tableau_energy_vis[i+u];
            u=u+1;
            }
        event_per_10kev[i]=(spectre_final_vis[i+u]+tableau_energy_vis[u])/2;
        energy_per_10kev[i]= (Emax+Emin)/2;
        i=u;
        }
    */
    FILE *gnuplotPipe4 = popen("gnuplot -persist", "w");
    if (gnuplotPipe4 == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture de Gnuplot.\n");
        return 1;
    }
    fprintf(gnuplotPipe4, "set title 'Juno 1 day data datking'\n");
    fprintf(gnuplotPipe4, "set xlabel 'visible energy [MeV]'\n");
    fprintf(gnuplotPipe4, "set ylabel 'Events'\n");
    //fprintf(gnuplotPipe4, "set yrange [0:0.5]\n");


    fprintf(gnuplotPipe4, "plot '-' with linespoints title 'IH'\n");

    for (int i = 0; i < size; i++) {    
        fprintf(gnuplotPipe4, "%g %g\n", tableau_energy_vis[i], Np*spectre_final_vis[i]*3600*24);
    }
    fprintf(gnuplotPipe4, "e\n");
    
    
    fflush(gnuplotPipe4);

    double flux_per_fission[size];
    for(int i=0;i<size;i++){
        flux_per_fission[i]=flux(tableau_energy[i]);
    }

    double flux_total[size];
    for(int i=0;i<size;i++){
        flux_total[i]=total_reactor_flux(flux_per_fission[i],36e9);
    }
    double spectre_final[size];
    for(int i=0;i<size;i++){
        //total_spectre_final[i]=total_flux[i]*sigma(tableau_energy[i])*probability(tableau_energy[i], 'I', 1);     
        spectre_final[i]=calcul_spectre(flux_total[i],tableau_energy[i]);
    }

    FILE *gnuplotPipe5 = popen("gnuplot -persist", "w");
    if (gnuplotPipe5 == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture de Gnuplot.\n");
        return 1;
    }
    fprintf(gnuplotPipe5, "set title 'Spectre des antineutrino detectés en fonction de l énergie'\n");
    fprintf(gnuplotPipe5, "set xlabel 'Antineutrino energy [MeV]'\n");
    fprintf(gnuplotPipe5, "set ylabel 'Nombre d antineutrinos detectés par jour'\n");

    fprintf(gnuplotPipe5, "plot '-' with linespoints title 'IH'\n");

    for (int i = 0; i < size; i++) {    
        fprintf(gnuplotPipe5, "%g %g\n", tableau_energy[i],spectre_final[i]*Np*3600*24);
    }
    fprintf(gnuplotPipe5, "e\n");
    fflush(gnuplotPipe5);
    /*
    double *E_converted = convert_Evis(tableau_energy, size);
    printf("Énergies converties :\n");
    for (int i = 0; i < size/1000; i++) {
        printf("%f\n", E_converted[i]);
    }*/
    double IBD_initial=0;//par seconde
    double h=tableau_energy[10]-tableau_energy[9];
    for(int i=0;i<size-1;i++){
        if(tableau_energy[i]>1.8&&tableau_energy[i]<12){//cf article 2
        IBD_initial += h*(spectre_final[i] + spectre_final[i+1])/2;}
        
    }
    double IBD_detected_per_d=Np*IBD_initial*3600*24;

    double new_spectre[size];
    for(int u=0;u<size;u++){
        new_spectre[u]=0;
        for(int i = 0; i < size; i++){
           //new_spectre[u] += spectre_final[i]*Np*3600*24*gauss_pdf(tableau_energy[u],tableau_energy[i]);
           double product=gauss_pdf(tableau_energy[u],tableau_energy[i])*spectre_final[i];
           if(product>0){
                new_spectre[u] +=product;
           }
           
    }printf("new_energy vaut %g\n",new_spectre[u]);
           }
        FILE *gnuplotPipe6 = popen("gnuplot -persist", "w");
    if (gnuplotPipe6 == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture de Gnuplot.\n");
        return 1;
    }
    fprintf(gnuplotPipe6, "set title 'Spectre des antineutrino detectés en fonction de l énergie visible'\n");
    fprintf(gnuplotPipe6, "set xlabel 'Visible energy [MeV]'\n");
    fprintf(gnuplotPipe6, "set ylabel 'Nombre d antineutrinos detectés par jour'\n");

    fprintf(gnuplotPipe6, "plot '-' with linespoints title 'IH'\n");

    for (int i = 0; i < size; i++) {    
        fprintf(gnuplotPipe6, "%g %g\n", tableau_energy[i],new_spectre[i]*Np*3600*24);
    }
    fprintf(gnuplotPipe6, "e\n");
    fflush(gnuplotPipe6);
    //recherche du 83
    //printf("il y a %g event par jour \n", Np*integrale_spectre(5.0,25.0,1000)*3600*24);
    printf("On peut s'attendre à %g détections par jour sans sélection \n", IBD_detected_per_d);
    printf("On peut s'attendre à %g détections par jour après la sélectivité Fiducial volume (91% ) \n", IBD_detected_per_d*0.91);
    printf("On peut s'attendre à %g détections par jour après la sélectivité Energy cut (97.8% ) \n", IBD_detected_per_d*0.91*0.978);
    printf("On peut s'attendre à %g détections par jour après la sélectivité Time cut (99.1% ) \n", IBD_detected_per_d*0.91*0.978*0.991);
    printf("On peut s'attendre à %g détections par jour après la sélectivité Vertex cut (98.7% ) \n", IBD_detected_per_d*0.91*0.978*0.991*0.987);
    printf("On peut s'attendre à %g détections par jour après la sélectivité Muon veto (83% ) \n", IBD_detected_per_d*0.91*0.978*0.991*0.987*0.83);
    printf("On peut s'attendre à %g détections par jour après la sélectivité Combined (73% ) \n", IBD_detected_per_d*0.73);

    // Attente de l'utilisateur avant de fermer la fenêtre Gnuplot
    printf("Appuyez sur Entrée pour fermer le graphique...\n");
    getchar();

    // Fermeture du pipe Gnuplot
    fclose(gnuplotPipe2);
    fclose(gnuplotPipe3);
    fclose(gnuplotPipe4);
    fclose(gnuplotPipe5);

    // Exécution des macros ROOT pour les plots
    printf("####################################################\n\n");
    printf("Current directory during the execution of main.cpp\n");
    printf("\t");
    system("pwd");
    printf("\n");
    system("root -q ../src/plotter_spectra.C");
    return 0;
}