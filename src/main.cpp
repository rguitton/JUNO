#include "fonctions.h"
#include "plotter.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstdlib>
#include <fstream>

using namespace std;

const char* dir = "../output/";


 double trouverMinimum(const double tableau[], int taille) {
    double minimum = tableau[0]; // Supposons que le premier élément est le minimum initial

    // Parcours du tableau pour trouver le minimum
    for (int i = 1; i < taille; ++i) {
        if (tableau[i] < minimum) {
            minimum = tableau[i]; // Mise à jour du minimum si on trouve un élément plus petit
        }
    }

    return minimum;
    }

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
    double tab_flux_cross_total[size];

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

    double sin2_teta_12=0.307;
    double sin2_teta_13=0.0218;   

    double teta_12=asin(sqrt(sin2_teta_12));
    double teta_13=asin(sqrt(sin2_teta_13));

    double delta2_m21=7.5e-5;//valeur absolue
    double delta2_m32=2.4e-3;//on peut aussi faire simplmement un changement de signe
    //double delta2_m31=2.53e-3;

    double delta2_m31_NH=delta2_m32+delta2_m21;
    double delta2_m31_IH=delta2_m32-delta2_m21;

    int N=100;

    double T_delta2_m21[N];
    double T_delta2_m32[N];
    double T_delta2_m31_NH[N];
    double T_delta2_m31_IH[N];

    double T_sin2_teta_12[N];
    double T_sin2_teta_13[N];

    for(int j=0;j<N;j++){
        T_delta2_m21[j]=delta2_m21+1e-5*(j-(N/2))/2083;//50/0.024=2083
        T_delta2_m32[j]=delta2_m32+1e-3*(j-(N/2))/10638;//50/0.0047=10638       //il faut balayer plus de valeur pour arriver jusqu'à dchi²=1
        //T_delta2_m31_NH[j]=delta2_m31_NH+1e-3*(j-(N/2))/10638;//50/0.0047=10638
        //T_delta2_m31_IH[j]=-(delta2_m31_IH+1e-3*(j-(N/2))/10638);//50/0.0047=10638
        T_delta2_m31_NH[j]=delta2_m31_NH+1e-3*(j-(N/2))/500;//50/0.0047=10638
        T_delta2_m31_IH[j]=-(delta2_m31_IH+1e-3*(j-(N/2))/800);//50/0.0047=10638
        T_sin2_teta_12[j]=pow(sin(teta_12),2)+1.0*(j-(N/2))/31250;//50/0.0016=31250
        T_sin2_teta_13[j]=pow(sin(teta_13),2)+1.0*(j-(N/2))/19231;//50/0.0026=19231
        printf("T_delta2_m31_IH vaut %g et T_delta2_m31_NH vaut %g \n", T_delta2_m31_IH[j],T_delta2_m31_NH[j]);
        }

    for (int i=0;i<size;i++){

        tableau_lenght[i]=5+double(i)*28/size;
        tableau_energy[i]=1.8+double(i)*10/size;
 
        //tableau_lenght[i]=5+double(i)*28/size;
        //tableau_energy[i]=1.81+double(i)*7.19/size;
    }

     for (int i=0;i<size;i++){
        double energy=L/tableau_lenght[i];
        spectra_initial[i]=flux(energy)*sigma(energy);
        spectra_p21[i]=flux(energy)*sigma(energy)*probability(energy, 'I', 0, sin2_teta_13, sin2_teta_12, delta2_m21, delta2_m32,delta2_m31_IH );
        spectra_NH[i]=flux(energy)*sigma(energy)*probability(energy, 'N', 1, sin2_teta_13, sin2_teta_12, delta2_m21, delta2_m32,delta2_m31_NH);
        spectra_IH[i]=flux(energy)*sigma(energy)*probability(energy, 'I', 1, sin2_teta_13, sin2_teta_12, delta2_m21, delta2_m32,delta2_m31_IH);

        tab_flux_U_235[i]=fission_fraction_U_235*flux_U_235(tableau_energy[i]);
        tab_flux_U_238[i]=fission_fraction_U_238*flux_U_238(tableau_energy[i]);
        tab_flux_PU_239[i]=fission_fraction_PU_239*flux_PU_239(tableau_energy[i]);
        tab_flux_PU_241[i]=fission_fraction_PU_241*flux_PU_241(tableau_energy[i]);

        tab_cross[i]=sigma(tableau_energy[i]);

        tab_flux_cross_U_235[i]=tab_cross[i]*tab_flux_U_235[i];
        tab_flux_cross_U_238[i]=tab_cross[i]*tab_flux_U_238[i];
        tab_flux_cross_PU_239[i]=tab_cross[i]*tab_flux_PU_239[i];
        tab_flux_cross_PU_241[i]=tab_cross[i]*tab_flux_PU_241[i];
        tab_flux_cross_total[i]= tab_flux_cross_U_235[i] + tab_flux_cross_U_238[i] + tab_flux_cross_PU_239[i] + tab_flux_cross_PU_241[i];
        tab_flux_cross_total[i]=tab_flux_cross_total[i]*20; //arbitrary unit (we just wan to plot the shape)
    }


// PLOTTER : spectra fig2.4, art[1]
    cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% spectra"  << endl;
    for (int i = 0; i < size; i++) {    
        cout << tableau_lenght[i] << " " <<  spectra_initial[i] << endl;
    }
    cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"  << endl;

    const double* all_yaxis1[4]={spectra_initial,spectra_p21,spectra_NH,spectra_IH};
    plotter_spectra(size,tableau_lenght,all_yaxis1,dir);


// PLOTTER : flux from reactor fig2.6, art[1]
    //const double* all_yaxis2[6]={tab_flux_U_235,tab_flux_U_238,tab_flux_PU_239,tab_flux_PU_241,tab_cross,tab_flux_cross_total};
    const double* all_yaxis2[6]={tab_flux_U_235,tab_flux_U_238,tab_flux_PU_239,tab_flux_PU_241,tab_cross,tab_flux_cross_total};
    plotter_flux(size, tableau_energy, all_yaxis2, dir);


    /*
    fprintf(gnuplotPipe2, "set title 'Donnée de l article'\n");
    fprintf(gnuplotPipe2, "set datafile separator ','\n plot '../data/fig_2_6.txt' using 1:2 with linespoints axes x1y1 title 'U235', '' using 3:4 with linespoints axes x1y1 title 'Pu239', '' using 5:6 with linespoints axes x1y1 title 'U238', '' using 7:8 with linespoints axes x1y1 title 'Pu241'\n");

    fprintf(gnuplotPipe2, "unset multiplot\n");
    fflush(gnuplotPipe2); 
    */

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
        spectre_final_dep[i]=calcul_spectre(flux_total_dep[i],tableau_energy_deposited[i], sin2_teta_13, sin2_teta_12, delta2_m21, delta2_m32, 'I',delta2_m31_IH);
        }
 

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
        spectre_final_vis[i]=calcul_spectre(flux_total_vis[i],tableau_energy_vis[i], sin2_teta_13, sin2_teta_12, delta2_m21, delta2_m32, 'I',delta2_m31_IH);
        }


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

//PLOTTER : visible energy spectrum: Figure 4, article 2
    cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% visible energy"  << endl;
    for (int i = 0; i < size; i++) {    
        cout << tableau_energy_vis[i] << " " <<  Np*spectre_final_vis[i]*3600*24 << endl;
    }
    cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"  << endl;
    double* tmp = new double[size];
	for (int i=0;i<size;i++){
		tmp[i]= spectre_final_vis[i]*Np*3600*24;
	}
    plotter_visible_energy_spectrum(size,tableau_energy_vis,spectre_final_vis,dir);
    delete[] tmp;
   

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
        spectre_final[i]=calcul_spectre(flux_total[i],tableau_energy[i], sin2_teta_13, sin2_teta_12, delta2_m21, delta2_m32, 'I',delta2_m31_IH);
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
           }}
        FILE *gnuplotPipe6 = popen("gnuplot -persist", "w");
    if (gnuplotPipe6 == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture de Gnuplot.\n");
        return 1;
    }
    fprintf(gnuplotPipe6, "set title 'Spectre des antineutrino detectés en fonction de l énergie visible'\n");
    fprintf(gnuplotPipe6, "set xlabel 'Visible energy [MeV]'\n");
    fprintf(gnuplotPipe6, "set ylabel 'Nombre d antineutrinos detectés par jour'\n");
    fprintf(gnuplotPipe6, "set xrange [1:12]\n");

    fprintf(gnuplotPipe6, "plot '-' with linespoints title 'IH'\n");

    for (int i = 0; i < size; i++) {    
        fprintf(gnuplotPipe6, "%g %g\n", tableau_energy[i]-0.8,new_spectre[i]*Np*3600*24/10000);
    }
    fprintf(gnuplotPipe6, "e\n");
    fflush(gnuplotPipe6);

        FILE *gnuplotPipe3 = popen("gnuplot -persist", "w");
    if (gnuplotPipe3 == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture de Gnuplot.\n");
        return 1;
    }
    fprintf(gnuplotPipe3, "set title 'Juno 6 years data datking'\n");
    fprintf(gnuplotPipe3, "set xlabel 'Deposited energy [MeV]'\n");
    fprintf(gnuplotPipe3, "set ylabel 'Events number'\n");
    fprintf(gnuplotPipe3, "set xrange [1:12]\n");

    fprintf(gnuplotPipe3, "plot '-' with linespoints title 'IH', '-' with linespoints title 'teta13_u=0' \n");

    for (int i = 0; i < size; i++) {    
        fprintf(gnuplotPipe3, "%g %g\n", tableau_energy[i]-0.8, Np*new_spectre[i]*3600*24*365*6/10000);
    }
    
    
    fprintf(gnuplotPipe3, "e\n");

    fflush(gnuplotPipe3);


    //recherche du 83
    //printf("il y a %g event par jour \n", Np*integrale_spectre(5.0,25.0,1000)*3600*24);
    printf("On peut s'attendre à %g détections par jour sans sélection \n", IBD_detected_per_d);
    printf("On peut s'attendre à %g détections par jour après la sélectivité Fiducial volume (91%) \n", IBD_detected_per_d*0.91);
    printf("On peut s'attendre à %g détections par jour après la sélectivité Energy cut (97.8%) \n", IBD_detected_per_d*0.91*0.978);
    printf("On peut s'attendre à %g détections par jour après la sélectivité Time cut (99.1%) \n", IBD_detected_per_d*0.91*0.978*0.991);
    printf("On peut s'attendre à %g détections par jour après la sélectivité Vertex cut (98.7%) \n", IBD_detected_per_d*0.91*0.978*0.991*0.987);
    printf("On peut s'attendre à %g détections par jour après la sélectivité Muon veto (83%) \n", IBD_detected_per_d*0.91*0.978*0.991*0.987*0.83);
    printf("On peut s'attendre à %g détections par jour après la sélectivité Combined (73%) \n", IBD_detected_per_d*0.73);


    /*
    int count=0;
    for(int i=0;i<size;i++){
        if(tableau_energy[i]>1.8&&tableau_energy[i]<8.0){//cf article 1 p43
            count=count+1;
    }}

    double M_neutrino[count];
    double T_neutrino[count];
    double chi=0;
    for(int i=0;i<size-1;i++){
        if(tableau_energy[i]>1.8&&tableau_energy[i]<8.0){//cf article 1 p43
            T_neutrino[i]=Np*h*(spectre_final[i] + spectre_final[i+1])/2*3600*24*365*6;
             M_neutrino[i]=Np*h*(new_spectre[i]+new_spectre[i+1])/2*3600*24*365*6;
            //printf("Mi vaut %g \n", M_neutrino[i]);
            //printf("Ti vaut %g \n", T_neutrino[i]);
            chi+=pow((M_neutrino[i]-T_neutrino[i]),2)/M_neutrino[i];
            }
    }

    printf("chi vaut %g \n", chi);
    */

    double Ti_spectre_final_teta_13[N][size];//N represente l'indice de ligne, et donc de la variation (à 50 on est à 0)
    double Ti_spectre_final_teta_12[N][size];
    double Ti_spectre_final_delta2_m21[N][size];
    double Ti_spectre_final_delta2_m31_IH[N][size];
    double Ti_spectre_final_delta2_m31_NH[N][size];

    for(int u=0;u<N;u++){
    for(int i=0;i<size;i++){//on peut se concentrer sur 1 plot
        //total_spectre_final[i]=total_flux[i]*sigma(tableau_energy[i])*probability(tableau_energy[i], 'I', 1);     
        Ti_spectre_final_teta_13[u][i]=calcul_spectre(flux_total[i],tableau_energy[i], T_sin2_teta_13[u], T_sin2_teta_12[50], T_delta2_m21[50], T_delta2_m32[50], 'I',delta2_m31_IH)*Np;
        Ti_spectre_final_teta_12[u][i]=calcul_spectre(flux_total[i],tableau_energy[i], T_sin2_teta_13[50], T_sin2_teta_12[u], T_delta2_m21[50], T_delta2_m32[50], 'I',delta2_m31_IH)*Np;
        //Ti_spectre_final_teta_12[u][i]=calcul_spectre(flux_total[i],tableau_energy[i], T_sin2_teta_13[50], T_sin2_teta_12[u], T_delta2_m21[50], T_delta2_m32[50], 'I',delta2_m31_IH)*Np;

        Ti_spectre_final_delta2_m21[u][i]=calcul_spectre(flux_total[i],tableau_energy[i], T_sin2_teta_13[50], T_sin2_teta_12[50], T_delta2_m21[u], T_delta2_m32[50], 'I',delta2_m31_IH)*Np;
        Ti_spectre_final_delta2_m31_IH[u][i]=calcul_spectre(flux_total[i],tableau_energy[i], T_sin2_teta_13[50], T_sin2_teta_12[50], T_delta2_m21[50], T_delta2_m32[50], 'I',T_delta2_m31_IH[u])*Np;
        Ti_spectre_final_delta2_m31_NH[u][i]=calcul_spectre(flux_total[i],tableau_energy[i], T_sin2_teta_13[50], T_sin2_teta_12[50], T_delta2_m21[50], T_delta2_m32[50], 'N',T_delta2_m31_NH[u])*Np;

    }}

    double Ti_new_spectre_teta_13[N][size];
    double Ti_new_spectre_teta_12[N][size];
    double Ti_new_spectre_delta2_m21[N][size];
    double Ti_new_spectre_delta2_m31_IH[N][size];
    double Ti_new_spectre_delta2_m31_NH[N][size];

    for(int k=0;k<N;k++){
    for(int u=0;u<size;u++){
        Ti_new_spectre_teta_13[k][u]=0;
        Ti_new_spectre_teta_12[k][u]=0;
        Ti_new_spectre_delta2_m21[k][u]=0;
        Ti_new_spectre_delta2_m31_IH[k][u]=0;
        Ti_new_spectre_delta2_m31_NH[k][u]=0;

        for(int i = 0; i < size; i++){
           //new_spectre[u] += spectre_final[i]*Np*3600*24*gauss_pdf(tableau_energy[u],tableau_energy[i]);
            double product_teta13=gauss_pdf(tableau_energy[u],tableau_energy[i])*Ti_spectre_final_teta_13[k][i];
            double product_teta12=gauss_pdf(tableau_energy[u],tableau_energy[i])*Ti_spectre_final_teta_12[k][i];
            double product_delta2_m21=gauss_pdf(tableau_energy[u],tableau_energy[i])*Ti_spectre_final_delta2_m21[k][i];
            double product_delta2_m31_IH=gauss_pdf(tableau_energy[u],tableau_energy[i])*Ti_spectre_final_delta2_m31_IH[k][i];
            double product_delta2_m31_NH=gauss_pdf(tableau_energy[u],tableau_energy[i])*Ti_spectre_final_delta2_m31_NH[k][i];

            if(product_teta13>0 || product_teta12>0 || product_delta2_m21>0||product_delta2_m31_IH>0||product_delta2_m31_NH>0){
            Ti_new_spectre_teta_13[k][u] +=product_teta13;
            Ti_new_spectre_teta_12[k][u] +=product_teta12;
            Ti_new_spectre_delta2_m21[k][u] +=product_delta2_m21;
            Ti_new_spectre_delta2_m31_IH[k][u] +=product_delta2_m31_IH;
            Ti_new_spectre_delta2_m31_NH[k][u]+=product_delta2_m31_NH;
            }
           }
           //printf("Ti_new_spectre_teta_13 vaut %g pour k vaut %d et u vaut %d \n", Ti_new_spectre_teta_13[k][u],k,u);
           }
    }
    double chi_teta13[N];
    double chi_teta12[N];
    double chi_delta2_m21[N];
    double chi_delta2_m31_IH[N];
    double chi_delta2_m31_NH[N];

    for(int u=0;u<N;u++){
        chi_teta13[u]=0;
        chi_teta12[u]=0;
        chi_delta2_m21[u]=0;
        chi_delta2_m31_IH[u]=0;
        chi_delta2_m31_NH[u]=0;

        for(int i=0;i<size;i++){
            if(tableau_energy[i]>1.8&&tableau_energy[i]<8.0){//cf article 1 p43
                chi_teta13[u]+=pow(Ti_new_spectre_teta_13[50][i]-Ti_new_spectre_teta_13[u][i],2)/(Ti_new_spectre_teta_13[50][i]);
                chi_teta12[u]+=pow(Ti_new_spectre_teta_12[50][i]-Ti_new_spectre_teta_12[u][i],2)/(Ti_new_spectre_teta_12[50][i]);
                chi_delta2_m21[u]+=pow(Ti_new_spectre_delta2_m21[50][i]-Ti_new_spectre_delta2_m21[u][i],2)/(Ti_new_spectre_delta2_m21[50][i]);
                chi_delta2_m31_IH[u]+=pow(Ti_new_spectre_delta2_m31_IH[50][i]-Ti_new_spectre_delta2_m31_IH[u][i],2)/(Ti_new_spectre_delta2_m31_IH[50][i]);
                chi_delta2_m31_NH[u]+=pow(Ti_new_spectre_delta2_m31_IH[50][i]-Ti_new_spectre_delta2_m31_NH[u][i],2)/(Ti_new_spectre_delta2_m31_IH[50][i]);
                }
        }
        

    }

     FILE *gnuplotPipe7 = popen("gnuplot -persist", "w");
    if (gnuplotPipe7 == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture de Gnuplot.\n");
        return 1;
    }
    fprintf(gnuplotPipe7, "set title 'évolution de delta chi'\n");
    fprintf(gnuplotPipe7, "set xlabel 'delta13'\n");
    fprintf(gnuplotPipe7, "set xrange [0.0190:0.0244]\n");

    fprintf(gnuplotPipe7, "set ylabel 'deltachi'\n");

    fprintf(gnuplotPipe7, "plot '-' with linespoints title 'test'\n");

    for (int i = 0; i < size; i++) {    
        fprintf(gnuplotPipe7, "%g %g\n", T_sin2_teta_13[i], chi_teta13[i]);
    }
    fprintf(gnuplotPipe7, "e\n");

    fflush(gnuplotPipe7);



      FILE *gnuplotPipe8 = popen("gnuplot -persist", "w");
    if (gnuplotPipe8 == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture de Gnuplot.\n");
        return 1;
    }
    fprintf(gnuplotPipe8, "set title 'évolution de delta chi'\n");
    fprintf(gnuplotPipe8, "set xlabel 'variation sin2(teta12)'\n");
    fprintf(gnuplotPipe8, "set xrange [0.3055:0.3085]\n");

    fprintf(gnuplotPipe8, "set ylabel 'deltachi'\n");

    fprintf(gnuplotPipe8, "plot '-' with linespoints title 'graph'\n");

    for (int i = 0; i < size; i++) {    
        fprintf(gnuplotPipe8, "%g %g\n", T_sin2_teta_12[i], chi_teta12[i]);
        //printf("T_sin2_teta_12[i] vaut à la fin %g et chi_teta12[i] %g \n", T_sin2_teta_12[i],chi_teta12[i]);
    }
    fprintf(gnuplotPipe8, "e\n");

    fflush(gnuplotPipe8);


      FILE *gnuplotPipe9 = popen("gnuplot -persist", "w");
    if (gnuplotPipe9 == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture de Gnuplot.\n");
        return 1;
    }
    fprintf(gnuplotPipe9, "set title 'évolution de delta chi'\n");
    fprintf(gnuplotPipe9, "set xlabel 'variation _delta2_m21'\n");
    fprintf(gnuplotPipe9, "set xrange [7.25e-5:7.75e-5]\n");

    fprintf(gnuplotPipe9, "set ylabel 'deltachi'\n");

    fprintf(gnuplotPipe9, "plot '-' with linespoints title 'graph'\n");

    for (int i = 0; i < size; i++) {    
        fprintf(gnuplotPipe9, "%g %g\n",T_delta2_m21[i],chi_delta2_m21[i]);
        //printf("T_sin2_teta_12[i] vaut à la fin %g et chi_teta12[i] %g \n", T_delta2_m21[i],chi_delta2_m21[i]);
    }
    fprintf(gnuplotPipe9, "e\n");

    fflush(gnuplotPipe9);

 for (int i = 0; i < size; i++) {   
        printf("chi_delta2_m31_NH vaut %g et i vaut %d et tableau energy %g\n", chi_delta2_m31_NH[i], i, tableau_energy[i]-0.8);
        }

    FILE *gnuplotPipe10 = popen("gnuplot -persist", "w");
    if (gnuplotPipe10 == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture de Gnuplot.\n");
        return 1;
    }
    fprintf(gnuplotPipe10, "set title 'évolution de delta chi inverted mass'\n");
    fprintf(gnuplotPipe10, "set xlabel 'variation delta2 m31'\n");
    //fprintf(gnuplotPipe10, "set xrange [2.38e-3:2.5e-3]\n");
    //fprintf(gnuplotPipe10, "set xrange [0:2.5e-3]\n");

    fprintf(gnuplotPipe10, "set ylabel 'deltachi'\n");

fprintf(gnuplotPipe10, "plot '-' with points title 'chi delta2 m31 IH'\n");

    for (int i = 100; i < size; i++) {   

        //fprintf(gnuplotPipe10, "%g %g\n",  T_delta2_m31_NH[i],chi_delta2_m31_IH[i]);
        fprintf(gnuplotPipe10, "%g %g\n",  tableau_energy[i]-0.8,Ti_spectre_final_delta2_m31_NH[50][i]);

        //printf("T_delta2_m32[i] vaut à la fin %g et chi_delta2_m32[i] %g \n", T_delta2_m32[i],chi_delta2_m32[i]);
    }
    fprintf(gnuplotPipe10, "e\n");
    fflush(gnuplotPipe10);

    FILE *gnuplotPipe11 = popen("gnuplot -persist", "w");
    if (gnuplotPipe11 == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture de Gnuplot.\n");
        return 1;
    }
    fprintf(gnuplotPipe11, "set title 'comparaison de delta chi'\n");
    fprintf(gnuplotPipe11, "set xlabel 'variation  delta2 m31'\n");
    //fprintf(gnuplotPipe11, "set xrange [2.38e-3:2.5e-3]\n");
    //fprintf(gnuplotPipe10, "set xrange [0:2.5e-3]\n");
    fprintf(gnuplotPipe11, "set xrange [2.38e-3:2.6e-3]\n");

    fprintf(gnuplotPipe11, "set ylabel 'deltachi'\n");

    fprintf(gnuplotPipe11, "plot '-' with linespoints title 'chi delta2 m31 NH', '-' with linespoints title 'chi delta2 m31 IH'\n");
    
    for (int i = 0; i < N; i++) {   
        fprintf(gnuplotPipe11, "%g %g\n",  T_delta2_m31_NH[i],chi_delta2_m31_NH[i]);
        //printf("T_delta2_m32[i] vaut à la fin %g et chi_delta2_m32[i] %g \n", T_delta2_m32[i],chi_delta2_m32[i]);
    }
    fprintf(gnuplotPipe11, "e\n");

    for (int i = 0; i < N; i++) {   
        fprintf(gnuplotPipe11, "%g %g\n",  T_delta2_m31_NH[i],chi_delta2_m31_IH[i] );
        //printf("T_delta2_m32[i] vaut à la fin %g et chi_delta2_m32[i] %g \n", T_delta2_m32[i],chi_delta2_m32[i]);
    }
    fprintf(gnuplotPipe11, "e\n");

    fflush(gnuplotPipe11);

        FILE *gnuplotPipe12 = popen("gnuplot -persist", "w");
    if (gnuplotPipe12 == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture de Gnuplot.\n");
        return 1;
    }
    fprintf(gnuplotPipe12, "set title 'Juno 6 years data datking test'\n");
    fprintf(gnuplotPipe12, "set xlabel 'Deposited energy [MeV]'\n");
    fprintf(gnuplotPipe12, "set ylabel 'Events number'\n");
    fprintf(gnuplotPipe12, "set xrange [1:12]\n");

    fprintf(gnuplotPipe12, "plot '-' using 1:2 with linespoints title '_teta_12=50'\n");

    for (int i = 0; i < size; i++) {    
        fprintf(gnuplotPipe12, "%g %g %g\n", tableau_energy[i]-0.8, Ti_spectre_final_teta_12[50][i]);
        //mauvaise forme de spectre après convolution pour teta12
        //manque de l'oscillation teta12
    }
    
    
    fprintf(gnuplotPipe12, "e\n");
    fflush(gnuplotPipe12);
    // Attente de l'utilisateur avant de fermer la fenêtre Gnuplot
    printf("Appuyez sur Entrée pour fermer le graphique...\n");
    getchar();

    // Fermeture du pipe Gnuplot
    //fclose(gnuplotPipe1);
    //fclose(gnuplotPipe2);
    fclose(gnuplotPipe3);
    fclose(gnuplotPipe4);
    fclose(gnuplotPipe5);
    fclose(gnuplotPipe6);
    fclose(gnuplotPipe7);
    fclose(gnuplotPipe8);
    fclose(gnuplotPipe9);
    fclose(gnuplotPipe10);
    fclose(gnuplotPipe11);
    fclose(gnuplotPipe12);

   
    printf("Le minimum de chi_delta2_m31_NH est %g \n", trouverMinimum(chi_delta2_m31_NH,N));
    
    
    double tab_min[50];
    double tableau_longueur[50];

    for(int i=0;i<50;i++){
        tableau_longueur[i]=i*3e3;
    }


//    double Ti_spectre_final_delta2_m31_IH[N][size];
//    double Ti_spectre_final_delta2_m31_NH[N][size];

for(int w=0;w<50;w++){

    for(int u=0;u<N;u++){
    for(int i=0;i<size;i++){//on peut se concentrer sur 1 plot
        Ti_spectre_final_delta2_m31_IH[u][i]=calcul_spectre_lenght(flux_total[i],tableau_energy[i], T_sin2_teta_13[50], T_sin2_teta_12[50], T_delta2_m21[50], T_delta2_m32[50], 'I',T_delta2_m31_IH[u],tableau_longueur[w] )*Np;
        Ti_spectre_final_delta2_m31_NH[u][i]=calcul_spectre_lenght(flux_total[i],tableau_energy[i], T_sin2_teta_13[50], T_sin2_teta_12[50], T_delta2_m21[50], T_delta2_m32[50], 'N',T_delta2_m31_NH[u],tableau_longueur[w])*Np;

    }}

//    double Ti_new_spectre_delta2_m31_IH[N][size];
//    double Ti_new_spectre_delta2_m31_NH[N][size];

    for(int k=0;k<N;k++){
    for(int u=0;u<size;u++){

        Ti_new_spectre_delta2_m31_IH[k][u]=0;
        Ti_new_spectre_delta2_m31_NH[k][u]=0;

        for(int i = 0; i < size; i++){

            double product_delta2_m31_IH=gauss_pdf(tableau_energy[u],tableau_energy[i])*Ti_spectre_final_delta2_m31_IH[k][i];
            double product_delta2_m31_NH=gauss_pdf(tableau_energy[u],tableau_energy[i])*Ti_spectre_final_delta2_m31_NH[k][i];

            if(product_delta2_m31_IH>0||product_delta2_m31_NH>0){
            Ti_new_spectre_delta2_m31_IH[k][u] +=product_delta2_m31_IH;
            Ti_new_spectre_delta2_m31_NH[k][u]+=product_delta2_m31_NH;
            }
           }
           //printf("Ti_new_spectre_teta_13 vaut %g pour k vaut %d et u vaut %d \n", Ti_new_spectre_teta_13[k][u],k,u);
           }
    }

//    double chi_delta2_m31_IH[N];
//    double chi_delta2_m31_NH[N];

    for(int u=0;u<N;u++){

        chi_delta2_m31_IH[u]=0;
        chi_delta2_m31_NH[u]=0;

        for(int i=0;i<size;i++){
            if(tableau_energy[i]>1.8&&tableau_energy[i]<8.0){//cf article 1 p43
                chi_delta2_m31_IH[u]+=pow(Ti_new_spectre_delta2_m31_IH[50][i]-Ti_new_spectre_delta2_m31_IH[u][i],2)/(Ti_new_spectre_delta2_m31_IH[50][i]);
                chi_delta2_m31_NH[u]+=pow(Ti_new_spectre_delta2_m31_IH[50][i]-Ti_new_spectre_delta2_m31_NH[u][i],2)/(Ti_new_spectre_delta2_m31_IH[50][i]);
                }
        }
        

    } 
    tab_min[w]=trouverMinimum(chi_delta2_m31_NH,N);
    printf("le min vaut %g pour i w vaut %d \n", tab_min[w],w);

}


    return 0;
}