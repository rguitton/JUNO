#include "fonctions.h"
#include "plotter.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cstdlib>
#include <fstream>

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

using namespace std;

const char* dir = "../output/";

int main()
{
    //Définition des variables 

    int size=1000;
    int N=100;

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

    double L=53;
    double energy;

    double delta2_m21=7.5e-5;//valeur absolue
    double sin2_teta_12=0.308;
    double teta_12=asin(sqrt(sin2_teta_12));

    //Normal ordering mass (Table 1.1 article [1])
    double sin2_teta_13_NH=0.0234;   

    double teta_13_NH=asin(sqrt(sin2_teta_13_NH));

    double delta2_m31_NH=2.47e-3;
    double delta2_m32_NH=delta2_m31_NH-delta2_m21;

    //Inverted ordering mass (Table 1.1 article [1])
    double sin2_teta_13_IH=0.0240;   

    double teta_13_IH=asin(sqrt(sin2_teta_13_IH));

    double delta2_m31_IH=2.42e-3;
    double delta2_m32_IH=delta2_m31_IH+delta2_m21;

    double Np=1.44e33; //Nombre de protons disponible pour la détection

    double T_delta2_m21[N];
    double T_delta2_m32_IH[N];
    double T_delta2_m32_NH[N];

    double T_delta2_m31_NH[N];
    double T_delta2_m31_IH[N];
    double T_delta2_m31_moy[N];

    double T_sin2_teta_12[N];
        
    double T_sin2_teta_13_NH[N];
    double T_sin2_teta_13_IH[N];

    for(int j=0;j<N;j++){
        T_delta2_m21[j]=delta2_m21+1e-5*(j-(N/2))/1000;

        T_delta2_m32_IH[j]=delta2_m32_IH+1e-3*(j-(N/2))/1000;
        T_delta2_m32_NH[j]=delta2_m32_NH+1e-3*(j-(N/2))/1000;

        T_delta2_m31_NH[j]=delta2_m31_NH+1e-3*(j-(N/2))/1000;
        T_delta2_m31_IH[j]=delta2_m31_IH+1e-3*(j-(N/2))/1000;

        T_sin2_teta_12[j]=pow(sin(teta_12),2)+1.0*(j-(N/2))/1000;
        T_sin2_teta_13_IH[j]=pow(sin(teta_13_IH),2)+1.0*(j-(N/2))/1000;
        T_sin2_teta_13_NH[j]=pow(sin(teta_13_NH),2)+1.0*(j-(N/2))/1000;
        }

    //definition des tableaux énergie et longueur 
    for (int i=0;i<size;i++){
        tableau_lenght[i]=5+double(i)*28/size;
        tableau_energy[i]=1.8+double(i)*7.2/size;
    }
    //on définit ensuite les tableaux associés à chaque 
     for (int i=0;i<size;i++){
        double energy=L/tableau_lenght[i];
        spectra_initial[i]=flux(energy)*sigma(energy);
        spectra_p21[i]=flux(energy)*sigma(energy)*probability(energy, 0, sin2_teta_13_IH, sin2_teta_12, delta2_m21, delta2_m32_IH,delta2_m31_IH );
        spectra_NH[i]=flux(energy)*sigma(energy)*probability(energy, 1, sin2_teta_13_NH, sin2_teta_12, delta2_m21, delta2_m32_NH,delta2_m31_NH);
        spectra_IH[i]=flux(energy)*sigma(energy)*probability(energy, 1, sin2_teta_13_IH, sin2_teta_12, delta2_m21, delta2_m32_IH,delta2_m31_IH);

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

    //On reproduit la figure 2.4 de l'article [1] 
    const double* all_yaxis1[4]={spectra_initial,spectra_p21,spectra_NH,spectra_IH};
    plotter_spectra(size,tableau_lenght,all_yaxis1,dir);

    //On reproduit la figure 2.6 de l'article [1]
    const double* all_yaxis2[6]={tab_flux_U_235,tab_flux_U_238,tab_flux_PU_239,tab_flux_PU_241,tab_cross,tab_flux_cross_total};
    plotter_flux(size, tableau_energy, all_yaxis2, dir);

    //on traces les données de l'article récuperées graphiquement pour pouvoir les comparer

        //fprintf(gnuplotPipe2, "set title 'Donnée de l article'\n");
        //fprintf(gnuplotPipe2, "set datafile separator ','\n plot '../data/fig_2_6.txt' using 1:2 with points axes x1y1 title 'U235', '' using 3:4 with points axes x1y1 title 'Pu239', '' using 5:6 with points axes x1y1 title 'U238', '' using 7:8 with points axes x1y1 title 'Pu241'\n");

        //fprintf(gnuplotPipe2, "unset multiplot\n");
        //fflush(gnuplotPipe2);


    //On veut obtenir le spectre des neutrinos détectés en fonction de leur énergie réelle sans prendre en compte l'incertitude sur l'énergie visible

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
        spectre_final[i]=calcul_spectre(flux_total[i],tableau_energy[i], sin2_teta_13_IH, sin2_teta_12, delta2_m21, delta2_m32_IH,delta2_m31_IH);
    }

    //On reproduit la figure 4 de l'article [2] avec en abscisse l'énergie de l'antineutrino
    double* tmp = new double[size];
	for (int i=0;i<size;i++){
		tmp[i]= spectre_final[i]*Np*3600*24;
	}
    plotter_energy_spectrum(size,tableau_energy,tmp,dir);
    delete[] tmp;

    //Cette fois on cherche à déterminer le nombre de neutrinos détectables par jour par JUNO 

    double IBD_initial=0; //par seconde
    double h=tableau_energy[10]-tableau_energy[9];
    for(int i=0;i<size-1;i++){
        if(tableau_energy[i]>1.8&&tableau_energy[i]<12){//cf article 2
            IBD_initial += h*(spectre_final[i] + spectre_final[i+1])/2;
            }
        
        }
    
    double IBD_detected_per_d=Np*IBD_initial*3600*24; //nombre de neutrinons détectés par jour
    
    //On peut donc estimer le nombre d'anti-neutrinos détectés par jour en fonction des différents critères de sélection, on retrouve bien la Table 2.1 de l'article [1]
    
    printf("On peut s'attendre à %g détections par jour sans sélection \n", IBD_detected_per_d);
    printf("On peut s'attendre à %g détections par jour après la sélectivité Fiducial volume (91%) \n", IBD_detected_per_d*0.91);
    printf("On peut s'attendre à %g détections par jour après la sélectivité Energy cut (97.8%) \n", IBD_detected_per_d*0.91*0.978);
    printf("On peut s'attendre à %g détections par jour après la sélectivité Time cut (99.1%) \n", IBD_detected_per_d*0.91*0.978*0.991);
    printf("On peut s'attendre à %g détections par jour après la sélectivité Vertex cut (98.7%) \n", IBD_detected_per_d*0.91*0.978*0.991*0.987);
    printf("On peut s'attendre à %g détections par jour après la sélectivité Muon veto (83%) \n", IBD_detected_per_d*0.91*0.978*0.991*0.987*0.83);
    printf("On peut s'attendre à %g détections par jour après la sélectivité Combined (73%) \n", IBD_detected_per_d*0.73);

    //Pour pouvoir tracé le spectre des anti-neutrinos en fonction de l'énergie vissible nous devons considérer l'incertitude liée à celle-ci (Equation (9) article [2])

    double new_spectre[size];
    for(int u=0;u<size;u++){
        new_spectre[u]=0;
        for(int i = 0; i < size; i++){
           double product=gauss_pdf(tableau_energy[u],tableau_energy[i])*spectre_final[i]; //on convolue par une gaussienne 
           if(product>0){
                new_spectre[u] +=product;
           }
           }}

    //on trace le spectre des anti-neutrinons en fonction de l'énergie visible  (Figure 4 article 2)

    double* tmpx = new double[size];
    double* tmpy = new double[size];
	for (int i=0;i<size;i++){
        tmpx[i]=tableau_energy[i]-0.8;
		tmpy[i]= new_spectre[i]*Np*3600*24/10000;
	}
    plotter_visible_energy_spectrum(size,tmpx,tmpy,dir);

    //On trace le même spectre mais en considérant une periode de 6 ans 

	for (int i=0;i<size;i++){
        tmpx[i]=tableau_energy[i]-0.8;
		tmpy[i]= Np*new_spectre[i]*3600*24*365*6/10000;
	}
    plotter_visible_energy_spectrum_six_years(size,tmpx,tmpy,dir);
    delete[] tmpx; delete[] tmpy;


    

    //recherche du 83

    //Dans cette partie nous allons essayer de déterminer les fonctions chi associés aux incertitudes de chaques paramètres 
    //On définit des matrices, la ligne correspond à la variation du paramètre associé et la colonne représente le spectre associé à chaque paramètre
    //L'indice de ligne 0 correspond au paramètre-sigma et l'indice N-1 au paramètre+sigma 

    double Ti_spectre_final_teta_13[N][size];
    double Ti_spectre_final_teta_12[N][size];
    double Ti_spectre_final_delta2_m21[N][size];
    double Ti_spectre_final_delta2_m31_IH[N][size];
    double Ti_spectre_final_delta2_m31_NH[N][size];

    for(int u=0;u<N;u++){
    for(int i=0;i<size;i++){//on peut se concentrer sur 1 plot
        Ti_spectre_final_teta_13[u][i]=calcul_spectre(flux_total[i],tableau_energy[i], T_sin2_teta_13_IH[u], T_sin2_teta_12[50], T_delta2_m21[50], T_delta2_m32_IH[50],T_delta2_m31_IH[50])*Np;
        Ti_spectre_final_teta_12[u][i]=calcul_spectre(flux_total[i],tableau_energy[i], T_sin2_teta_13_IH[50], T_sin2_teta_12[u], T_delta2_m21[50], T_delta2_m32_IH[50],T_delta2_m31_IH[50])*Np;

        Ti_spectre_final_delta2_m21[u][i]=calcul_spectre(flux_total[i],tableau_energy[i], T_sin2_teta_13_IH[50], T_sin2_teta_12[50], T_delta2_m21[u], T_delta2_m32_IH[50],T_delta2_m31_IH[50])*Np;

        Ti_spectre_final_delta2_m31_NH[u][i]=calcul_spectre(flux_total[i],tableau_energy[i], T_sin2_teta_13_NH[50], T_sin2_teta_12[50], T_delta2_m21[50], T_delta2_m32_NH[50], T_delta2_m31_NH[u])*Np;
        Ti_spectre_final_delta2_m31_IH[u][i]=calcul_spectre(flux_total[i],tableau_energy[i], T_sin2_teta_13_IH[50], T_sin2_teta_12[50], T_delta2_m21[50], T_delta2_m32_IH[50], T_delta2_m31_IH[u])*Np;

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
            if(tableau_energy[i]>1.8&&tableau_energy[i]<8.0){
                //On utilise la formule de Pearson pour le calcul de chi (équation (2.9) article [1])
                chi_teta13[u]+=pow(Ti_new_spectre_teta_13[50][i]-Ti_new_spectre_teta_13[u][i],2)/(Ti_new_spectre_teta_13[50][i]); //l'indice de ligne 50 correspond à la valeur moyenne du paramètre 
                chi_teta12[u]+=pow(Ti_new_spectre_teta_12[50][i]-Ti_new_spectre_teta_12[u][i],2)/(Ti_new_spectre_teta_12[50][i]);
                chi_delta2_m21[u]+=pow(Ti_new_spectre_delta2_m21[50][i]-Ti_new_spectre_delta2_m21[u][i],2)/(Ti_new_spectre_delta2_m21[50][i]);
                chi_delta2_m31_IH[u]+=pow(Ti_new_spectre_delta2_m31_IH[50][i]-Ti_new_spectre_delta2_m31_IH[u][i],2)/(Ti_new_spectre_delta2_m31_IH[50][i]);
                chi_delta2_m31_NH[u]+=pow(Ti_new_spectre_delta2_m31_IH[50][i]-Ti_new_spectre_delta2_m31_NH[u][i],2)/(Ti_new_spectre_delta2_m31_IH[50][i]);
                
            }
        }
    }    
    //On peut ensuite tracer chacune des fonctions chi pour chaques paramètres 

    //Pour SIN^2(θ_{13})
 
    FILE *gnuplotPipe6 = popen("gnuplot -persist", "w");
    if (gnuplotPipe6 == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture de Gnuplot.\n");
        return 1;
    }
    fprintf(gnuplotPipe6, "set multiplot layout 2, 2 title 'Variation de Δχ^2'\n");

    // Affichage de SIN^2(θ_{13})
    fprintf(gnuplotPipe6, "set xlabel 'SIN^2(θ_{13})'\n");
    fprintf(gnuplotPipe6, "set ylabel 'Δχ^2'\n");
    fprintf(gnuplotPipe6, "plot '-' with points title 'Δχ^2_{SIN^2(θ_{13})}'\n");
    for (int i = 0; i < N; i++) {    
        fprintf(gnuplotPipe6, "%g %g\n", T_sin2_teta_13_IH[i], chi_teta13[i]);
    }
    fprintf(gnuplotPipe6, "e\n");

    // Affichage de SIN^2(θ_{12})
    fprintf(gnuplotPipe6, "set xlabel 'SIN^2(θ_{12})'\n");
    fprintf(gnuplotPipe6, "set ylabel 'Δχ^2'\n");
    fprintf(gnuplotPipe6, "plot '-' with points title 'Δχ^2_{SIN^2(θ_{12})}'\n");
    for (int i = 0; i < N; i++) {    
        fprintf(gnuplotPipe6, "%g %g\n", T_sin2_teta_12[i], chi_teta12[i]);
    }
    fprintf(gnuplotPipe6, "e\n");

    // Affichage de Δm^2_{21}
    fprintf(gnuplotPipe6, "set xlabel 'Δm^2_{21}'\n");
    fprintf(gnuplotPipe6, "set ylabel 'Δχ^2'\n");
    fprintf(gnuplotPipe6, "plot '-' with points title 'Δχ^2_{Δm^2_{21}}'\n");
    for (int i = 0; i < N; i++) {    
        fprintf(gnuplotPipe6, "%g %g\n", T_delta2_m21[i], chi_delta2_m21[i]);
    }
    fprintf(gnuplotPipe6, "e\n");

    // Affichage de Δm^2_{31}
    fprintf(gnuplotPipe6, "set xlabel 'Δm^2_{31}'\n");
    fprintf(gnuplotPipe6, "set ylabel 'Δχ^2'\n");
    fprintf(gnuplotPipe6, "plot '-' with points title 'Δχ^2_{Δm^2_{31}}'\n");
    for (int i = 0; i < N; i++) {   
        fprintf(gnuplotPipe6, "%g %g\n",  T_delta2_m31_IH[i],chi_delta2_m31_IH[i]);
    }
    fprintf(gnuplotPipe6, "e\n");

    fprintf(gnuplotPipe6, "unset multiplot\n");
    fflush(gnuplotPipe6);

    //#######
    const double* tab_xaxis[4] ={T_sin2_teta_13_IH, T_sin2_teta_12, T_delta2_m21, T_delta2_m31_IH};
    const double* tab_yaxis[4] ={chi_teta13, chi_teta12,chi_delta2_m21, chi_delta2_m31_IH};
    plotter_chi2(N,2,2,tab_xaxis,tab_yaxis,dir);




    //On cherche ensuite à comparer la valeur de Δχ² en fonction de si l'on considère une certaine hierarchie de masse et si on en observe une autre:
    //On reproduit la figure 2.8 de l'article [1]

    FILE *gnuplotPipe7 = popen("gnuplot -persist", "w");
    if (gnuplotPipe7 == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture de Gnuplot.\n");
        return 1;
    }

    fprintf(gnuplotPipe7, "set title 'Comparaison des Δχ^2 '\n");
    fprintf(gnuplotPipe7, "set xlabel 'Δm^2_{31}'\n");

    fprintf(gnuplotPipe7, "set ylabel 'Δχ^2'\n");

    fprintf(gnuplotPipe7, "plot '-' with points title 'Δχ^2_{Δm^2_{31}} False MH', '-' with points title 'Δχ^2_{Δm^2_{31}} True MH'\n");
    
    for (int i = 0; i < N; i++) {   
        fprintf(gnuplotPipe7, "%g %g\n",  T_delta2_m31_IH[i],chi_delta2_m31_NH[i]);
    }
    fprintf(gnuplotPipe7, "e\n");

    for (int i = 0; i < N; i++) {   
        fprintf(gnuplotPipe7, "%g %g\n",  T_delta2_m31_IH[i],chi_delta2_m31_IH[i]);
    }
    fprintf(gnuplotPipe7, "e\n");
    fflush(gnuplotPipe7);

    
    //Pour finir nous allons essayer de retouver la figure 2.7 de l'article [1] qui permet de déterminer la distance optimale pour maximiser la différence des Δχ^2 True MH et False MH 
    //Cette étape est la plus longue en terme de temps de calcul car elle nécessite de recalculer le Δχ^2 à chaque longueur 
 
    double tab_min[50];
    double tableau_longueur[50];

    for(int i=0;i<50;i++){
        tableau_longueur[i]=i*3e3; //on va de 0 à 150km par pas de 3km
    }

for(int w=1;w<50;w++){//boucle sur la distance

    for(int u=0;u<N;u++){
    for(int i=0;i<size;i++){//on peut se concentrer sur 1 plot
 
        Ti_spectre_final_delta2_m31_NH[u][i]=calcul_spectre_lenght(flux_total[i],tableau_energy[i], T_sin2_teta_13_NH[50], T_sin2_teta_12[50], T_sin2_teta_12[50], T_delta2_m32_NH[50], T_delta2_m31_NH[u],tableau_longueur[w] )*Np;
        Ti_spectre_final_delta2_m31_IH[u][i]=calcul_spectre_lenght(flux_total[i],tableau_energy[i], T_sin2_teta_13_NH[50], T_sin2_teta_12[50], T_sin2_teta_12[50], T_delta2_m32_IH[50],T_delta2_m31_IH[u],tableau_longueur[w] )*Np;
    }}

    for(int k=0;k<N;k++){//boucle sur la variation du paramètre
        for(int u=0;u<size;u++){//boucle sur l'énergie

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
        }
    }
    for(int u=0;u<N;u++){

        //chi_delta2_m31_IH[u]=0;
        chi_delta2_m31_NH[u]=0;

        for(int i=0;i<size;i++){
            if(tableau_energy[i]>1.8&&tableau_energy[i]<8.0){//cf article 1 p43
                //chi_delta2_m31_IH[u]+=pow(Ti_new_spectre_delta2_m31_IH[50][i]-Ti_new_spectre_delta2_m31_IH[u][i],2)/(Ti_new_spectre_delta2_m31_IH[50][i]);
                chi_delta2_m31_NH[u]+=pow(Ti_new_spectre_delta2_m31_IH[50][i]-Ti_new_spectre_delta2_m31_NH[u][i],2)/(Ti_new_spectre_delta2_m31_IH[50][i]);
                }
        }
        

    } 
    tab_min[w]=trouverMinimum(chi_delta2_m31_NH,N); //-trouverMinimum(chi_delta2_m31_IH,N);
    //tab_min[w]=chi_delta2_m31_NH[N/2];
    printf("le min vaut %g pour L vaut %d km \n", tab_min[w],3*w);

}
    FILE *gnuplotPipe8 = popen("gnuplot -persist", "w");
    if (gnuplotPipe8 == NULL) {
        fprintf(stderr, "Erreur lors de l'ouverture de Gnuplot.\n");
        return 1;
    }

    fprintf(gnuplotPipe8, "set title 'Δχ^2 en fonction de la distance'\n");
    fprintf(gnuplotPipe8, "set xlabel 'Distance des réacteurs (km)'\n");

    fprintf(gnuplotPipe8, "set ylabel 'Δχ^2'\n");

    fprintf(gnuplotPipe8, "plot '-' with points title 'Δχ^2 en fonction de la distance'\n");
    
    for (int i = 0; i < 50; i++) {   
        fprintf(gnuplotPipe8, "%g %g\n",  tableau_longueur[i]*1e-3,tab_min[i]);
    }

    fprintf(gnuplotPipe8, "e\n");
    fflush(gnuplotPipe8);

    // Attente de l'utilisateur avant de fermer la fenêtre Gnuplot
    printf("Appuyez sur Entrée pour fermer le graphique...\n");
    getchar();

    // Fermeture du pipe Gnuplot


    fclose(gnuplotPipe6);
    fclose(gnuplotPipe7);
    fclose(gnuplotPipe8);
    return 0;
}