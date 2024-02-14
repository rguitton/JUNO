/*******************************************
plotter.C

Macro ROOT qui sert à tracer le spectre en fonction de l'énergie du neutrino. 
Prend en entrée les fichiers textes :

	spectra_initial.txt
	spectra_NH.txt
	spectra_IH.txt
	spectra_p21.txt

générés par la fonction main.cpp

Ces fichiers sont écrits dans le format %lf %lf  %lf %lf.

*******************************************/

#include "../include/fonctions.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <filesystem>
#include <string>

using namespace std;
// Please specify the directory by (un)commenting one of the following line
const TString dir = "../results/";


void plotter_spectra(){

	/*
	std::filesystem::path currentPath = std::filesystem::current_path();
	const TString dir = TString::Format("%s/results/",currentPath.c_str());
	cout << "Dans plotter.C \n" << "\t pwd : " << dir.Data() << endl;
	*/
	TCanvas *c = new TCanvas("spectra","spectra",1366,768);
	c->SetTickx();
	c->SetTicky();
	// spectra_initial
	TGraphErrors *g_initial = new TGraphErrors("../results/spectra_initial.txt");
	g_initial->SetTitle("Reactor antineutrino spectrum;L/E (km/MeV);Arbitrary Unit");
	g_initial->SetMarkerStyle(1);
	g_initial->SetLineColor(kBlack);
	g_initial->SetLineStyle(2);
	g_initial->GetXaxis()->SetTitleSize(0.045);
	g_initial->GetYaxis()->SetTitleSize(0.05);
	g_initial->GetXaxis()->SetTickSize(0.02);
	g_initial->GetYaxis()->SetTickSize(0.015);
	g_initial->GetXaxis()->SetDecimals();
	g_initial->Draw("APL");
	// spectra_p21
	TGraphErrors *g_p21 = new TGraphErrors("../results/spectra_p21.txt");
	//g_p21->SetTitle(";L/E (km/MeV);Arbitrary Unit");
	g_p21->SetMarkerStyle(1);
	g_p21->SetLineColor(kBlack);
	g_p21->SetLineStyle(1);
	g_p21->GetXaxis()->SetTitleSize(0.045);
	g_p21->GetYaxis()->SetTitleSize(0.05);
	g_p21->GetXaxis()->SetTickSize(0.02);
	g_p21->GetYaxis()->SetTickSize(0.015);
	g_p21->GetXaxis()->SetDecimals();
	g_p21->Draw("SAME");
	// spectra_NH
	TGraphErrors *g_NH = new TGraphErrors("../results/spectra_NH.txt");
	//g_NH->SetTitle(";L/E (km/MeV);Arbitrary Unit");
	g_NH->SetMarkerStyle(1);
	g_NH->SetLineColor(kBlue);
	g_NH->SetLineStyle(1);
	g_NH->GetXaxis()->SetTitleSize(0.045);
	g_NH->GetYaxis()->SetTitleSize(0.05);
	g_NH->GetXaxis()->SetTickSize(0.02);
	g_NH->GetYaxis()->SetTickSize(0.015);
	g_NH->GetXaxis()->SetDecimals();
	g_NH->Draw("SAME");
	// spectra_IH
	TGraphErrors *g_IH = new TGraphErrors("../results/spectra_IH.txt");
	//g_IH->SetTitle(";L/E (km/MeV);Arbitrary Unit");
	g_IH->SetMarkerStyle(1);
	g_IH->SetLineColor(kRed);
	g_IH->SetLineStyle(1);
	g_IH->GetXaxis()->SetTitleSize(0.045);
	g_IH->GetYaxis()->SetTitleSize(0.05);
	g_IH->GetXaxis()->SetTickSize(0.02);
	g_IH->GetYaxis()->SetTickSize(0.015);
	g_IH->GetXaxis()->SetDecimals();
	g_IH->Draw("SAME");
	// legend
	auto legend = new TLegend(0.6,0.6,0.8,0.8);
	//legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	legend->AddEntry(g_initial,"Non oscillation","lpf");
	legend->AddEntry(g_p21,"#theta_{12} oscillation","lpf");
	legend->AddEntry(g_NH,"Normal hierarchy","lpf");
	legend->AddEntry(g_IH,"Inverted hierarchy","lpf");
	legend->Draw("SAME");


	
	
// Save plots
	c->Print(TString::Format("%sspectra.pdf",dir.Data()));
	//c->Print(TString::Format("%sspectra.png",dir.Data()));
}

