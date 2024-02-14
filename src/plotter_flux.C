/*******************************************
plotter.C

Macro ROOT qui sert à tracer le spectre en fonction de l'énergie du neutrino. 
Prend en entrée les fichiers textes :

	flux_U_235.txt
	flux_NH.txt
	flux_IH.txt
	flux_U_238.txt

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


void plotter_flux(){

	/*
	std::filesystem::path currentPath = std::filesystem::current_path();
	const TString dir = TString::Format("%s/results/",currentPath.c_str());
	cout << "Dans plotter.C \n" << "\t pwd : " << dir.Data() << endl;
	*/

	TCanvas *c = new TCanvas("flux","flux",1366,768);
	c->SetTickx();
	c->SetTicky();

	TPad *p1 = new TPad("p1", "", 0, 0, 1, 1);
	//p1->SetGrid();
	TPad *p2 = new TPad("p2", "", 0, 0, 1, 1);
	p2->SetFillStyle(4000); // will be transparent
	// p2->SetGrid();
	//TPad *p0 = new TPad("p0", "", 0, 0, 1, 1);

	p1->Draw();
	p1->cd();

	// flux_U_235
	TGraphErrors *g_U_235 = new TGraphErrors("../results/flux_U_235.txt");
	g_U_235->SetTitle("Reactor antineutrino flux;Antineutrino Energy (MeV);Antineutrinos / MeV / Fission");
	g_U_235->SetMarkerStyle(1);
	g_U_235->SetLineColor(kBlack);
	g_U_235->SetLineStyle(1);
	g_U_235->GetXaxis()->SetTitleSize(0.045);
	g_U_235->GetYaxis()->SetTitleSize(0.05);
	g_U_235->GetXaxis()->SetTickSize(0.02);
	g_U_235->GetYaxis()->SetTickSize(0.015);
	g_U_235->GetXaxis()->SetDecimals();
	g_U_235->GetYaxis()->SetNdivisions(510, kTRUE);
	g_U_235->Draw("APL"); gPad->Update();
	// flux_U_238
	TGraphErrors *g_U_238 = new TGraphErrors("../results/flux_U_238.txt");
	//g_U_238->SetTitle("^{238}U ");
	g_U_238->SetMarkerStyle(1);
	g_U_238->SetLineColor(kGreen);
	g_U_238->SetLineStyle(1);
	g_U_238->GetXaxis()->SetTitleSize(0.045);
	g_U_238->GetYaxis()->SetTitleSize(0.05);
	g_U_238->GetXaxis()->SetTickSize(0.02);
	g_U_238->GetYaxis()->SetTickSize(0.015);
	g_U_238->GetXaxis()->SetDecimals();
	g_U_238->Draw("SAME"); gPad->Update();
	// flux_NH
	TGraphErrors *g_PU_239 = new TGraphErrors("../results/flux_PU_239.txt");
	//g_PU_239->SetTitle("^{239}PU ");
	g_PU_239->SetMarkerStyle(1);
	g_PU_239->SetLineColor(kMagenta);
	g_PU_239->SetLineStyle(1);
	g_PU_239->GetXaxis()->SetTitleSize(0.045);
	g_PU_239->GetYaxis()->SetTitleSize(0.05);
	g_PU_239->GetXaxis()->SetTickSize(0.02);
	g_PU_239->GetYaxis()->SetTickSize(0.015);
	g_PU_239->GetXaxis()->SetDecimals();
	g_PU_239->Draw("SAME"); gPad->Update();
	// flux_IH
	TGraphErrors *g_PU_241 = new TGraphErrors("../results/flux_PU_241.txt");
	//g_PU_241->SetTitle("^{241}PU ");
	g_PU_241->SetMarkerStyle(1);
	g_PU_241->SetLineColor(kOrange);
	g_PU_241->SetLineStyle(1);
	g_PU_241->GetXaxis()->SetTitleSize(0.045);
	g_PU_241->GetYaxis()->SetTitleSize(0.05);
	g_PU_241->GetXaxis()->SetTickSize(0.02);
	g_PU_241->GetYaxis()->SetTickSize(0.015);
	g_PU_241->GetXaxis()->SetDecimals();
	g_PU_241->Draw("SAME"); gPad->Update();
	// legend
	auto legend = new TLegend(0.11,0.15,0.2,0.4);
	legend->SetFillColor(0);
	//legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	legend->AddEntry(g_U_235,"^{235}U ","lpf");
	legend->AddEntry(g_U_238,"^{238}U ","lpf");
	legend->AddEntry(g_PU_239,"^{239}PU ","lpf");
	legend->AddEntry(g_PU_241,"^{241}PU ","lpf");
	legend->Draw("SAME"); gPad->Update();
	
	/*
	//cross-section
	TGraphErrors *g_cross= new TGraphErrors("../results/cross.txt");
	//g_cross->SetTitle("Antineutrino flux");
	g_cross->SetMarkerStyle(1);
	g_cross->SetLineColor(kBlue);
	g_cross->SetLineStyle(1);
	//g_cross->Draw("APL");
	*/

	Double_t xmin = p1->GetUxmin(); cout << "xmin = " << xmin << endl;
	Double_t xmax = p1->GetUxmax(); cout << "xmax = " << xmax << endl;
	Double_t dx = (xmax - xmin) / 0.8; // 10 percent margins left and right

	//p0->Draw();
	//p0->cd();
	//g_cross->Draw("SAME");
	Double_t ymin = 3.63287e-45;
	//Double_t ymin = g_cross->GetYaxis()->GetXmin(); cout << "ymin = " << ymin << endl;
	//Double_t ymin = p2->GetUymin(); cout << "ymin = " << ymin << endl;
	Double_t ymax = 6e-42;
	//Double_t ymax = g_cross->GetYaxis()->GetXmax(); cout << "ymax = " << ymax << endl;
	//Double_t ymax = p2->GetUymax(); cout << "ymax = " << ymax << endl;
	Double_t dy = (ymax - ymin) / 0.8; // 10 percent margins top and bottom
	p2->Range(xmin-0.1*dx, ymin-0.1*dy, xmax+0.1*dx, ymax+0.1*dy);
	p2->Draw();
	p2->cd();

	//cross-section
	TGraphErrors *g_cross= new TGraphErrors("../results/cross.txt");
	//g_cross->SetTitle("Antineutrino flux");
	g_cross->SetMarkerStyle(1);
	g_cross->SetLineColor(kBlue);
	g_cross->SetLineStyle(1);
	g_cross->Draw("L");
	gPad->Update();

	//spectre : cross*flux
	TGraphErrors *g_cross_flux= new TGraphErrors("../results/flux_cross_U_235.txt");
	g_cross_flux->SetMarkerStyle(1);
	g_cross_flux->SetLineColor(kRed);
	g_cross_flux->SetLineStyle(1);
	g_cross_flux->Draw("L");
	gPad->Update();
	


	TGaxis *axis = new TGaxis(xmax, ymin, xmax, ymax, ymin, ymax, 510, "+L");
	axis->SetLineColor(kBlue);
	axis->SetLabelColor(kBlue);
	axis->Draw();
	gPad->Update();

	

	c->cd();

// Save plots
	c->Print(TString::Format("%sflux.pdf",dir.Data()));
	//c->Print(TString::Format("%sflux.png",dir.Data()));
}

