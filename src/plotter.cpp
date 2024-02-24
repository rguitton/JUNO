#include "plotter.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <filesystem>
#include <string>
//ROOT Libraries
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TAxis.h>


using namespace std;
//const TString dir = "../results/";


void plotter_spectra(int size, const double* xaxis, const double** all_yaxis, const char* dir){
    
	TCanvas *c = new TCanvas("spectra","spectra",1366,768);
	c->SetTickx();
	c->SetTicky();
	// spectra_initial
	TGraphErrors *g_initial = new TGraphErrors(size,xaxis,all_yaxis[0]);
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
	TGraphErrors *g_p21 = new TGraphErrors(size,xaxis,all_yaxis[1]);
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
	TGraphErrors *g_NH = new TGraphErrors(size,xaxis,all_yaxis[2]);
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
	TGraphErrors *g_IH = new TGraphErrors(size,xaxis,all_yaxis[3]);
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
	c->Print(TString::Format("%sspectra.pdf",dir));
	//c->Print(TString::Format("%sspectra.png",dir.Data()));
}


void plotter_flux(int size, const double* xaxis, const double** all_yaxis, const char* dir){
	
}