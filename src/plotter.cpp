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
#include <TGaxis.h>
#include <TText.h>
#include <TH1F.h>


using namespace std;
//const TString dir = "../results/";


void plotter_spectra(int size, const double* xaxis, const double** all_yaxis, const char* dir){
    
	TCanvas *c = new TCanvas("spectra","spectra",1366,768);
	c->SetTickx();
	c->SetTicky(); 

	// spectra_initial
	//TGraphErrors *g_initial = new TGraphErrors(size,xaxis,all_yaxis[0]);
	TGraphErrors g_tmp0;
	fillGraphIgnoringNaN(g_tmp0,xaxis,all_yaxis[0],0,0,size);
	TGraphErrors *g_initial = new TGraphErrors(g_tmp0);
	g_initial->SetTitle("Reactor antineutrino spectrum;L/E (km/MeV);Arbitrary Unit");
	g_initial->SetMarkerStyle(1);
	g_initial->SetLineColor(kBlack);
	g_initial->SetLineStyle(2);
	g_initial->GetXaxis()->SetTitleSize(0.045);
	g_initial->GetYaxis()->SetTitleSize(0.05);
	g_initial->GetXaxis()->SetTickSize(0.02);
	g_initial->GetYaxis()->SetTickSize(0.015);
	g_initial->GetXaxis()->SetDecimals();
	g_initial->Draw("APL"); c->Update();
	// spectra_p21
	//TGraphErrors *g_p21 = new TGraphErrors(size,xaxis,all_yaxis[1]);
	//g_p21->SetTitle(";L/E (km/MeV);Arbitrary Unit");
	TGraphErrors g_tmp1;
	fillGraphIgnoringNaN(g_tmp1,xaxis,all_yaxis[1],0,0,size);
	TGraphErrors *g_p21 = new TGraphErrors(g_tmp1);
	g_p21->SetMarkerStyle(1);
	g_p21->SetLineColor(kBlack);
	g_p21->SetLineStyle(1);
	g_p21->GetXaxis()->SetTitleSize(0.045);
	g_p21->GetYaxis()->SetTitleSize(0.05);
	g_p21->GetXaxis()->SetTickSize(0.02);
	g_p21->GetYaxis()->SetTickSize(0.015);
	g_p21->GetXaxis()->SetDecimals();
	g_p21->Draw("SAME"); c->Update();
	// spectra_NH
	//TGraphErrors *g_NH = new TGraphErrors(size,xaxis,all_yaxis[2]);
	//g_NH->SetTitle(";L/E (km/MeV);Arbitrary Unit");
	TGraphErrors g_tmp2;
	fillGraphIgnoringNaN(g_tmp2,xaxis,all_yaxis[2],0,0,size);
	TGraphErrors *g_NH = new TGraphErrors(g_tmp2);
	g_NH->SetMarkerStyle(1);
	g_NH->SetLineColor(kBlue);
	g_NH->SetLineStyle(1);
	g_NH->GetXaxis()->SetTitleSize(0.045);
	g_NH->GetYaxis()->SetTitleSize(0.05);
	g_NH->GetXaxis()->SetTickSize(0.02);
	g_NH->GetYaxis()->SetTickSize(0.015);
	g_NH->GetXaxis()->SetDecimals();
	g_NH->Draw("SAME"); c->Update();
	// spectra_IH
	//TGraphErrors *g_IH = new TGraphErrors(size,xaxis,all_yaxis[3]);
	//g_IH->SetTitle(";L/E (km/MeV);Arbitrary Unit");
	TGraphErrors g_tmp3;
	fillGraphIgnoringNaN(g_tmp3,xaxis,all_yaxis[3],0,0,size);
	TGraphErrors *g_IH = new TGraphErrors(g_tmp3);
	g_IH->SetMarkerStyle(1);
	g_IH->SetLineColor(kRed);
	g_IH->SetLineStyle(1);
	g_IH->GetXaxis()->SetTitleSize(0.045);
	g_IH->GetYaxis()->SetTitleSize(0.05);
	g_IH->GetXaxis()->SetTickSize(0.02);
	g_IH->GetYaxis()->SetTickSize(0.015);
	g_IH->GetXaxis()->SetDecimals();
	g_IH->Draw("SAME"); c->Update();
	// legend
	auto legend = new TLegend(0.6,0.6,0.8,0.8);
	//legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
	legend->AddEntry(g_initial,"Non oscillation","lpf");
	legend->AddEntry(g_p21,"#theta_{12} oscillation","lpf");
	legend->AddEntry(g_NH,"Normal hierarchy","lpf");
	legend->AddEntry(g_IH,"Inverted hierarchy","lpf");
	legend->Draw("SAME"); c->Update();

// Save plots
	c->Print(TString::Format("%sspectra.pdf",dir));
	//c->Print(TString::Format("%sspectra.png",dir));
	c->Close();
}


void plotter_flux(int size, const double* xaxis, const double** all_yaxis, const char* dir){
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
	//TGraphErrors *g_U_235 = new TGraphErrors(size,xaxis,all_yaxis[0]);
	TGraphErrors g_tmp0;
	fillGraphIgnoringNaN(g_tmp0,xaxis,all_yaxis[0],0,0,size);
	TGraphErrors *g_U_235 = new TGraphErrors(g_tmp0);
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
	TText *t_U_235 = new TText(2.5, 0.48,"antineutrino flux");
	t_U_235->SetTextAlign(22);
	t_U_235->SetTextColor(kBlack);
	t_U_235->SetTextFont(43);
	t_U_235->SetTextSize(25);
	t_U_235->SetTextAngle(-60);
	t_U_235->Draw(); gPad->Update();
	// flux_U_238
	//TGraphErrors *g_U_238 = new TGraphErrors(size,xaxis,all_yaxis[1]);
	//g_U_238->SetTitle("^{238}U ");
	TGraphErrors g_tmp1;
	fillGraphIgnoringNaN(g_tmp1,xaxis,all_yaxis[1],0,0,size);
	TGraphErrors *g_U_238 = new TGraphErrors(g_tmp1);
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
	//TGraphErrors *g_PU_239 = new TGraphErrors(size,xaxis,all_yaxis[2]);
	//g_PU_239->SetTitle("^{239}PU ");
	TGraphErrors g_tmp2;
	fillGraphIgnoringNaN(g_tmp2,xaxis,all_yaxis[2],0,0,size);
	TGraphErrors *g_PU_239 = new TGraphErrors(g_tmp2);
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
	//TGraphErrors *g_PU_241 = new TGraphErrors(size,xaxis,all_yaxis[3]);
	//g_PU_241->SetTitle("^{241}PU ");
	TGraphErrors g_tmp3;
	fillGraphIgnoringNaN(g_tmp3,xaxis,all_yaxis[3],0,0,size);
	TGraphErrors *g_PU_241 = new TGraphErrors(g_tmp3);
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
	legend->AddEntry(g_PU_239,"^{239}PU ","lpf");
	legend->AddEntry(g_U_238,"^{238}U ","lpf");
	legend->AddEntry(g_PU_241,"^{241}PU ","lpf");
	legend->Draw("SAME"); gPad->Update();
	

	Double_t xmin = p1->GetUxmin(); cout << "xmin = " << xmin << endl;
	Double_t xmax = p1->GetUxmax(); cout << "xmax = " << xmax << endl;
	Double_t dx = (xmax - xmin) / 0.8; // 10 percent margins left and right

	//hand written, would be better to get them automatically
	Double_t ymin = 3.63287e-45;
	Double_t ymax = 6e-42;

	Double_t dy = (ymax - ymin) / 0.8; // 10 percent margins top and bottom
	p2->Range(xmin-0.1*dx, ymin-0.1*dy, xmax+0.1*dx, ymax+0.1*dy);
	p2->Draw();
	p2->cd();

	//cross-section
	//TGraphErrors *g_cross= new TGraphErrors(size,xaxis,all_yaxis[4]);
	TGraphErrors g_tmp4;
	fillGraphIgnoringNaN(g_tmp4,xaxis,all_yaxis[4],0,0,size);
	TGraphErrors *g_cross = new TGraphErrors(g_tmp4);
	//g_cross->SetTitle("Antineutrino flux");
	g_cross->SetMarkerStyle(1);
	g_cross->SetLineColor(kBlue);
	g_cross->SetLineStyle(1);
	g_cross->Draw("L"); gPad->Update();
	TText *t_cross = new TText(7.5, 3.2e-42,"cross section");
	t_cross->SetTextAlign(22);
	t_cross->SetTextColor(kBlue);
	t_cross->SetTextFont(43);
	t_cross->SetTextSize(25);
	t_cross->SetTextAngle(45);
	t_cross->Draw(); gPad->Update();


	//spectre : flux*cross
	//TGraphErrors *g_flux_cross_total= new TGraphErrors(size,xaxis,all_yaxis[5]);
	TGraphErrors g_tmp5;
	fillGraphIgnoringNaN(g_tmp5,xaxis,all_yaxis[5],0,0,size);
	TGraphErrors *g_flux_cross_total = new TGraphErrors(g_tmp5);
	g_flux_cross_total->SetMarkerStyle(1);
	g_flux_cross_total->SetLineColor(kRed);
	g_flux_cross_total->SetLineStyle(1);
	g_flux_cross_total->Draw("L"); gPad->Update(); 
	TText *t_flux_cross_total = new TText(4, 4.6e-42,"measured spectrum");
	t_flux_cross_total->SetTextAlign(22);
	t_flux_cross_total->SetTextColor(kRed);
	t_flux_cross_total->SetTextFont(43);
	t_flux_cross_total->SetTextSize(25);
	t_flux_cross_total->SetTextAngle(0);
	t_flux_cross_total->Draw(); gPad->Update();
	


	TGaxis *axis = new TGaxis(xmax, ymin, xmax, ymax, ymin, ymax, 510, "+L");
	axis->SetTitle("Inverse Beta Decay Cross Section (cm^{2})");
	axis->SetTitleColor(kBlue);
	axis->SetLineColor(kBlue);
	axis->SetLabelColor(kBlue);
	axis->Draw();
	gPad->Update();

	

	c->cd();

// Save plots
	c->Print(TString::Format("%sflux.pdf",dir));
	//c->Print(TString::Format("%sflux.png",dir));
	c->Close();
}

void plotter_visible_energy_spectrum(int size, const double* xaxis, const double* yaxis, const char* dir){
	
	TCanvas *c = new TCanvas("visible_energy_spectrum","visible_energy_spectrum",1366,768);
	c->SetTickx();
	c->SetTicky();

	//TGraphErrors *g = new TGraphErrors(size,xaxis,yaxis);
	TGraphErrors g_tmp;
	fillGraphIgnoringNaN(g_tmp,xaxis,yaxis,0,0,size);
	TGraphErrors *g = new TGraphErrors(g_tmp);
	g->SetTitle(";Visible Energy [MeV];Events/0.02 [MeV^{-1}day^{-1}]");
	g->SetMarkerStyle(1);
	g->SetLineColor(kBlack);
	g->SetLineStyle(1);
	g->GetXaxis()->SetTitleSize(0.045);
	g->GetYaxis()->SetTitleSize(0.05);
	g->GetXaxis()->SetTickSize(0.02);
	g->GetYaxis()->SetTickSize(0.015);
	g->GetXaxis()->SetDecimals();
	g->Draw("APL"); gPad->Update();

	

	c->Print(TString::Format("%svisible_energy_spectrum.pdf",dir));
	//c->Print(TString::Format("%svisible_energy_spectrum.png",dir));
	c->Close();
}

void fillGraphIgnoringNaN(TGraphErrors &graph, const double *x, const double *y, const double *ex, const double *ey, int nPoints) {
    for (int i = 0; i < nPoints; ++i) {
        // Check if both x and y values are not NaN
        if (!std::isnan(x[i]) && !std::isnan(y[i])) {
            graph.SetPoint(graph.GetN(), x[i], y[i]);
            //graph.SetPointError(graph.GetN() - 1, ex[i], ey[i]);
			graph.SetPointError(graph.GetN() - 1, 0,0);
        }
    }
}