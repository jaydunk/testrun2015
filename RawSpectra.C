#ifndef __CINT__
#include <fstream.h>
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TObject.h"
#include "TTree.h"
#include <iostream.h>
#include <fstream.h>
#include <stdio.h>
#include <cmath>
#include <math.h>
#include <string.h>
#endif

const int ECalNChannel = 16;
const int HodNChannel  = 16;

const float HodoscopeCutoff[16] = {250, 250, 250, 250,
							250, 250, 250, 250,
	250, 250, 250, 250, 
	250, 250, 250, 250};

const float Ce1Cutoff = 1600;

RawSpectra() {
	
	gStyle->SetOptStat("nem");
	gStyle->SetStatX(.90);
	gStyle->SetStatY(.88);
	gStyle->SetStatW(.18);
	gStyle->SetStatH(.20);
	
	char fname[200];
	char pname[200];
	int eventno = 0;
	
	TChain *chain = new TChain("T");
	cout << "Enter multiple .root files. Type 'end' when finished."<< endl;
	while(true) {
		cout << "Enter .root file name, or 'end' if done adding files: " << endl;
		cin >> fname;
		if(strcmp(fname,"end")==0) break;
		chain->Add(fname);
	}
	
	cout << "chain entries: " << chain->GetEntries() << endl;
	int nentries = chain->GetEntries();
	
	float PbGThreshold = 120;
	
	TH2D *CerenkovCor = new TH2D("CerenkovCor", "Ce1 Ce2 Ped Sub Correlation", 100, .5, 600.5, 100, .5, 600.5);
	
	TH1D *PbGPedSubtracted = new TH1D("PbGPedSubtracted", "Lead Glass Pedestal Subtracted ADC", 4096, .5, 4096.5);
	TH2D *PbGEXYCor = new TH2D("PbGEXYCor", "PbG Position Energy Correlation", 8, -19.2, 19.2, 8, -19.2, 19.2);
	
	char buffer1[100];
	char buffer2[100];
	
	/****Raw Spectra Histograms****/
	TH1D* HCalRawADCSpectrum[16];
	for(int i=0; i<16; i++) {
		sprintf(buffer1, "HCalADCHist%d", i);
		sprintf(buffer2, "HCal ADC Channel %d", i);
		HCalRawADCSpectrum[i] = new TH1D(buffer1, buffer2, 4096, .5, 4096.5);
	}
	TH1D* ECalRawADCSpectrum[16];
	for(int i=0; i<16; i++) {
		sprintf(buffer1, "ECalADCHist%d", i);
		sprintf(buffer2, "ECal ADC Channel %d", i);
		ECalRawADCSpectrum[i] = new TH1D(buffer1, buffer2, 4096, .5, 4096.5);
	}
	
	TH1D* Sc1RawADCSpectrum = new TH1D("Sc1ADCHist", "Sc1 ADC", 4096, .5, 4096.5);
	TH1D* MonRawADCSpectrum = new TH1D("MonADCHist", "Mon ADC", 4096, .5, 4096.5);
	TH1D* Ce1RawADCSpectrum = new TH1D("Ce1ADCHist", "Ce1 ADC", 4096, .5, 4096.5);
	TH1D* Ce2RawADCSpectrum = new TH1D("Ce2ADCHist", "Ce2 ADC", 4096, .5, 4096.5);
	TH1D* PbGRawADCSpectrum = new TH1D("PbGADCHist", "PbG ADC", 4096, .5, 4096.5);
	
	TH1D* HodRawADCSpectrum[16];
	for(int i=0; i<16; i++) {
		sprintf(buffer1, "HodADCHist%d", i);
		sprintf(buffer2, "Hod ADC Channel %d", i);
		HodRawADCSpectrum[i] = new TH1D(buffer1, buffer2, 4096, .5, 4096.5);
	}
	
	
	TH1D* XMult = new TH1D("XMult", "Multiplicity in X direction", 9, -.5, 8.5);
	TH1D* YMult = new TH1D("YMult", "Multiplicity in Y direction", 9, -.5, 8.5);
	TH1D* BeamProfileX = new TH1D("BeamProfileX", "Beam Profile X", 8, -19.2, 19.2);
	TH1D* BeamProfileY = new TH1D("BeamProfileY", "Beam Profile Y", 8, -19.2, 19.2);
	TH2D* BeamXY = new TH2D("BeamXY", "Beam Profile XY", 8, -19.2, 19.2, 8, -19.2, 19.2);
	
	int eventindex = 0;
	int goodevents = 0;
	
	float ECalADC[16];
	float Sc1ADC;
	float MonADC;
	float Ce1ADC;
	float Ce2ADC;
	float PbGADC;
	float HodADC[16];
	
	bool HighHodHit = false;
	float xposition = 0;
	float yposition = 0;
	
	float PbGPedestalSubtracted = 0;
	float Ce1PedestalSubtracted = 0;
	float Ce2PedestalSubtracted = 0;
	
	float PbGHitsAbovePed = 0;
	float Ce1HitsAbovePed = 0;
	float Ce2HitsAbovePed = 0;
	
	chain->SetBranchAddress("ECalRawADC", ECalADC);
	chain->SetBranchAddress("HodRawADC", HodADC);
	chain->SetBranchAddress("Sc1RawADC", &Sc1ADC);
	chain->SetBranchAddress("MonRawADC", &MonADC);
	chain->SetBranchAddress("Ce1RawADC", &Ce1ADC);
	chain->SetBranchAddress("Ce2RawADC", &Ce2ADC);
	chain->SetBranchAddress("PbGRawADC", &PbGADC);
	
	for (eventindex=0; eventindex<nentries; eventindex++) {
		
		chain->GetEntry(eventindex);
		
		if(eventindex%10000==0) cout << "Processing event: " << eventindex << endl;
		
		/*Calculate Multiplicities*/
		
		
		int xmult = 0;
		int ymult = 0;
		HighHodHit = false;
		for(int i=0; i<8; i++) {
			if(HodADC[i]>=HodoscopeCutoff[i]) {
				ymult++;
			}
			if(HodADC[i]>2275) {
				HighHodHit = true;
			}
		}
		
		for(int i=0; i<8; i++) {
			if(HodADC[8+i]>=HodoscopeCutoff[8+i]) {
				xmult++;
			}
			if(HodADC[8+i]>2275) {
				HighHodHit = true;
			}
		}
		
		
		
		for(int i=0; i<8; i++) {
			if(HodADC[i]>=HodoscopeCutoff[i]) {
				for(int j=0; j<8; j++) {
					if(HodADC[8+j]>=HodoscopeCutoff[8+j]) {
						xposition = 16.8 - 4.8*j;
						yposition = -16.8 + 4.8*i;
					}
				}
			}
		}
		
		XMult->Fill(xmult);
		YMult->Fill(ymult);
		
		/****************Cuts******************/
		//if(!(xmult==1&&ymult==1)) {eventindex++;continue;} //xmult ymult both must be 1
		//if(!(Sc1ADC<1550)) {eventindex++;continue;}
		//if(HighHodHit) {eventindex++;continue;}
		//if(!(Ce1ADC>1470)) {eventindex++;continue;} // cut on cerenkov1 signal
		//if(!(Ce1ADC<1590)) {eventindex++;continue;}
	    //if(Ce2HitsAbovePed>0.0) {eventindex++;continue;} // cut on ce2
		//if(PbGHitsAbovePed > 0.0) {eventindex++; continue;} //cut on lead glass signal
		//if(xposition>5.0||xposition<-5.0||yposition>5.0||yposition<-5.0) {eventindex++;continue;} //include these positions
		//if(!(xposition>5.0||xposition<0.0||yposition>0.0||yposition<-5.0)) {eventindex++;continue;} //exclude these positions
		//if(!(yposition<10.0&&yposition>5.0&&xposition<-10.0&&xposition>-15.0)) {eventindex++;continue;}
		/*************************************/

		for(int i=0; i<16; i++) {
			ECalRawADCSpectrum[i]->Fill(ECalADC[i]);
		}
		for(int i=0; i<16; i++) {
			HodRawADCSpectrum[i]->Fill(HodADC[i]);
		}
		
		Sc1RawADCSpectrum->Fill(Sc1ADC);
		MonRawADCSpectrum->Fill(MonADC);
		Ce1RawADCSpectrum->Fill(Ce1ADC);
		Ce2RawADCSpectrum->Fill(Ce2ADC);
		PbGRawADCSpectrum->Fill(PbGADC);

		eventindex++;
		goodevents++;
	}
		
	/********************Plots*********************/
	
	/*******************************/
	/*********Raw Spectra***********/
	/*******************************/
	
	TCanvas *cECal = new TCanvas("cECal", "ECal ADC", 1100, 700);
	TCanvas *cHod = new TCanvas("cHod", "Hod ADC", 1100, 700);
	TCanvas *cSc1 = new TCanvas("cSc1", "Sc1 ADC", 700, 500);
	TCanvas *cMon = new TCanvas("cMon", "Mon ADC", 700, 500);
	TCanvas *cCe1 = new TCanvas("cCe1", "Ce1 ADC", 700, 500);
	TCanvas *cCe2 = new TCanvas("cCe2", "Ce2 ADC", 700, 500);
	TCanvas *cPbG = new TCanvas("cPbG", "PbG ADC", 700, 500);
	
	cECal->Divide(4,4);
	for(int i=0; i<16; i++) {
		cECal->cd(i+1);
		ECalRawADCSpectrum[i]->Draw();
		ECalRawADCSpectrum[i]->GetXaxis()->SetRangeUser(0,1000);
	}
	cECal->Update();
	
	cHod->Divide(4,4);
	for(int i=0; i<16; i++) {
		cHod->cd(i+1);
		HodRawADCSpectrum[i]->Draw();
		HodRawADCSpectrum[i]->GetXaxis()->SetRangeUser(0,1000);
		gPad->SetLogy();
	}
	cHod->Update();
	
	cSc1->cd();
	Sc1RawADCSpectrum->Draw();
	cSc1->Update();
	
	cMon->cd();
	MonRawADCSpectrum->Draw();
	cMon->Update();
	
	cCe1->cd();
	Ce1RawADCSpectrum->Draw();
	cCe1->Update();
	
	cCe2->cd();
	Ce2RawADCSpectrum->Draw();
	cCe2->Update();
	
	cPbG->cd();
	PbGRawADCSpectrum->Draw();
	cPbG->Update();
	
	/****************************/
	/********Hodoscope***********/
	/****************************/
	
	/*
	TCanvas *cxm = new TCanvas("cxm", "X Multiplicity", 700, 500);
	TCanvas *cym = new TCanvas("cym", "Y Multiplicity", 700, 500);
	TCanvas *cx = new TCanvas("cx", "X Beam Profile", 700, 500);
	TCanvas *cy = new TCanvas("cy", "Y Beam Profile", 700, 500);
	TCanvas *cxy = new TCanvas("cxy", "2D Beam Profile", 1000, 700);
	
	cxm->cd();
	XMult->Draw();
	gPad->SetLogy();
	cxm->Update();
	
	cym->cd();
	YMult->Draw();
	gPad->SetLogy();
	cym->Update();
	
	cx->cd();
	BeamProfileX->Draw();
	BeamProfileX->GetXaxis()->SetTitle("Position in mm East->West");
	cx->Update();
	
	cy->cd();
	BeamProfileY->Draw();
	BeamProfileY->GetXaxis()->SetTitle("Position in mm Down->Up");
	cy->Update();
	
	cxy->cd();
	BeamXY->Draw("zcol text");
	BeamXY->GetXaxis()->SetTitle("Position in mm East->West");
	BeamXY->GetYaxis()->SetTitle("Position in mm Down->Up");
	cy->Update();
	*/
	 
	/*************************/
	/***Ckov and Lead Glass***/
	/*************************/
	
	/*
	TCanvas *cCeCor = new TCanvas("cCeCor", "Cerenkov Correlation", 1000, 700);
	TCanvas *cPbGPedSub = new TCanvas("cPbGPedSub", "PbG Ped Sub", 1000, 700);
	TCanvas *cPbGEXY = new TCanvas("cPbGEXY", "PbG E-XY Cor", 1000, 700);
	
	cCeCor->cd();
	CerenkovCor->Draw("zcol");
	CerenkovCor->GetXaxis()->SetTitle("Ce1 Pedestal Subtracted ADC");
	CerenkovCor->GetYaxis()->SetTitle("Ce2 Pedestal Subtracted ADC");
	cCeCor->Update();
	
	cPbGPedSub->cd();
	PbGPedSubtracted->Draw();
	cPbGPedSub->Update();
	
	cPbGEXY->cd();
	PbGEXYCor->Draw("zcol");
	cPbGEXY->Update();
	*/
}


