#ifndef __CINT__
#include <fstream.h>
#include "TFile.h"
#include "TCanvas.h"
#include "TPostScript.h"
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TObject.h"
#include "TKey.h"
#include "TTree.h"
#include "TObject.h"
#include "TIterator.h"
#include <iostream.h>
#include <fstream.h>
#include <stdio.h>
#include <cmath>
#include <math.h>
#include <string.h>
#endif

const int ECalNChannel = 16;
const int HodNChannel  = 16;

const double ECal_Pedestal_Width = 1.7;

//Manually set thresholds

const float HodoscopeCutoff[HodNChannel] = {250, 250, 250, 250,
											250, 250, 250, 250,
											250, 250, 250, 250, 
											250, 250, 250, 250};


float ECalPedestal[ECalNChannel];
float ECalGainFactor[ECalNChannel];

void LabelAxes(TH1* hist, char* xtitle, char* ytitle) {

	hist->GetXaxis()->SetTitleFont(52);
	hist->GetXaxis()->SetTitle(xtitle);
	hist->GetYaxis()->SetTitleFont(52);
	hist->GetYaxis()->SetTitle(ytitle);
	return;

}

void SetPedestals() {
	for (int i=0; i<ECalNChannel; i++) {
		ECalPedestal[i] = 100;
	}

	//set towers manually below
	ECalPedestal[0] = 100.7,ECalPedestal[1] = 100.3,ECalPedestal[2] = 99.6,ECalPedestal[3] = 101.2,
	ECalPedestal[4] = 102.1,ECalPedestal[5] = 101.7,ECalPedestal[6] = 101.5,ECalPedestal[7] = 102.5,
	ECalPedestal[8] = 102.4,ECalPedestal[9] = 103.2,ECalPedestal[10] = 103.0,ECalPedestal[11] = 103.2,
	ECalPedestal[12] = 103.2,ECalPedestal[13] = 103.5,ECalPedestal[14] = 103.7,ECalPedestal[15] = 105.3;
}

void SetGainFactors() {
	for (int i=0; i<ECalNChannel; i++) {
		ECalGainFactor[i] = 1.0;
	}

	//set towers manually below
	ECalGainFactor[0] = .977, ECalGainFactor[1] = 1.061, ECalGainFactor[2] = .985, ECalGainFactor[3] = 1.022,
	ECalGainFactor[4] = .987, ECalGainFactor[5] = .996, ECalGainFactor[6] = .996, ECalGainFactor[7] = 1.015,
	ECalGainFactor[8] = .986, ECalGainFactor[9] = .990, ECalGainFactor[10] = .991, ECalGainFactor[11] = .985,
	ECalGainFactor[12] = 1.024, ECalGainFactor[13] = 1.020, ECalGainFactor[14] = .981, ECalGainFactor[15] = .973;
}

void ComputeOneFile(char* file, TH1D *AttLength, int offset) {
	
	gStyle->SetOptStat("nemr");
	
	int eventno = 0;
	
	TChain *chain = new TChain("T");
	chain->Add(file);
	cout << "chain entries: " << chain->GetEntries() << endl;
	int nentries = chain->GetEntries();
	
	//initialize pedestal values;
	float PbGThreshold = 120;
	SetPedestals();
	SetGainFactors();

	/*ECal histograms*/
	TH2D *HodX_vs_ECalSumElectron = new TH2D("HodX_vs_ECalSumElectron", "Hodoscope X vs ECal elec signal", 8, -19.2, 19.2, 6000, .5, 6000.5);
	TH1D *HodXProjections[8];

	int eventindex = 0;
	int goodevents = 0;
	
	float ECalADC[ECalNChannel];
	float HodADC[HodNChannel];
	float Sc1ADC;
	float MonADC;
	float Ce1ADC;
	float Ce2ADC;
	float PbGADC;
	
	bool HighHodHit = false;
	float xposition = 0;
	float yposition = 0;
	
	float PbGPedestalSubtracted = 0;
	float Ce1PedestalSubtracted = 0;
	float Ce2PedestalSubtracted = 0;
	
	float PbGHitsAbovePed = 0;
	float Ce1HitsAbovePed = 0;
	float Ce2HitsAbovePed = 0;
	
	int iEl, iEx, iEy;

	int ECalMultiplicity = 0;
	float ECalHitsAbovePed[ECalNChannel];
	for(int i=0;i<ECalNChannel;i++) ECalHitsAbovePed[i]=0;
	
	chain->SetBranchAddress("ECalRawADC", ECalADC);
 	chain->SetBranchAddress("HodRawADC", HodADC);
	chain->SetBranchAddress("Sc1RawADC", &Sc1ADC);
	chain->SetBranchAddress("MonRawADC", &MonADC);
	chain->SetBranchAddress("Ce1RawADC", &Ce1ADC);
	chain->SetBranchAddress("Ce2RawADC", &Ce2ADC);
	chain->SetBranchAddress("PbGRawADC", &PbGADC);
	
	int consec_pulser_events = 0, nonpulser_events_counter = 0;
	int spill = 1;

	for(eventindex=0; eventindex<nentries; eventindex++) {	

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
		
		if (Sc1ADC<130 && xmult==0 && ymult==0) {
			consec_pulser_events++;
		}
		else {
			consec_pulser_events = 0;
			nonpulser_events_counter++;
		}

		if (consec_pulser_events >= 10) nonpulser_events_counter = 0;
		if (nonpulser_events_counter == 50) {
			spill++;
		}//threshold value for new spill
		 
		float ECalPedestalSubtracted;
		for(int i=0; i<ECalNChannel; i++) {
			ECalPedestalSubtracted = ECalADC[i] - ECalPedestal[i];
			if(ECalPedestalSubtracted>2*ECal_Pedestal_Width) {
				ECalHitsAbovePed[i] = ECalGainFactor[i]*ECalPedestalSubtracted;
			} 
			else {
				ECalHitsAbovePed[i] = 0.0;
			}
			
		}
		
		/****************Cuts******************/
		if(!(xmult==1&&ymult==1)) {continue;} //xmult ymult both must be 1
		//if(!(Sc1ADC>200)) {continue;}
		//if(HighHodHit) {continue;}
		//if(!(MonADC<300)) {continue;} //Monitor Cut
		//if(!(Sc1ADC>700)) {continue;} //Trigger counter cut
		//if(!(yposition>0&&xposition>-15&&xposition<-5)) {continue;} //include these positions
		//if(!(xposition>5.0||xposition<0.0||yposition>0.0||yposition<-5.0)) {continue;} //exclude these positions
		bool ecut = (xmult==1&&ymult==1)&&(Ce1ADC>100);
		/*************************************/
		
		//ECal Analysis
		float ECalSum = 0;
		int maxCh = -1;
		for(int i=0; i<ECalNChannel; i++) {
			ECalSum += ECalHitsAbovePed[i];
			if(maxCh==-1 || ECalHitsAbovePed[i] > ECalHitsAbovePed[maxCh]) maxCh = i;
		}

		
		if (ecut) {
			//HodX_vs_ECalSumElectron->Fill(xposition, ECalSum);
			HodX_vs_ECalSumElectron->Fill(xposition, ECalHitsAbovePed[6]);
		}

		goodevents++;
	}
	
	for (int i=0; i<8; i++) {
		HodXProjections[i] = HodX_vs_ECalSumElectron->ProjectionY(Form("_py_%d_%d", offset,i), i+1,i+1, "");
	}

	/********************Plots*********************/
	
	gStyle->SetOptFit(0001);
	
	/**************************/
	/***********ECAL***********/
	/**************************/

	TF1 *fit_hod_bin[8];
	TCanvas *cHodDependence = new TCanvas("cHodDependence", "Hodoscope X dependence", 1000, 700); 
	cHodDependence->Divide(4,2);
	for (int i=0; i<8; i++) {
		cHodDependence->cd(i+1);
		fit_hod_bin[i] = new TF1(Form("fit_hod_bin_%d", i), "gaus", 600, 1000);
		HodXProjections[i]->Rebin(60);
		HodXProjections[i]->Fit(fit_hod_bin[i], "R");
		HodXProjections[i]->Draw();
		AttLength->SetBinContent(offset+(i+1), fit_hod_bin[i]->GetParameter(1));
		AttLength->SetBinError(offset+(i+1), fit_hod_bin[i]->GetParError(1));
	}
	cHodDependence->Update();
	delete[] fit_hod_bin;
	delete[] HodXProjections;
	delete HodX_vs_ECalSumElectron;
	delete chain;
}

AttenuationLength() {

	TH1D *AttLength = new TH1D("AttLength", "Attenuation Length", 24, 0.0, .5*24.5);
	ComputeOneFile("/Users/jaydunkelberger/testrun2015/data/2015_5_25/data50.txt.root", AttLength, 0);
	ComputeOneFile("/Users/jaydunkelberger/testrun2015/data/2015_5_25/data49.txt.root", AttLength, 8);
	ComputeOneFile("/Users/jaydunkelberger/testrun2015/data/2015_5_25/data51.txt.root", AttLength, 16);

	TCanvas *cAttLength = new TCanvas("cAttLength", "Attenuation in Fibers", 900, 700);
	cAttLength->cd();
	AttLength->Fit("expo");
	cout << "att length: " << -1.0/AttLength->GetFunction("expo")->GetParameter(1) << endl;
	cout << "att length error: " << AttLength->GetFunction("expo")->GetParError(1) << endl;
	AttLength->Draw("][");
	cAttLength->Update();
}
