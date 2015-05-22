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

const bool ecal_plots = false;
const bool pbg_plots = true;

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

Analysis() {
	
	gStyle->SetOptStat("nemr");
	
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
	
	//initialize pedestal values;
	float PbGThreshold = 120;
	SetPedestals();
	SetGainFactors();

	char buffer1[100];
	char buffer2[100];
	
	/****Raw Spectra Histograms****/
	TH1D* ECalRawADCSpectrum[16];
	for(int i=0; i<ECalNChannel; i++) {
		sprintf(buffer1, "ECalADCHist%d", i);
		sprintf(buffer2, "ECal ADC Channel %d", i);
		ECalRawADCSpectrum[i] = new TH1D(buffer1, buffer2, 4096, .5, 4096.5);
	}
	TH1D* HodRawADCSpectrum[16];
	for(int i=0; i<16; i++) {
		sprintf(buffer1, "HodADCHist%d", i);
		sprintf(buffer2, "Hod ADC Channel %d", i);
		HodRawADCSpectrum[i] = new TH1D(buffer1, buffer2, 4096, .5, 4096.5);
	}

	TH1D* Sc1RawADCSpectrum = new TH1D("Sc1ADCHist", "Sc1 ADC", 4096, .5, 4096.5);
	TH1D* MonRawADCSpectrum = new TH1D("MonADCHist", "Mon ADC", 4096, .5, 4096.5);
	TH1D* Ce1RawADCSpectrum = new TH1D("Ce1ADCHist", "Ce1 ADC", 4096, .5, 4096.5);
	TH1D* Ce2RawADCSpectrum = new TH1D("Ce2ADCHist", "Ce2 ADC", 4096, .5, 4096.5);
	TH1D* PbGRawADCSpectrum = new TH1D("PbGADCHist", "PbG ADC", 4096, .5, 4096.5);
	
	TH1D* XMult = new TH1D("XMult", "Multiplicity in X direction", 9, -.5, 8.5);
	TH1D* YMult = new TH1D("YMult", "Multiplicity in Y direction", 9, -.5, 8.5);
	TH1D* BeamProfileX = new TH1D("BeamProfileX", "Beam Profile X", 8, -19.2, 19.2);
	TH1D* BeamProfileY = new TH1D("BeamProfileY", "Beam Profile Y", 8, -19.2, 19.2);
	TH2D* BeamXY = new TH2D("BeamXY", "Beam Profile XY", 8, -19.2, 19.2, 8, -19.2, 19.2);
		
	/*ECal histograms*/
	TH1D *ECalSumHist = new TH1D("ECalSumHist", "ECal Tower Sum", 6000, .5, 6000.5);
	TH1D *ECalSumElectron = new TH1D("ECalSumElectron", "ECal Electron Sum", 6000, .5, 6000.5);
	TH1D *ECalMultiplicityHist = new TH1D("ECalMultiplicityHist", "ECal Multiplicity", 17, -0.5, 16.5);
	TH1D *ECalMultiplicityElectron = new TH1D("ECalMultiplicityElectron", "ECal Multiplicity Electron", 17, -0.5, 16.5);
	TH2D *HodX_vs_ECalSumElectron = new TH2D("HodX_vs_ECalSumElectron", "Hodoscope X vs ECal elec signal", 8, -19.2, 19.2, 6000, .5, 6000.5);
	TH2D *HodY_vs_ECalSumElectron = new TH2D("HodY_vs_ECalSumElectron", "Hodoscope Y vs ECal elec signal", 8, -19.2, 19.2, 6000, .5, 6000.5);
	TH1D *ECalLocX = new TH1D("ECalLocX", "ECal Loc X", 200, -30.0, 30.0);
	TH1D *ECalLocY = new TH1D("ECalLocY", "ECal Loc Y", 200, -30.0, 30.0);
	TH2D *ECalLocX_vs_Sum = new TH2D("ECalLocX_vs_Sum", "ECal Loc X vs Sum", 100, -15.0, 15.0, 6000, .5, 6000.5);
	TH2D *ECalLocY_vs_Sum = new TH2D("ECalLocY_vs_Sum", "ECal Loc Y vs Sum", 100, -15.0, 15.0, 6000, .5, 6000.5);

	/*PbG Resolution*/
	TH1D *ElectroninPbG = new TH1D("ElectroninPbG", "Electron Signal in PbG", 4096, .5, 4096.5);

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
		
		XMult->Fill(xmult);
		YMult->Fill(ymult);
		
		 
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
		//if(!(Sc1ADC>1200)) {continue;} //Trigger counter cut
		//if(!(yposition>-5&&yposition<5&&xposition>-5&&xposition<5)) {continue;} //include these positions
		//if(!(xposition>5.0||xposition<0.0||yposition>0.0||yposition<-5.0)) {continue;} //exclude these positions
		bool ecut = (xmult==1&&ymult==1)&&(Ce1ADC>100)&&(Ce2ADC<100);
		/*************************************/
		
		//PbG
		if (ecut) {
			ElectroninPbG->Fill(PbGADC);
		}

		//ECal Analysis
		float ECalSum = 0;
		float ECalSumCorrected = 0.;
		int maxCh = -1;
		ECalMultiplicity = 0;
		for(int i=0; i<ECalNChannel; i++) {
			ECalSum += ECalHitsAbovePed[i];
			if(maxCh==-1 || ECalHitsAbovePed[i] > ECalHitsAbovePed[maxCh]) maxCh = i;
			if(ECalHitsAbovePed[i] > 0.0) ECalMultiplicity++;
		}

		ECalSumHist->Fill(ECalSum);
		ECalMultiplicityHist->Fill(ECalMultiplicity);

		//ECal Tower coordinates
		bool CalOK = false;
		float ECalX = 0.0; 
		float ECalY = 0.0;
		float ECalXloc = 0.0;
		float ECalYloc = 0.0;
		
		if(ECalSum > 0.0) {
			float wtSum = 0.0;
			
			for(int i=0; i<ECalNChannel; i++) {
				float weight = ECalHitsAbovePed[i] > 0. ? 3.8 + log(ECalHitsAbovePed[i]/ECalSum) : 0;
				if(weight>0.) {
					wtSum += weight;
					ECalX += weight*(i%4);
					ECalY += weight*(i/4);
				}
			}
			
			if(wtSum) {
				ECalX /= wtSum;
				ECalY /= wtSum;
				
				ECalXloc = 25.*(ECalX - maxCh%4);
				ECalYloc = 25.*(ECalY - maxCh/4);
				
				CalOK = true;
			}
		}
		
		if(maxCh==-1) {
			iEx = iEy = iEl = 0;
		}
		else {
			iEl = maxCh;
			iEx = maxCh%4;
			iEy = maxCh/4;
		}

		if (ecut) {
			ECalSumElectron->Fill(ECalSum);
			ECalMultiplicityElectron->Fill(ECalMultiplicity);	
			HodX_vs_ECalSumElectron->Fill(xposition, ECalSum);
			HodY_vs_ECalSumElectron->Fill(yposition, ECalSum);
			ECalLocX->Fill(ECalXloc);
			ECalLocY->Fill(ECalYloc);
			ECalLocX_vs_Sum->Fill(ECalXloc, ECalSum);
			ECalLocY_vs_Sum->Fill(ECalYloc, ECalSum);
		}

		//Hodoscope
		for(int i=0; i<8; i++) {
			if(HodADC[i]>=HodoscopeCutoff[i]) {
				for(int j=0; j<8; j++) {
					if(HodADC[8+j]>=HodoscopeCutoff[8+j]) {
						BeamXY->Fill(16.8 - 4.8*j, -16.8 + 4.8*i);
					}
				}
			}
		}

		for(int i=0; i<ECalNChannel; i++) {
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

		goodevents++;
	}
	
	BeamProfileX->Add(BeamXY->ProjectionX());
	BeamProfileY->Add(BeamXY->ProjectionY());
	
	/********************Plots*********************/
	
	gStyle->SetOptFit(0001);
	
	/*******************************/
	/*********Raw Spectra***********/
	/*******************************/
	
	TCanvas *cECal = new TCanvas("cECal", "ECal ADC", 1050, 850);
	TCanvas *cHod = new TCanvas("cHod", "Hod ADC", 1050, 850);
	TCanvas *cSc1 = new TCanvas("cSc1", "Sc1 ADC", 700, 500);
	TCanvas *cMon = new TCanvas("cMon", "Mon ADC", 700, 500);
	TCanvas *cCe1 = new TCanvas("cCe1", "Ce1 ADC", 700, 500);
	TCanvas *cCe2 = new TCanvas("cCe2", "Ce2 ADC", 700, 500);
	TCanvas *cPbG = new TCanvas("cPbG", "PbG ADC", 700, 500);
	
	cECal->Divide(4,4);
	for(int i=0; i<ECalNChannel; i++) {
		cECal->cd(i+1);
		ECalRawADCSpectrum[i]->Draw();
		ECalRawADCSpectrum[i]->GetXaxis()->SetRangeUser(0,3000);
		ECalRawADCSpectrum[i]->GetXaxis()->SetLabelSize(.05);
		ECalRawADCSpectrum[i]->GetYaxis()->SetLabelSize(.05);
		LabelAxes(ECalRawADCSpectrum[i], "ADC", "counts");
		gPad->SetLogy();
	}
	cECal->Update();
	
	cHod->Divide(4,4);
	for(int i=0; i<HodNChannel; i++) {
		cHod->cd(i+1);
		HodRawADCSpectrum[i]->Draw();
		HodRawADCSpectrum[i]->GetXaxis()->SetRangeUser(0,2000);
		HodRawADCSpectrum[i]->GetXaxis()->SetLabelSize(.05);
		HodRawADCSpectrum[i]->GetYaxis()->SetLabelSize(.05);
		LabelAxes(HodRawADCSpectrum[i], "ADC", "counts");
		gPad->SetLogy();
	}
	cHod->Update();
	
	cSc1->cd();
	Sc1RawADCSpectrum->Draw();
	LabelAxes(Sc1RawADCSpectrum, "ADC", "counts");
	cSc1->Update();
	
	cMon->cd();
	MonRawADCSpectrum->Draw();
	LabelAxes(MonRawADCSpectrum, "ADC", "counts");
	cMon->Update();
	
	cCe1->cd();
	Ce1RawADCSpectrum->Draw();
	LabelAxes(Ce1RawADCSpectrum, "ADC", "counts");
	cCe1->Update();
	
	cCe2->cd();
	Ce2RawADCSpectrum->Draw();
	LabelAxes(Ce2RawADCSpectrum, "ADC", "counts");
	cCe2->Update();
	
	cPbG->cd();
	PbGRawADCSpectrum->Draw();
	LabelAxes(PbGRawADCSpectrum, "ADC", "counts");
	cPbG->Update();
	
	/****************************/
	/********Hodoscope***********/
	/****************************/
	
	
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
	LabelAxes(BeamProfileX, "Position in mm West #rightarrow East", "counts");
	cx->Update();
	
	cy->cd();
	BeamProfileY->Draw();
	LabelAxes(BeamProfileY, "Position in mm Down #rightarrow Up", "counts");
	cy->Update();
	
	cxy->cd();
	BeamXY->Draw("zcol text");
	LabelAxes(BeamXY, "Position in mm West #rightarrow East", "Position in mm Down #rightarrow Up");
	cy->Update();
	
	/**************************/
	/********Lead Glass********/
	/**************************/
	if (pbg_plots) {
	TCanvas *cePbG = new TCanvas("cePbG", "Electron in PbG", 700, 500);
	cePbG->cd();
	LabelAxes(ElectroninPbG, "ADC", "counts");
	TF1 *fefit = new TF1("fefit", "gaus", 200, 600);
	ElectroninPbG->Fit(fefit, "R");
	ElectroninPbG->Draw();
	cePbG->Update();
	}
	/**************************/
	/***********ECAL***********/
	/**************************/
	
	if (ecal_plots) {
	TCanvas *cECalSum = new TCanvas("cECalSum", "Sum in ECal", 700, 500);
	cECalSum->cd();
	LabelAxes(ECalSumHist, "ADC Sum", "counts");
	ECalSumHist->Draw();
	cECalSum->Update();
	
	TCanvas *cECalElectronSum = new TCanvas("cECalElectronSum", "Sum for Electrons in ECal", 700, 500);
	cECalElectronSum->cd();
	LabelAxes(ECalSumElectron, "ADC Sum", "counts");
	ECalSumElectron->Draw();
	cECalElectronSum->Update();

	TCanvas *cECalMult = new TCanvas("cECalMult", "Mult in ECal", 700, 500);
	cECalMult->cd();
	LabelAxes(ECalMultiplicityHist, "Multiplicity", "counts");
	ECalMultiplicityHist->Draw();
	cECalMult->Update();
	
	TCanvas *cECalElectronMult = new TCanvas("cECalElectronMult", "Mult for Electrons in ECal", 700, 500);
	cECalElectronMult->cd();
	LabelAxes(ECalMultiplicityElectron, "Multiplicity", "counts");
	ECalMultiplicityElectron->Draw();
	cECalElectronSum->Update();

	TCanvas *cHodX_vs_ECalElectron = new TCanvas("cHodX_vs_ECalElectron", "Hodoscope X vs Electron Signal", 1000, 800);
	cHodX_vs_ECalElectron->cd();
	LabelAxes(HodX_vs_ECalSumElectron, "X (mm)", "ECal ADC Sum");
	HodX_vs_ECalSumElectron->Draw("zcol");
	cHodX_vs_ECalElectron->Update();
	
	TCanvas *cHodY_vs_ECalElectron = new TCanvas("cHodY_vs_ECalElectron", "Hodoscope Y vs Electron Signal", 1000, 800);
	cHodY_vs_ECalElectron->cd();
	LabelAxes(HodY_vs_ECalSumElectron, "Y (mm)", "ECal ADC Sum");
	HodY_vs_ECalSumElectron->Draw("zcol");
	cHodY_vs_ECalElectron->Update();
	}
	/**************************/
	/*********Profiles*********/
	/**************************/

	if (ecal_plots) {
	TCanvas *cECalLocX = new TCanvas("cECalLocX", "ECal Loc X", 700, 500);
	cECalLocX->cd();
	LabelAxes(ECalLocX, "Loc X", "counts");
	ECalLocX->Draw();
	cECalLocX->Update();
	
	TCanvas *cECalLocY = new TCanvas("cECalLocY", "ECal Loc Y", 700, 500);
	cECalLocY->cd();
	LabelAxes(ECalLocY, "Loc Y", "counts");
	ECalLocY->Draw();
	cECalLocY->Update();
	
	TCanvas *cECalLocX_vs_Sum = new TCanvas("cECalLocX_vs_Sum", "ECal Loc X vs ECal Sum", 900, 700);
	cECalLocX_vs_Sum->cd();
	LabelAxes(ECalLocX_vs_Sum, "Loc X", "ECalSum");
	ECalLocX_vs_Sum->Draw("zcol");
	cECalLocX_vs_Sum->Update();
	
	TCanvas *cECalLocY_vs_Sum = new TCanvas("cECalLocY_vs_Sum", "ECal Loc Y vs ECal Sum", 900, 700);
	cECalLocY_vs_Sum->cd();
	LabelAxes(ECalLocY_vs_Sum, "Loc Y", "ECalSum");
	ECalLocY_vs_Sum->Draw("zcol");
	cECalLocY_vs_Sum->Update();
	
	//Profile fits for corrections
	TH1D *LocX_Sum_Profile = ECalLocX_vs_Sum->ProfileX("LocX_Sum_Profile");
	TH1D *LocY_Sum_Profile = ECalLocY_vs_Sum->ProfileX("LocY_Sum_Profile");

	TF1 *fquadx = new TF1("fquadx", "[0]*(1 - [2]*x*x + [1]*x)", -10, 10);
	TF1 *fquady = new TF1("fquady", "[0]*(1 - [2]*x*x + [1]*x)", -10, 10);
	
	TCanvas *cECalLoc_Profiles = new TCanvas("cECalLoc_Profiles", "ECal Loc Variable Profiles", 1300, 700);
	cECalLoc_Profiles->Divide(2,1);
	cECalLoc_Profiles->cd(1);
	LocX_Sum_Profile->Draw();
	LocX_Sum_Profile->Fit(fquadx, "R");
	cECalLoc_Profiles->cd(2);
	LocY_Sum_Profile->Draw();
	LocY_Sum_Profile->Fit(fquady, "R");
	cECalLoc_Profiles->Update();
	}
}


