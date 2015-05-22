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

//alternative setups

//event struct

const int ECalNChannel = 16;
const int HodNChannel = 16;



typedef struct {
	
	int eventno;
	float ECalRawADC[ECalNChannel];
	float HodRawADC[HodNChannel];
	float Sc1RawADC;
	float MonRawADC;
	float Ce1RawADC;
	float Ce2RawADC;
	float PbGRawADC;
		
} Event;

ReadoutHistograms() {
	
	float SingleEventRawADC[64];
	Event evn;
	
	int eventno;

	float ECalRawADC[ECalNChannel];
	float HodRawADC[HodNChannel];
	for(int i=0; i<ECalNChannel; i++) ECalRawADC[i] = 0;
	for(int i=0; i<HodNChannel; i++) HodRawADC[i] = 0;
	float Sc1RawADC = 0;
	float MonRawADC = 0;
	float Ce1RawADC = 0;
	float Ce2RawADC = 0;
	float PbGRawADC = 0;
	
	//Channel lists modify here for configuration changes
	
	int ECalChannelList[ECalNChannel];
	for (int i=0; i<ECalNChannel; i++) ECalChannelList[i] = i;
	
	int HodChannelList[HodNChannel];
	for(int i=0; i<HodNChannel; i++) HodChannelList[i] = i+16;

	int Sc1ChannelList = 32;
	int MonChannelList = 33;
	int Ce1ChannelList = 34; //Inner
	int Ce2ChannelList = 35; //Outer
	int PbGChannelList = 36;

	char buffer1[100];
	char buffer2[100];
	TH1D* ChannelHist[64];
	for (int i=0; i<64; i++) {
		sprintf(buffer1, "ChannelHist%d", i);
		sprintf(buffer2, "Channel %d", i);
		ChannelHist[i] = new TH1D(buffer1, buffer2, 4096, .5, 4096.5);
	}
	int index = 0;
	long int number;
	int thisevent;
	char fname[200];
	char outfilename[200];
	char linebuffer[1000];
	int ch1, ch2, ch3;
	cout << "Enter file name: " <<endl;
	cin >> fname;
	
	strncpy(outfilename, fname, 200);
	strcat(outfilename, ".root");
	
	TFile f(outfilename, "recreate");
	TTree T("T", "event data");
	
	//branches
	T.Branch("eventno", &evn.eventno, "eventno/I");
	T.Branch("ECalRawADC", evn.ECalRawADC, "ECalRawADC[16]/F");
	T.Branch("HodRawADC", evn.HodRawADC, "HodRawADC[16]/F");
	T.Branch("Sc1RawADC", &evn.Sc1RawADC, "Sc1RawADC/F");
	T.Branch("MonRawADC", &evn.MonRawADC, "MonRawADC/F");
	T.Branch("Ce1RawADC", &evn.Ce1RawADC, "Ce1RawADC/F");
	T.Branch("Ce2RawADC", &evn.Ce2RawADC, "Ce2RawADC/F");
	T.Branch("PbGRawADC", &evn.PbGRawADC, "PbGRawADC/F");

	bool endofevent = false;
	ifstream indata;
	indata.open(fname);
	indata.clear();
	
	while(!indata.eof()) {

		indata >> number;
		
		if(index%65==0) {
			thisevent = number;
			if((index/65)%10000==0) cout << "Processing event: " << index/65 << endl;
		}
		else {
			SingleEventRawADC[index%65-1] = number;
			if(index%65==64) endofevent = true;
		}
		
		//end of single event in loop add to appropriate arrays
		if(endofevent) {
			eventno = thisevent;
			
			for(int i=0; i<ECalNChannel; i++) {
				ECalRawADC[i] = SingleEventRawADC[ECalChannelList[i]];
			}
			
			Sc1RawADC = SingleEventRawADC[Sc1ChannelList];
			MonRawADC = SingleEventRawADC[MonChannelList];
			Ce1RawADC = SingleEventRawADC[Ce1ChannelList];
			Ce2RawADC = SingleEventRawADC[Ce2ChannelList];
			PbGRawADC = SingleEventRawADC[PbGChannelList];
			
			for(int i=0; i<HodNChannel; i++) {
				HodRawADC[i] = SingleEventRawADC[HodChannelList[i]];
			}
			
			//Write to event struct
			evn.eventno = eventno;

			for (int i=0; i<ECalNChannel; i++) {
				evn.ECalRawADC[i] = ECalRawADC[i];
			}
			for (int i=0; i<HodNChannel; i++) {
				evn.HodRawADC[i] = HodRawADC[i];
			}
			evn.Sc1RawADC = Sc1RawADC;
			evn.MonRawADC = MonRawADC;
			evn.Ce1RawADC = Ce1RawADC;
			evn.Ce2RawADC = Ce2RawADC;
			evn.PbGRawADC = PbGRawADC;
			
			
			T.Fill();
			
			for(int i=0; i<64; i++) SingleEventRawADC[i] = 0;
			endofevent = false;
		}
		
		index++;

	}
	indata.close();	
	
	f.Write();

}
