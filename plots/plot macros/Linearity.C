/*Linearity*/

Linearity() {
	gStyle->SetOptFit();
	
	float Ebeam[] = {2.0,3.0,4.0,6.0,8.0};
	float Ebeam_err[] = {0.,0.,0.,0.,0.};
	float beam_spread[] = {0.,0.,0.,0.,0.};

	float Ebeam_no3[] = {2.0,4.0,6.0,8.0};
	float Ebeam_no3_err[] = {0.,0.,0.,0.,0.};
	float beam_no3_spread[] = {.027,.027,0.023,0.023};

	float ECal_mean[] = {270.5, 475.9, 586.8, 921.5, 1323.};
	float PbG_mean[] = {775.3, 1320., 1613., 2556., 3546.};

	float OldPbG_mean[] = {601.5, 1105, 1622, 2074};

	TF1 *fLinFitEcal = new TF1("fLinFitEcal", "pol1", 0, 10);
	TF1 *fLinFitPbG = new TF1("fLinFitPbG", "pol1", 0, 10);
	TF1 *fLinFitOldPbG = new TF1("fLinFitOldPbG", "pol1", 0, 10);

	float ECal_mean_err[] = {0.,0.,0.,0.,0.};
	float PbG_mean_err[] = {0.,0.,0.,0.,0.};
	float OldPbG_mean_err[] = {0.,0.,0.,0.};
	TGraphErrors *gLinECal = new TGraphErrors(5, Ebeam, ECal_mean, Ebeam_err, ECal_mean_err);
	TGraphErrors *gLinPbG = new TGraphErrors(5, Ebeam, PbG_mean, Ebeam_err, PbG_mean_err);
	TGraphErrors *gLinOldPbG = new TGraphErrors(4, Ebeam_no3, OldPbG_mean, Ebeam_no3_err, OldPbG_mean_err);

	gLinECal->Fit(fLinFitEcal, "NR");
	gLinPbG->Fit(fLinFitPbG, "NR");
	gLinOldPbG->Fit(fLinFitOldPbG, "NR");

	fLinFitEcal->SetLineStyle(2);
	fLinFitEcal->SetLineColor(kGreen+2);
	fLinFitPbG->SetLineStyle(2);
	fLinFitPbG->SetLineColor(kRed+1);
	fLinFitOldPbG->SetLineStyle(2);
	fLinFitOldPbG->SetLineColor(kMagenta+1);

	gLinECal->SetMarkerStyle(20);
	gLinECal->SetMarkerSize(2.0);
	gLinECal->SetMarkerColor(kGreen+2);
	gLinPbG->SetMarkerStyle(21);
	gLinPbG->SetMarkerSize(2.0);
	gLinPbG->SetMarkerColor(kRed+1);
	gLinOldPbG->SetMarkerStyle(22);
	gLinOldPbG->SetMarkerSize(2.0);
	gLinOldPbG->SetMarkerColor(kMagenta+1);

	TLegend *llin = new TLegend(.18, .66, .61, .83);
	llin->AddEntry(fLinFitEcal, Form("ECAL: %.2f*E + %.2f", fLinFitEcal->GetParameter(1), fLinFitEcal->GetParameter(0)), "L");
	llin->AddEntry(fLinFitPbG, Form("2015 Lead Glass: %.2f*E + %.2f", fLinFitPbG->GetParameter(1), fLinFitPbG->GetParameter(0)), "L");
	llin->AddEntry(fLinFitOldPbG, Form("2014 Lead Glass: %.2f*E + %.2f", fLinFitOldPbG->GetParameter(1), fLinFitOldPbG->GetParameter(0)), "L");
	llin->SetFillStyle(0);
	llin->SetTextSize(.03);
	llin->SetTextFont(42);

	TCanvas *cLin = new TCanvas("cLin", "Linearity", 900, 700);
	cLin->cd();
	gLinPbG->Draw("AP");
	gLinECal->Draw("Psame");
	gLinOldPbG->Draw("Psame");
	gLinPbG->GetXaxis()->SetTitle("Beam Energy (GeV)");
	gLinPbG->GetYaxis()->SetTitle("Detector Response");
	fLinFitEcal->Draw("same");
	fLinFitPbG->Draw("same");
	fLinFitOldPbG->Draw("same");
	llin->Draw("same");
	gLinPbG->GetXaxis()->SetRangeUser(0, 10);
	gLinPbG->GetYaxis()->SetRangeUser(0, 4000);
	cLin->Update();

	float ECalDev[5] = {0.,0.,0.,0.,0.}, ECalDev_err[5] = {0.,0.,0.,0.,0.};
	float PbGDev[5] = {0.,0.,0.,0.,0.}, PbGDev_err[5] = {0.,0.,0.,0.,0.};
	float OldPbGDev[4] = {0.,0.,0.,0.}, OldPbGDev_err[4] = {0.,0.,0.,0.};
	for (int i=0; i<5; i++) {
		ECalDev[i] = (ECal_mean[i] - fLinFitEcal->Eval(Ebeam[i]))/ECal_mean[i];
		PbGDev[i] = (PbG_mean[i] - fLinFitPbG->Eval(Ebeam[i]))/PbG_mean[i];
		if (i<4) OldPbGDev[i] = (OldPbG_mean[i] - fLinFitOldPbG->Eval(Ebeam_no3[i]))/OldPbG_mean[i];
	}

	TGraphErrors *gECalDev = new TGraphErrors(5, Ebeam, ECalDev, Ebeam_err, ECalDev_err);
	TGraphErrors *gPbGDev = new TGraphErrors(5, Ebeam, PbGDev, Ebeam_err, PbGDev_err);
	TGraphErrors *gOldPbGDev = new TGraphErrors(4, Ebeam_no3, OldPbGDev, Ebeam_no3_err, OldPbGDev_err);
	TH1D *zeroline = new TH1D("zeroline", "zeroline", 1, 0, 10.0);
	zeroline->SetBinContent(1, 0.0);
	zeroline->SetLineStyle(2);

	gECalDev->SetMarkerStyle(20);
	gECalDev->SetMarkerSize(2.0);
	gECalDev->SetMarkerColor(kGreen+2);
	gPbGDev->SetMarkerStyle(21);
	gPbGDev->SetMarkerSize(2.0);
	gPbGDev->SetMarkerColor(kRed+1);
	gOldPbGDev->SetMarkerStyle(22);
	gOldPbGDev->SetMarkerSize(2.0);
	gOldPbGDev->SetMarkerColor(kMagenta+1);

	TCanvas *cDev = new TCanvas("cDev", "Deviation from Linearity", 900, 700);
	cDev->cd();
	zeroline->Draw("][");
	gECalDev->Draw("sameP");
	gPbGDev->Draw("sameP");
	gOldPbGDev->Draw("sameP");
	zeroline->GetXaxis()->SetTitle("Beam Energy (GeV)");
	zeroline->GetYaxis()->SetTitle("Deviation from Linearity");
	zeroline->GetXaxis()->SetRangeUser(0.0, 10.0);
	zeroline->GetYaxis()->SetRangeUser(-.10, .10);
	cDev->Update();

}