/*Ecal and PbGlass resolution*/

Resolutions() {

	float Ebeam[] = {2.0,3.0,4.0,6.0,8.0,12.0};
	float Ebeam_err[] = {0.,0.,0.,0.,0.,0.};
	float beam_spread[] = {.027,.027,0.023,0.023,0.023,0.023};
	//float beam_spread[] = {0.,0.,0.,0.,0.,0.};


	float ECal_mean[] = {270.5, 475.9, 586.8, 921.5, 1323., 3022.};
	float ECal_sigma[] = {24.69, 31.54, 35.64, 47.51, 55.71, 133.};

	float PbG_mean[] = {775.3, 1320., 1613., 2556., 3546., 1846.};
	float PbG_sigma[] = {49.63, 71.96, 81., 125.3, 137., 77.99};
	
	float ECal_res[6], ECal_res_err[6];
	float PbG_res[6], PbG_res_err[6];
	for (int i=0; i<6; i++) {
		ECal_res[i] = sqrt(ECal_sigma[i]*ECal_sigma[i]/(ECal_mean[i]*ECal_mean[i]) - beam_spread[i]*beam_spread[i]);
		PbG_res[i] = sqrt(PbG_sigma[i]*PbG_sigma[i]/(PbG_mean[i]*PbG_mean[i]) - beam_spread[i]*beam_spread[i]);
		ECal_res_err[i] = 0.;
		PbG_res_err[i] = 0.;
	}

	TGraphErrors *gECalRes = new TGraphErrors(6, Ebeam, ECal_res, Ebeam_err, ECal_res_err);
	TGraphErrors *gPbGRes = new TGraphErrors(6, Ebeam, PbG_res, Ebeam_err, PbG_res_err);
	TF1 *fResECal = new TF1("fResECal", "[1]/sqrt(x) + [0]", 2, 14);
	TF1 *fResPbG = new TF1("fResPbG", "[1]/sqrt(x) + [0]", 2, 14);
	fResECal->SetLineStyle(2);
	fResECal->SetLineColor(kGreen+2);
	fResPbG->SetLineStyle(2);
	fResPbG->SetLineColor(kRed+1);

	gECalRes->Fit(fResECal, "NR");
	gPbGRes->Fit(fResPbG, "NR");

	gECalRes->SetMarkerStyle(20);
	gECalRes->SetMarkerSize(2.0);
	gECalRes->SetMarkerColor(kGreen+2);
	gECalRes->SetLineColor(kGreen+2);
	gPbGRes->SetMarkerStyle(21);
	gPbGRes->SetMarkerSize(2.0);
	gPbGRes->SetMarkerColor(kRed+1);
	gPbGRes->SetLineColor(kRed+1);

	TLegend *lres = new TLegend(.45, .66, .88, .83);
	lres->AddEntry(fResECal, Form("ECAL: #frac{%.2f %%}{#sqrt{E}} + %.2f", fResECal->GetParameter(1)*100.0, fResECal->GetParameter(0)*100.0), "L");
	lres->AddEntry(fResPbG, Form("Lead Glass: #frac{%.2f %%}{#sqrt{E}} + %.2f", fResPbG->GetParameter(1)*100.0, fResPbG->GetParameter(0)*100.0), "L");
	lres->SetFillStyle(0);
	lres->SetTextSize(.035);
	lres->SetTextFont(42);

	TCanvas *cRes = new TCanvas("cRes", "Resolutions", 900, 700);
	cRes->cd();
	gECalRes->Draw("AP");
	gECalRes->GetXaxis()->SetTitle("Beam Energy (GeV)");
	gECalRes->GetYaxis()->SetTitle("Resolution");
	gPbGRes->Draw("Psame");
	fResECal->Draw("same");
	fResPbG->Draw("same");
	lres->Draw("same");
	cRes->Update();

}