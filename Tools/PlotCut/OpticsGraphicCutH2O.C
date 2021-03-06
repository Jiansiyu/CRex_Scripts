/*
 * OpticsGraphicCutH2O.C
 *
 *  Created on: Apr 10, 2020
 *      Author: newdriver
 */

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TChain.h>
#include <TCut.h>
#include <TCutG.h>
#include <TPad.h>
#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1.h>
#include <TF1.h>
#include <TMath.h>
#include <TF1NormSum.h>
#include <TPaveText.h>
#include <map>
#include <vector>
#include <random>
#include <iostream>
#include <sys/stat.h>

#include <TComplex.h>
#include <TVirtualPad.h>

#include <TSpectrum2.h>
#include <TF2.h>
#include <TObject.h>
#include "TMinuit.h"
#include <TFile.h>
#include <fstream>
#include <TLatex.h>
#include <TSystem.h>
#include "TPaveText.h"
#include "TText.h"
#include <TApplication.h>

int FoilID=0;

int col=3;
int row=1;
int row_min=0;
int row_count=10;
const UInt_t NSieveCol = 13;
const UInt_t NSieveRow = 7;

//////////////////////////////////////////////////////////////////////////////
// Work Directory
// cut options
// Need to change
//////////////////////////////////////////////////////////////////////////////
TString prepcut;
TString generalcut;
TString generalcutR="R.tr.n==1 && R.vdc.u1.nclust==1&& R.vdc.v1.nclust==1 && R.vdc.u2.nclust==1 && R.vdc.v2.nclust==1 && fEvtHdr.fEvtType==1 && R.gold.p > 2.14 && R.gold.p < 2.2";
TString generalcutL="L.tr.n==1 && L.vdc.u1.nclust==1&& L.vdc.v1.nclust==1 && L.vdc.u2.nclust==1 && L.vdc.v2.nclust==1 && fEvtHdr.fEvtType==1 && L.gold.p > 2.14 && L.gold.p < 2.19";



inline TF1 *SpectrumFit(TH1F *spectrum){

	// start the fit
	TH1F *momentum=(TH1F *) spectrum->Clone("spectrum");

	// search for the highest peak, take it to be the ground states peak
	auto CGroundp=momentum->GetXaxis()->GetBinCenter(momentum->GetMaximumBin());

	// probably need to use a file to config the peak
	double_t fgroudGausPar[3];
	double_t ffirstGuasPar[3];
	TF1 *fgroudGaus=new TF1("groudstatesgaus","gaus",CGroundp-0.0005,CGroundp+0.0005);
	momentum->Fit("groudstatesgaus","R","ep",fgroudGaus->GetXmin(),fgroudGaus->GetXmax());
	fgroudGaus->GetParameters(fgroudGausPar);

	// do a crystal ball fit function for the ground state peak
	double_t fgroundCrystalballPar[5];
	TF1 *fgroundCrystalball=new TF1("fgroundCrystal","crystalball",fgroudGausPar[1]-0.0030,fgroudGaus->GetXmax()+0.0003);
	fgroundCrystalball->SetParameters(fgroudGausPar[0],fgroudGausPar[1],fgroudGausPar[2],1.64,1.1615);
	momentum->Fit("fgroundCrystal","R","same",fgroundCrystalball->GetXmin(),fgroundCrystalball->GetXmax());
	fgroundCrystalball->GetParameters(fgroundCrystalballPar);

	// fit the first excited states
	auto C1stp=fgroudGausPar[1]-0.0062;
	TF1 *ffirstGuas=new TF1 ("firststatesgaus","gaus",C1stp-0.001,C1stp+0.00155);
	momentum->Fit("firststatesgaus","R","ep",ffirstGuas->GetXmin(),ffirstGuas->GetXmax());
	//ffirstGuas->Draw("same");
	ffirstGuas->GetParameters(ffirstGuasPar);

	double fsecondGausPar[3];
	fsecondGausPar[1]=2.1615;
	TF1 *fsecondGaus=new TF1("secondCgaus","gaus",fsecondGausPar[1]-0.001,fsecondGausPar[1]+0.001);
	momentum->Fit("secondCgaus","R","ep",fsecondGaus->GetXmin(),fsecondGaus->GetXmax());
//	fsecondGaus->Draw("same");
	fsecondGaus->GetParameters(fsecondGausPar);

	double fthirdGausPar[3];
	fthirdGausPar[1]=2.1593;
	TF1 *fthirdGaus=new TF1("secondCgaus","gaus",fthirdGausPar[1]-0.0005,fthirdGausPar[1]+0.0005);
	momentum->Fit("secondCgaus","R","ep",fthirdGaus->GetXmin(),fthirdGaus->GetXmax());
//	fthirdGaus->Draw("same");
	fthirdGaus->GetParameters(fthirdGausPar);

	//
	double_t ffirstCrystalPar[5];
	TF1 *ffirstCrystal=new TF1("ffirstCrystal","crystalball",ffirstGuasPar[1]-0.0025,ffirstGuas->GetXmax());
	ffirstCrystal->SetParameters(ffirstGuasPar[0],ffirstGuasPar[1],ffirstGuasPar[2],1.64,1.1615);
	momentum->Fit("ffirstCrystal","R","ep",ffirstCrystal->GetXmin(),ffirstCrystal->GetXmax());
//	ffirstCrystal->Draw("same");
	ffirstCrystal->GetParameters(ffirstCrystalPar);


	double fhGausPar[3];
	fhGausPar[1]=2.1565;
	TF1 *fhGaus=new TF1("secondCgaus","gaus",fhGausPar[1]-0.002,fhGausPar[1]+0.0015);
	momentum->Fit("secondCgaus","R","ep",fhGaus->GetXmin(),fhGaus->GetXmax());
	fhGaus->GetParameters(fhGausPar);

		//
	double_t fhCrystalPar[5];
	TF1 *fhCrystal=new TF1("ffirstCrystal","crystalball",fhGausPar[1]-0.0025,fhGausPar[1]+0.0025);
	fhCrystal->SetParameters(fhGausPar[0],fhGausPar[1],fhGausPar[2],1.64,1.1615);
	momentum->Fit("ffirstCrystal","R","ep",fhCrystal->GetXmin(),fhCrystal->GetXmax());
//	fhCrystal->Draw("same");
	fhCrystal->GetParameters(fhCrystalPar);




	double_t fCrystalMomentumPar[15];
	TF1 *fCrystalMomentum=new TF1("fCrystalMomentum","crystalball(0)+crystalball(5)+crystalball(10)",fhCrystal->GetXmin(),fgroundCrystalball->GetXmax());
	std::copy(fgroundCrystalballPar,fgroundCrystalballPar+5,fCrystalMomentumPar);
	std::copy(ffirstCrystalPar,ffirstCrystalPar+5,fCrystalMomentumPar+5);
	std::copy(fhCrystalPar,fhCrystalPar+5,fCrystalMomentumPar+10);

	fCrystalMomentum->SetParameters(fCrystalMomentumPar);
	momentum->Fit("fCrystalMomentum","","",fCrystalMomentum->GetXmin(),fCrystalMomentum->GetXmax());
	//fCrystalMomentum->Draw("same");
	fCrystalMomentum->GetParameters(fCrystalMomentumPar);

	delete fgroundCrystalball;
	delete ffirstGuas;
	return fCrystalMomentum;


}

inline double GetPointingAngle(double DeltaE, double BeamE=2.17568){
	//In pointing measurement,
	// calculate the HRS angle according to H and O seperation

	Double_t amu = 931.494028 * 1.e-3;  // amu to GeV
	Double_t H2O = 18 * amu;
	Double_t TargetH = amu;
	double_t TargetO = 16 * amu;
	Double_t mass_tg = H2O;
	double beamE = BeamE;         // in GeV
	// get the pointing measurement
	double deltaE = DeltaE; // in GeV

	double a = 4.0 * beamE * deltaE;
	double b = 2.0 * deltaE * (TargetH + TargetO)+ 2.0 * beamE * (TargetH - TargetO);
	double c = deltaE * TargetH * TargetO / beamE;
	double sin2theta = ((0. - b) - TMath::Sqrt(b * b - 4.0 * a * c)) / (2.0 * a);
	double HRSAngleReal=TMath::ASin(TMath::Sqrt(sin2theta)) * 2.0 * 180.0 / TMath::Pi();
	return HRSAngleReal;
}

inline Bool_t IsFileExist (const std::string& name) {
	  struct stat buffer;
	  return (stat (name.c_str(), &buffer) == 0);
}

Int_t OpticsGraphicCutH2O(UInt_t runID,TString folder="/home/newdriver/Storage/Research/CRex_Experiment/RasterReplay/Replay/Result") {
	// prepare the data
	TChain *chain=new TChain("T");
	TString rootDir(folder.Data());
	TString HRS="R";
	if(runID>20000){ //RHRS
		if(IsFileExist(Form("%s/prexRHRS_%d_-1.root",rootDir.Data(),runID))){
			std::cout<<"Add File::"<<Form("%s/prexRHRS_%d_-1.root",rootDir.Data(),runID)<<std::endl;
			chain->Add(Form("%s/prexRHRS_%d_-1.root",rootDir.Data(),runID));

			TString filename;
			int16_t split=1;
			filename=Form("%s/prexRHRS_%d_-1_%d.root",rootDir.Data(),runID,split);
			while (IsFileExist(filename.Data())){
				std::cout<<"Add File::"<<filename.Data()<<std::endl;
				chain->Add(filename.Data());
				split++;
				filename=Form("%s/prexRHRS_%d_-1_%d.root",rootDir.Data(),runID,split);
			}
		}else{
			std::cout<<"Looking file :"<<Form("%s/prexRHRS_%d_-1.root",rootDir.Data(),runID)<<std::endl;
		}
	}else{
		HRS="L";
		if(IsFileExist(Form("%s/prexLHRS_%d_-1.root",rootDir.Data(),runID))){
			std::cout<<"Add File::"<<Form("%s/prexLHRS_%d_-1.root",rootDir.Data(),runID)<<std::endl;
			chain->Add(Form("%s/prexLHRS_%d_-1.root",rootDir.Data(),runID));

			TString filename;
			int16_t split=1;
			filename=Form("%s/prexLHRS_%d_-1_%d.root",rootDir.Data(),runID,split);
			while (IsFileExist(filename.Data())){
				std::cout<<"Add File::"<<filename.Data()<<std::endl;
				chain->Add(filename.Data());
				split++;
				filename=Form("%s/prexLHRS_%d_-1_%d.root",rootDir.Data(),runID,split);
			}
		}else{
			std::cout<<"Looking file :"<<Form("%s/prexLHRS_%d_-1.root",rootDir.Data(),runID)<<std::endl;
		}
	}

	if(HRS=="L"){
		generalcut=generalcutL;
	}else{
		generalcut=generalcutR;
	}

	TCanvas *mainPatternCanvas=(TCanvas *)gROOT->GetListOfCanvases()->FindObject("cutPro");
	if(!mainPatternCanvas){
		mainPatternCanvas=new TCanvas("cutPro","cutPro",600,600);
	}else{
		mainPatternCanvas->Clear();
	}
//	TCanvas *mainPatternCanvas=new TCanvas("cut","cut",600,600);
	mainPatternCanvas->Draw();
	TH2F *TargetThPhHH=(TH2F *)gROOT->FindObject("th_vs_ph");
	if(TargetThPhHH) TargetThPhHH->Delete();
	TargetThPhHH=new TH2F("th_vs_ph","th_vs_ph",1000,-0.03,0.03,1000,-0.045,0.045);

	chain->Project(TargetThPhHH->GetName(),Form("%s.gold.th:%s.gold.ph",HRS.Data(),HRS.Data()),generalcut.Data());
	TargetThPhHH->Draw("zcol");

	mainPatternCanvas->Update();
	mainPatternCanvas->ToggleEventStatus();
//	mainPatternCanvas->AddExec("ex", "DynamicCoordinates()");
	mainPatternCanvas->AddExec("ex", "DynamicCanvas()");
	std::cout<<"This is an test point"<<std::endl;
	return 1;
}

void DynamicCanvas(){
	//check which button is clicked
	//if the S button clicked, save the current  cut
	//if the the d button clicked, skip the current hole and continue with the next one

	int event = gPad->GetEvent();
	if (event == kNoEvent)
		return;

	TObject *select = gPad->GetSelected();
	if (!select)
		return;
	if (!select->InheritsFrom(TH2::Class())) {
		gPad->SetUniqueID(0);
		return;
	}
	if (event!=kButton1Down) return;


	TFile *f1=new TFile("test_temp.root","UPDATE");
		assert(f1);

	// link the root tree and check which HRS we are working on
	TChain *chain = (TChain *) gROOT->FindObject("T");
	TString HRS("R");
	TString filename(chain->GetFile()->GetName());
	if (filename.Contains("RHRS")) {
	} else if (filename.Contains("LHRS")) {
		HRS = "L";
	}

	TH2 *h = (TH2*) select;
	gPad->GetCanvas()->FeedbackMode(kTRUE);

	// if the button is clicked
	// get the mouse click position in histogram
	double_t x = (gPad->PadtoX(gPad->AbsPixeltoX(gPad->GetEventX())));
	double_t y = (gPad->PadtoY(gPad->AbsPixeltoY(gPad->GetEventY())));

	// create new canvas
	TCanvas *SieveRecCanvas = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(
			"SieveRecCanvas");
	if (SieveRecCanvas) {
		SieveRecCanvas->Clear();
//		delete SieveRecCanvas->GetPrimitive("Projection");
	} else
		SieveRecCanvas = new TCanvas("SieveRecCanvas", "Projection Canvas",
				1000, 1000);

	SieveRecCanvas->Divide(1, 2);
	SieveRecCanvas->cd(2)->Divide(4, 1);
	//get the hsitogram and start rec
	SieveRecCanvas->cd(2)->cd(2);

	TH2F *selectedSievehh = (TH2F *) gROOT->FindObject("Sieve_Selected_th_ph");
	if (selectedSievehh) {
		selectedSievehh->Clear();
	}
	selectedSievehh = new TH2F("Sieve_Selected_th_ph", "Sieve_Selected_th_ph",
			100, h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(), 100,
			h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax());
	chain->Project(selectedSievehh->GetName(),
			Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()),
			Form("sqrt((%s.gold.th-%f)^2+ (%s.gold.ph-%f)^2)<0.003 && %s",
					HRS.Data(), y, HRS.Data(), x, generalcut.Data()));
	selectedSievehh->SetContour(10);
	selectedSievehh->GetXaxis()->SetTitle(Form("%s.gold.ph", HRS.Data()));
	selectedSievehh->GetYaxis()->SetTitle(Form("%s.gold.th", HRS.Data()));
	selectedSievehh->Draw("CONT LIST");

	SieveRecCanvas->Update(); // update the canvas to let the pattern buffer in root

	// extract the contour
	TObjArray *conts = (TObjArray*) gROOT->GetListOfSpecials()->FindObject(
			"contours");
	if (!conts)
		return;
	TList *lcontour1 = (TList*) conts->At(0);
	if (!lcontour1)
		return;
	TGraph *gc1 = (TGraph*) lcontour1->First();
	if (!gc1)
		return;
	if (gc1->GetN() < 10)
		return;

	//TODO need to change the name of
	TCutG *cutg = new TCutG(Form("hcut_R_%ld",random()),
			gc1->GetN(), gc1->GetX(), gc1->GetY());
	cutg->SetLineWidth(2);
	cutg->SetLineColor(kRed);
	cutg->SetVarX(Form("%s.gold.ph", HRS.Data()));
	cutg->SetVarY(Form("%s.gold.th", HRS.Data()));
	cutg->Draw("same");

	// plot the cut on the canvas
	SieveRecCanvas->cd(2)->cd(1);

	TH2F *patternCheck = (TH2F *) gROOT->FindObject("Sieve_Pattern_Check");
	if (patternCheck) {
		patternCheck->Clear();
	}
	patternCheck = new TH2F("Sieve_Pattern_Check", "Sieve_Pattern_Check",
			h->GetXaxis()->GetNbins(), h->GetXaxis()->GetXmin(),
			h->GetXaxis()->GetXmax(), h->GetYaxis()->GetNbins(),
			h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax());
	chain->Project(patternCheck->GetName(),
			Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()),generalcut.Data());
	patternCheck->Draw("zcol");
	cutg->Draw("same");

	SieveRecCanvas->cd(1);
	SieveRecCanvas->cd(1)->SetLogy();
	// plot the dp and fit
	TH1F *momentum=new TH1F("H2O gold.p",Form("H2O gold.p "),500,2.1,2.2);
	chain->Project(momentum->GetName(),Form("%s.gold.p",HRS.Data()),Form("%s && %s",generalcut.Data(),cutg->GetName()));
	momentum->Draw();

	auto fit = SpectrumFit(momentum);
	fit->Draw("same");

	double FitParameters[15];
	fit->GetParameters(FitParameters);
	double H20DeltaE=FitParameters[1]-FitParameters[11];

	TPaveText *pt = new TPaveText(0.1,0.8,0.3,0.9,"NDC");
	pt->AddText(Form("%1.3f MeV (%1.4f Degree)",1000.0*H20DeltaE,GetPointingAngle(H20DeltaE)));
	pt->Draw("same");

	/*// get the maximum bin, this should be the first excited states
	auto CGroundp=momentum->GetXaxis()->GetBinCenter(momentum->GetMaximumBin());

	momentum->GetXaxis()->SetRangeUser(CGroundp-0.02,CGroundp+0.0044*2);
	momentum->GetXaxis()->SetTitle("gold.p");
	momentum->GetYaxis()->SetTitle("#");

	momentum->Draw();
	double_t fgroudGausPar[3];
	double_t ffirstGuasPar[3];
	TF1 *fgroudGaus=new TF1("groudstatesgaus","gaus",CGroundp-0.0005,CGroundp+0.0005);
	momentum->Fit("groudstatesgaus","R","ep",fgroudGaus->GetXmin(),fgroudGaus->GetXmax());
	fgroudGaus->Draw("same");
	fgroudGaus->GetParameters(fgroudGausPar);

	auto C1stp=2.1565;//CGroundp-0.016504;

	TF1 *ffirstGuas=new TF1 ("firststatesgaus","gaus",C1stp-0.001,C1stp+0.00155);
	momentum->Fit("firststatesgaus","R","ep",ffirstGuas->GetXmin(),ffirstGuas->GetXmax());
	ffirstGuas->Draw("same");
	ffirstGuas->GetParameters(ffirstGuasPar);


	// change the gause fit to cristal ball
	double_t fgroundCrystalballPar[5];
	TF1 *fgroundCrystalball=new TF1("fgroundCrystal","crystalball",fgroudGausPar[1]-0.0030,fgroudGaus->GetXmax()+0.0003);
	fgroundCrystalball->SetParameters(fgroudGausPar[0],fgroudGausPar[1],fgroudGausPar[2],1.64,1.1615);
	momentum->Fit("fgroundCrystal","R","same",fgroundCrystalball->GetXmin(),fgroundCrystalball->GetXmax());
	fgroundCrystalball->GetParameters(fgroundCrystalballPar);
	//fgroundCrystalball->Draw("same");

	double_t ffirstCrystalPar[5];
	TF1 *ffirstCrystal=new TF1("ffirstCrystal","crystalball",ffirstGuasPar[1]-0.0025,ffirstGuas->GetXmax());
	ffirstCrystal->SetParameters(ffirstGuasPar[0],ffirstGuasPar[1],ffirstGuasPar[2],1.64,1.1615);
	momentum->Fit("ffirstCrystal","R","ep",ffirstCrystal->GetXmin(),ffirstCrystal->GetXmax());
	ffirstCrystal->GetParameters(ffirstCrystalPar);
	//	ffirstCrystal->Draw("same");
	// fit together
	double_t fCrystalMomentumPar[10];
	TF1 *fCrystalMomentum=new TF1("fCrystalMomentum","crystalball(0)+crystalball(5)",ffirstCrystal->GetXmin(),fgroundCrystalball->GetXmax());
	std::copy(fgroundCrystalballPar,fgroundCrystalballPar+5,fCrystalMomentumPar);
	std::copy(ffirstCrystalPar,ffirstCrystalPar+5,fCrystalMomentumPar+5);
	fCrystalMomentum->SetParameters(fCrystalMomentumPar);
	momentum->Fit("fCrystalMomentum","","",fCrystalMomentum->GetXmin(),fCrystalMomentum->GetXmax());
	fCrystalMomentum->Draw("same");
	fCrystalMomentum->GetParameters(fCrystalMomentumPar);

	SieveRecCanvas->Update();
	// plot the reconstrcution peak
	TLine *groudposLine=new TLine(fCrystalMomentumPar[1],0,fCrystalMomentumPar[1],fgroudGausPar[0]*1.1);
	groudposLine->SetLineColor(3);
	groudposLine->SetLineWidth(2);
	groudposLine->Draw("same");

	TLine *firstposLine=new TLine(fCrystalMomentumPar[6],0,fCrystalMomentumPar[6],ffirstGuasPar[0]*2.0);
	firstposLine->SetLineColor(3);
	firstposLine->SetLineWidth(2);
	firstposLine->Draw("same");

	TPaveText *pt = new TPaveText(0.1,0.8,0.3,0.9,"NDC");
	double_t deltaE=fCrystalMomentumPar[1]-fCrystalMomentumPar[6];
	pt->AddText(Form("%1.3f MeV (%1.4f Degree)",1000.0*deltaE,GetPointingAngle(deltaE)));
	pt->Draw("same");

	TLatex *t1 = new TLatex(fgroudGausPar[1] + 2 * fgroudGausPar[2],fgroudGausPar[0], Form("P=%2.5fGeV #pm %1.5f", fCrystalMomentumPar[1],fCrystalMomentumPar[2]));
	t1->SetTextSize(0.055);
	t1->SetTextAlign(12);
	t1->SetTextColor(2);
	t1->Draw("same");

	TLatex *t2 = new TLatex(ffirstGuasPar[1]+ffirstGuasPar[2]*0.5,fCrystalMomentumPar[5], Form("P=%2.5fGeV #pm %1.5f", fCrystalMomentumPar[6],fCrystalMomentumPar[7]));
	t2->SetTextSize(0.055);
	t2->SetTextAlign(12);
	t2->SetTextColor(2);
	t2->Draw("same");

	//plot the bigger plot for first excited states
	SieveRecCanvas->cd(2)->cd(3);
	TH1F *groundStats=(TH1F *)momentum->Clone("H2O p.ground");
	groundStats->GetXaxis()->SetRangeUser(fCrystalMomentumPar[1]-0.002,fCrystalMomentumPar[1]+0.002);
	groundStats->Draw();
	//groudposLine->Draw("same");

	SieveRecCanvas->cd(2)->cd(4);
	TH1F *firstStats=(TH1F *)momentum->Clone("H2O p.first");
	firstStats->GetXaxis()->SetRangeUser(fCrystalMomentumPar[6]-0.002,fCrystalMomentumPar[6]+0.002);
	firstStats->Draw();
	//firstposLine->Draw("same");

	SieveRecCanvas->Update();

	// Get Main Canvas and plot the pos on the canvas
	TH2F *hSieveHole = (TH2F *) gROOT->FindObject("sieveholeh");
	if (patternCheck) {
		patternCheck->Clear();
	}
	hSieveHole=new TH2F("sieveholeh", "sieveholeh", 1000, h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(), 1000,
			h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax());
	chain->Project(hSieveHole->GetName(),
				Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()),
				Form("%s && %s",generalcut.Data(),cutg->GetName()));
	TCanvas *SieveMainCanvas = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(
				"cutPro");


	if(SieveMainCanvas){
		SieveMainCanvas->cd();
		cutg->Draw("same");
		TLatex *t1 = new TLatex(hSieveHole->GetMean(1)-0.001,hSieveHole->GetMean(2)+0.001, Form("%1.4fMeV", fCrystalMomentumPar[1]));
		t1->SetTextSize(0.02);
		t1->Draw("same");
		TLatex *t2 = new TLatex(hSieveHole->GetMean(1)-0.001,hSieveHole->GetMean(2)-0.001, Form("%1.3fMeV", 1000.0*(fCrystalMomentumPar[1]-fCrystalMomentumPar[6])));
		t2->SetTextSize(0.02);
		t2->Draw("same");
	}



	hSieveHole->Delete();
	cutg->Delete();*/
}

