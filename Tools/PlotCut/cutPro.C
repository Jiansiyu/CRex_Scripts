/*
 * cutPro.C
 *
 *  Created on: Dec 12, 2019
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

struct SievePos{
	SievePos(int16_t row=0, int16_t col=3){
		sieve_col=col;
		sieve_row=row;
	};
	int16_t sieve_row;
	int16_t sieve_col;
	SievePos GetNext(){
		return this;
	}
	SievePos GetCurrent(){
		return this;
	}
private:
	const int16_t sieve_row_min;
	const int16_t sieve_row_max;
	const int16_t sieve_col_min;
	const int16_t sieve_col_max;
};

SievePos sievePosition;

inline Bool_t IsFileExist (const std::string& name) {
	  struct stat buffer;
	  return (stat (name.c_str(), &buffer) == 0);
}

Int_t cutPro(UInt_t runID,TString folder="/home/newdriver/Storage/Research/CRex_Experiment/optReplay/Result/AfterCorrection") {
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


	TCanvas *mainPatternCanvas=new TCanvas("cut","cut",600,600);
	mainPatternCanvas->ToggleEventStatus();

	TH2F *TargetThPhHH=new TH2F("th vs ph","th vs ph",1000,-0.027,0.02,1000,-0.043,0.043);
	chain->Project(TargetThPhHH->GetName(),Form("%s.gold.th:%s.gold.ph",HRS.Data(),HRS.Data()));
	TargetThPhHH->Draw("zcol");

	auto box_save=new TBox(0.2,0.2,0.8,0.3);
	box_save->SetFillColor(5);
	box_save->Draw("same");

	mainPatternCanvas->AddExec("ex", "DynamicCoordinates()");

	std::cout<<"This is an test point"<<std::endl;
	return 1;
}

//Recognize the save patter
void HoleContourRec(int Sieve_row, int Sieve_col){

};


void DynamicCoordinates()
{
	// check button clicked
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

	// check the input keys
	if(event == kKeyPress) {

		std::cout<<"Key Pressed"<<std::endl;
		switch ((int) (gPad->GetEventX())) {
		case 'd':
			std::cout << "Delete Button Clicked" << std::endl;
			break;
		case 's':
			std::cout << "Save Button Clicked" << std::endl;
			break;
		default:
			std::cout<<(char)(gPad->GetEventX())<<std::endl;
		}
	}else if (event==kButton1Down) {
		std::cout<<"Button 1 Down"<<std::endl;
	}


	/*
	TH2 *h = (TH2*) select;
	gPad->GetCanvas()->FeedbackMode(kTRUE);

	// get the mouse click position in histogram
	double_t x = (gPad->PadtoX(gPad->AbsPixeltoX(gPad->GetEventX())));
	double_t y = (gPad->PadtoY(gPad->AbsPixeltoY(gPad->GetEventY())));

	// link the root tree and check which HRS we are working on
	TChain *chan = (TChain *) gROOT->FindObject("T");
	TString HRS("R");
	TString filename(chan->GetFile()->GetName());
	if (filename.Contains("RHRS")) {
	} else if(filename.Contains("LHRS")){
		HRS = "L";
	}

	//create or set the new canvas c2
	TCanvas *c2 = (TCanvas*) gROOT->GetListOfCanvases()->FindObject("c2");
	if (c2) {
		c2->Clear();
		delete c2->GetPrimitive("Projection");
	} else
		c2 = new TCanvas("c2", "Projection Canvas", 710, 10, 700, 500);

	c2->Divide(1, 2);
	c2->cd(2)->Divide(3,1);


	c2->cd(2)->cd(1);
	// project the data with cut, and check the counter
	TH2F *patternCut = (TH2F *) gROOT->FindObject("Target_th_ph");
	if (patternCut)
		patternCut->Clear();
	patternCut = new TH2F("Target_th_ph", "Target_th_ph", 100,
			h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(), 100,
			h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax());
	chan->Project(patternCut->GetName(),
			Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()),
			Form("sqrt((%s.gold.th-%f)^2+ (%s.gold.ph-%f)^2)<0.003 ",
					HRS.Data(), y, HRS.Data(), x));

	patternCut->SetContour(10);
	patternCut->GetXaxis()->SetTitle("R.gold.ph");
	patternCut->GetYaxis()->SetTitle("R.gold.th");
	patternCut->Draw("CONT LIST");

	c2->Update();

	TObjArray *conts = (TObjArray*) gROOT->GetListOfSpecials()->FindObject(
			"contours");
	if (!conts)
		return;
	TList *lcontour1 = (TList*) conts->At(2);
	if (!lcontour1)
		return;
	TGraph *gc1 = (TGraph*) lcontour1->First();
	if (!gc1)
		return;
	if (gc1->GetN() < 10)
		return;

	TCutG *cutg = new TCutG("cutg", gc1->GetN(), gc1->GetX(), gc1->GetY());
	cutg->SetLineColor(kRed);
	cutg->Draw("same");
	cutg->SetName("test");
	cutg->SetVarX("R.gold.ph");
	cutg->SetVarY("R.gold.th");

	c2->cd(2)->cd(2);
	auto projectx = patternCut->ProjectionX();
	projectx->Draw();
	projectx->Fit("gaus");

	c2->cd(2)->cd(3);
	auto projecty = patternCut->ProjectionY();
	projecty->Draw();
	projecty->Fit("gaus");

	c2->cd(1);
	TH2F *patternCheck = new TH2F("Sieve_Pattern_Check", "Sieve_Pattern_Check",
			h->GetXaxis()->GetNbins(), h->GetXaxis()->GetXmin(),
			h->GetXaxis()->GetXmax(), h->GetYaxis()->GetNbins(),
			h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax());
	chan->Project(patternCheck->GetName(),Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()));
	patternCheck->Draw("zcol");
	cutg->Draw("same");
	c2->Update();
	*/
}
