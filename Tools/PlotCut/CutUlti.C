/*
 * CutUlti.C
 *
 *  Created on: Dec 21, 2019
 *      Author: newdriver
 *      new code used cut the sievehole
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
#include <TSystem.h>
#include <TApplication.h>
#include <map>
#include <TVector2.h>
#include <set>
#include <algorithm>
#include <functional>
int FoilID=0;
int col=3;
int row=1;
int row_min=0;
int row_count=10;
const UInt_t NSieveCol = 13;
const UInt_t NSieveRow = 7;
UInt_t SieveUID=0;

//////////////////////////////////////////////////////////////////////////////
// Work Directory
// cut options
// Need to change
//////////////////////////////////////////////////////////////////////////////

TString generalcut;
TString generalcutR="R.tr.n==1 && R.vdc.u1.nclust==1&& R.vdc.v1.nclust==1 && R.vdc.u2.nclust==1 && R.vdc.v2.nclust==1 && fEvtHdr.fEvtType==1";
TString generalcutL="L.tr.n==1 && L.vdc.u1.nclust==1&& L.vdc.v1.nclust==1 && L.vdc.u2.nclust==1 && L.vdc.v2.nclust==1 && fEvtHdr.fEvtType==1";

//////////////////////////////////////////////////////////////////////////////
// Work Directory
//////////////////////////////////////////////////////////////////////////////
TString HRS="R";
TString WorkDir = "Result/";
TString CutSuf = ".FullCut.root";
TString CutDescFileSufVertex = ".VertexCut.cut";
TString CutDescFileSufDp = ".DpCut.cut";
TString CutDescFileSufSieve = ".SieveCut.%d_%d.cut";
TString RootFileName;     //

inline Bool_t IsFileExist (const std::string& name) {
	  struct stat buffer;
	  return (stat (name.c_str(), &buffer) == 0);
}

Int_t CutUlti(UInt_t runID,UInt_t current_col,TString folder="/home/newdriver/Storage/Research/CRex_Experiment/optReplay/Result/LHRS_AfterThetaPhi") {
	// prepare the data
	TChain *chain=new TChain("T");
	TString rootDir(folder.Data());

	if(runID>20000){ //RHRS
		generalcut=generalcutR;
		if(IsFileExist(Form("%s/prexRHRS_%d_-1.root",rootDir.Data(),runID))){
			std::cout<<"Add File::"<<Form("%s/prexRHRS_%d_-1.root",rootDir.Data(),runID)<<std::endl;
			RootFileName=Form("prexRHRS_%d_-1.root",runID);
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
		generalcut=generalcutL;
		HRS="L";
		if(IsFileExist(Form("%s/prexLHRS_%d_-1.root",rootDir.Data(),runID))){
			RootFileName=Form("prexLHRS_%d_-1.root",runID);
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

	TCanvas *mainPatternCanvas=(TCanvas *)gROOT->GetListOfCanvases()->FindObject("cutPro");
	if(!mainPatternCanvas){
		mainPatternCanvas=new TCanvas("cutPro","cutPro",800,800);
	}else{
		mainPatternCanvas->Clear();
	}
	mainPatternCanvas->ToggleEventStatus();

	mainPatternCanvas->Draw();
	TH2F *TargetThPhHH=(TH2F *)gROOT->FindObject("th_vs_ph");
	if(TargetThPhHH) TargetThPhHH->Delete();
	TargetThPhHH=new TH2F("th_vs_ph","th_vs_ph",1000,-0.03,0.04,1000,-0.06,0.06);

	chain->Project(TargetThPhHH->GetName(),Form("%s.gold.th:%s.gold.ph",HRS.Data(),HRS.Data()),generalcut.Data());
	TargetThPhHH->Draw("zcol");
	// input how start row and how many holes in this row
	col=current_col;
    int nhol = 0;
    std::cout << "How many holes in this No." << col << " column?" << std::endl;
    std::cin >> nhol;
    row_count=nhol;
    if(nhol < 0)return 0;
    std::cout << "min hole id : ";
    int rmin = -1;
    std::cin >> rmin;
    row_min=rmin;
    row=row_min;

    if(rmin < 0)return 0;
    mainPatternCanvas->Update();

	mainPatternCanvas->AddExec("ex", "DynamicCanvas()");
	std::cout<<"This is an test point"<<std::endl;
	return 1;
}


//Recognize the save patter
// save the rec hole to the folder
void SavePatternHole(){
	//search all the holes in this col and save in the folder
	std::cout<<std::endl<<std::endl;
	std::cout<<"*******Save process start ........*******"<<std::endl;
	TChain *chain = (TChain *) gROOT->FindObject("T");
	if(!chain) std::cout<<"[ERROR] CAN NOT FIND TREE"<<std::endl;
	TH2F *h=(TH2F*)gROOT->FindObject("th_vs_ph");
	if(!h)
		h=new TH2F("th_vs_ph1","th_vs_ph1",1000,-0.027,0.03,1000,-0.06,0.06);
	// get the list of the cut file
	//locate the cut file and get the number
	TString CutFileName = WorkDir + RootFileName + CutSuf;
	TFile *cutIO = new TFile(CutFileName, "UPDATE");
	// get list of canvas
	std::map<UInt_t, Int_t> sieveEventCount;
	std::map<UInt_t,TVector2> SieveIDMap;      // map between the UID and the real Col Row ID
	std::map<UInt_t,TVector2> SieveCenterPos; // cenrtrol position of each hole,
	std::map<UInt_t,TCutG *> SieveCut;
	TCutG *cutg;
	TH2F *sievehole_temp;
	for(UInt_t sieve_iter=0 ; sieve_iter<NSieveCol * NSieveRow;sieve_iter++){
		// check the existance of the cut
		cutg = (TCutG *) gROOT->FindObject(
				Form("ahcut_%d", sieve_iter));
		if (cutg) {
			sievehole_temp=new TH2F(Form("sieve_th_ph_cut%d",sieve_iter), Form("sieve_th_ph_cut%d",sieve_iter),
					h->GetXaxis()->GetNbins(), h->GetXaxis()->GetXmin(),
					h->GetXaxis()->GetXmax(), h->GetYaxis()->GetNbins(),
					h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax());
			chain->Project(sievehole_temp->GetName(),Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()),Form("%s && %s",cutg->GetName(),generalcut.Data()));
			sieveEventCount[sieve_iter]=sievehole_temp->GetEntries();
			TVector2 pos(sievehole_temp->GetMean(1),sievehole_temp->GetMean(2)); //get the position, this is used for plot the ID
			SieveIDMap[sieve_iter]=pos;
		}
	}

	// analysis real pos according to relative pos
	// find the three  largest event witch is A B C hole
	// need to check the  difference of LHRS RHRS
	// Declaring the type of Predicate that accepts 2 pairs and return a bool
	typedef std::function<
			bool(std::pair<UInt_t, Int_t>, std::pair<UInt_t, Int_t>)> Comparator;
	// Defining a lambda function to compare two pairs. It will compare two pairs using second field
	Comparator compFunctor =
			[](std::pair<UInt_t, Int_t> elem1 ,std::pair<UInt_t, Int_t> elem2)
			{
				return elem1.second > elem2.second;
			};
	// Declaring a set that will store the pairs using above comparision logic
	std::set<std::pair<UInt_t, Int_t>, Comparator> setOfWords(
			sieveEventCount.begin(), sieveEventCount.end(), compFunctor);

	// finish sort the UIDs according to the event number
	// get the Largeest three, and compare the pos to find A

	// find the first tree elemnts, and sort ythe theta angle to seperate ABC

	UInt_t counter_temp=0;
	for(auto pairElement = setOfWords.begin(); pairElement!= setOfWords.end(); pairElement++){
		if(counter_temp>=3)break;
		std::cout<<"ID  :"<<(pairElement->first)<<std::endl;
		counter_temp++;

	}
//	for (std::pair<UInt_t, Int_t> element : setOfWords){
//
//	}







	cutIO->Close();
}

// dynamic canvas
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

	// need to create  an file used for save the file
	TString CutFileName = WorkDir + RootFileName + CutSuf;
	TFile *cutIO = new TFile(CutFileName, "UPDATE");
	assert(cutIO);
	// link the root tree and check which HRS we are working on
	TChain *chain = (TChain *) gROOT->FindObject("T");
	// check the input keys
	if(event == kKeyPress) {

		std::cout<<"Key Pressed"<<std::endl;
		switch ((int) (gPad->GetEventX())) {
		case '+':
			row++;
//			CurrentStatus();
			break;
		case '-':
			row--;
//			CurrentStatus();
			break;
		case 's':
			std::cout << "Save Button Clicked" << std::endl;
			SavePatternHole();
			break;
		case 'q':
			{std::cout << "Quit Button Clicked" << std::endl;
//			gApplication->Clear();
//			gApplication->Terminate();
			break;}
		default:
			std::cout<<(char)(gPad->GetEventX())<<std::endl;
		}
	} else if (event == kButton1Down) {
		TH2 *h = (TH2*) select;
		gPad->GetCanvas()->FeedbackMode(kTRUE);

		// if the button is clicked
		//Rec the sieve pattern
		// get the mouse click position in histogram
		double_t x = (gPad->PadtoX(gPad->AbsPixeltoX(gPad->GetEventX())));
		double_t y = (gPad->PadtoY(gPad->AbsPixeltoY(gPad->GetEventY())));

		// create new canvas
		TCanvas *SieveRecCanvas = (TCanvas*) gROOT->GetListOfCanvases()->FindObject("SieveRecCanvas");
		if(SieveRecCanvas){
			SieveRecCanvas->Clear();
			delete SieveRecCanvas->GetPrimitive("Projection");
		}else {
			SieveRecCanvas = new TCanvas("SieveRecCanvas","Projection Canvas", 1000,1000);
		}
		SieveRecCanvas->Divide(1,2);
		SieveRecCanvas->cd(2)->Divide(3,1);

		//get the hsitogram and start rec
		SieveRecCanvas->cd(2)->cd(1);

		TH2F *selectedSievehh=(TH2F *)  gROOT->FindObject(Form("Sieve_Selected_th_ph_%d",SieveUID));
		if(selectedSievehh){
			selectedSievehh->Clear();
		}
		selectedSievehh = new TH2F(Form("Sieve_Selected_th_ph_%d",SieveUID),
				Form("Sieve_Selected_th_ph_%d",SieveUID), 100, h->GetXaxis()->GetXmin(),
				h->GetXaxis()->GetXmax(), 100, h->GetYaxis()->GetXmin(),
				h->GetYaxis()->GetXmax());
		chain->Project(selectedSievehh->GetName(),
				Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()),
				Form("sqrt((%s.gold.th-%f)^2+ (%s.gold.ph-%f)^2)<0.003 ",
						HRS.Data(), y, HRS.Data(), x), generalcut.Data());

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
		TList *lcontour1 = (TList*) conts->At(2);
		if (!lcontour1)
			return;
		TGraph *gc1 = (TGraph*) lcontour1->First();
		if (!gc1)
			return;
		if (gc1->GetN() < 10)
			return;

		//TODO need to change the name of the cut file, put an universal ID
		TCutG *cutg = new TCutG(Form("ahcut_%d", SieveUID),gc1->GetN(), gc1->GetX(), gc1->GetY());
		cutg->SetLineColor(kRed);
		cutg->SetVarX(Form("%s.gold.ph", HRS.Data()));
		cutg->SetVarY(Form("%s.gold.th", HRS.Data()));
		cutg->Draw("same");
		cutg->Write();

		SieveRecCanvas->cd(2)->cd(2);
		auto projectx = selectedSievehh->ProjectionX();
		projectx->SetName(Form("project_%d_x",SieveUID));
		projectx->Draw();
		projectx->Fit("gaus");

		SieveRecCanvas->cd(2)->cd(3);
		auto projecty = selectedSievehh->ProjectionY();
		projecty->SetName(Form("project_%d_y",SieveUID));
		projecty->Draw();
		projecty->Fit("gaus");

		// plot the cut on the canvas
		SieveRecCanvas->cd(1);

		TH2F *patternCheck = (TH2F *) gROOT->FindObject("Sieve_Pattern_Check");
		if (patternCheck) {
			patternCheck->Clear();
		}
		patternCheck = new TH2F("Sieve_Pattern_Check", "Sieve_Pattern_Check",
				h->GetXaxis()->GetNbins(), h->GetXaxis()->GetXmin(),
				h->GetXaxis()->GetXmax(), h->GetYaxis()->GetNbins(),
				h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax());
		chain->Project(patternCheck->GetName(),
				Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()));
		patternCheck->Draw("zcol");
		cutg->Draw("same");

		TPaveText *label = new TPaveText(x, y, x + 0.003, y + 0.005, "NB");
		label->AddText(Form("(%d)", SieveUID));
		label->Draw("same");
		row++;
		SieveUID++;
		SieveRecCanvas->Update();
		SieveRecCanvas->Update();

	}
	cutIO->Write();
	cutIO->Close();

}
