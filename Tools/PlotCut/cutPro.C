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
#include <TFile.h>
#include <fstream>
#include <TSystem.h>
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
//TString generalcut="R.tr.n==1 && R.vdc.u1.nclust==1&& R.vdc.v1.nclust==1 && R.vdc.u2.nclust==1 && R.vdc.v2.nclust==1 && fEvtHdr.fEvtType==1";
TString generalcut="L.tr.n==1 && L.vdc.u1.nclust==1&& L.vdc.v1.nclust==1 && L.vdc.u2.nclust==1 && L.vdc.v2.nclust==1 && fEvtHdr.fEvtType==1";

//////////////////////////////////////////////////////////////////////////////
// Work Directory
//////////////////////////////////////////////////////////////////////////////
TString WorkDir = "Result/";

TString CutSuf = ".FullCut.root";
TString CutDescFileSufVertex = ".VertexCut.cut";
TString CutDescFileSufDp = ".DpCut.cut";
TString CutDescFileSufSieve = ".SieveCut.%d_%d.cut";
TString RootFileName;

inline Bool_t IsFileExist (const std::string& name) {
	  struct stat buffer;
	  return (stat (name.c_str(), &buffer) == 0);
}

Int_t cutPro(UInt_t runID,UInt_t current_col,TString folder="/home/newdriver/Storage/Research/CRex_Experiment/optReplay/Result/LHRS_AfterThetaPhi") {
	// prepare the data
	TChain *chain=new TChain("T");
	TString rootDir(folder.Data());
	TString HRS="R";
	if(runID>20000){ //RHRS
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
		mainPatternCanvas=new TCanvas("cutPro","cutPro",600,600);
	}else{
		mainPatternCanvas->Clear();
	}
//	TCanvas *mainPatternCanvas=new TCanvas("cut","cut",600,600);
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
	mainPatternCanvas->ToggleEventStatus();
//	mainPatternCanvas->AddExec("ex", "DynamicCoordinates()");
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
	std::cout<<"	Searching for holes in col ("<<col<<")"<<std::endl;

	TCanvas *SaveCheckCanvas=(TCanvas *) gROOT->GetListOfCanvases()->FindObject("SieveSaveCheck");
	if(!SaveCheckCanvas){
		SaveCheckCanvas =new TCanvas("SieveSaveCheck","SieveSaveCheck",1000,1000);
	}else{
		SaveCheckCanvas->Clear();
	}
	SaveCheckCanvas->Divide(1,2);
	SaveCheckCanvas->cd(1)->Divide(2,1);
	SaveCheckCanvas->cd(2)->Divide(row_count,1);
	TString CutFileName = WorkDir + RootFileName + CutSuf;
	TString TempString(Form(CutDescFileSufSieve.Data(), FoilID, col));
	TString PlotDir(RootFileName + Form(".hcut_R_%d_%d/", FoilID, col));
	TString CutDescName = WorkDir + RootFileName + TempString;
	// prepare the filename and folder
	//TString CutDescName = WorkDir + RootFileName + CutDescFileSufDp;
	//TString CutFileName = WorkDir + RootFileName + CutSuf; // used for save the root file and cut file
	TFile *f1=new TFile(CutFileName,"UPDATE");
	assert(f1);

	std::fstream cutdesc(CutDescName, std::ios_base::out);
	assert(cutdesc.is_open());

	SaveCheckCanvas->cd(1)->cd(1);
	// attached the root file and plot the canvas used for check
	TChain *chain = (TChain *) gROOT->FindObject("T");
	if(!chain) std::cout<<"[ERROR] CAN NOT FIND TREE"<<std::endl;
	TH2F *h=(TH2F*)gROOT->FindObject("th_vs_ph");
	if(!h)
		h=new TH2F("th_vs_ph1","th_vs_ph1",1000,-0.027,0.03,1000,-0.06,0.06);

	TH2F *SavesieveCheck = new TH2F("SaveSieveCheck", "SaveSieveCheck",
			h->GetXaxis()->GetNbins(), h->GetXaxis()->GetXmin(),
			h->GetXaxis()->GetXmax(), h->GetYaxis()->GetNbins(),
			h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax());

	TString HRS("R");
	TString filename(chain->GetFile()->GetName());
	if (filename.Contains("RHRS")) {
	} else if (filename.Contains("LHRS")) {
		HRS = "L";
	}
	chain->Project(SavesieveCheck->GetName(),Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()),Form("%s",generalcut.Data()));
	//chain->Draw(Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()),Form("%s",generalcut.Data()),"zcol");
	SavesieveCheck->Draw("zcol");
	// search how many cut file exist in the root buffer
	for (int i = 0 ; i < row_min; i ++)cutdesc<<"fEvtHdr.fRun==0"<< std::endl;
	TH2F *sievehole[row_count];
	TH1F *sieveholemomentum[row_count];
	TF1  *sieveholemomentumGausFit[row_count];
	TCutG *cutg;
	for (int row_iter = row_min; row_iter<row_min+row_count;row_iter++){
		cutg=(TCutG *) gROOT->FindObject(Form("hcut_R_%d_%d_%d",FoilID,col,row_iter));
		if(cutg){

			SaveCheckCanvas->cd(2)->cd(row_iter - row_min + 1);
			sievehole[row_iter-row_min] = new TH2F(
					Form("hcut_R_%d_%d_%d_hh", FoilID, col, row_iter),
					Form("hcut_R_%d_%d_%d_hh", FoilID, col, row_iter),
					h->GetXaxis()->GetNbins(), h->GetXaxis()->GetXmin(),
					h->GetXaxis()->GetXmax(), h->GetYaxis()->GetNbins(),
					h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax());
			chain->Project(sievehole[row_iter-row_min]->GetName(), Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()),Form("%s && %s",cutg->GetName(),generalcut.Data()));

			sieveholemomentum[row_iter-row_min]=new TH1F(Form("hcut_R_%d_%d_%d_h_momentum", FoilID, col, row_iter),Form("hcut_R_%d_%d_%d_momentum", FoilID, col, row_iter),400,2.1,2.2);
			chain->Project(sieveholemomentum[row_iter-row_min]->GetName(),Form("%s.gold.p",HRS.Data()),Form("%s && %s",cutg->GetName(),generalcut.Data()));
			sieveholemomentumGausFit[row_iter-row_min]=new TF1(Form("1ststatesDpgaushcut_R_%d_%d_%d", FoilID, col, row_iter),"gaus",
					sieveholemomentum[row_iter-row_min]->GetXaxis()->GetBinCenter(sieveholemomentum[row_iter-row_min]->GetMaximumBin())-0.002,
					sieveholemomentum[row_iter-row_min]->GetXaxis()->GetBinCenter(sieveholemomentum[row_iter-row_min]->GetMaximumBin())+0.002);
			sieveholemomentumGausFit[row_iter-row_min]->SetParameter(1,sieveholemomentum[row_iter-row_min]->GetXaxis()->GetBinCenter(sieveholemomentum[row_iter-row_min]->GetMaximumBin()));
			sieveholemomentum[row_iter - row_min]->Fit(
					Form("1ststatesDpgaushcut_R_%d_%d_%d", FoilID, col,
							row_iter), "R", "ep",
					sieveholemomentumGausFit[row_iter-row_min]->GetXmin(),
					sieveholemomentumGausFit[row_iter-row_min]->GetXmax());

			sieveholemomentum[row_iter - row_min]->Draw();
			sieveholemomentumGausFit[row_iter-row_min]->Draw("same");

			auto groudpcenter=sieveholemomentumGausFit[row_iter-row_min]->GetParameter(1);
			auto groudpsigma=sieveholemomentumGausFit[row_iter-row_min]->GetParameter(2);
			TPaveText *pt = new TPaveText(0.1,0.6,0.9,0.9,"NDC");
			    pt->AddText(Form("Mean %2.5f \n sigma %2.5f  (rec peak %2.5f)", groudpcenter, groudpsigma,sieveholemomentum[row_iter-row_min]->GetXaxis()->GetBinCenter(sieveholemomentum[row_iter-row_min]->GetMaximumBin())));
			    pt->Draw("same");

//			chain->Project(sievehole[row_iter-row_min]->GetName(), Form("L.gold.th:L.gold.ph", HRS.Data(), HRS.Data()));

			SaveCheckCanvas->cd(1)->cd(1);
			cutg->SetName(Form("hcut_R_%d_%d_%d", FoilID, col, row_iter));
			cutg->SetVarX(Form("%s.gold.ph",HRS.Data()));
			cutg->SetVarY(Form("%s.gold.th",HRS.Data()));
			cutg->SetLineColor(kMagenta);
			cutg->SetLineWidth(2);
			cutg->Draw("PL same");

			//plot the momentum and apply cut on the momentum



			cutg->Write("", TObject::kOverwrite); // Overwrite old cut
			if(groudpsigma>0.0008)groudpsigma=0.0008;
			cutdesc << Form("hcut_R_%d_%d_%d", FoilID, col, row_iter) << " && ";
			cutdesc << (const char*)generalcut <<" && "
					<<Form("abs(%s.gold.p-%f)<3*%f",HRS.Data(),groudpcenter,groudpsigma)
					<< std::endl;

			SaveCheckCanvas->cd(1)->cd(2);
			sievehole[row_iter-row_min]->Draw("same");
			//sievehole[row_iter-row_min]->Delete();
			cutg->Draw("same");

		}else{
			//if the cut does not exist,then write the cut
			cutdesc << "fEvtHdr.fRun==0" << std::endl;
		}

	}

	for(int i = row_min+row_count; i < NSieveRow; i++)
		cutdesc << "fEvtHdr.fRun==0" << std::endl;


	f1->Write();
	f1->ls();
//	SavesieveCheck->Clear();
	//SavesieveCheck->Delete();
	cutdesc.close();
}


void CurrentStatus(){
	std::cout<<"***  's' save the current cut to file *****"<<std::endl;
	std::cout<<"***  '-' decrease one hole ID         *****"<<std::endl;
	std::cout<<"***  '+' increase one hole ID         *****"<<std::endl;
	std::cout<<"Currently working on ::"<<std::endl;
	std::cout<<"	hole col :"<<col<<std::endl;
	std::cout<<"    hole row :"<<row<< "  minID:"<<row_min<<"  count:"<<row_count<<std::endl;
	std::cout<<">>>Please click the center of the hole ("<<col<<" ,"<<row<<" ) ...."<<std::endl;
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


	// link the root tree and check which HRS we are working on
	TChain *chain = (TChain *) gROOT->FindObject("T");
	TString HRS("R");
	TString filename(chain->GetFile()->GetName());
	if (filename.Contains("RHRS")) {
	} else if (filename.Contains("LHRS")) {
		HRS = "L";
	}

	// check the input keys
	if(event == kKeyPress) {

		std::cout<<"Key Pressed"<<std::endl;
		switch ((int) (gPad->GetEventX())) {
		case '+':
			row++;
			CurrentStatus();
			break;
		case '-':
			row--;
			CurrentStatus();
			break;
		case 's':
			std::cout << "Save Button Clicked" << std::endl;
			SavePatternHole();
			break;
		case 'q':
			{std::cout << "Quit Button Clicked" << std::endl;
			gApplication->Clear();
			gApplication->Terminate();
			break;}
		default:
			std::cout<<(char)(gPad->GetEventX())<<std::endl;
		}
	}else
	if (event==kButton1Down) {
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
		}else
			SieveRecCanvas = new TCanvas("SieveRecCanvas","Projection Canvas", 1000,1000);

			SieveRecCanvas->Divide(1,2);
			SieveRecCanvas->cd(2)->Divide(3,1);

			//get the hsitogram and start rec
			SieveRecCanvas->cd(2)->cd(1);

			TH2F *selectedSievehh=(TH2F *)  gROOT->FindObject("Sieve_Selected_th_ph");
			if(selectedSievehh){
				selectedSievehh->Clear();
			}
			selectedSievehh = new TH2F("Sieve_Selected_th_ph",
					"Sieve_Selected_th_ph",
					100,
					h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(),
					100,
					h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
			chain->Project(selectedSievehh->GetName(),
					Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()),
					Form("sqrt((%s.gold.th-%f)^2+ (%s.gold.ph-%f)^2)<0.003 ",
							HRS.Data(), y, HRS.Data(), x),generalcut.Data());
			selectedSievehh->SetContour(10);
			selectedSievehh->GetXaxis()->SetTitle(Form("%s.gold.ph",HRS.Data()));
			selectedSievehh->GetYaxis()->SetTitle(Form("%s.gold.th",HRS.Data()));
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

			//TODO need to change the name of
			TCutG *cutg = new TCutG(Form("hcut_R_%d_%d_%d", FoilID, col, row), gc1->GetN(), gc1->GetX(), gc1->GetY());
			cutg->SetLineColor(kRed);
			cutg->SetName(Form("hcut_R_%d_%d_%d", FoilID, col, row));
			cutg->SetVarX(Form("%s.gold.ph",HRS.Data()));
			cutg->SetVarY(Form("%s.gold.th",HRS.Data()));
			cutg->Draw("same");

			SieveRecCanvas->cd(2)->cd(2);
			auto projectx = selectedSievehh->ProjectionX();
			projectx->Draw();
			projectx->Fit("gaus");

			SieveRecCanvas->cd(2)->cd(3);
			auto projecty = selectedSievehh->ProjectionY();
			projecty->Draw();
			projecty->Fit("gaus");

			// plot the cut on the canvas
			SieveRecCanvas->cd(1);

			TH2F *patternCheck=(TH2F *)  gROOT->FindObject("Sieve_Pattern_Check");
			if(patternCheck){
				patternCheck->Clear();
			}
			patternCheck = new TH2F("Sieve_Pattern_Check",
					"Sieve_Pattern_Check", h->GetXaxis()->GetNbins(),
					h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(),
					h->GetYaxis()->GetNbins(), h->GetYaxis()->GetXmin(),
					h->GetYaxis()->GetXmax());
			chain->Project(patternCheck->GetName(),
					Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()));
			patternCheck->Draw("zcol");
			cutg->Draw("same");

			TPaveText *label = new TPaveText(x,y,x+0.003,y+0.005,"NB");;
			label->AddText(Form("(%d %d)",col,row));
			label->Draw("same");

			row++;
			SieveRecCanvas->Update();
			SieveRecCanvas->Update();

	}
}

/*
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

 if (event==kButton1Down) {
		std::cout<<"Button 1 Down"<<std::endl;




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
 }
}*/
