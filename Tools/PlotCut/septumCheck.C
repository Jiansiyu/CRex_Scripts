/*
 * septumCheck.C
 *
 *  Created on: Jun 30, 2020
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
#include <TApplication.h>
#include <TLegend.h>

TString prepcut;
TString generalcut;
TString generalcutR="R.tr.n==1";
//TString generalcutR="R.tr.n==1 && R.vdc.u1.nclust==1&& R.vdc.v1.nclust==1 && R.vdc.u2.nclust==1 && R.vdc.v2.nclust==1 && R.gold.p<0.98 && R.gold.p > 0.93 && fEvtHdr.fEvtType==1";
TString generalcutL="L.tr.n==1 && L.vdc.u1.nclust==1&& L.vdc.v1.nclust==1 && L.vdc.u2.nclust==1 && L.vdc.v2.nclust==1 && fEvtHdr.fEvtType==1 && L.gold.p > 2.14 && L.gold.p < 2.2";

std::map<int,int> septumADJrunList;
std::map<int,int> septumNONADJrunList;


std::string OutfName("./SeptumCheck.root");

void loadRunList(){

	septumNONADJrunList[0]=21371;
	septumNONADJrunList[-1]=21381;
	septumNONADJrunList[1]=21380;

	septumADJrunList[0]=20827;
	septumADJrunList[-1]=20826;
	septumADJrunList[1]=20825;
}

inline Bool_t IsFileExist (const std::string& name) {
	  struct stat buffer;
	  return (stat (name.c_str(), &buffer) == 0);
}

TChain *LoadRootFiles(UInt_t runID,
		TString folder =
				"/home/newdriver/Storage/Research/PRex_Experiment/PRex_Replay/replay/Result"){
	TChain *chain=new TChain("T");
	// if the folder itself is and root file
	TString HRS="R";
	if(runID<20000){HRS="L";};

	if(folder.EndsWith(".root")){
		chain->Add(folder.Data());
	}else{
		TString rootDir(folder.Data());
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
	}
	return chain;
}

// get the run list and extract the information on target

int getfRun(TChain *chain){
	TH1F *eventIDhh = new TH1F("eventID", "eventID", 800000, 0, 800000);
	chain->Project(eventIDhh->GetName(), "fEvtHdr.fRun");
	int eventID = int(eventIDhh->GetMean());
	eventIDhh->Delete();
	return eventID;
}

void getGraphicCut(int runID){
	TString HRS="R";
	if(runID<20000){HRS="L";};

	if(HRS=="L"){
		generalcut=generalcutL;
	}else{
		generalcut=generalcutR;
	}

	auto chain=LoadRootFiles(runID);

	TCanvas *mainPatternCanvas=(TCanvas *)gROOT->GetListOfCanvases()->FindObject("cutPro");
	if(!mainPatternCanvas){
		mainPatternCanvas=new TCanvas("cutPro","cutPro",600,600);
	}else{
		mainPatternCanvas->Clear();
	}
	mainPatternCanvas->Draw();
	TH2F *TargetThPhHH=(TH2F *)gROOT->FindObject("th_vs_ph");
	if(TargetThPhHH) TargetThPhHH->Delete();
	TargetThPhHH=new TH2F("th_vs_ph","th_vs_ph",1000,-0.03,0.03,1000,-0.045,0.045);

	chain->Project(TargetThPhHH->GetName(),Form("%s.gold.th:%s.gold.ph",HRS.Data(),HRS.Data()),generalcut.Data());
	TargetThPhHH->Draw("zcol");

	mainPatternCanvas->Update();
	mainPatternCanvas->ToggleEventStatus();
	mainPatternCanvas->AddExec("ex", "DynamicCanvas()");
	std::cout<<"This is an test point"<<std::endl;
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
	// link the root tree and check which HRS we are working on

	TFile *file;
	if(IsFileExist(OutfName.c_str())){
		file=new TFile(OutfName.c_str(),"update");
	}else{
		file=new TFile(OutfName.c_str(),"recreate");
	}

	TChain *chain = (TChain *) gROOT->FindObject("T");
	TString HRS("R");
	TString filename(chain->GetFile()->GetName());
	if (filename.Contains("RHRS")) {
	} else if (filename.Contains("LHRS")) {
		HRS = "L";
	}

	int RunNumber=getfRun(chain);

	TH2 *h = (TH2*) select;
	gPad->GetCanvas()->FeedbackMode(kTRUE);
	// if the button is clicked
	// get the mouse click position in histogram
	double_t x = (gPad->PadtoX(gPad->AbsPixeltoX(gPad->GetEventX())));
	double_t y = (gPad->PadtoY(gPad->AbsPixeltoY(gPad->GetEventY())));

	// create new canvas
	TCanvas *SieveRecCanvas = (TCanvas*) gROOT->GetListOfCanvases()->FindObject(
			Form("SieveRecCanvas_run%d",RunNumber));
	if (SieveRecCanvas) {
		SieveRecCanvas->Clear();
	} else
		SieveRecCanvas = new TCanvas(Form("SieveRecCanvas_run%d",RunNumber),Form("SieveRecCanvas_run%d",RunNumber),
				1000, 1000);

	SieveRecCanvas->Divide(1, 2);
	SieveRecCanvas->cd(2)->Divide(4, 1);
	//get the hsitogram and start rec
	SieveRecCanvas->cd(1);

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
	TCutG *cutg = new TCutG(Form("hcut_R_%d",RunNumber),
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
	cutg->Write();
	SieveRecCanvas->Write();
	sleep(2);

	patternCheck->Delete();
	selectedSievehh->Delete();

	file->Write();
	file->Close();

	h->Delete();
	lcontour1->Delete();
	gc1->Delete();
	cutg->Delete();

}

void septumCheck(){

	loadRunList();
	TFile *cutFileIO=new TFile(OutfName.c_str(),"READ");
	if(cutFileIO->IsZombie()){
		std::cout<<"[ERROR]:: CAN NOT FIND CUT FILE \" "<<OutfName.c_str()<<"\""<<std::endl;
	}


	// read in the file and appy the cut

	std::map<int,TH1F *>sepADJFocalX;
	std::map<int,TH1F *>sepADJFocalY;
	std::map<int,TH1F *>sepADJFocalTheta;
	std::map<int,TH1F *>sepADJFocalPhi;

	std::map<int,TH1F *>sepNonADJFocalX;
	std::map<int,TH1F *>sepNonADJFocalY;
	std::map<int,TH1F *>sepNonADJFocalTheta;
	std::map<int,TH1F *>sepNonADJFocalPhi;

	int maximumX=0;
	int maximumY=0;
	int maximumTheta=0;
	int maximumPhi=0;


	TLegend *lgend=new TLegend(0.1,0.7,0.48,0.9);
	lgend->SetLineWidth(2);

	for (auto item = septumADJrunList.begin();item!=septumADJrunList.end(); item++){
		TString HRS="R";
		if(item->second<20000){HRS="L";};
		sepADJFocalX[item->second]=new TH1F(Form("Run%d_%s.tr.x",item->second,HRS.Data()),Form("Run%d_%s.tr.x",item->first,HRS.Data()),2000,-0.4,0.4);
		sepADJFocalY[item->second]=new TH1F(Form("Run%d_%s.tr.y",item->second,HRS.Data()),Form("Run%d_%s.tr.y",item->second,HRS.Data()),1000,-0.01,0.01);
		sepADJFocalTheta[item->second]=new TH1F(Form("Run%d_%s.tr.th",item->second,HRS.Data()),Form("Run%d_%s.tr.th",item->second,HRS.Data()),2000,-0.06,0.06);
		sepADJFocalPhi[item->second]=new TH1F(Form("Run%d_%s.tr.ph",item->second,HRS.Data()),Form("Run%d_%s.tr.ph",item->second,HRS.Data()),1000,-0.01,0.01);



		sepADJFocalX[item->second]->SetLineColor(8+item->first);
		sepADJFocalY[item->second]->SetLineColor(8+item->first);
		sepADJFocalTheta[item->second]->SetLineColor(8+item->first);
		sepADJFocalPhi[item->second]->SetLineColor(8+item->first);

		lgend->AddEntry(sepADJFocalX[item->second],Form("Septum Adjusted run%d Dp%d", item->second,item->first));
	}
	for (auto item = septumNONADJrunList.begin();item!=septumNONADJrunList.end(); item++){
		TString HRS="R";
		if(item->second<20000){HRS="L";};
		sepNonADJFocalX[item->second]=new TH1F(Form("Run%d_%s.tr.x",item->second,HRS.Data()),Form("Run%d_%s.tr.x",item->first,HRS.Data()),2000,-0.4,0.4);
		sepNonADJFocalY[item->second]=new TH1F(Form("Run%d_%s.tr.y",item->second,HRS.Data()),Form("Run%d_%s.tr.y",item->second,HRS.Data()),1000,-0.01,0.01);
		sepNonADJFocalTheta[item->second]=new TH1F(Form("Run%d_%s.tr.th",item->second,HRS.Data()),Form("Run%d_%s.tr.th",item->second,HRS.Data()),2000,-0.06,0.06);
		sepNonADJFocalPhi[item->second]=new TH1F(Form("Run%d_%s.tr.ph",item->second,HRS.Data()),Form("Run%d_%s.tr.ph",item->second,HRS.Data()),1000,-0.01,0.01);


		sepNonADJFocalX[item->second]->SetLineColor(2+item->first);
		sepNonADJFocalY[item->second]->SetLineColor(2+item->first);
		sepNonADJFocalTheta[item->second]->SetLineColor(2+item->first);
		sepNonADJFocalPhi[item->second]->SetLineColor(2+item->first);

		lgend->AddEntry(sepNonADJFocalX[item->second],Form("Septum Constant run%d Dp%d", item->second,item->first));
		}


	for (auto item = septumADJrunList.begin();item!=septumADJrunList.end(); item++){
		// check the existance of the cut
		int RunNumber=item->second;

		TString HRS="R";
		if(item->second<20000){HRS="L";};

		auto cutg=(TCutG *) gROOT->FindObject(Form("hcut_R_%d",RunNumber));
		if (cutg){
			auto chain = LoadRootFiles(item->second);
			//get the momentum and cut on the peak of the momentum

			double SieveMomFitPar[3];
			TH1F *sieveholemomentum=new TH1F(Form("hcut_R_%d_mom", RunNumber),Form("hcut_R_%d_mom", RunNumber),600,0.94,1.0);
			chain->Project(sieveholemomentum->GetName(),Form("%s.gold.p",HRS.Data()),Form("%s",cutg->GetName()));
			sieveholemomentum->Fit("gaus","","",sieveholemomentum->GetBinCenter(sieveholemomentum->GetMaximumBin())-0.001,sieveholemomentum->GetBinCenter(sieveholemomentum->GetMaximumBin())+0.001);
			sieveholemomentum->GetFunction("gaus")->GetParameters(SieveMomFitPar);
			std::string momCut(Form("abs(%s.gold.p-%f)<%f*%f",HRS.Data(),SieveMomFitPar[1],3.0,SieveMomFitPar[2]));
			sieveholemomentum->Delete();

			chain->Project(sepADJFocalX[item->second]->GetName(),Form("%s.tr.x",HRS.Data()),Form("%s && %s",cutg->GetName(),momCut.c_str()));
			chain->Project(sepADJFocalY[item->second]->GetName(),Form("%s.tr.y",HRS.Data()),Form("%s && %s",cutg->GetName(),momCut.c_str()));
			chain->Project(sepADJFocalTheta[item->second]->GetName(),Form("%s.tr.th",HRS.Data()),Form("%s && %s",cutg->GetName(),momCut.c_str()));
			chain->Project(sepADJFocalPhi[item->second]->GetName(),Form("%s.tr.ph",HRS.Data()),Form("%s && %s",cutg->GetName(),momCut.c_str()));

			if (maximumX < sepADJFocalX[item->second]->GetMaximumBin())maximumX = sepADJFocalX[item->second]->GetMaximumBin();
			if (maximumY < sepADJFocalY[item->second]->GetMaximumBin())maximumY = sepADJFocalY[item->second]->GetMaximumBin();
			if (maximumTheta < sepADJFocalTheta[item->second]->GetMaximumBin())maximumTheta = sepADJFocalTheta[item->second]->GetMaximumBin();
			if (maximumPhi < sepADJFocalPhi[item->second]->GetMaximumBin())maximumPhi = sepADJFocalPhi[item->second]->GetMaximumBin();

		}
	}

	for (auto item = septumNONADJrunList.begin();item!=septumNONADJrunList.end(); item++){
		// check the existance of the cut
		int RunNumber=item->second;

		TString HRS="R";
		if(item->second<20000){HRS="L";};

		auto cutg=(TCutG *) gROOT->FindObject(Form("hcut_R_%d",RunNumber));
		if (cutg){

			auto chain = LoadRootFiles(item->second);
			double SieveMomFitPar[3];
			TH1F *sieveholemomentum=new TH1F(Form("hcut_R_%d_mom", RunNumber),Form("hcut_R_%d_mom", RunNumber),600,0.94,1.0);
			chain->Project(sieveholemomentum->GetName(),Form("%s.gold.p",HRS.Data()),Form("%s",cutg->GetName()));
			sieveholemomentum->Fit("gaus","","",sieveholemomentum->GetBinCenter(sieveholemomentum->GetMaximumBin())-0.001,sieveholemomentum->GetBinCenter(sieveholemomentum->GetMaximumBin())+0.001);
			sieveholemomentum->GetFunction("gaus")->GetParameters(SieveMomFitPar);
			std::string momCut(Form("abs(%s.gold.p-%f)<%f*%f",HRS.Data(),SieveMomFitPar[1],3.0,SieveMomFitPar[2]));
			sieveholemomentum->Delete();

			chain->Project(sepNonADJFocalX[item->second]->GetName(),Form("%s.tr.x",HRS.Data()),Form("%s && %s",cutg->GetName(),momCut.c_str()));
			chain->Project(sepNonADJFocalY[item->second]->GetName(),Form("%s.tr.y",HRS.Data()),Form("%s && %s",cutg->GetName(),momCut.c_str()));
			chain->Project(sepNonADJFocalTheta[item->second]->GetName(),Form("%s.tr.th",HRS.Data()),Form("%s && %s",cutg->GetName(),momCut.c_str()));
			chain->Project(sepNonADJFocalPhi[item->second]->GetName(),Form("%s.tr.ph",HRS.Data()),Form("%s && %s",cutg->GetName(),momCut.c_str()));

			if (maximumX < sepNonADJFocalX[item->second]->GetMaximumBin())maximumX = sepNonADJFocalX[item->second]->GetMaximumBin();
			if (maximumY < sepNonADJFocalY[item->second]->GetMaximumBin())maximumY = sepNonADJFocalY[item->second]->GetMaximumBin();
			if (maximumTheta < sepNonADJFocalTheta[item->second]->GetMaximumBin())maximumTheta = sepNonADJFocalTheta[item->second]->GetMaximumBin();
			if (maximumPhi < sepNonADJFocalPhi[item->second]->GetMaximumBin())maximumPhi = sepNonADJFocalPhi[item->second]->GetMaximumBin();


		}
	}

	TCanvas *canv=new TCanvas("Septum_Focal_Variable_Check","Septum_Focal_Variable_Check",1960,1280);
	canv->Divide(2,2);
	canv->Draw();

	canv->cd(1);
	for (auto item = sepNonADJFocalX.begin(); item!=sepNonADJFocalX.end();item++){
		item->second->GetYaxis()->SetRangeUser(0,20000);
		if (item == sepNonADJFocalX.begin()){
			item->second->Draw();
		}else{

			item->second->Draw("same");
		}
	}

	for (auto item = sepADJFocalX.begin(); item!=sepADJFocalX.end();item++){
//			item->second->GetYaxis()->SetRangeUser(0,maximumX*2);
			item->second->Draw("same");
	}

	lgend->Draw("same");

	canv->cd(2);
	for (auto item = sepNonADJFocalY.begin(); item!=sepNonADJFocalY.end();item++){
		item->second->GetYaxis()->SetRangeUser(0,5000);
		if (item == sepNonADJFocalY.begin()){
			item->second->Draw();
		}else{

			item->second->Draw("same");
		}
	}

	for (auto item = sepADJFocalY.begin(); item!=sepADJFocalY.end();item++){
//			item->second->GetYaxis()->SetRangeUser(0,10000);
			item->second->Draw("same");
	}


	canv->cd(3);
	for (auto item = sepNonADJFocalTheta.begin(); item!=sepNonADJFocalTheta.end();item++){
		item->second->GetYaxis()->SetRangeUser(0,5000);
		if (item == sepNonADJFocalTheta.begin()){
			item->second->Draw();
		}else{

			item->second->Draw("same");
		}
	}

	for (auto item = sepADJFocalTheta.begin(); item!=sepADJFocalTheta.end();item++){
//		    item->second->GetYaxis()->SetRangeUser(0,maximumTheta*2);
			item->second->Draw("same");
	}

	canv->cd(4);
	for (auto item = sepNonADJFocalPhi.begin(); item!=sepNonADJFocalPhi.end();item++){
		item->second->GetYaxis()->SetRangeUser(0,1400);
		if (item == sepNonADJFocalPhi.begin()){
			item->second->Draw();
		}else{

			item->second->Draw("same");
		}
	}

	for (auto item = sepADJFocalPhi.begin(); item!=sepADJFocalPhi.end();item++){
//			item->second->GetYaxis()->SetRangeUser(0,maximumPhi*2);
			item->second->Draw("same");
	}

	canv->Update();
}
