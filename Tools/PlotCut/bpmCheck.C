/*
 * bpmCheck.C
 *
 *  Created on: Jun 29, 2020
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

std::map<int,int> runListArray;




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

void bpmCheck(UInt_t runID,TString folder){

	runListArray[1694]=20825;
	runListArray[1695]=20826;
	runListArray[1696]=20827;

	runListArray[2566]=21642;
	runListArray[2565]=21641;
	runListArray[2556]=21632;
	runListArray[2550]=21626;

	TCanvas *canv=new TCanvas("canvas","canvas",1960,1080);
	canv->Divide(2,2);

	canv->Draw();
	auto chainLHRS=LoadRootFiles(runID,folder);
	auto chainRHRS=LoadRootFiles(runListArray[runID],folder);

	TH1F *LHRS_x=new TH1F(Form("LHRS%d_x",runID),Form("LHRS%d_x",runID),1000,-2,3);
	TH1F *LHRS_y=new TH1F(Form("LHRS%d_y",runID),Form("LHRS%d_y",runID),1000,-2,3);
	TH1F *RHRS_x=new TH1F(Form("RHRS%d_x",runListArray[runID]),Form("RHRS%d_x",runListArray[runID]),1000,-2,2);
	TH1F *RHRS_y=new TH1F(Form("RHRS%d_y",runListArray[runID]),Form("RHRS%d_y",runListArray[runID]),1000,-2,2);

	chainLHRS->Project(LHRS_x->GetName(),"targx");
	chainLHRS->Project(LHRS_y->GetName(),"targy");

	chainRHRS->Project(RHRS_x->GetName(),"targx");
	chainRHRS->Project(RHRS_y->GetName(),"targy");


	canv->cd(1);
	LHRS_x->Draw();

	canv->cd(2);
	RHRS_x->Draw();

	canv->cd(3);
	LHRS_y->Draw();

	canv->cd(4);
	RHRS_y->Draw();

	canv->Update();
}

void Getbpm(UInt_t runID, TString folder="/home/newdriver/Storage/Research/PRex_Workspace/PREX-MPDGEM/PRexScripts/Tools/PlotCut/Result/PRex/Cut20200719/rootfiles"){

	// project the parameters
	TH1F *target_x=new TH1F(Form("target_%d_x",runID),Form("target_%d_x",runID),1000,-5,5);
	TH1F *target_y=new TH1F(Form("target_%d_y",runID),Form("target_%d_y",runID),1000,-5,5);

	auto chain= LoadRootFiles(runID,folder);
	chain->Project(target_x->GetName(),"targx");
	chain->Project(target_y->GetName(),"targy");

	TCanvas *canv=new TCanvas("targx","targy",1960,1080);
	canv->Divide(1,2);

	canv->cd(1);
	target_x->Draw();
	target_x->Fit("gaus");
	double meanx=target_x->GetFunction("gaus")->GetParameter(1);
	double sigmax=target_x->GetFunction("gaus")->GetParameter(2);
	target_x->GetXaxis()->SetRangeUser(meanx-8*sigmax,meanx+8*sigmax);
	TLatex *text1=new TLatex(target_x->GetFunction("gaus")->GetParameter(1),target_x->GetFunction("gaus")->GetParameter(0),Form("%f",target_x->GetFunction("gaus")->GetParameter(1)));
	text1->Draw("same");

	canv->cd(2);
	target_y->Draw();
	target_y->Fit("gaus");
	TLatex *text2=new TLatex(target_y->GetFunction("gaus")->GetParameter(1),target_y->GetFunction("gaus")->GetParameter(0),Form("%f",target_y->GetFunction("gaus")->GetParameter(1)));
	target_y->GetXaxis()->SetRangeUser(target_y->GetFunction("gaus")->GetParameter(1)-8*target_y->GetFunction("gaus")->GetParameter(2),target_y->GetFunction("gaus")->GetParameter(1)+8*target_y->GetFunction("gaus")->GetParameter(2));
	text2->Draw("same");
	canv->Update();
	canv->SaveAs(Form("bpmInform%d.jpg",runID));





}
