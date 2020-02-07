/*
 * plotTemp.C
 *
 *  Created on: Jan 11, 2020
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

inline Bool_t IsFileExist (const std::string& name) {
	  struct stat buffer;
	  return (stat (name.c_str(), &buffer) == 0);
}

// plot cut
void plotTemp(UInt_t runID,TString folder="/home/newdriver/Storage/Research/CRex_Experiment/optReplay/Result", double_t groundp=0.9476, double_t firstp=0.94325){

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

	//plot check the plot
	TH1F *bpmAX=new TH1F("BeamPosAX","BeamPosAX",100,-0.8,0.8);  // beam position  on beam A
	TH1F *bpmAY=new TH1F("BeamPosAY","BeamPosAY",100,-0.8,0.8);

	TH1F *bpmBX = new TH1F("BeamPosBX", "BeamPosBX", 100, -0.8, 0.8); // beam position  on beam A
	TH1F *bpmBY = new TH1F("BeamPosBY", "BeamPosBY", 100, -0.8, 0.8);

	TH1F *bpmCX = new TH1F("BeamPosCX", "BeamPosCX", 100, -0.8, 0.8); // beam position  on beam A
	TH1F *bpmCY = new TH1F("BeamPosCY", "BeamPosCY", 100, -0.8, 0.8);

	TH1F *bpmDX = new TH1F("BeamPosDX", "BeamPosDX", 100, -0.8, 0.8); // beam position  on beam A
	TH1F *bpmDY = new TH1F("BeamPosDY", "BeamPosDY", 100, -0.8, 0.8);


	TH1F *bpmEX=new TH1F("BeamPosEX","BeamPosEX",100,-0.8,0.8);  // beam position  on beam E
	TH1F *bpmEY=new TH1F("BeamPosEY","BeamPosEY",100,-0.8,0.8);

	TH1F *bpmTargetX=new TH1F("BeamTargetX","BeamTargetX",100,-0.5,0.5);  // beam position  on beam E
//	TH1F *bpmTargetY=new TH1F("BeamTargetEY","BeamPosEY",100,-0.5,0.5);


	chain->Project(bpmAX->GetName(),"IPM1H04A.XPOS");
	chain->Project(bpmAY->GetName(),"IPM1H04A.YPOS");

	chain->Project(bpmBX->GetName(),"IPM1H04B.XPOS");

	chain->Project(bpmCX->GetName(),"IPM1H04C.XPOS");

	chain->Project(bpmDX->GetName(),"IPM1H04D.XPOS");


	chain->Project(bpmEX->GetName(),"IPM1H04E.XPOS");

	chain->Project(bpmEY->GetName(),"IPM1H04E.YPOS");

	chain->Project(bpmTargetX->GetName(), "IPM1H04A.XPOS-IPM1H04E.XPOS");

	bpmAX->SetLineColor(2);
	bpmBX->SetLineColor(3);
	bpmCX->SetLineColor(4);
	bpmDX->SetLineColor(5);
	bpmEX->SetLineColor(6);
//	bpmAX->SetLineColor(2);
//	bpmAY->SetLineColor(3);
//	bpmEX->SetLineColor(4);
//	bpmEY->SetLineColor(6);

	TCanvas *beamPosition = new TCanvas("Beam Pos check", "Beam Pos check", 1100,1000);
	beamPosition->cd();
	bpmAX->Draw();
	bpmBX->Draw("same");
	bpmCX->Draw("same");
	bpmDX->Draw("same");
	bpmEX->Draw("same");

//	bpmAX->Draw("");
//	bpmAY->Draw("SAME");
//	bpmEX->Draw("SAME");
//	bpmEY->Draw("SAME");

//	bpmTargetX->Draw();
//	bpmAX->Draw("PLC PMC");
//	bpmAY->Draw("SAME PLC PMC");
//	bpmEX->Draw("SAME PLC PMC");
//	bpmEY->Draw("SAME PLC PMC");
	beamPosition->Draw();
	 gPad->BuildLegend();




}


