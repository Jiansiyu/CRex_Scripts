/*
 * OpticsEnergySpectrum.cpp
 *
 *  Created on: Nov 13, 2019
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
std::vector <TString> cutNames;

inline Bool_t IsFileExist (const std::string& name) {
	  struct stat buffer;
	  return (stat (name.c_str(), &buffer) == 0);
}
UInt_t OpticsEnergySpectrum(UInt_t runID,TString folder="/home/newdriver/Storage/Research/PRex_Experiment/prex_analyzer/optReplay/Result", double_t groundp=0.9476, double_t firstp=0.94325){
	TChain *chain=new TChain("T");
	TString rootDir(folder.Data());

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

	TCanvas *a=new TCanvas("cut","cut",600,600);
	a->Draw();
	a->cd(0);
	TH2F *HHistThPh=new TH2F("th vs ph","th vs ph",1000,-0.027,0.02,1000,-0.043,0.043);
	chain->Project(HHistThPh->GetName(),Form("R.gold.th:R.gold.ph"));
	HHistThPh->Draw("zcol");
	TCutG *cutg=new TCutG("gcut");
	cutg=(TCutG*)(TVirtualPad::Pad()->WaitPrimitive("CUTG", "CutG")); // making cut, store to
	cutg->SetName(Form("fcut")); //
	cutg->SetVarX("R.gold.ph");
	cutg->SetVarY("R.gold.th");

	cutg->SetLineColor(kMagenta);
    cutg->SetLineWidth(2);
    cutg->Draw("PL");
    a->Update();


    TCanvas *canvasp=new TCanvas("a","q",1000,1000);
    canvasp->cd();
    // plot the momentum spectrum and fit the
    TH1F *TargetMomentumH=new TH1F("gold.p","gold.p", 500,0.93,0.96);
    chain->Project(TargetMomentumH->GetName(),"R.gold.p",Form("%s",cutg->GetName()));


    TargetMomentumH->Draw();

    // fit the ground states momentum peak
    double_t fgroundGausPar[3];
    TF1 *fgroudGaus=new TF1("groudstatesgaus","gaus",groundp-0.0004,groundp+0.0004);
    TargetMomentumH->Fit("groudstatesgaus","R","ep",fgroudGaus->GetXmin(),fgroudGaus->GetXmax());
    fgroudGaus->GetParameters(fgroundGausPar);

    //add the single crystal ball function for the ground states
    TF1 *fgroundCrystalball=new TF1("fgroundCrystal","crystalball",groundp-0.0035,fgroudGaus->GetXmax());
    fgroundCrystalball->SetParameters(fgroundGausPar[0],fgroundGausPar[1],fgroundGausPar[2],1.64,1.1615);
    TargetMomentumH->Fit("fgroundCrystal","R","same",fgroundCrystalball->GetXmin(),fgroundCrystalball->GetXmax());

    double_t fgroundCrystalballPar[5];
    fgroundCrystalball->GetParameters(fgroundCrystalballPar);
    fgroundCrystalball->Draw("same");
    // fit the first exited states
    // fit with the gaus function
    double_t ffirstGuasPar[3];
    TF1 *ffirstGuas=new TF1 ("firststatesgaus","gaus",firstp-0.00025,firstp+0.0005);
    TargetMomentumH->Fit("firststatesgaus","R","ep",ffirstGuas->GetXmin(),ffirstGuas->GetXmax());
    ffirstGuas->GetParameters(ffirstGuasPar);


    // fit with crystall functions
    double_t ffirstCrystalPar[5];
    TF1 *ffirstCrystal=new TF1("ffirstCrystal","crystalball",firstp-0.0025,ffirstGuas->GetXmax());
    ffirstCrystal->SetParameters(ffirstGuasPar[0],ffirstGuasPar[1],ffirstGuasPar[2],1.64,1.1615);
    TargetMomentumH->Fit("ffirstCrystal","R","ep",ffirstCrystal->GetXmin(),ffirstCrystal->GetXmax());
    ffirstCrystal->GetParameters(ffirstCrystalPar);

    ffirstCrystal->Draw("same");

    TPaveText *pt = new TPaveText(0.1,0.8,0.2,0.9,"NDC");
    pt->AddText(Form("%f (%2.2f\%)",1000.0*(fgroundCrystalballPar[1]-ffirstCrystalPar[1]),100.0*abs(abs(fgroundCrystalballPar[1]-ffirstCrystalPar[1])-0.00443891)/0.00443891));
    pt->Draw("same");
    		// fit the two crystall ball functions


//    double_t fCrystalMomentumPar[10];
//    TF1 *fCrystalMomentum=new TF1("fCrystalMomentum","crystalball(0)+crystalball(5)",ffirstCrystal->GetXmin(),fgroundCrystalball->GetXmax());
//    std::copy(fgroundCrystalballPar,fgroundCrystalballPar+5,fCrystalMomentumPar);
//    std::copy(ffirstCrystalPar,ffirstCrystalPar+5,fCrystalMomentumPar+5);
//    fCrystalMomentum->SetParameters(fCrystalMomentumPar);
//    TargetMomentumH->Fit("fCrystalMomentum","","",fCrystalMomentum->GetXmin(),fCrystalMomentum->GetXmax());
//

    gPad->SetLogy(1);

    TCanvas *canvasdp=new TCanvas("dp","dp",1000,1000);
    canvasdp->cd();
    TH1F *TargetDpH=new TH1F("gold.dp","gold.dp", 1000,-0.1,0.1);
    chain->Project(TargetDpH->GetName(),"R.gold.dp",Form("%s",cutg->GetName()));
    TargetDpH->Draw("same");


	return 1;
}


