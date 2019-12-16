/*
 * OpticsGraphicCut.C
 *
 *  Created on: Dec 11, 2019
 *      Author: newdriver
 */



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
UInt_t OpticsGraphicCut(UInt_t runID,TString folder="/home/newdriver/Storage/Research/CRex_Experiment/optReplay/Result", double_t groundp=0.9476, double_t firstp=0.94325){

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


	TCanvas *a=new TCanvas("cut","cut",600,600);
	a->Draw();
	a->cd(0);
	TH2F *HHistThPh=new TH2F("th vs ph","th vs ph",1000,-0.027,0.02,1000,-0.043,0.043);
	chain->Project(HHistThPh->GetName(),Form("%s.gold.th:%s.gold.ph",HRS.Data(),HRS.Data()));
	HHistThPh->Draw("zcol");
	TCutG *cutg=new TCutG("gcut");
	cutg=(TCutG*)(TVirtualPad::Pad()->WaitPrimitive("CUTG", "CutG")); // making cut, store to
	cutg->SetName(Form("fcut")); //
	cutg->SetVarX(Form("%s.gold.ph",HRS.Data()));
	cutg->SetVarY(Form("%s.gold.th",HRS.Data()));

	cutg->SetLineColor(kMagenta);
    cutg->SetLineWidth(2);
    cutg->Draw("PL");
    a->Update();


    // plot all the plot that plot that with the cut
    TCanvas *canvasp=new TCanvas("a","q",1000,1000);
    canvasp->cd();
    // plot the momentum spectrum and fit the
    TH1F *TargetMomentumH=new TH1F("gold.p","gold.p", 500,0.93,0.96);
    chain->Project(TargetMomentumH->GetName(),Form("%s.gold.p",HRS.Data()),Form("%s",cutg->GetName()));


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

    gPad->SetLogy(1);

    TCanvas *canvasdp=new TCanvas("dp","dp",1000,1000);
    canvasdp->cd();
    TH1F *TargetDpH=new TH1F("gold.dp","gold.dp", 1000,-0.1,0.1);
    chain->Project(TargetDpH->GetName(),"R.gold.dp",Form("%s",cutg->GetName()));
    TargetDpH->Draw("same");


	return 1;
}

UInt_t OpticsGraphicCutGeneral(UInt_t runID,TString plotString,TString folder="/home/newdriver/Storage/Research/CRex_Experiment/optReplay/Result"){

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


	TCanvas *a=new TCanvas("cut","cut",600,600);
	a->Draw();
	a->cd(0);
	TH2F *HHistThPh=new TH2F("th vs ph","th vs ph",1000,-0.027,0.02,1000,-0.043,0.043);
	chain->Project(HHistThPh->GetName(),Form("%s.gold.th:%s.gold.ph",HRS.Data(),HRS.Data()));
	HHistThPh->Draw("zcol");
	TCutG *cutg=new TCutG("gcut");
	cutg=(TCutG*)(TVirtualPad::Pad()->WaitPrimitive("CUTG", "CutG")); // making cut, store to
	cutg->SetName(Form("fcut")); //
	cutg->SetVarX(Form("%s.gold.ph",HRS.Data()));
	cutg->SetVarY(Form("%s.gold.th",HRS.Data()));

	cutg->SetLineColor(kMagenta);
    cutg->SetLineWidth(2);
    cutg->Draw("PL");
    a->Update();
//    TCanvas	*canvasdpx=new TCanvas("dp-x","dp-x",600,600);
//    canvasdpx->Divide(1,2);
//
//    TH1F *dp=new TH1F("dp","dp",1000,-0.04,0.02);
    chain->Draw(plotString.Data(),Form("%s",cutg->GetName()),"zcol");

	return 1;
}


UInt_t OpticsGraphicCutDp(UInt_t runID,TString folder="/home/newdriver/Storage/Research/CRex_Experiment/optReplay/Result/AfterCorrection"){

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


	TCanvas *a=new TCanvas("cut","cut",600,600);
	a->Draw();
	a->cd(0);
	TH2F *HHistThPh=new TH2F("th vs ph","th vs ph",1000,-0.027,0.02,1000,-0.043,0.043);
	chain->Project(HHistThPh->GetName(),Form("%s.gold.th:%s.gold.ph",HRS.Data(),HRS.Data()));
	HHistThPh->Draw("zcol");
	TCutG *cutg=new TCutG("gcut");
	cutg=(TCutG*)(TVirtualPad::Pad()->WaitPrimitive("CUTG", "CutG")); // making cut, store to
	cutg->SetName(Form("fcut")); //
	cutg->SetVarX(Form("%s.gold.ph",HRS.Data()));
	cutg->SetVarY(Form("%s.gold.th",HRS.Data()));

	cutg->SetLineColor(kMagenta);
    cutg->SetLineWidth(2);
    cutg->Draw("PL");
    a->Update();
    TCanvas	*canvasdpx=new TCanvas("dp-x","dp-x",600,600);
    canvasdpx->Divide(2,2);



    TH1F *dp_temp=new TH1F("temp","temp",1000,-0.5,0.5);
    chain->Project(dp_temp->GetName(),Form("%s.gold.dp",HRS.Data()),Form("%s",cutg->GetName()));
    auto dp_peak=dp_temp->GetXaxis()->GetBinCenter(dp_temp->GetMaximumBin());

    canvasdpx->cd(1);
    TH1F *dp=new TH1F("gold.dp","gold.dp",1000,dp_peak-0.005*2,dp_peak+0.005);
    chain->Project(dp->GetName(),Form("%s.gold.dp",HRS.Data()),Form("%s",cutg->GetName()));
    dp->GetXaxis()->SetTitle("gold.dp");
    dp->GetYaxis()->SetTitle("#");
    dp->Draw();

    double fgrounddpPar[3];
    double fg1stdpPar[3];
    TF1 *fgroudDpGaus=new TF1("groudstatesDpgaus","gaus",dp->GetXaxis()->GetBinCenter(dp->GetMaximumBin())-0.0005,dp->GetXaxis()->GetBinCenter(dp->GetMaximumBin())+0.0005);
    dp->Fit("groudstatesDpgaus","R","ep",fgroudDpGaus->GetXmin(),fgroudDpGaus->GetXmax());
    fgroudDpGaus->Draw("same");

    TF1 *f1stDpGaus=new TF1("1ststatesDpgaus","gaus",dp->GetXaxis()->GetBinCenter(dp->GetMaximumBin())-0.0025-0.001,dp->GetXaxis()->GetBinCenter(dp->GetMaximumBin())-0.0025+0.001);
    dp->Fit("1ststatesDpgaus","R","ep",f1stDpGaus->GetXmin(),f1stDpGaus->GetXmax());
    f1stDpGaus->Draw("same");

    TPaveText *ptdp = new TPaveText(0.1,0.8,0.2,0.9,"NDC");
    ptdp->AddText(Form("%f -%f",fgroudDpGaus->GetParameter(1),f1stDpGaus->GetParError(1)));
     ptdp->Draw("same");




    canvasdpx->cd(2);
    TH1F *focal_x=new TH1F ("x-focal","x-focal",1000,-0.3,0.3);
    chain->Project(focal_x->GetName(),Form("%s.tr.r_x",HRS.Data()),Form("%s",cutg->GetName()));
    focal_x->GetXaxis()->SetRangeUser(focal_x->GetXaxis()->GetBinCenter(focal_x->GetMaximumBin())-0.06,focal_x->GetXaxis()->GetBinCenter(focal_x->GetMaximumBin())+0.02);
    focal_x->GetXaxis()->SetTitle("tr.r_x");
    focal_x->GetYaxis()->SetTitle("#");
    focal_x->Draw();

    canvasdpx->cd(3);
    TH2F *focalx_dp=new TH2F("dp vs focal_x","dp vs focal_x",1000,-0.6,0.3,1000,-0.05,0.015);
//    TH2F *focalx_dp=new TH2F("dp vs focal_x","dp vs focal_x",1000,focal_x->GetXaxis()->GetXmin(),focal_x->GetXaxis()->GetXmax(),1000,dp->GetXaxis()->GetXmin(),dp->GetXaxis()->GetXmax());
    chain->Project(focalx_dp->GetName(),Form("%s.gold.dp:%s.tr.r_x",HRS.Data(),HRS.Data()),Form("%s",cutg->GetName()));
    focalx_dp->GetXaxis()->SetTitle("tr.r_x");
    focalx_dp->GetYaxis()->SetTitle("gold.dp");
    focalx_dp->Draw();
    focalx_dp->Fit("pol1");
    TPaveText *pt3 = new TPaveText(0.1,0.8,0.2,0.9,"NDC");
    pt3->AddText(Form("dp=%f * focal_x + %f ",focalx_dp->GetFunction("pol1")->GetParameter(1),focalx_dp->GetFunction("pol1")->GetParameter(0)));
    pt3->Draw("same");



    canvasdpx->cd(4);
    TH1F *momentum=new TH1F("C-12 gold.p","C-12 gold.p",1500,2.1,2.2);
    chain->Project(momentum->GetName(),Form("%s.gold.p",HRS.Data()),Form("%s",cutg->GetName()));
    // get the maximum bin, this should be the first excited states
    auto CGroundp=momentum->GetXaxis()->GetBinCenter(momentum->GetMaximumBin());
    auto C1stp=CGroundp-0.0044;
    momentum->GetXaxis()->SetRangeUser(CGroundp-0.0044*3,CGroundp+0.0044*2);
    momentum->GetXaxis()->SetTitle("gold.p");
    momentum->GetYaxis()->SetTitle("#");

	momentum->Draw();
    double_t fgroudGausPar[3];
    double_t ffirstGuasPar[3];
    TF1 *fgroudGaus=new TF1("groudstatesgaus","gaus",CGroundp-0.001,CGroundp+0.001);
    momentum->Fit("groudstatesgaus","R","ep",fgroudGaus->GetXmin(),fgroudGaus->GetXmax());
    fgroudGaus->Draw("same");
    fgroudGaus->GetParameters(fgroudGausPar);
    TF1 *ffirstGuas=new TF1 ("firststatesgaus","gaus",C1stp-0.0005,C1stp+0.0005);
    momentum->Fit("firststatesgaus","R","ep",ffirstGuas->GetXmin(),ffirstGuas->GetXmax());
    ffirstGuas->Draw("same");
    ffirstGuas->GetParameters(ffirstGuasPar);

    TPaveText *pt = new TPaveText(0.1,0.8,0.2,0.9,"NDC");
    pt->AddText(Form("%1.3f MeV (%2.2f\%)",1000.0*(fgroudGausPar[1]-ffirstGuasPar[1]),100.0*abs(abs(fgroudGausPar[1]-ffirstGuasPar[1])-0.00443891)/0.00443891));
    pt->Draw("same");


	return 1;
}

