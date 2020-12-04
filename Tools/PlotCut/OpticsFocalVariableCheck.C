/*
 * OpticsFocalVariableCheck.C
 *
 *  Created on: May 19, 2020
 *      Author: newdriver
 */




/*
 * OpticsPointingDpCheck.C
 *
 *  Created on: May 7, 2020
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
#include <TLegend.h>
#include <TApplication.h>
#include <TArrow.h>
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
TString generalcutR="R.tr.n==1 && R.vdc.u1.nclust==1&& R.vdc.v1.nclust==1 && R.vdc.u2.nclust==1 && R.vdc.v2.nclust==1 && R.gold.p > 2.14 && R.gold.p < 2.2  ";
TString generalcutL="L.tr.n==1 && L.vdc.u1.nclust==1&& L.vdc.v1.nclust==1 && L.vdc.u2.nclust==1 && L.vdc.v2.nclust==1  && L.gold.p > 2.14 && L.gold.p < 2.2";

inline Bool_t IsFileExist (const std::string& name) {
	  struct stat buffer;
	  return (stat (name.c_str(), &buffer) == 0);
}


// does it needed to add another function to predict the position of each peak
// add an global fit function used for the fit
TF1 * SpectroCrystalFit_C12(TH1F *momentumSpectro){
	TF1 *globalFit;

	TH1F *FitInitialh=(TH1F *)momentumSpectro->Clone("initialh");
	// fit the highest peak, this should be the ground states peak
	auto CGroundp=FitInitialh->GetXaxis()->GetBinCenter(FitInitialh->GetMaximumBin());
	auto C1stp=CGroundp-0.00443891;
	return globalFit;
}

TF1 *SpectroCrystalFitDp_C12(TH1F*momentumSpectro){


	auto CGroundDp=momentumSpectro->GetXaxis()->GetBinCenter(momentumSpectro->GetMaximumBin());
	//start the fit and get the mean ans sigma
	momentumSpectro->Fit("gaus","RQ0","ep",CGroundDp-0.0003,CGroundDp+0.0003);

	double_t fgroundCrystalballPar[5];

	TF1 *fgroundCrystalball = new TF1("fgroundCrystal", "crystalball",
			momentumSpectro->GetFunction("gaus")->GetParameter(1)
					- 5 * momentumSpectro->GetFunction("gaus")->GetParameter(2),
			momentumSpectro->GetFunction("gaus")->GetParameter(1)
					+ 5 * momentumSpectro->GetFunction("gaus")->GetParameter(2));
	fgroundCrystalball->SetParameters(
			momentumSpectro->GetFunction("gaus")->GetParameter(0),
			momentumSpectro->GetFunction("gaus")->GetParameter(1),
			momentumSpectro->GetFunction("gaus")->GetParameter(2), 1.64, 1.1615);

	momentumSpectro->Fit("fgroundCrystal","RQ0","ep",fgroundCrystalball->GetXmin(),fgroundCrystalball->GetXmax());
	fgroundCrystalball->GetParameters(fgroundCrystalballPar);


	TH1F *test=(TH1F *)momentumSpectro->Clone("fitTest");
	test->GetXaxis()->SetRangeUser(momentumSpectro->GetXaxis()->GetXmin(),fgroundCrystalballPar[1]-5*fgroundCrystalballPar[2]);

	double_t ffirstGuasPar[3];
	auto C1stp=test->GetXaxis()->GetBinCenter(test->GetMaximumBin());
	test->Delete();
	TF1 *ffirstGuas=new TF1 ("firststatesgaus","gaus",C1stp-3*fgroundCrystalballPar[2],C1stp+3*fgroundCrystalballPar[2]);
	momentumSpectro->Fit("firststatesgaus","R0Q","ep",ffirstGuas->GetXmin(),ffirstGuas->GetXmax());
	ffirstGuas->GetParameters(ffirstGuasPar);

	double_t ffirstCrystalPar[5];
	TF1 *ffirstCrystal=new TF1("ffirstCrystal","crystalball",ffirstGuasPar[1]-0.0025,ffirstGuas->GetXmax());
	ffirstCrystal->SetParameters(ffirstGuasPar[0],ffirstGuasPar[1],ffirstGuasPar[2],1.64,1.1615);
	momentumSpectro->Fit("ffirstCrystal","R","ep",ffirstCrystal->GetXmin(),ffirstCrystal->GetXmax());
	ffirstCrystal->GetParameters(ffirstCrystalPar);

	double_t fCrystalMomentumPar[10];
	TF1 *fCrystalMomentum=new TF1("fCrystalMomentum","crystalball(0)+crystalball(5)",ffirstCrystal->GetXmin(),fgroundCrystalball->GetXmax());
	std::copy(fgroundCrystalballPar,fgroundCrystalballPar+5,fCrystalMomentumPar);
	std::copy(ffirstCrystalPar,ffirstCrystalPar+5,fCrystalMomentumPar+5);
	fCrystalMomentum->SetParameters(fCrystalMomentumPar);
	momentumSpectro->Fit("fCrystalMomentum","","",fCrystalMomentum->GetXmin(),fCrystalMomentum->GetXmax());


	return fCrystalMomentum;
}

TF1 *SpectroCrystalFitDp_H2O(TH1F*momentumSpectro){
	auto CGroundDp=momentumSpectro->GetXaxis()->GetBinCenter(momentumSpectro->GetMaximumBin());
	momentumSpectro->Fit("gaus","RQ0","ep",CGroundDp-0.0003,CGroundDp+0.0003);
	double_t fgroundCrystalballPar[5];

	TF1 *fgroundCrystalball =
			new TF1("fgroundCrystal", "crystalball",
					momentumSpectro->GetFunction("gaus")->GetParameter(1)
							- 5
									* momentumSpectro->GetFunction("gaus")->GetParameter(
											2),
					momentumSpectro->GetFunction("gaus")->GetParameter(1)
							+ 5
									* momentumSpectro->GetFunction("gaus")->GetParameter(
											2));
	fgroundCrystalball->SetParameters(
			momentumSpectro->GetFunction("gaus")->GetParameter(0),
			momentumSpectro->GetFunction("gaus")->GetParameter(1),
			momentumSpectro->GetFunction("gaus")->GetParameter(2), 1.64,
			1.1615);

	momentumSpectro->Fit("fgroundCrystal", "RQ0", "ep",
			fgroundCrystalball->GetXmin(), fgroundCrystalball->GetXmax());
	fgroundCrystalball->GetParameters(fgroundCrystalballPar);

	double_t ffirstGuasPar[3];
	auto C1stp=fgroundCrystalballPar[1]-16/21800;//momentumSpectro->GetXaxis()->GetBinCenter(momentumSpectro->GetMaximumBin(momentumSpectro->GetXaxis()->GetXmin()),fgroundCrystalballPar[1]-6*fgroundCrystalballPar[2]);
	TF1 *ffirstGuas=new TF1 ("firststatesgaus","gaus",C1stp-3*fgroundCrystalballPar[2],C1stp+3*fgroundCrystalballPar[2]);
	momentumSpectro->Fit("firststatesgaus","R0Q","ep",ffirstGuas->GetXmin(),ffirstGuas->GetXmax());
	ffirstGuas->GetParameters(ffirstGuasPar);

	double_t ffirstCrystalPar[5];
	TF1 *ffirstCrystal=new TF1("ffirstCrystal","crystalball",ffirstGuasPar[1]-0.0025,ffirstGuas->GetXmax());
	ffirstCrystal->SetParameters(ffirstGuasPar[0],ffirstGuasPar[1],ffirstGuasPar[2],1.64,1.1615);
	momentumSpectro->Fit("ffirstCrystal","R","ep",ffirstCrystal->GetXmin(),ffirstCrystal->GetXmax());
	ffirstCrystal->GetParameters(ffirstCrystalPar);

	double_t fCrystalMomentumPar[10];
	TF1 *fCrystalMomentum=new TF1("fCrystalMomentum","crystalball(0)+crystalball(5)",ffirstCrystal->GetXmin(),fgroundCrystalball->GetXmax());
	std::copy(fgroundCrystalballPar,fgroundCrystalballPar+5,fCrystalMomentumPar);
	std::copy(ffirstCrystalPar,ffirstCrystalPar+5,fCrystalMomentumPar+5);
	fCrystalMomentum->SetParameters(fCrystalMomentumPar);
	momentumSpectro->Fit("fCrystalMomentum","","",fCrystalMomentum->GetXmin(),fCrystalMomentum->GetXmax());

	return fCrystalMomentum;

}


TChain *LoadrootFile(UInt_t runID,TString folder="/home/newdriver/Storage/Research/CRex_Experiment/RasterReplay/Replay/Result"){
	TChain *chain=new TChain("T");
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

Int_t OpticsPointingDpCheck(UInt_t runID,TString folder="/home/newdriver/Storage/Research/CRex_Experiment/RasterReplay/Replay/Result/") {
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
	gStyle->SetTimeOffset(0);
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
	SieveRecCanvas->cd(2)->Divide(5, 1);
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
	TCutG *cutg = new TCutG(Form("hcut_R_%ld",random()),
			gc1->GetN(), gc1->GetX(), gc1->GetY());
	cutg->SetLineWidth(2);
	cutg->SetLineColor(kRed);
	cutg->SetVarX(Form("%s.gold.ph", HRS.Data()));
	cutg->SetVarY(Form("%s.gold.th", HRS.Data()));
	cutg->Draw("same");

	// plot the cut on the canvas
	SieveRecCanvas->Update();

	//Load ALl the root files
	std::map<int,TChain *> chainArray;

	if (HRS == 'L') {
		std::cout<<"Working on LHRS"<<std::endl;
		chainArray[-2] = LoadrootFile(2566);
		chainArray[-1] = LoadrootFile(2565);
		chainArray[0] = LoadrootFile(2550);
		chainArray[1] = LoadrootFile(2556);
//		chainArray[999] = LoadrootFile(2555);
	} else {
		std::cout<<"Working on RHRS"<<std::endl;
		chainArray[-2] = LoadrootFile(21642);
		chainArray[-1] = LoadrootFile(21641);
		chainArray[0] = LoadrootFile(21626);
		chainArray[1] = LoadrootFile(21632);
//		chainArray[999] = LoadrootFile(21789);
	}

	std::map<int,TH2F *>SieveThetaPhiList;

	int canvas_counter_temp=1;
	for (auto item = chainArray.begin(); item!= chainArray.end();item++){
		if((item->first)<10){
			SieveThetaPhiList[item->first]=new TH2F(Form("th_vs_ph_Dp%1d",item->first),Form("th_vs_ph_Dp%1d",item->first),1000,-0.03,0.03,1000,-0.045,0.045);
			item->second->Project(SieveThetaPhiList[item->first]->GetName(),Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()),
				Form("%s", generalcut.Data()));

			SieveRecCanvas->cd(2)->cd(canvas_counter_temp);
			SieveThetaPhiList[item->first]->Draw("zcol");
			cutg->Draw("same");
		}else{
			SieveThetaPhiList[item->first]=new TH2F(Form("th_vs_ph_H20%1d",item->first),Form("th_vs_ph_H20%1d",item->first),1500,-0.03,0.03,1000,-0.045,0.045);
			item->second->Project(SieveThetaPhiList[item->first]->GetName(),Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()),
							Form("%s", generalcut.Data()));
			SieveRecCanvas->cd(2)->cd(canvas_counter_temp);
			SieveThetaPhiList[item->first]->Draw("zcol");
			cutg->Draw("same");
		}
		canvas_counter_temp++;
	}

	SieveRecCanvas->Update();

	// put the Dp in same canvas
	// create an new canvas
	std::map<int,TH1F *>OptDpArrayH;
	int maximumPeakHight=0;
	for (auto item = chainArray.begin(); item!= chainArray.end();item++){
		if((item->first)<10){
			OptDpArrayH[item->first]=new TH1F(Form("Dp_hist%d",item->first),Form("Dp_hist%d",item->first),1000,-0.03,0.02);
			OptDpArrayH[item->first]->GetYaxis()->SetRange(0,10000);
			item->second->Project(OptDpArrayH[item->first]->GetName(),Form("%s.gold.dp",HRS.Data()),
			Form("%s && %s", generalcut.Data(),cutg->GetName()));
		}else{
			OptDpArrayH[item->first]=new TH1F(Form("Dp_hist:H2O%d",item->first),Form("Dp_hist:H_{2}O%d",item->first),1000,-0.1,0.1);

			item->second->Project(OptDpArrayH[item->first]->GetName(),Form("%s.gold.dp",HRS.Data()),
		    Form("%s && %s", generalcut.Data(),cutg->GetName()));
		}
		if(maximumPeakHight<OptDpArrayH[item->first]->GetMaximumBin()){
			maximumPeakHight=OptDpArrayH[item->first]->GetMaximumBin();
		}
	}
	TCanvas *DpCanvas = (TCanvas*) gROOT->GetListOfCanvases()->FindObject("DpCanvas");
	if (DpCanvas) {
		DpCanvas->Clear();

	} else {
		DpCanvas = new TCanvas("DpCanvas", "DpCanvas", 600, 600);
	}
	DpCanvas->Draw();
	DpCanvas->cd();
	auto legend = new TLegend(0.1,0.7,0.48,0.9);
	legend->SetHeader("Dp Scan","C"); //option "C" allows to center the header
	std::map<int,TF1 *> fitFunctionsList;

	std::map<int, double *>FitPars;

	//0.0117515,0.0015453,-0.00859458,-0.0180413
	std::map<int, double> DpTheoreticalMap;
	if (HRS == "L") {
		DpTheoreticalMap[-2] = 0.0154806;
		DpTheoreticalMap[-1] = 0.00542595;
		DpTheoreticalMap[0] = -0.00496891;
		DpTheoreticalMap[1] = -0.010256;
	} else {//0.0141027,0.00387559,-0.00628469,-0.0157513,0.012034,0.00182792,-0.008312,-0.0177587
		DpTheoreticalMap[-2] = 0.0141027;
		DpTheoreticalMap[-1] = 0.00387559;
		DpTheoreticalMap[0] = -0.00628469;
		DpTheoreticalMap[1] = -0.0157513;
	}

	std::map<int, double >CentralPArray;
	for (auto item = OptDpArrayH.begin(); item!= OptDpArrayH.end();item++){

		if((item->first)<10){
			item->second->SetLineColor(8+item->first);
			legend->AddEntry((item->second),Form("C_{12} Dp:%2d%% scan",item->first));
		}else{
			item->second->SetLineColor(kRed);
			legend->AddEntry((item->second),Form("H_{2}O Dp"));
		}

		item->second->GetYaxis()->SetRangeUser(0,10000);
		item->second->SetLineWidth(2);

		if(item==OptDpArrayH.begin()){
			item->second->Draw();
		}else{
			item->second->Draw("same");
		}

		fitFunctionsList[item->first]=SpectroCrystalFitDp_C12(item->second);
		fitFunctionsList[item->first]->SetLineColor(42);
		fitFunctionsList[item->first]->Draw("same");
		FitPars[item->first]=new double[10];
		fitFunctionsList[item->first]->GetParameters(FitPars[item->first]);


		if(HRS=="L"){
			TH1F *HallProbHH = new TH1F("HallLProb", "HallLProb", 1000, -1, 0);
			chainArray[item->first]->Project(HallProbHH->GetName(),
					"HacL_D_LS450_FLD_DATA", generalcut.Data());
			if (HallProbHH->GetEntries() != 0) {
				double CentralP = std::abs(
						(HallProbHH->GetMean()) * 0.95282 / 0.33930);
				std::cout << "CentralMomentum is ::" << (CentralP) << std::endl;
				CentralPArray[item->first]=CentralP;
			}else {
				std::cout << "\033[1;33m [Warning]\033[0m Missing HallLProb:"
						<< std::endl;
			}
			HallProbHH->Delete();
		}else{
			//HacR_D1_NMR_SIG
			double CentralP;
			TH1F *HallR_NMR = new TH1F("HallR_NMR", "HallR_NMR", 1000, 0.7, 0.9);
			chainArray[item->first]->Project(HallR_NMR->GetName(), "HacR_D1_NMR_SIG",
					generalcut.Data());
			if (HallR_NMR->GetEntries()) {
				double Mag = HallR_NMR->GetMean();
				CentralP = 2.702 * (Mag) - 1.6e-03 * (Mag) * (Mag) * (Mag);
				CentralPArray[item->first]=CentralP;

			} else {
				std::cout << "\033[1;33m [Warning]\033[0m Missing HallR_NMR:"
						<< std::endl;
			}
			HallR_NMR->Delete();
		}
	}

	// start searching for the focal plane variables
	// general cut + central sieve hole cut + ground states cut
	double groundPeakSigma=3.0;

	std::map<int,TCut> dpCutArray;   // list of the cut
	for (auto fitpar_iter= FitPars.begin(); fitpar_iter!=FitPars.end();fitpar_iter++){
		auto groundPeakMean=fitpar_iter->second[1];
		auto groundPeakRMS =fitpar_iter->second[2];
		dpCutArray[fitpar_iter->first]=TCut(Form("abs(%s.gold.dp-%f)<%f*%f",HRS.Data(),groundPeakMean,groundPeakSigma,groundPeakRMS));
	}



	std::map<int, TH1F *> sieveFocalXArray;
	std::map<int, TH1F *> sieveFocalYArray;
	std::map<int, TH1F *> sieveFocalThetaArray;
	std::map<int, TH1F *> sieveFocalPhiArray;


	//project the focal plane variables
	for (auto item = chainArray.begin(); item!= chainArray.end();item++){
		sieveFocalXArray[item->first]=new TH1F(Form("Dp%d_%s.tr.x",item->first,HRS.Data()),Form("Dp%d_%s.tr.x",item->first,HRS.Data()),2000,-0.4,0.4);
		sieveFocalYArray[item->first]=new TH1F(Form("Dp%d_%s.tr.y",item->first,HRS.Data()),Form("Dp%d_%s.tr.y",item->first,HRS.Data()),1000,-0.01,0.01);
		sieveFocalThetaArray[item->first]=new TH1F(Form("Dp%d_%s.tr.th",item->first,HRS.Data()),Form("Dp%d_%s.tr.th",item->first,HRS.Data()),2000,-0.06,0.06);
		sieveFocalPhiArray[item->first]=new TH1F(Form("Dp%d_%s.tr.ph",item->first,HRS.Data()),Form("Dp%d_%s.tr.ph",item->first,HRS.Data()),1000,-0.01,0.01);

		TCut focalCut=TCut(generalcut.Data())+dpCutArray[item->first]+TCut(cutg->GetName());
		item->second->Project(sieveFocalXArray[item->first]->GetName(),Form("%s.tr.x",HRS.Data()),focalCut);
		item->second->Project(sieveFocalYArray[item->first]->GetName(),Form("%s.tr.y",HRS.Data()),focalCut);
		item->second->Project(sieveFocalThetaArray[item->first]->GetName(),Form("%s.tr.th",HRS.Data()),focalCut);
		item->second->Project(sieveFocalPhiArray[item->first]->GetName(),Form("%s.tr.ph",HRS.Data()),focalCut);

		//focalLegend[item->first]=new TLegend(0.1,0.7,0.48,0.9);

	}


	//focal plan variables
	std::map<int, double *> SieveFocalXGuasFitPar;
	std::map<int, double *> SieveFocalYGuasFitPar;

	TCanvas *focalPlaneDiag=new TCanvas("Focal Plane Variables","Folca Plane Variables",1000,1000);
	focalPlaneDiag->Draw();
	focalPlaneDiag->Divide(2,2);

	focalPlaneDiag->cd(1);
	TLegend * focalLegend_x;
	focalLegend_x=new TLegend(0.1,0.7,0.48,0.9);
	for(auto item=sieveFocalXArray.begin();item!=sieveFocalXArray.end(); item++){
		item->second->SetLineColor(8+item->first);
		if(item->first < 100){
			focalLegend_x->AddEntry((item->second),Form("C_{12} Dp%d%%",item->first));
		}else{
			focalLegend_x->AddEntry((item->second),Form("H_{2}O %d",item->first));
		}
		focalLegend_x->SetLineWidth(2);
		if (item == sieveFocalXArray.begin()){
			item->second->Draw();
		}else{
			item->second->Draw("same");
		}

		// get peak value and put in the canvas
		// fit the peak with gaus
		double GuasFitPars[3];
		SieveFocalXGuasFitPar[item->first]=new double[3];
		auto maximumPeak=item->second->GetXaxis()->GetBinCenter(item->second->GetMaximumBin());
		item->second->Fit("gaus","R0Q","ep",maximumPeak-0.005,maximumPeak+0.005);
		item->second->GetFunction("gaus")->GetParameters(GuasFitPars);
		item->second->GetFunction("gaus")->GetParameters(SieveFocalXGuasFitPar[item->first]);

		TLatex *peakInfor=new TLatex(GuasFitPars[1],GuasFitPars[0],Form("%f",GuasFitPars[1]));
		peakInfor->SetLineWidth(2);
		peakInfor->SetTextSize(0.04);
		peakInfor->Draw("same");
	}
	focalLegend_x->Draw("same");

	focalPlaneDiag->cd(2);
	TLegend * focalLegend_y;
	focalLegend_y=new TLegend(0.1,0.7,0.48,0.9);
	for(auto item=sieveFocalYArray.begin();item!=sieveFocalYArray.end(); item++){
		item->second->SetLineColor(8+item->first);
		item->second->SetLineWidth(2);
		if(item->first < 100){
			focalLegend_y->AddEntry((item->second),Form("C_{12} Dp%d%%",item->first));
		}else{
			focalLegend_y->AddEntry((item->second),Form("H_{2}O %d",item->first));
		}
		focalLegend_y->SetLineWidth(2);
		if (item == sieveFocalYArray.begin()){
			item->second->Draw();
		}else{
			item->second->Draw("same");
		}
		// get peak value and put in the canvas
		// fit the peak with gaus
		SieveFocalYGuasFitPar[item->first]=new double[3];
		double GuasFitPars[3];
		auto maximumPeak=item->second->GetXaxis()->GetBinCenter(item->second->GetMaximumBin());
		item->second->Fit("gaus","R0Q","ep",maximumPeak-0.005,maximumPeak+0.005);
		item->second->GetFunction("gaus")->GetParameters(GuasFitPars);
		item->second->GetFunction("gaus")->GetParameters(SieveFocalYGuasFitPar[item->first]);
		TLatex *peakInfor=new TLatex(GuasFitPars[1],GuasFitPars[0],Form("%f",GuasFitPars[1]));
		peakInfor->SetLineWidth(2);
		peakInfor->SetTextSize(0.04);
		peakInfor->Draw("same");
	}
	focalLegend_y->Draw("same");

	focalPlaneDiag->cd(3);
	TLegend * focalLegend_th;
	focalLegend_th=new TLegend(0.1,0.7,0.48,0.9);
	for(auto item=sieveFocalThetaArray.begin();item!=sieveFocalThetaArray.end(); item++){
		item->second->SetLineColor(8+item->first);
		item->second->SetLineWidth(2);
		if(item->first < 100){
			focalLegend_th->AddEntry((item->second),Form("C_{12} Dp%d%%",item->first));
		}else{
			focalLegend_th->AddEntry((item->second),Form("H_{2}O %d",item->first));
		}
		focalLegend_th->SetLineWidth(2);
		if (item == sieveFocalThetaArray.begin()){
			item->second->Draw();
		}else{
			item->second->Draw("same");
		}

		// get peak value and put in the canvas
		// fit the peak with gaus
		double GuasFitPars[3];
		auto maximumPeak=item->second->GetXaxis()->GetBinCenter(item->second->GetMaximumBin());
		item->second->Fit("gaus","R0Q","ep",maximumPeak-0.005,maximumPeak+0.005);
		item->second->GetFunction("gaus")->GetParameters(GuasFitPars);
		TLatex *peakInfor=new TLatex(GuasFitPars[1],GuasFitPars[0],Form("%f",GuasFitPars[1]));
		peakInfor->SetLineWidth(2);
		peakInfor->SetTextSize(0.04);
		peakInfor->Draw("same");
	}
	focalLegend_th->Draw("same");

	focalPlaneDiag->cd(4);
	TLegend * focalLegend_phi;
	focalLegend_phi=new TLegend(0.1,0.7,0.48,0.9);
	for(auto item=sieveFocalPhiArray.begin();item!=sieveFocalPhiArray.end(); item++){
		item->second->SetLineColor(8+item->first);
//		item->second->SetLineWidth(2);
		if(item->first < 100){
			focalLegend_phi->AddEntry((item->second),Form("C_{12} Dp%d%%",item->first));
		}else{
			focalLegend_phi->AddEntry((item->second),Form("H_{2}O %d",item->first));
		}
		focalLegend_phi->SetLineWidth(2);
		if (item == sieveFocalPhiArray.begin()){
			item->second->Draw();
		}else{
			item->second->Draw("same");
		}
		// get peak value and put in the canvas
		// fit the peak with gaus
		double GuasFitPars[3];
		auto maximumPeak=item->second->GetXaxis()->GetBinCenter(item->second->GetMaximumBin());
		item->second->Fit("gaus","R0Q","ep",maximumPeak-0.005,maximumPeak+0.005);
		item->second->GetFunction("gaus")->GetParameters(GuasFitPars);
		TLatex *peakInfor=new TLatex(GuasFitPars[1],GuasFitPars[0],Form("%f",GuasFitPars[1]));
		peakInfor->SetLineWidth(2);
		peakInfor->SetTextSize(0.04);
		peakInfor->Draw("same");
	}
	focalLegend_phi->Draw("same");
	focalPlaneDiag->Update();

	// check the cut, check the sieve patter and check Dp plot
	// check the goodness of cut used in the
	TCanvas *SieveThetaPhiPCanv=new TCanvas("Check Sieve Theta Phi P","Check Sieve Theta Phi P",1000,1000);
	SieveThetaPhiPCanv->Divide(SieveThetaPhiList.size(),3);
	// draw the initial canvas theta and phi
	canvas_counter_temp=1;
	for(auto item = SieveThetaPhiList.begin(); item != SieveThetaPhiList.end();item++){
		SieveThetaPhiPCanv->cd(canvas_counter_temp);
		SieveThetaPhiList[item->first]->Draw("zcol");
		canvas_counter_temp++;
	}

	std::map<int,TH2F *>SieveThetaPhiWithFocalCutList;
	canvas_counter_temp=1;
	for (auto item = chainArray.begin(); item!= chainArray.end();item++){
		TCut focalCut=TCut(generalcut.Data())+dpCutArray[item->first]+TCut(cutg->GetName());
		if((item->first)<10){
			SieveThetaPhiWithFocalCutList[item->first]=new TH2F(Form("th_vs_ph_Dp%1d_focalcut",item->first),Form("th_vs_ph_Dp%1d_focalcu",item->first),1000,-0.03,0.03,1000,-0.045,0.045);
			item->second->Project(SieveThetaPhiWithFocalCutList[item->first]->GetName(),Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()),
				focalCut);

			SieveThetaPhiPCanv->cd(SieveThetaPhiList.size()+canvas_counter_temp);
			SieveThetaPhiWithFocalCutList[item->first]->Draw("zcol");
			cutg->Draw("same");
		}else{
			SieveThetaPhiWithFocalCutList[item->first]=new TH2F(Form("th_vs_ph_H20%1d_focalcu",item->first),Form("th_vs_ph_H20%1d_focalcu",item->first),1500,-0.03,0.03,1000,-0.045,0.045);
			item->second->Project(SieveThetaPhiWithFocalCutList[item->first]->GetName(),Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()),
							focalCut);
			SieveThetaPhiPCanv->cd(SieveThetaPhiList.size()+canvas_counter_temp);
			SieveThetaPhiWithFocalCutList[item->first]->Draw("zcol");
			cutg->Draw("same");
		}
		canvas_counter_temp++;
	}

	// plot the dp informations
	std::map<int,TH1F *>OptDpArrayWithFocalCutH;
	canvas_counter_temp=1;
	for (auto item = chainArray.begin(); item!= chainArray.end();item++){
		TCut focalCut=TCut(generalcut.Data())+dpCutArray[item->first]+TCut(cutg->GetName());
		if((item->first)<10){
			OptDpArrayWithFocalCutH[item->first]=new TH1F(Form("Dp_hist:%d",item->first),Form("Dp_hist:%d",item->first),1000,-0.03,0.02);
			OptDpArrayWithFocalCutH[item->first]->GetYaxis()->SetRange(0,10000);
			item->second->Project(OptDpArrayWithFocalCutH[item->first]->GetName(),Form("%s.gold.dp",HRS.Data()),
			focalCut);
		}else{
			OptDpArrayWithFocalCutH[item->first]=new TH1F(Form("Dp_hist:H2O%d",item->first),Form("Dp_hist:H_{2}O%d",item->first),1000,-0.1,0.1);

			item->second->Project(OptDpArrayWithFocalCutH[item->first]->GetName(),Form("%s.gold.dp",HRS.Data()),
		    focalCut);
		}
		SieveThetaPhiPCanv->cd(SieveThetaPhiList.size()*2+canvas_counter_temp);
		OptDpArrayWithFocalCutH[item->first]->Draw();
		canvas_counter_temp++;
	}
	SieveThetaPhiPCanv->Update();


	// plot the detail of the x and theta

	TCanvas *focalPlaneXCanvas=new TCanvas("Focal Plane Variables x","Focal Plane Variables x",1000,1000);
	focalPlaneXCanvas->Divide(1,2);
	focalPlaneXCanvas->cd(2)->Divide(4,1);

	for(auto item=sieveFocalXArray.begin();item!=sieveFocalXArray.end(); item++){
			item->second->SetLineColor(8+item->first);
			focalPlaneXCanvas->cd(1);
			if (item == sieveFocalXArray.begin()){
				item->second->Draw();
			}else{
				item->second->Draw("same");
			}

			// get peak value and put in the canvas
			// fit the peak with gaus
			double GuasFitPars[3];
			auto maximumPeak=item->second->GetXaxis()->GetBinCenter(item->second->GetMaximumBin());
			item->second->Fit("gaus","R0Q","ep",maximumPeak-0.005,maximumPeak+0.005);
			item->second->GetFunction("gaus")->GetParameters(GuasFitPars);

			TLatex *peakInfor=new TLatex(GuasFitPars[1],GuasFitPars[0],Form("%f",GuasFitPars[1]));
			peakInfor->SetLineWidth(2);
			peakInfor->SetTextSize(0.04);
			peakInfor->Draw("same");


			// locate each canvas set the range
			focalPlaneXCanvas->cd(2)->cd(item->first+3);
			auto plot=(TH1F *)item->second->Clone(Form("%s_d",item->second->GetName()));
			plot->GetXaxis()->SetRangeUser(GuasFitPars[1]-8*GuasFitPars[2],GuasFitPars[1]+8*GuasFitPars[2]);
			plot->SetLineWidth(3);
			plot->Draw();
			peakInfor->Draw("same");

		}
	focalPlaneXCanvas->cd(1);
	focalLegend_x->Draw("same");
	focalPlaneXCanvas->Update();

	//draw the detail of the theta part
	TCanvas *focalPlaneThetaCanvas=new TCanvas("Focal Plane Variables theta","Focal Plane Variables theta",1000,1000);
	focalPlaneThetaCanvas->Divide(1,2);
	focalPlaneThetaCanvas->cd(2)->Divide(4,1);
	for(auto item=sieveFocalThetaArray.begin();item!=sieveFocalThetaArray.end(); item++){
			item->second->SetLineColor(8+item->first);
			item->second->SetLineWidth(2);

			focalPlaneThetaCanvas->cd(1);
			if (item == sieveFocalThetaArray.begin()){
				item->second->Draw();
			}else{
				item->second->Draw("same");
			}

			// get peak value and put in the canvas
			// fit the peak with gaus
			double GuasFitPars[3];
			auto maximumPeak=item->second->GetXaxis()->GetBinCenter(item->second->GetMaximumBin());
			item->second->Fit("gaus","R0Q","ep",maximumPeak-0.005,maximumPeak+0.005);
			item->second->GetFunction("gaus")->GetParameters(GuasFitPars);
			TLatex *peakInfor=new TLatex(GuasFitPars[1],GuasFitPars[0],Form("%f",GuasFitPars[1]));
			peakInfor->SetLineWidth(2);
			peakInfor->SetTextSize(0.04);
			peakInfor->Draw("same");

			focalPlaneThetaCanvas->cd(2)->cd(item->first+3);
			auto plot=(TH1F *)item->second->Clone(Form("%s_d",item->second->GetName()));
			plot->GetXaxis()->SetRangeUser(GuasFitPars[1]-8*GuasFitPars[2],GuasFitPars[1]+8*GuasFitPars[2]);
			plot->SetLineWidth(3);
			plot->Draw();
			peakInfor->Draw("same");
		}
	focalPlaneThetaCanvas->cd(1);
	focalLegend_th->Draw("same");
	focalPlaneThetaCanvas->Update();

	//plot the 2-D focal plane x-y
	// used for plot the focal plane x.y in the central sieve cut, in ground states
	std::map<int, TH2F *>SieveFocalGroundXY;
	for (auto item = chainArray.begin(); item!= chainArray.end();item++){
		if(item->first < 10){
			SieveFocalGroundXY[item->first]=new TH2F(Form("C_{12} Dp%d%% focal plane x vs. y",item->first),Form("C_{12} Dp%d%% focal plane x vs. y",item->first),1000,-0.2,0.3,1000,-0.01,0.01);
			// project the data
			TCut focalCut=TCut(generalcut.Data())+dpCutArray[item->first]+TCut(cutg->GetName());
			item->second->Project(SieveFocalGroundXY[item->first]->GetName(),Form("%s.tr.y:%s.tr.x",HRS.Data(),HRS.Data()),focalCut);
			SieveFocalGroundXY[item->first]->GetXaxis()->SetTitle(Form("%s.tr.x",HRS.Data()));
			SieveFocalGroundXY[item->first]->GetYaxis()->SetTitle(Form("%s.tr.y",HRS.Data()));
		}else{
			SieveFocalGroundXY[item->first]=new TH2F(Form("H_{2}O Dp%d%% focal plane x vs. y",item->first),Form("C_{12} Dp%d%% focal plane x vs. y",item->first),1000,-0.4,0.4,1000,-0.01,0.01);
			// project the data
			TCut focalCut=TCut(generalcut.Data())+dpCutArray[item->first]+TCut(cutg->GetName());
			item->second->Project(SieveFocalGroundXY[item->first]->GetName(),Form("%s.tr.y:%s.tr.x",HRS.Data(),HRS.Data()),focalCut);
		}
	}

	//draw the plot on canvas
	TCanvas *sieveFocalXYCanv=new TCanvas("Focal Canvas X. Y","Focal Canvas X. Y",1000,1000);
	sieveFocalXYCanv->Divide(1,2);
	sieveFocalXYCanv->cd(2)->Divide(4,1);
	sieveFocalXYCanv->Draw();
	// plot the canvas holes
	for (auto item = SieveFocalGroundXY.begin(); item!= SieveFocalGroundXY.end();item++){
		if (item->first < 10){
			sieveFocalXYCanv->cd(1);
			if(item==SieveFocalGroundXY.begin()){
				item->second->Draw("zcol");
			}else{

				item->second->Draw("zcol same");
			}
			// get the fit parameter and set the range the set the cut
			std::string str("");
			TLatex *text=new TLatex(SieveFocalXGuasFitPar[item->first][1],SieveFocalYGuasFitPar[item->first][1]+0.002,Form("Dp%d%%(%1.4f,%1.4f)",item->first,SieveFocalXGuasFitPar[item->first][1],SieveFocalYGuasFitPar[item->first][1]));
			text->SetTextColor(kRed);
//			text->SetTextSize(0.01);
			text->Draw("same");

			sieveFocalXYCanv->cd(2)->cd(item->first+3);
			auto plot=(TH2F *)item->second->Clone(Form("focal x vs. y dp%d%%_d",item->first));
			plot->SetTitle(Form("focal x vs. y dp%d%%_d",item->first));
			//plot->SetAxisRange(SieveFocalXGuasFitPar[item->first][1]-4*SieveFocalXGuasFitPar[item->first][0],SieveFocalXGuasFitPar[item->first][1]+4*SieveFocalXGuasFitPar[item->first][0],"X");
			//plot->SetAxisRange(SieveFocalYGuasFitPar[item->first][1]-4*SieveFocalYGuasFitPar[item->first][0],SieveFocalYGuasFitPar[item->first][1]+4*SieveFocalYGuasFitPar[item->first][0],"Y");
			plot->Draw("zcol");
		}
	}
	sieveFocalXYCanv->Update();
}

