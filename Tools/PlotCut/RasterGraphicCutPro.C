/*
 * OpticsGraphicCutPro.C
 *
 *  Created on: Jan 9, 2020
 *      Author: newdriver
 *
 *      How to do the correction :
 *
 *      > chop the rang, and according to the range, calculate the energy difference between the central,
 *      > calculate the momentum change
 *      > subtrack all the event with this ammount of momentum change to all the event in this enrgy range
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
TString generalcutR="R.gold.p > 2.14 && R.gold.p < 2.2 && R.tr.n==1 && R.vdc.u1.nclust==1&& R.vdc.v1.nclust==1 && R.vdc.u2.nclust==1 && R.vdc.v2.nclust==1 ";//&& fEvtHdr.fEvtType==1 ";
TString generalcutL="L.tr.n==1 && L.vdc.u1.nclust==1&& L.vdc.v1.nclust==1 && L.vdc.u2.nclust==1 && L.vdc.v2.nclust==1 && fEvtHdr.fEvtType==1 && L.gold.p > 2.14 && L.gold.p < 2.2";

// Fit the Carbon-12 with Multi-Crystalball function combination
//
TF1 * SpectroCrystalFit_C12(TH1F *momentumSpectro){
	TF1 *globalFit;

	TH1F *FitInitialh=(TH1F *)momentumSpectro->Clone("initialh");
	// fit the highest peak, this should be the ground states peak
	auto CGroundp=FitInitialh->GetXaxis()->GetBinCenter(FitInitialh->GetMaximumBin());
	auto C1stp=CGroundp-0.00443891;
	return globalFit;
}

inline Bool_t IsFileExist (const std::string& name) {
	  struct stat buffer;
	  return (stat (name.c_str(), &buffer) == 0);
}

std::vector<TString> chopRange(const double_t E_low, const double_t E_high,const  int nbin,const TString epicName){
	double numBin=nbin;
	double_t bin=(E_high-E_low)/numBin;
	std::vector<TString> cutList;
	for(double currentE= E_low;currentE<E_high; currentE+=bin){
		double_t lower_boundary=currentE;
		double_t Upper_boundary=currentE+bin;
		TString cut=Form("(%s>=%f)&&(%s<%f)",epicName.Data(),lower_boundary,epicName.Data(),Upper_boundary);
		cutList.push_back(cut);
	}
	return cutList;
}

Int_t OpticsGraphicCutPro(UInt_t runID,TString folder="/home/newdriver/Storage/Research/CRex_Experiment/optReplay/Result") {
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

	SieveRecCanvas->Divide(1, 3);
	SieveRecCanvas->cd(2)->Divide(4, 1);
	SieveRecCanvas->cd(3)->Divide(4, 1);
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
	TH1F *momentum=new TH1F("C-12 gold.p","C-12 gold.p",1500,2.1,2.2);
	chain->Project(momentum->GetName(),Form("%s.gold.p",HRS.Data()),Form("%s && %s",generalcut.Data(),cutg->GetName()));
	// get the maximum bin, this should be the first excited states
	auto CGroundp=momentum->GetXaxis()->GetBinCenter(momentum->GetMaximumBin());
	auto C1stp=CGroundp-0.00443891;
	momentum->GetXaxis()->SetRangeUser(CGroundp-0.0044*3,CGroundp+0.0044*2);
	momentum->GetXaxis()->SetTitle("gold.p");
	momentum->GetYaxis()->SetTitle("#");

	momentum->Draw();
	double_t fgroudGausPar[3];
	double_t ffirstGuasPar[3];
	TF1 *fgroudGaus=new TF1("groudstatesgaus","gaus",CGroundp-0.0005,CGroundp+0.0005);
	momentum->Fit("groudstatesgaus","R","ep",fgroudGaus->GetXmin(),fgroudGaus->GetXmax());
	//fgroudGaus->Draw("same");
	fgroudGaus->GetParameters(fgroudGausPar);

	{
		// check the first excited states
		TH1F *firstexcited_temp=new TH1F("C-12 gold.p_temp","C-12 gold.p_temp",2000,2.1,2.2);
		chain->Project(firstexcited_temp->GetName(),Form("%s.gold.p",HRS.Data()),Form("%s && %s && %s.gold.p> %f && %s.gold.p < %f ",generalcut.Data(),cutg->GetName(), HRS.Data(),CGroundp-0.0044*3,HRS.Data(),fgroudGausPar[1]-0.003));
		C1stp=firstexcited_temp->GetXaxis()->GetBinCenter(firstexcited_temp->GetMaximumBin());
		firstexcited_temp->Delete();
	}

	TF1 *ffirstGuas=new TF1 ("firststatesgaus","gaus",C1stp-0.0006,C1stp+0.0004);
	momentum->Fit("firststatesgaus","R","ep",ffirstGuas->GetXmin(),ffirstGuas->GetXmax());
	//ffirstGuas->Draw("same");
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

	TLine *firstposLine=new TLine(fCrystalMomentumPar[6],0,fCrystalMomentumPar[6],ffirstGuasPar[0]*1.1);
	firstposLine->SetLineColor(3);
	firstposLine->SetLineWidth(2);
	firstposLine->Draw("same");

	TPaveText *pt = new TPaveText(0.1,0.8,0.3,0.9,"NDC");
	pt->AddText(Form("%1.3f MeV (%2.2f\%%)",1000.0*(fCrystalMomentumPar[1]-fCrystalMomentumPar[6]),100.0*abs(abs(fCrystalMomentumPar[1]-fCrystalMomentumPar[6])-0.00443891)/0.00443891));
	pt->Draw("same");

	TLatex *t1 = new TLatex(fgroudGausPar[1] + 2 * fgroudGausPar[2],fgroudGausPar[0], Form("P=%2.5fGeV #sigma=%1.2f x 10^{-3}", fCrystalMomentumPar[1],fCrystalMomentumPar[2]*1000));
	t1->SetTextSize(0.055);
	t1->SetTextAlign(12);
	t1->SetTextColor(2);
	t1->Draw("same");

	TLatex *t2 = new TLatex(ffirstGuasPar[1] + 2 * ffirstGuasPar[2],ffirstGuasPar[0], Form("P=%2.5fGeV #sigma=%1.2f x 10^{-3}", fCrystalMomentumPar[6],fCrystalMomentumPar[7]*1000));
	t2->SetTextSize(0.055);
	t2->SetTextAlign(12);
	t2->SetTextColor(2);
	t2->Draw("same");

	//plot the bigger plot for first excited states
	SieveRecCanvas->cd(2)->cd(3);
	TH1F *groundStats=(TH1F *)momentum->Clone("C-12 p.ground");
	groundStats->GetXaxis()->SetRangeUser(fCrystalMomentumPar[1]-0.002,fCrystalMomentumPar[1]+0.002);
	groundStats->Draw();
	groudposLine->Draw("same");

	SieveRecCanvas->cd(2)->cd(4);
	TH1F *firstStats=(TH1F *)momentum->Clone("C-12 p.first");
	firstStats->GetXaxis()->SetRangeUser(fCrystalMomentumPar[6]-0.002,fCrystalMomentumPar[6]+0.002);
	firstStats->Draw();
	firstposLine->Draw("same");

	SieveRecCanvas->Update();


	// plot the target variables plot the canvas on x and y
	SieveRecCanvas->cd(3)->cd(1);
	TH1F *CutSieveY=(TH1F *) gROOT->FindObject("CutSieveY");
	if(CutSieveY){
		CutSieveY->Clear();
	}else{
		CutSieveY=new TH1F("CutSieveY","CutSieveY",100,-0.0001,0.0001);
	}
	chain->Project(CutSieveY->GetName(), Form("%s.gold.y",HRS.Data()),Form("%s && %s ", generalcut.Data(), cutg->GetName()));
	CutSieveY->Fit("gaus");
	auto CutSieveYGausFit=CutSieveY->GetFunction("gaus");
	double CutSieveYGausFitPar[3];
	CutSieveYGausFit->GetParameters(CutSieveYGausFitPar);

	CutSieveY->Draw("same");
	TLatex *CutSieveYt1 = new TLatex(CutSieveYGausFitPar[1] - 10 * CutSieveYGausFitPar[2],CutSieveYGausFitPar[0], Form("%f x 10^{-3}",1000.0*CutSieveYGausFitPar[1]));
	CutSieveYt1->SetTextSize(0.055);
	CutSieveYt1->SetTextAlign(12);
	CutSieveYt1->SetTextColor(2);
	CutSieveYt1->Draw("same");
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

	// start to add the Raster information
    TCanvas *RasterDiagnoseCanv=(TCanvas *) gROOT->GetListOfCanvases()->FindObject("Raster Diagnose Canvas");
	if(! RasterDiagnoseCanv){
		RasterDiagnoseCanv=new TCanvas("Raster Diagnose","Raster Disgnose",1200,1800);
	}else{
		RasterDiagnoseCanv->Clear();
	}
	RasterDiagnoseCanv->Divide(4,3);
	
	// re-orgnize the plot
	TH2F * RasterTgXhh=(TH2F *)gROOT->FindObject("RasterXTgXhh");
	if(RasterTgXhh){
		RasterTgXhh->Clear();
	}else{
		// instead of using hard coded boundary, need to find a better way to define the boundary
		RasterTgXhh= new TH2F("RasterXTgXhh","Raster Current X vs. Target X",1000,-0.001,0.001,1000,20000,70000);
	}
	TText * CoeffTxt;	
	std::string CarbonGroundCut=Form("abs(%s.gold.p-%f)<3.0*%f",HRS.Data(),fCrystalMomentumPar[1],fCrystalMomentumPar[2]);
    chain->Project(RasterTgXhh->GetName(),Form("%srb.Raster2.rawcur.x:%s.gold.x",HRS.Data(),HRS.Data()),Form("%s && %s", generalcut.Data(),cutg->GetName()));
	
	RasterDiagnoseCanv->cd(1);
	RasterTgXhh->Draw("zcol");	
	RasterDiagnoseCanv->cd(1)->Update();	
	
	CoeffTxt=new TText(-0.0002,25000,Form("coefficient:%1.4f",RasterTgXhh->GetCorrelationFactor()));
	CoeffTxt->Draw("same");

	TH2F * RasterTgXSievePCuthh=(TH2F *)gROOT->FindObject("RasterXTgX_SP_hh");
        if(RasterTgXSievePCuthh){
                RasterTgXSievePCuthh->Clear();
        }else{
                // instead of using hard coded boundary, need to find a better way to define the boundary
                RasterTgXSievePCuthh= new TH2F("RasterXTgX_SP_hh","Raster Current X vs. Target X{Sieve P cut}",1000,-0.001,0.001,1000,20000,70000);
        }
        
       chain->Project(RasterTgXSievePCuthh->GetName(),Form("%srb.Raster2.rawcur.x:%s.gold.x",HRS.Data(),HRS.Data()),Form("%s && %s && %s", generalcut.Data(),cutg->GetName(), CarbonGroundCut.c_str()));
       RasterDiagnoseCanv->cd(2);
	RasterTgXSievePCuthh->Draw("zcol");
	CoeffTxt=new TText(-0.0002,25000,Form("coefficient:%1.4f",RasterTgXSievePCuthh->GetCorrelationFactor()));
        CoeffTxt->Draw("same");
	
	//momentum vs. raster current plot
	TH2F *RasterXTgPhh=(TH2F *) gROOT->FindObject("RasterTgPhh");
	if(RasterXTgPhh){
		RasterXTgPhh->Clear();
	}else{
		RasterXTgPhh=new TH2F("RasterTgPhh","Raster Current X vs. P",200,fCrystalMomentumPar[1]-4*fCrystalMomentumPar[2],fCrystalMomentumPar[1]+4*fCrystalMomentumPar[2],1000,20000,70000);
	}      
	chain->Project(RasterXTgPhh->GetName(),Form("%srb.Raster2.rawcur.x:%s.gold.p",HRS.Data(),HRS.Data()),Form("%s && %s && %s", generalcut.Data(), cutg->GetName(),CarbonGroundCut.c_str()));
	RasterDiagnoseCanv->cd(3);
	RasterXTgPhh->Draw("zcol");
	CoeffTxt=new TText(2.1745,25000,Form("coefficient:%1.4f",RasterXTgPhh->GetCorrelationFactor()));
        CoeffTxt->Draw("same");
	// project the same plot for Y dimension 
	

	TH2F * RasterTgYhh = (TH2F *) gROOT->FindObject("RasterYTgYhh");
	if (RasterTgYhh) {
		RasterTgYhh->Clear();
	} else {
		RasterTgYhh = new TH2F("RasterYTgYhh", "Raster Current Y vs. Target Y",
				1000, -0.001, 0.001, 1000, 20000, 70000);
	}
	chain->Project(RasterTgYhh->GetName(),Form("%srb.Raster2.rawcur.y:%s.gold.y",HRS.Data(),HRS.Data()),Form("%s && %s", generalcut.Data(),cutg->GetName()));
	RasterDiagnoseCanv->cd(5);
	RasterTgYhh->Draw("zcol");
	
	CoeffTxt=new TText(-0.0002,25000,Form("coefficient:%1.4f",RasterTgYhh->GetCorrelationFactor()));
    CoeffTxt->Draw("same");


	TH2F * RasterTgYSievePCuthh=(TH2F *)gROOT->FindObject("RasterYTgY_SP_hh");
    if(RasterTgYSievePCuthh){
        RasterTgYSievePCuthh->Clear();
    }else{
		RasterTgYSievePCuthh= new TH2F("RasterYTgY_SP_hh","Raster Current X vs. Target X{Sieve P cut}",1000,-0.001,0.001,1000,20000,70000);
    }
    chain->Project(RasterTgYSievePCuthh->GetName(),Form("%srb.Raster2.rawcur.y:%s.gold.y",HRS.Data(),HRS.Data()),Form("%s && %s && %s", generalcut.Data(),cutg->GetName(), CarbonGroundCut.c_str()));
	RasterDiagnoseCanv->cd(6);
	RasterTgYSievePCuthh->Draw("zcol");
	
	CoeffTxt=new TText(-0.0002,25000,Form("coefficient:%1.4f",RasterTgYSievePCuthh->GetCorrelationFactor()));
    CoeffTxt->Draw("same");

	TH2F *RasterYTgPhh=(TH2F *) gROOT->FindObject("RasterYTgPhh");
        if(RasterYTgPhh){
                RasterYTgPhh->Clear();
        }else{
                RasterYTgPhh=new TH2F("RasterYTgPhh","Raster_Current_Y_vs_P",80,fCrystalMomentumPar[1]-4*fCrystalMomentumPar[2],fCrystalMomentumPar[1]+4*fCrystalMomentumPar[2],1000,20000,70000);
        }

    chain->Project(RasterYTgPhh->GetName(),Form("%srb.Raster2.rawcur.y:%s.gold.p",HRS.Data(),HRS.Data()),Form("%s && %s && %s", generalcut.Data(), cutg->GetName(),CarbonGroundCut.c_str()));
    RasterDiagnoseCanv->cd(7);
    RasterYTgPhh->Draw("zcol");
	CoeffTxt=new TText(2.1743,25000,Form("coefficient:%1.4f",RasterYTgPhh->GetCorrelationFactor()));
    CoeffTxt->Draw("same");

    // fit the minimum bin center
	TH2F *RasterYTgPhhBinFiltered=(TH2F *)RasterYTgPhh->Clone("Filtered");
     for(auto binx=0; binx<(RasterYTgPhh->GetXaxis()->GetNbins());binx++){
    	 for(auto biny=0; biny < (RasterYTgPhh->GetYaxis()->GetNbins());biny++){
    		// apply the cut on the bin content size
    		double_t binValue=RasterYTgPhh->GetBinContent(binx,biny);
    		if(binValue<7.0){
    			RasterYTgPhhBinFiltered->SetBinContent(binx,biny,0);
    		}
    	 }
     }
     RasterDiagnoseCanv->cd(8);
     RasterYTgPhhBinFiltered->Fit("pol1");
     TF1 *RasterYTgPhhBinFilteredFit=RasterYTgPhhBinFiltered->GetFunction("pol1");
     RasterYTgPhhBinFiltered->Draw("colz");
     RasterDiagnoseCanv->cd(7);
     RasterYTgPhhBinFilteredFit->Draw("same");


    // generate X-Y corrolation
    // current X vs. target Y
    TH2F *RasterXTgYhh=(TH2F *) gROOT->FindObject("RasterXTgYhh");
    if(!RasterXTgYhh){
    	RasterXTgYhh=new TH2F("RasterXTgYhh", "Raster Current X vs. Target Y{Sieve Cut}",
    			1000, 0.00004, 0.00012, 1000, 20000, 70000);
    }else{
    	RasterXTgYhh->Clear();
    }
    chain->Project(RasterXTgYhh->GetName(),Form("%srb.Raster2.rawcur.x:%s.gold.y",HRS.Data(),HRS.Data()),Form("%s && %s", generalcut.Data(),cutg->GetName()));
    RasterDiagnoseCanv->cd(9);
    RasterXTgYhh->Draw("zcol");
    CoeffTxt=new TText(RasterXTgYhh->ProjectionX()->GetMean(),25000,Form("coefficient:%1.4f",RasterXTgYhh->GetCorrelationFactor()));
    CoeffTxt->Draw("same");
    TH2F *RasterXTgYCuthh=(TH2F *) gROOT->FindObject("RasterXTgYCuthh");
    if(!RasterXTgYCuthh){
    	RasterXTgYCuthh=new TH2F("RasterXTgYCuthh", "Raster Current X vs. Target Y{Sieve P Cut}",
    			1000, 0.00004, 0.00012, 1000, 20000, 70000);
    }else{
    	RasterXTgYCuthh->Clear();
    }
    chain->Project(RasterXTgYCuthh->GetName(),Form("%srb.Raster2.rawcur.x:%s.gold.y",HRS.Data(),HRS.Data()),Form("%s && %s && %s", generalcut.Data(),cutg->GetName(), CarbonGroundCut.c_str()));
    RasterDiagnoseCanv->cd(10);
    RasterXTgYCuthh->Draw("zcol");
    CoeffTxt=new TText(RasterXTgYCuthh->ProjectionX()->GetMean(),25000,Form("coefficient:%1.4f",RasterXTgYCuthh->GetCorrelationFactor()));
    CoeffTxt->Draw("same");


    // draw the plot on cuyrrent Y vs. Target X
    TH2F *RasterYTgXhh=(TH2F *) gROOT->FindObject("RasterYTgXhh");
    if(!RasterYTgXhh){
    	RasterYTgXhh=new TH2F("RasterYTgXhh", "Raster Current Y vs. Target X {Sieve Cut}",
    			1000, -0.0001, 0.0001, 1000, 20000, 70000);
    }else{
    	RasterYTgXhh->Clear();
    }
    chain->Project(RasterYTgXhh->GetName(),Form("%srb.Raster2.rawcur.y:%s.gold.x",HRS.Data(),HRS.Data()),Form("%s && %s", generalcut.Data(),cutg->GetName()));
    RasterDiagnoseCanv->cd(11);
    RasterYTgXhh->Draw("colz");

    // draw the plot with cut
    TH2F *RasterYTgXCuthh=(TH2F *) gROOT->FindObject("RasterYTgXCuthh");
    if(!RasterYTgXCuthh){
    	RasterYTgXCuthh=new TH2F("RasterYTgXCuthh", "Raster Current Y vs. Target X{Sieve P Cut}",
    			1000, -0.0001, 0.0001, 1000, 20000, 70000);
    }else{
    	RasterYTgXCuthh->Clear();
    }
    chain->Project(RasterYTgXCuthh->GetName(),Form("%srb.Raster2.rawcur.y:%s.gold.x",HRS.Data(),HRS.Data()),Form("%s && %s && %s", generalcut.Data(),cutg->GetName(), CarbonGroundCut.c_str()));
    RasterDiagnoseCanv->cd(12);
    RasterYTgXCuthh->Draw("colz");

    RasterDiagnoseCanv->Update();



    // first do the correction
    // chop the Y dimension, for x-dimension for need another way to correct it
    // re-bin and get the center momentum of the 50th bin as reference
    TCanvas *RasterCorrDiagnoseCanv= new TCanvas ("RasterCorrDiagnoseCanv", "Raster Correction Diagnose Canv", 1800,1200);
    RasterCorrDiagnoseCanv->Divide(1,3);
    RasterCorrDiagnoseCanv->cd(1)->Divide(4,1);
    RasterCorrDiagnoseCanv->cd(2)->Divide(4,1);
    RasterCorrDiagnoseCanv->cd(3)->Divide(2,1);
    RasterCorrDiagnoseCanv->Draw();

    // get the Y of the reference bin, and the other correction term will take this as reference
    TH1F* RasterCurrentXh=(TH1F *)gROOT->FindObject("RasterCurrentXh");
    if(RasterCurrentXh){
    	RasterCurrentXh->Clear();

    }else{

    	RasterCurrentXh=new TH1F("RasterCurrentXh","Raster Raw Current X ",1000, 25000, 70000);
    }
   chain->Project(RasterCurrentXh->GetName(),Form("%srb.Raster2.rawcur.x",HRS.Data()),Form("%s && %s && %s", generalcut.Data(),cutg->GetName(), CarbonGroundCut.c_str()));

   RasterCorrDiagnoseCanv->cd(1)->cd(1);
   RasterCurrentXh->Draw();


   TH2F *RasterCorrectedYTgPhh = (TH2F *) gROOT->FindObject("RasterCorrectedYTgPhh");
	if (RasterCorrectedYTgPhh) {
		RasterCorrectedYTgPhh->Clear();
	} else {
		RasterCorrectedYTgPhh = new TH2F("RasterCorrectedYTgPhh", "Raster Corrected raw current Y vs. P", 80,
				fCrystalMomentumPar[1] - 4 * fCrystalMomentumPar[2],
				fCrystalMomentumPar[1] + 4 * fCrystalMomentumPar[2], 1000,
				20000, 70000);
	}
	//take the center if the bin tobe the reference

	TH1F *MomhRefenceChop=(TH1F *)gROOT->FindObject("MomhRefenceChop");
	if(MomhRefenceChop){
		MomhRefenceChop->Clear();
	}else{
		MomhRefenceChop=new TH1F("MomhRefenceChop","MomhRefenceChop",RasterCorrectedYTgPhh->GetXaxis()->GetNbins(),RasterCorrectedYTgPhh->GetXaxis()->GetXmin(),RasterCorrectedYTgPhh->GetXaxis()->GetXmax());
		// project the momentum
	}
	auto cut_temp=chopRange(RasterCurrentXh->GetXaxis()->GetXmin(),RasterCurrentXh->GetXaxis()->GetXmax(), 100,Form("%srb.Raster2.rawcur.y", HRS.Data())).at(50);

	chain->Project(MomhRefenceChop->GetName(),
			Form("%s.gold.p", HRS.Data()),
			Form("%s && %s && %s && %s", generalcut.Data(), cutg->GetName(),CarbonGroundCut.c_str(),cut_temp.Data()));


	TH1F *RasterCorrectedmomentumhh=(TH1F*)gROOT->FindObject("Corrected C-12 gold.p");
	if(RasterCorrectedmomentumhh){
		RasterCorrectedmomentumhh->Clear();
	}else{
		RasterCorrectedmomentumhh=new TH1F("Corrected C-12 gold.p","Corrected  C-12 gold.p",1500,2.1,2.2);
	}

	double_t ReferenceP=2.17466;//MomhRefenceChop->GetMean();
	//draw the reference point on the canvas
	std::cout<<"referenece:"<<ReferenceP<<std::endl;

	// plot the reference point on the canvas
	RasterCorrDiagnoseCanv->cd(1)->cd(2);
	RasterYTgPhh->Draw("zcol");
	TH2F *RasterYTgPhhMarker = (TH2F *) gROOT->FindObject("RasterYTgPhhMarker");
	if (RasterYTgPhhMarker) {
		RasterYTgPhhMarker->Clear();
	} else {
		RasterYTgPhhMarker = new TH2F("RasterYTgPhhMarker",
				"RasterYTgPhhMarker", 80,
				fCrystalMomentumPar[1] - 4 * fCrystalMomentumPar[2],
				fCrystalMomentumPar[1] + 4 * fCrystalMomentumPar[2], 1000,
				20000, 70000);
	}
	RasterYTgPhhMarker->Fill(ReferenceP,RasterCurrentXh->GetXaxis()->GetXmin()+(RasterCurrentXh->GetXaxis()->GetXmax()-RasterCurrentXh->GetXaxis()->GetXmin())/2.0);
	RasterYTgPhhMarker->SetMarkerStyle(43);
	RasterYTgPhhMarker->SetMarkerColor(2);
	RasterYTgPhhMarker->SetMarkerSize(3);
	RasterYTgPhhMarker->Draw("same");
	RasterYTgPhhBinFiltered->Draw("same");
	RasterCorrDiagnoseCanv->Update();

	TH1F *RasterCorrectedmomentumhh_temp;
	TH2F *RasterCorrectedYTgPhh_temp;
	for (auto chopElement: chopRange(RasterCurrentXh->GetXaxis()->GetXmin(),RasterCurrentXh->GetXaxis()->GetXmax(),100, Form("%srb.Raster2.rawcur.y",HRS.Data()))){
		//take the element and took the reference
		//take the difference, and take take the standard divation
		TH1F *temp=(TH1F *)gROOT->FindObject("MomYhChop");
		if(temp){
			temp->Clear();
		}else{
			temp=new TH1F("MomYhChop","MomYhChop",RasterCorrectedYTgPhh->GetXaxis()->GetNbins(),RasterCorrectedYTgPhh->GetXaxis()->GetXmin(),RasterCorrectedYTgPhh->GetXaxis()->GetXmax());
			// project the momentum
		}
		chain->Project(temp->GetName(),Form("%s.gold.p", HRS.Data()),Form("%s && %s && %s && %s", generalcut.Data(),cutg->GetName(), CarbonGroundCut.c_str(),chopElement.Data()));


		RasterCorrectedmomentumhh_temp=new TH1F("RasterCorrectedmomentumhh_temp","RasterCorrectedmomentumhh_temp",1500,2.1,2.2);

		RasterCorrectedYTgPhh_temp = new TH2F("RasterCorrectedYTgPhh_temp", "RasterCorrectedYTgPhh_temp", 80,
					fCrystalMomentumPar[1] - 4 * fCrystalMomentumPar[2],
					fCrystalMomentumPar[1] + 4 * fCrystalMomentumPar[2], 1000,
					20000, 70000);

		//chop the momentum and correct the momentum
		double_t MomCorrection=ReferenceP-temp->GetMean();
		std::cout<<"Cut Range::"<<chopElement.Data()<<std::endl;
		chain->Project(RasterCorrectedYTgPhh_temp->GetName(),Form("%srb.Raster2.rawcur.y:(%s.gold.p+%f)",HRS.Data(),HRS.Data(),MomCorrection),Form("%s && %s && %s && %s", generalcut.Data(),cutg->GetName(), CarbonGroundCut.c_str(),chopElement.Data()));
		chain->Project(RasterCorrectedmomentumhh_temp->GetName(),Form("%s.gold.p+%f",HRS.Data(),MomCorrection),Form("%s && %s && %s", generalcut.Data(),cutg->GetName(),chopElement.Data()));

		RasterCorrectedYTgPhh->Add(RasterCorrectedYTgPhh_temp);
		RasterCorrectedmomentumhh->Add(RasterCorrectedmomentumhh_temp);

		RasterCorrectedYTgPhh_temp->Delete();
		RasterCorrectedmomentumhh_temp->Delete();
		temp->Delete();

		RasterCorrDiagnoseCanv->cd(1)->cd(3);
		RasterCorrectedYTgPhh->Draw("zcol");
		RasterCorrDiagnoseCanv->cd(1)->cd(4);
		RasterCorrectedmomentumhh->Draw("zcol");
		RasterCorrDiagnoseCanv->Update();

		// get the number of the event in each chop and plot the
	}

	RasterCorrDiagnoseCanv->cd(2)->cd(1);
	RasterCorrectedYTgPhh->Draw("zcol");

	RasterCorrDiagnoseCanv->cd(3)->cd(1);
	RasterCorrDiagnoseCanv->cd(3)->cd(1)->SetLogy();
	RasterCorrectedmomentumhh->GetXaxis()->SetRangeUser(CGroundp-0.0044*3,CGroundp+0.0044*2);
	RasterCorrectedmomentumhh->Fit("fCrystalMomentum","","",fCrystalMomentum->GetXmin(),fCrystalMomentum->GetXmax());
	RasterCorrectedmomentumhh->Draw("zcol");
	TH1F *RasterMomemtum=(TH1F *)momentum->Clone("InitialMom");
	RasterMomemtum->SetLineColor(3);
	RasterMomemtum->Draw("same");
	// add labels
	double_t fCrystalRasterCorrectedmomentumPar[10];
	RasterCorrectedmomentumhh->GetFunction("fCrystalMomentum")->GetParameters(fCrystalRasterCorrectedmomentumPar);
	TLatex *RasterMomemtumLabel1 = new TLatex(fgroudGausPar[1] + 2 * fCrystalRasterCorrectedmomentumPar[2],fCrystalRasterCorrectedmomentumPar[0], Form("P=%2.5fGeV #sigma=%1.2f x 10^{-3}", fCrystalRasterCorrectedmomentumPar[1],fCrystalRasterCorrectedmomentumPar[2]*1000));
	RasterMomemtumLabel1->SetTextSize(0.055);
	RasterMomemtumLabel1->SetTextAlign(12);
	RasterMomemtumLabel1->SetTextColor(2);
	RasterMomemtumLabel1->Draw("same");

	TLatex *RasterMomemtumLabel2 = new TLatex(fgroudGausPar[1] + 2 * fCrystalRasterCorrectedmomentumPar[2],fCrystalRasterCorrectedmomentumPar[0], Form("P=%2.5fGeV #sigma=%1.2f x 10^{-3}", fCrystalRasterCorrectedmomentumPar[6],fCrystalRasterCorrectedmomentumPar[7]*1000));
	RasterMomemtumLabel2->SetTextSize(0.055);
	RasterMomemtumLabel2->SetTextAlign(12);
	RasterMomemtumLabel2->SetTextColor(2);
	RasterMomemtumLabel2->Draw("same");

	RasterCorrDiagnoseCanv->Update();


	// plot the canvas, using the fit functions
	//Method two
	double_t RasterYTgPhhBinFilteredFitPar[2];
	RasterYTgPhhBinFilteredFit->GetParameters(RasterYTgPhhBinFilteredFitPar);
    TH2F *RasterCorrected2YTgPhh = (TH2F *) gROOT->FindObject("RasterCorrected2YTgPhh");
	if (RasterCorrected2YTgPhh ) {
		RasterCorrected2YTgPhh ->Clear();
	} else {
		RasterCorrected2YTgPhh  = new TH2F("RasterCorrected2YTgPhh", "Raster Corrected raw current Y vs. P", 80,
				fCrystalMomentumPar[1] - 4 * fCrystalMomentumPar[2],
				fCrystalMomentumPar[1] + 4 * fCrystalMomentumPar[2], 1000,
				20000, 70000);
	}

	TH1F *RasterCorrected2momentumhh=(TH1F*)gROOT->FindObject("RasterCorrected2momentumhh");
	if(RasterCorrected2momentumhh){
		RasterCorrected2momentumhh->Clear();
	}else{
		RasterCorrected2momentumhh=new TH1F("RasterCorrected2momentumhh","Corrected  C-12 gold.p",1500,2.1,2.2);
	}

	TString MomCorrection_temp=Form("(%f-(%srb.Raster2.rawcur.y-%f)/%f)",ReferenceP,HRS.Data(),RasterYTgPhhBinFilteredFitPar[0],RasterYTgPhhBinFilteredFitPar[1]);
	chain->Project(RasterCorrected2YTgPhh->GetName(),Form("%srb.Raster2.rawcur.y:(%s.gold.p+%s)",HRS.Data(),HRS.Data(),MomCorrection_temp.Data()),Form("%s && %s && %s", generalcut.Data(),cutg->GetName(),CarbonGroundCut.c_str()));
	chain->Project(RasterCorrected2momentumhh->GetName(),Form("(%s.gold.p+%s)",HRS.Data(),MomCorrection_temp.Data()),Form("%s && %s", generalcut.Data(),cutg->GetName()));
	RasterCorrDiagnoseCanv->cd(2)->cd(2);
	RasterCorrected2YTgPhh->Draw("zcol");

	RasterCorrDiagnoseCanv->cd(3)->cd(2);
	RasterCorrDiagnoseCanv->cd(3)->cd(2)->SetLogy();
	Double_t fCrystalRasterCorrected2momentumPar[10];
	std::copy(fCrystalMomentumPar,fCrystalMomentumPar+10,fCrystalRasterCorrected2momentumPar);
	TF1 *fCrystalRasterCorrected2momentumhh=new TF1("fCrystalRasterCorrected2momentumhh","crystalball(0)+crystalball(5)",ffirstCrystal->GetXmin(),fgroundCrystalball->GetXmax());
	RasterCorrectedmomentumhh->GetFunction("fCrystalMomentum")->GetParameters(fCrystalRasterCorrected2momentumPar);
	fCrystalRasterCorrected2momentumhh->SetParameters(fCrystalRasterCorrected2momentumPar);
	RasterCorrected2momentumhh->Fit("fCrystalRasterCorrected2momentumhh","","",fCrystalRasterCorrected2momentumhh->GetXmin(),fCrystalRasterCorrected2momentumhh->GetXmax());
	fCrystalRasterCorrected2momentumhh->GetParameters(fCrystalRasterCorrected2momentumPar);

	TH1F *RasterMomemtum1=(TH1F *)momentum->Clone("InitialMom1");
	RasterCorrected2momentumhh->GetXaxis()->SetRangeUser(CGroundp-0.0044*3,CGroundp+0.0044*2);
	RasterCorrected2momentumhh->Draw();
	RasterMomemtum1->SetLineColor(3);
	RasterMomemtum1->Draw("same");

	TLatex *RasterMomemtum1Label1 = new TLatex(fgroudGausPar[1] + 2 * fCrystalRasterCorrectedmomentumPar[2],fCrystalRasterCorrectedmomentumPar[0], Form("P=%2.5fGeV #sigma=%1.2f x 10^{-3}", fCrystalRasterCorrected2momentumPar[1],fCrystalRasterCorrected2momentumPar[2]*1000));
	RasterMomemtum1Label1->SetTextSize(0.055);
	RasterMomemtum1Label1->SetTextAlign(12);
	RasterMomemtum1Label1->SetTextColor(2);
	RasterMomemtum1Label1->Draw("same");

	TLatex *RasterMomemtum2Label2 = new TLatex(fgroudGausPar[1] + 2 * fCrystalRasterCorrectedmomentumPar[2],fCrystalRasterCorrectedmomentumPar[0], Form("P=%2.5fGeV #sigma=%1.2f x 10^{-3}", fCrystalRasterCorrected2momentumPar[6],fCrystalRasterCorrected2momentumPar[7]*1000));
	RasterMomemtum2Label2->SetTextSize(0.055);
	RasterMomemtum2Label2->SetTextAlign(12);
	RasterMomemtum2Label2->SetTextColor(2);
	RasterMomemtum2Label2->Draw("same");

	RasterCorrDiagnoseCanv->Update();

}

