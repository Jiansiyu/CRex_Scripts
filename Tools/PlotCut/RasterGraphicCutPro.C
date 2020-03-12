/*
 * OpticsGraphicCutPro.C
 *
 *  Created on: Jan 9, 2020
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

inline Bool_t IsFileExist (const std::string& name) {
	  struct stat buffer;
	  return (stat (name.c_str(), &buffer) == 0);
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

	SieveRecCanvas->Divide(1, 2);
	SieveRecCanvas->cd(2)->Divide(4, 1);
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


}

