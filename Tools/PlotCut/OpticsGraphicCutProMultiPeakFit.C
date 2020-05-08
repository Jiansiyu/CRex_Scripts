/*
 * OpticsGraphicCutPro.C
 *
 *  Created on: Jan 9, 2020
 *      Author: newdriver
 */




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
TString generalcutR="R.gold.p > 2.14 && R.gold.p < 2.2";;//"R.tr.n==1 && R.vdc.u1.nclust==1&& R.vdc.v1.nclust==1 && R.vdc.u2.nclust==1 && R.vdc.v2.nclust==1 && fEvtHdr.fEvtType==1 && R.gold.p > 2.14 && R.gold.p < 2.2";
TString generalcutL="L.tr.n==1 && L.vdc.u1.nclust==1&& L.vdc.v1.nclust==1 && L.vdc.u2.nclust==1 && L.vdc.v2.nclust==1 && fEvtHdr.fEvtType==1 && L.gold.p > 2.14 && L.gold.p < 2.2";


struct PeakFitInfor{
public:
	PeakFitInfor(double_t mean, double_t sigma,double_t height=10){
		Fit_mean=mean;
		Fit_sigma=sigma;
		Fit_height=height;
	};
	~PeakFitInfor(){};

	double_t GetMean(){
		return Fit_mean;
	}

	double_t GetSigma(){
		return Fit_sigma;
	}

	double_t GetHeight(){
		return Fit_height;
	}

private:
	double_t Fit_mean;
	double_t Fit_sigma;
	double_t Fit_height;

};

inline Bool_t IsFileExist (const std::string& name) {
	  struct stat buffer;
	  return (stat (name.c_str(), &buffer) == 0);
}

Int_t OpticsGraphicCutPro(UInt_t runID,TString folder="/home/newdriver/Storage/Research/CRex_Experiment/RasterReplay/Replay/Result") {
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
//	fCrystalMomentum->Draw("same");
	fCrystalMomentum->GetParameters(fCrystalMomentumPar);

	SieveRecCanvas->Update();
	

	// fit the second and the third peak in the historgram
	double C2stp=fCrystalMomentumPar[6]-((7.65407-4.43982)*0.001);
	TF1 *fSecondGuas=new TF1 ("secondstatesgaus","gaus",C2stp-0.0006,C2stp+0.0004);
	double fSecondGuasPar[3];
	momentum->Fit("secondstatesgaus","R","ep",fSecondGuas->GetXmin(),fSecondGuas->GetXmax());
	fSecondGuas->GetParameters(fSecondGuasPar);
	
	// Fit the third Peak
	double C3stp=fCrystalMomentumPar[6]-(9.641-4.43982)*0.001;
	TF1 *fThirdGuas=new TF1 ("Thirdstatesgaus","gaus",C3stp-0.0008,C3stp+0.0004);
	double fThirdGuasPar[3];
	momentum->Fit("Thirdstatesgaus","R","ep",fThirdGuas->GetXmin(),fThirdGuas->GetXmax());
	fThirdGuas->GetParameters(fThirdGuasPar);
	

	// put all the peaks on the same fit plot 5 + 5 + 3 + 3 parameters

	double MomCombinedFitPar[16];
	std::cout<<"Start the Global Fit with: "<<sizeof(MomCombinedFitPar)/sizeof(double)<<"  Parameters!!"<<std::endl;
	std::copy(fCrystalMomentumPar,fCrystalMomentumPar+10,MomCombinedFitPar);
	std::copy(fSecondGuasPar,fSecondGuasPar+3,MomCombinedFitPar+10);
	std::copy(fThirdGuasPar,fThirdGuasPar+3,MomCombinedFitPar+13);
	TF1 *fMomentumCombinedFit=new TF1("fMomentumCombinedFit","crystalball(0)+crystalball(5)+gaus(10)+gaus(13)",fThirdGuas->GetXmin(),fgroundCrystalball->GetXmax());
	fMomentumCombinedFit->SetParameters(MomCombinedFitPar);
	momentum->Fit("fMomentumCombinedFit","","",fMomentumCombinedFit->GetXmin(),fMomentumCombinedFit->GetXmax());
	fMomentumCombinedFit->Draw("same");

	std::map<unsigned int, PeakFitInfor*>SpectrumFitResult;
	SpectrumFitResult[0]=new PeakFitInfor(MomCombinedFitPar[1],MomCombinedFitPar[2],MomCombinedFitPar[0]);
	SpectrumFitResult[1]=new PeakFitInfor(MomCombinedFitPar[6],MomCombinedFitPar[7],MomCombinedFitPar[5]);
	SpectrumFitResult[2]=new PeakFitInfor(MomCombinedFitPar[11],MomCombinedFitPar[12],MomCombinedFitPar[10]);
	SpectrumFitResult[3]=new PeakFitInfor(MomCombinedFitPar[14],MomCombinedFitPar[15],MomCombinedFitPar[13]);

	for (auto peak_indix=SpectrumFitResult.begin(); peak_indix!=SpectrumFitResult.end(); peak_indix++){
		TLatex *PeakMomInfor=new TLatex(peak_indix->second->GetMean(),peak_indix->second->GetHeight(),Form("P=%2.4f #pm %1.2fx10^{-3}GeV",peak_indix->second->GetMean(),peak_indix->second->GetSigma()*1000.0));
		PeakMomInfor->SetTextSize(0.035);
		PeakMomInfor->SetTextAlign(12);
		PeakMomInfor->SetTextColor(2);
		//PeakMomInfor->SetTextAngle(30);
		PeakMomInfor->Draw("same");
	}


	//create the peak separation
	for(auto Spectro_index=SpectrumFitResult.size()-1;Spectro_index>=1;Spectro_index--){
		//get the energy separation and the error propagation
		double_t SeparationE=SpectrumFitResult.at(0)->GetMean()-SpectrumFitResult.at(Spectro_index)->GetMean();
		double_t SeparationError=sqrt(pow(SpectrumFitResult.at(0)->GetSigma(),2)+pow(SpectrumFitResult.at(Spectro_index)->GetSigma(),2)); //
		TLatex *deltaEtxt=new TLatex(SpectrumFitResult.at(Spectro_index)->GetMean(), SpectrumFitResult.at(Spectro_index)->GetHeight()/10,Form("#DeltaP=%1.3f",SeparationE*1000.0));
		deltaEtxt->SetTextSize(0.035);
		deltaEtxt->SetTextAlign(12);
		deltaEtxt->SetTextColor(2);

		deltaEtxt->Draw("same");
	}

	// create the line of each peak
	for (auto peak_indix=SpectrumFitResult.begin(); peak_indix!=SpectrumFitResult.end(); peak_indix++){
		TLine *peakline=new TLine(peak_indix->second->GetMean(),0,peak_indix->second->GetMean(),peak_indix->second->GetHeight()*1.1);
		peakline->SetLineColor(3);
		peakline->SetLineWidth(2);
		peakline->Draw("same");
	}

	// create the peak separation for each individual
	TPaveText *pt = new TPaveText(0.1,0.7,0.3,0.9,"NDC");
	pt->AddText("Energy levels of ^{12}C  (MeV #pm keV)");
	pt->AddText("E_1=4.43982 #pm 0.21");
	pt->AddText("E_2=7.65407 #pm 0.19");
	pt->AddText("E_3=9.641   #pm 5   ");
	pt->Draw("same");




	//plot the bigger plot for first excited states
	SieveRecCanvas->cd(2)->cd(3);
	TH1F *groundStats=(TH1F *)momentum->Clone("C-12 p.ground");
	groundStats->GetXaxis()->SetRangeUser(fCrystalMomentumPar[1]-0.002,fCrystalMomentumPar[1]+0.002);
	groundStats->Draw();

	SieveRecCanvas->cd(2)->cd(4);
	TH1F *firstStats=(TH1F *)momentum->Clone("C-12 p.first");
	firstStats->GetXaxis()->SetRangeUser(fCrystalMomentumPar[6]-0.002,fCrystalMomentumPar[6]+0.002);
	firstStats->Draw();

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

	SieveRecCanvas->Update();
	hSieveHole->Delete();
	
	
}

