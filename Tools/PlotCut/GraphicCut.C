/*
 * Graphic Cut Function that used for apply Graphic Cut
 *
 * author: Siyu Jian
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


TString generalcut;
TString generalcutR="R.tr.n==1 && R.vdc.u1.nclust==1&& R.vdc.v1.nclust==1 && R.vdc.u2.nclust==1 && R.vdc.v2.nclust==1 && fEvtHdr.fEvtType==1 && R.gold.p > 2.14 && R.gold.p < 2.2";
TString generalcutL="L.tr.n==1 && L.vdc.u1.nclust==1&& L.vdc.v1.nclust==1 && L.vdc.u2.nclust==1 && L.vdc.v2.nclust==1 && fEvtHdr.fEvtType==1 && L.gold.p > 2.14 && L.gold.p < 2.2";

std::vector <TString> cutNames;

inline Bool_t IsFileExist (const std::string& name) {
	  struct stat buffer;
	  return (stat (name.c_str(), &buffer) == 0);
}
int GraphicCut(UInt_t runID,TString folder="/home/newdriver/Storage/Research/CRex_Experiment/optReplay/Result") {

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

	TCanvas *a=new TCanvas("cut","cut",600,600);

	a->Draw();
	a->cd(0);

	TH2F *HHistThPh=new TH2F("th vs ph","th vs ph",1000,-0.045,0.045,1000,-0.045,0.045);
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

//    TCanvas *b=new TCanvas("t","t",600,600);
//
//
//    TH1F *HistP=new TH1F("gold.p","gold.p",1000,0.93,0.95);
//    TString generalCut="fEvtHdr.fEvtType==1";
//    chain->Project(HistP->GetName(),"R.gold.p",Form("%s && %s",generalCut.Data(),cutg->GetName()));
//    HistP->Draw();
//
//
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
	SieveRecCanvas->cd(2)->cd(2);


	TH2F *h =(TH2F *) HHistThPh->Clone();
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
		TH1F *momentum=new TH1F("C-12 gold.p","C-12 gold.p",1000,2.1,2.2);
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
			TH1F *firstexcited_temp=new TH1F("C-12 gold.p_temp","C-12 gold.p_temp",500,2.1,2.2);
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
		pt->AddText(Form("%1.3f MeV (%2.2f\%)",1000.0*(fCrystalMomentumPar[1]-fCrystalMomentumPar[6]),100.0*abs(abs(fCrystalMomentumPar[1]-fCrystalMomentumPar[6])-0.00443891)/0.00443891));
		pt->Draw("same");

		TLatex *t1 = new TLatex(fgroudGausPar[1] + 2 * fgroudGausPar[2],fgroudGausPar[0], Form("P=%2.5fMeV", fCrystalMomentumPar[1]));
		t1->SetTextSize(0.055);
		t1->SetTextAlign(12);
		t1->SetTextColor(2);
		t1->Draw("same");

		TLatex *t2 = new TLatex(ffirstGuasPar[1] + 2 * ffirstGuasPar[2],ffirstGuasPar[0], Form("P=%2.5fMeV", fCrystalMomentumPar[6]));
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
	return 1;
}
