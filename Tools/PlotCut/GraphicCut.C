/*
 * Graphic Cut Function that used for apply Graphic Cut
 *
 * author: Siyu Jian
 */

#include <TCut.h>
#include <TCutG.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TString.h>
#include <TChain.h>
#include <TTree.h>
#include <map>
#include <vector>
#include <random>
#include <TROOT.h>
#include <iostream>
#include <sys/stat.h>
std::vector <TString> cutNames;

inline Bool_t IsFileExist (const std::string& name) {
	  struct stat buffer;
	  return (stat (name.c_str(), &buffer) == 0);
}
int GraphicCut(Int_t runID) {

	TChain *chain=new TChain("T");

	TString rootDir="/home/newdriver/Storage/Research/PRex_Experiment/prex_analyzer/optReplay/Result";
			//"/home/newdriver/Storage/Research/PRex_Workspace/PREX-MPDGEM/PRexScripts/Tools/ClumAlign/replay/Result";

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

    TCanvas *b=new TCanvas("t","t",600,600);

    
    TH1F *HistP=new TH1F("gold.p","gold.p",1000,0.93,0.95);
    TString generalCut="fEvtHdr.fEvtType==1";
    chain->Project(HistP->GetName(),"R.gold.p",Form("%s && %s",generalCut.Data(),cutg->GetName()));
    HistP->Draw();
//
//    TH2F *HistP=new TH2F("gold.p","gold.p",500,-0.1,0.1,500,-0.1,0.1);
//    //TString generalCut="fEvtHdr.fEvtType==1";//"R.vdc.u1.nclust==1&& R.vdc.v1.nclust==1 && R.vdc.u2.nclust==1 && R.vdc.v2.nclust==1 && R.gold.dp<1 && R.gold.dp > -0.1 && fEvtHdr.fEvtType==1";
//    TString generalCut="R.vdc.u1.nclust==1&& R.vdc.v1.nclust==1 && R.vdc.u2.nclust==1 && R.vdc.v2.nclust==1 && R.gold.dp<1 && R.gold.dp > -0.1";
//
//    chain->Project(HistP->GetName(),"RGEM.tr.th:RGEM.tr.ph",Form("%s && %s",generalCut.Data(),cutg->GetName()));
//    HistP->Draw("ZCOL");
//
//    b->cd(2);
//    TH2F *HistAll=new TH2F("GEM TH vs GEM ph all","GEM TH vs GEM ph all",500,-0.1,0.1,500,-0.1,0.1);
//    //TString generalCut="fEvtHdr.fEvtType==1";//"R.vdc.u1.nclust==1&& R.vdc.v1.nclust==1 && R.vdc.u2.nclust==1 && R.vdc.v2.nclust==1 && R.gold.dp<1 && R.gold.dp > -0.1 && fEvtHdr.fEvtType==1";
////    TString generalCut="R.vdc.u1.nclust==1&& R.vdc.v1.nclust==1 && R.vdc.u2.nclust==1 && R.vdc.v2.nclust==1 && R.gold.dp<1 && R.gold.dp > -0.1";
//
//    chain->Project(HistAll->GetName(),"RGEM.tr.th:RGEM.tr.ph",Form("%s",generalCut.Data()));
//    HistAll->Draw("ZCOL");
    
    b->Update();




	return 1;
}
