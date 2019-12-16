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
int GraphicCutCheck(Int_t runID) {

	TChain *chain=new TChain("T");

	TString rootDir="/home/newdriver/Storage/Research/PRex_Workspace/PREX-MPDGEM/PRexScripts/Tools/ClumAlign/replay/Result";

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

    TCanvas *ZalignYCanv=new TCanvas("Y","Y",600,600);
    TCanvas *ZalignXCanv=new TCanvas("X","X",600,600);
    ZalignYCanv->Divide(3,2);
    ZalignXCanv->Divide(3,2);


    TCanvas *ZalignYCanvProf=new TCanvas("YProf","YProf",600,600);
    TCanvas *ZalignXCanvProf=new TCanvas("XProf","XProf",600,600);
    ZalignYCanvProf->Divide(3,2);
    ZalignXCanvProf->Divide(3,2);

    TString generalCut="R.vdc.u1.nclust==1&& R.vdc.v1.nclust==1 && R.vdc.u2.nclust==1 && R.vdc.v2.nclust==1 && R.gold.dp<1 && R.gold.dp > -0.1";


    std::map<Int_t, TH2F *>ZalignY_ph;
    std::map<Int_t, TH2F *>ZalignX_th;



    for (int i =1; i <7; i ++){
    	ZalignY_ph[i]=new TH2F(Form("RGEM.rgems.y%d.coord.resid:RGEM.tr.ph",i),Form("RGEM.rgem.y%d.coord,resid:RGEM.tr.ph",i),1000,-0.03,0.03,1000,-0.03,0.03);
    	ZalignX_th[i]=new TH2F(Form("RGEM.rgem.x%d.coord.resid:RGEM.tr.th",i),Form("RGEM.rgem.x%d.coord,resid:RGEM.tr.th",i),1000,-0.03,0.03,1000,-0.03,0.03);

    	chain->Project(ZalignY_ph[i]->GetName(),Form("RGEM.rgems.y%d.coord.resid:RGEM.tr.ph",i),Form("%s && %s",generalCut.Data(),cutg->GetName()));
    	chain->Project(ZalignX_th[i]->GetName(),Form("RGEM.rgems.x%d.coord.resid:RGEM.tr.th",i),Form("%s && %s",generalCut.Data(),cutg->GetName()));
    	ZalignYCanv->cd(i);
    	ZalignY_ph[i]->Draw("zcol");


    	ZalignXCanv->cd(i);
    	ZalignX_th[i]->Draw("zcol");

    	ZalignYCanvProf->cd(i);
    	ZalignY_ph[i]->Draw("profileY");

    	ZalignXCanvProf->cd(i);
    	ZalignX_th[i]->Draw("profileY");


    }

	return 1;
}
