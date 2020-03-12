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
#include <sstream>
#include <TComplex.h>
#include <TVirtualPad.h>

#include <TSpectrum2.h>
#include <TF2.h>
#include <TObject.h>
#include "TMinuit.h"
#include <TFile.h>
#include <fstream>
#include <TSystem.h>
#include <TApplication.h>
#include <TLatex.h>
#include <TGApplication.h>

#include <boost/filesystem.hpp>
R__LOAD_LIBRARY(/usr/lib/x86_64-linux-gnu/libboost_filesystem.so)


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
TString generalcutR="R.tr.n==1 && R.vdc.u1.nclust==1&& R.vdc.v1.nclust==1 && R.vdc.u2.nclust==1 && R.vdc.v2.nclust==1 && R.gold.dp<1 && R.gold.dp > -0.1 && fEvtHdr.fEvtType==1";
TString generalcutL="L.tr.n==1 && L.vdc.u1.nclust==1&& L.vdc.v1.nclust==1 && L.vdc.u2.nclust==1 && L.vdc.v2.nclust==1 && fEvtHdr.fEvtType==1 && L.gold.p > 2.14 && L.gold.p < 2.2";
//////////////////////////////////////////////////////////////////////////////
// Work Directory
//////////////////////////////////////////////////////////////////////////////
//TString WorkDir = "Result/Test/";
//TString WorkDir = "/home/newdriver/Storage/Research/CRex_Experiment/optReplay/Result/RHRS_Feb292020/";
//TString WorkDir = "/home/newdriver/Storage/Research/PRex_Workspace/PREX-MPDGEM/PRexScripts/Tools/PlotCut/Result/RHRS_20200311/";
TString WorkDir = "/home/newdriver/Storage/Research/PRex_Workspace/PREX-MPDGEM/PRexScripts/Tools/PlotCut/Result/Cut20200311/RHRS/";

TString CutSuf = ".FullCut.root";
TString CutDescFileSufVertex = ".VertexCut.cut";
TString CutDescFileSufDp = ".DpCut.cut";
TString CutDescFileSufSieve = ".SieveCut.%d_%d.cut";
TString RootFileName;

//LHRS
int numberofSieveHoles[13]={0,0,0,5,6,5,5,6,5,5,4,3,2};
int minSieveHoles[13]=     {0,0,0,1,0,1,1,0,1,1,1,2,2};
//RHRS
//int numberofSieveHoles[13]={0,0,0,6,6,5,5,6,5,5,4,5,3};
//int minSieveHoles[13]=     {0,0,0,0,0,1,1,0,1,1,1,1,1};

inline Bool_t IsFileExist (const std::string& name) {
	  struct stat buffer;
	  return (stat (name.c_str(), &buffer) == 0);
}

Int_t cutPro(UInt_t runID,UInt_t current_col=3,TString folder="/home/newdriver/Storage/Research/CRex_Experiment/RasterReplay/Replay/Result/") {
	// need to check the folder
	std::string bufferedWorkFolder;
	std::string bufferedSourceDir;
	int bufferedCol=-1;
	int bufferedRunID=-1;


	if(!boost::filesystem::is_regular_file("logfile.txt")){
		bufferedSourceDir=folder;
		bufferedWorkFolder=WorkDir;
		bufferedCol=3;
		bufferedRunID=runID;

		// create the folder and save the infor
	}else{
		std::ifstream textinfile("logfile.txt");
		textinfile>>bufferedSourceDir>>bufferedWorkFolder>>bufferedRunID>>bufferedCol;

		if((bufferedSourceDir==folder)&&(bufferedWorkFolder==WorkDir)&&(bufferedRunID==runID)){
			bufferedCol++;
		}else{
			bufferedSourceDir=folder;
			bufferedWorkFolder=WorkDir;
			bufferedRunID=runID;
			bufferedCol=3;
		}
		std::cout<<"source dir:"<<bufferedSourceDir.c_str()<<std::endl;
		std::cout<<"current work dir: " << bufferedWorkFolder.c_str()<<"\n  runID:"<<bufferedRunID<<"\n  current col: "<< bufferedCol<<std::endl;
		// update the run infor
	}

	std::ofstream textoutfile;
	textoutfile.open("logfile.txt", std::ios::trunc);
	textoutfile <<bufferedSourceDir.c_str()<<" "<<bufferedWorkFolder.c_str()<< " "<<bufferedRunID<<" "<<bufferedCol << std::endl;

	current_col=bufferedCol;

	gStyle->SetOptStat(0);
	// prepare the data
	TChain *chain=new TChain("T");
	TString rootDir(folder.Data());
	TString HRS="R";
	if(runID>20000){ //RHRS
		if(IsFileExist(Form("%s/prexRHRS_%d_-1.root",rootDir.Data(),runID))){
			std::cout<<"Add File::"<<Form("%s/prexRHRS_%d_-1.root",rootDir.Data(),runID)<<std::endl;
			RootFileName=Form("prexRHRS_%d_-1.root",runID);
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
			std::cout<<"Cannot find file :"<<Form("%s/prexRHRS_%d_-1.root",rootDir.Data(),runID)<<std::endl;
		}
	}else{
		HRS="L";
		if(IsFileExist(Form("%s/prexLHRS_%d_-1.root",rootDir.Data(),runID))){
			RootFileName=Form("prexLHRS_%d_-1.root",runID);
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
		mainPatternCanvas=new TCanvas("cutPro","cutPro",1000,1200);
	}else{
		mainPatternCanvas->Clear();
	}
//	TCanvas *mainPatternCanvas=new TCanvas("cut","cut",600,600);
	mainPatternCanvas->Draw();
	TH2F *TargetThPhHH=(TH2F *)gROOT->FindObject("th_vs_ph");
	if(TargetThPhHH) TargetThPhHH->Delete();
	TargetThPhHH=new TH2F("th_vs_ph","th_vs_ph",1000,-0.025,0.025,1000,-0.047,0.045);

	chain->Project(TargetThPhHH->GetName(),Form("%s.gold.th:%s.gold.ph",HRS.Data(),HRS.Data()),generalcut.Data());
	TargetThPhHH->Draw("zcol");
	mainPatternCanvas->SetGridx(10);
	mainPatternCanvas->SetGridy(10);

//	mainPatternCanvas->Update();
	// input how start row and how many holes in this row
	col=current_col;
    int nhol = 0;
    std::cout << "How many holes in this No." << col << " column?" << std::endl;
//    std::cin >> nhol;
    nhol=numberofSieveHoles[col];
    row_count=nhol;
    std::cout<<numberofSieveHoles[col]<<std::endl;
//    row_count=numberofSieveHoles[col];
    if(nhol < 0)return 0;
    std::cout << "min hole id : ";
    int rmin = -1;
//    std::cin >> rmin;
    rmin=minSieveHoles[col];
    row_min=rmin;
    std::cout<<minSieveHoles[col]<<std::endl;
    row=row_min;

    if(rmin < 0)return 0;
    TLatex *p=new TLatex(-0.025, 0.045,Form("Sieve %d, #frac{%d}{%d}",col,rmin,nhol));
    p->Draw("same");
    mainPatternCanvas->Update();
	mainPatternCanvas->ToggleEventStatus();
	mainPatternCanvas->AddExec("ex", "DynamicCanvas()");
	return 1;
}

//Recognize the save patter
// save the rec hole to the folder
// ground state
void SavePatternHole(double momentumSigmaCut=3.0){
	//search all the holes in this col and save in the folder
	std::cout<<std::endl<<std::endl;
	std::cout<<"*******Save process start ........*******"<<std::endl;
	std::cout<<"	Searching for holes in col ("<<col<<")"<<std::endl;

	TCanvas *SaveCheckCanvas=(TCanvas *) gROOT->GetListOfCanvases()->FindObject("SieveSaveCheck");
	if(!SaveCheckCanvas){
		SaveCheckCanvas =new TCanvas("SieveSaveCheck","SieveSaveCheck",1600,1080);
	}else{
		SaveCheckCanvas->Clear();
	}
	SaveCheckCanvas->Divide(1,2);
	SaveCheckCanvas->cd(1)->Divide(2,1);
	SaveCheckCanvas->cd(2)->Divide(row_count,1);

	TString workdir_temp=WorkDir;
	if(momentumSigmaCut>10.0){
		workdir_temp+="/WithOutMomCut/";
	}else{
		workdir_temp+="/GroundMomCut/";
	}
	// check the existance of the folder, if not create the folder
	if(!boost::filesystem::is_directory(workdir_temp.Data())){
		std::cout<<"Folder :"<< workdir_temp.Data()<<"  Does Not EXIST\n Trying to Create the folder"<<std::endl;
		if(!boost::filesystem::create_directories(workdir_temp.Data())){
			std::cout<<"ERROR:: cannot creat folde, please check the permission!!!!"<<std::endl;
			exit(-1);
		}
	}

	TString CutFileName = workdir_temp + RootFileName + CutSuf;
	TString TempString(Form(CutDescFileSufSieve.Data(), FoilID, col));
	TString PlotDir(RootFileName + Form(".hcut_R_%d_%d/", FoilID, col));
	TString CutDescName = workdir_temp + RootFileName + TempString;
	// prepare the filename and folder
	//TString CutDescName = WorkDir + RootFileName + CutDescFileSufDp;
	//TString CutFileName = WorkDir + RootFileName + CutSuf; // used for save the root file and cut file
	TFile *f1=new TFile(CutFileName,"UPDATE");
	assert(f1);

	std::fstream cutdesc(CutDescName, std::ios_base::out);
	assert(cutdesc.is_open());

	SaveCheckCanvas->cd(1)->cd(1);
	// attached the root file and plot the canvas used for check
	TChain *chain = (TChain *) gROOT->FindObject("T");
	if(!chain) std::cout<<"[ERROR] CAN NOT FIND TREE"<<std::endl;
	TH2F *h=(TH2F*)gROOT->FindObject("th_vs_ph");
	if(!h)
		h=new TH2F("th_vs_ph1","th_vs_ph1",1000,-0.03,0.03,1000,-0.045,0.045);

	TH2F *SavesieveCheck = new TH2F("SaveSieveCheck", "SaveSieveCheck",
			h->GetXaxis()->GetNbins(), h->GetXaxis()->GetXmin(),
			h->GetXaxis()->GetXmax(), h->GetYaxis()->GetNbins(),
			h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax());

	TString HRS("R");
	TString filename(chain->GetFile()->GetName());
	if (filename.Contains("RHRS")) {
	} else if (filename.Contains("LHRS")) {
		HRS = "L";
	}
	chain->Project(SavesieveCheck->GetName(),Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()),Form("%s",generalcut.Data()));
	//chain->Draw(Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()),Form("%s",generalcut.Data()),"zcol");
	SavesieveCheck->Draw("zcol");
	// search how many cut file exist in the root buffer
	for (int i = 0 ; i < row_min; i ++)cutdesc<<"fEvtHdr.fRun==0"<< std::endl;
	std::vector<int> sieveIDList;
	TH2F *sievehole[row_count];
	TH1F *sieveholemomentum[row_count];    // ground states momentum
	TF1  *sieveholemomentumGausFit[row_count];
	TLatex * momentumInfor1[row_count];
	TCutG *cutg;
	for (int row_iter = row_min; row_iter<row_min+row_count;row_iter++){

		cutg=(TCutG *) gROOT->FindObject(Form("hcut_R_%d_%d_%d",FoilID,col,row_iter));
		if(cutg){
			sieveIDList.push_back(row_iter-row_min);
			sievehole[row_iter-row_min] = new TH2F(
					Form("hcut_R_%d_%d_%d_hh", FoilID, col, row_iter),
					Form("hcut_R_%d_%d_%d_hh", FoilID, col, row_iter),
					h->GetXaxis()->GetNbins(), h->GetXaxis()->GetXmin(),
					h->GetXaxis()->GetXmax(), h->GetYaxis()->GetNbins(),
					h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax());
			SaveCheckCanvas->cd(2)->cd(row_iter - row_min + 1);
			chain->Project(sievehole[row_iter-row_min]->GetName(), Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()),Form("%s && %s",cutg->GetName(),generalcut.Data()));

			sieveholemomentum[row_iter-row_min]=new TH1F(Form("hcut_R_%d_%d_%d_h_momentum", FoilID, col, row_iter),Form("hcut_R_%d_%d_%d_momentum", FoilID, col, row_iter),600,2.1,2.25);
			chain->Project(sieveholemomentum[row_iter-row_min]->GetName(),Form("%s.gold.p",HRS.Data()),Form("%s && %s",cutg->GetName(),generalcut.Data()));
			sieveholemomentumGausFit[row_iter-row_min]=new TF1(Form("1ststatesDpgaushcut_R_%d_%d_%d", FoilID, col, row_iter),"gaus",
					sieveholemomentum[row_iter-row_min]->GetXaxis()->GetBinCenter(sieveholemomentum[row_iter-row_min]->GetMaximumBin())-0.002,
					sieveholemomentum[row_iter-row_min]->GetXaxis()->GetBinCenter(sieveholemomentum[row_iter-row_min]->GetMaximumBin())+0.002);
			sieveholemomentumGausFit[row_iter-row_min]->SetParameter(1,sieveholemomentum[row_iter-row_min]->GetXaxis()->GetBinCenter(sieveholemomentum[row_iter-row_min]->GetMaximumBin()));
			sieveholemomentum[row_iter - row_min]->Fit(
					Form("1ststatesDpgaushcut_R_%d_%d_%d", FoilID, col,
							row_iter), "R", "ep",
					sieveholemomentumGausFit[row_iter-row_min]->GetXmin(),
					sieveholemomentumGausFit[row_iter-row_min]->GetXmax());

			sieveholemomentum[row_iter - row_min]->GetXaxis()->SetRangeUser(
					sieveholemomentum[row_iter - row_min]->GetXaxis()->GetBinCenter(
							sieveholemomentum[row_iter - row_min]->GetMaximumBin())
							- 0.009,
					sieveholemomentum[row_iter - row_min]->GetXaxis()->GetBinCenter(
							sieveholemomentum[row_iter - row_min]->GetMaximumBin())
							+ 0.004);
			sieveholemomentum[row_iter - row_min]->Draw();
			sieveholemomentumGausFit[row_iter-row_min]->Draw("same");

			auto groudpcenter=sieveholemomentumGausFit[row_iter-row_min]->GetParameter(1);
			auto groudpsigma=sieveholemomentumGausFit[row_iter-row_min]->GetParameter(2);

			momentumInfor1[row_iter - row_min] =
					new TLatex(groudpcenter + 2 * groudpsigma-0.007,
							sieveholemomentumGausFit[row_iter - row_min]->GetParameter(0),
							Form("P_{0}= %2.5f",
									groudpcenter));//, groudpsigma
			momentumInfor1[row_iter - row_min]->SetTextSize(0.055);
			momentumInfor1[row_iter - row_min]->SetTextAlign(12);
			momentumInfor1[row_iter - row_min]->SetTextColor(2);
			momentumInfor1[row_iter - row_min]->Draw("same");

			SaveCheckCanvas->cd(2)->cd(row_iter - row_min + 1)->Update();
			if(groudpsigma>0.0008)groudpsigma=0.0008;
			// plot the boundary of the cut
			TLine *leftboundary=new TLine(groudpcenter-momentumSigmaCut*groudpsigma,0,groudpcenter-momentumSigmaCut*groudpsigma,(0.9*SaveCheckCanvas->cd(2)->cd(row_iter - row_min + 1)->GetUymax()));
			leftboundary->SetLineColor(3);
			leftboundary->SetLineWidth(2);
			leftboundary->Draw("same");

			TLine *rightboundary=new TLine(groudpcenter+momentumSigmaCut*groudpsigma,0,groudpcenter+momentumSigmaCut*groudpsigma,(0.9*SaveCheckCanvas->cd(2)->cd(row_iter - row_min + 1)->GetUymax()));
			rightboundary->SetLineColor(3);
			rightboundary->SetLineWidth(2);
			rightboundary->Draw("same");

			SaveCheckCanvas->cd(1)->cd(1);
			cutg->SetName(Form("hcut_R_%d_%d_%d", FoilID, col, row_iter));
			cutg->SetVarX(Form("%s.gold.ph",HRS.Data()));
			cutg->SetVarY(Form("%s.gold.th",HRS.Data()));
			cutg->SetLineColor(kMagenta);
			cutg->SetLineWidth(2);
			cutg->Draw("PL same");
			SaveCheckCanvas->cd(1)->cd(1)->SetGridx(20);
			SaveCheckCanvas->cd(1)->cd(1)->SetGridy(20);
			//plot the momentum and apply cut on the momentum

			cutg->Write("", TObject::kOverwrite); // Overwrite old cut
			if(groudpsigma>0.0008)groudpsigma=0.0008;
			cutdesc << Form("hcut_R_%d_%d_%d", FoilID, col, row_iter) << " && ";

			if(momentumSigmaCut>10.0){
				cutdesc << (const char*)generalcut << std::endl;
			}else{
			cutdesc << (const char*)generalcut <<" && "
					<<Form("abs(%s.gold.p-%f)<%f*%f",HRS.Data(),groudpcenter,momentumSigmaCut,groudpsigma)
					<< std::endl;
			}

			SaveCheckCanvas->cd(1)->cd(2);
			sievehole[row_iter-row_min]->Draw("same");
			TLatex *eventCountLable=new TLatex(sievehole[row_iter-row_min]->GetMean(1) + 0.005,sievehole[row_iter-row_min]->GetMean(2), Form("Entries(%d,%d): %2.0f",col,row_iter, (sievehole[row_iter-row_min]->GetEntries())));
			eventCountLable->SetTextSize(0.03);
			eventCountLable->SetTextAlign(12);
			eventCountLable->SetTextColor(2);
			eventCountLable->Draw("same");
			SaveCheckCanvas->cd(1)->cd(2)->SetGridx(20);
			SaveCheckCanvas->cd(1)->cd(2)->SetGridy(20);
			//sievehole[row_iter-row_min]->Delete();
			cutg->Draw("same");
			SaveCheckCanvas->Update();

		}else{
			//if the cut does not exist,then write the cut
			cutdesc << "fEvtHdr.fRun==0" << std::endl;
		}

	}

	for(int i = row_min+row_count; i < NSieveRow; i++)
		cutdesc << "fEvtHdr.fRun==0" << std::endl;
	SaveCheckCanvas->SetName(Form("CutProfcut_R_%d_%d",FoilID, col));
	SaveCheckCanvas->Write("", TObject::kOverwrite);
	SaveCheckCanvas->SaveAs(Form("%s/%s.hcut_R_%d_%d.jpg",workdir_temp.Data(),RootFileName.Data(),FoilID, col));

	for(auto i : sieveIDList){
//	for (unsigned int i = 0; i < row_count; i++) {
		if (!sievehole[i]->IsZombie()) {
			sievehole[i]->Delete();
		}
		if (!sieveholemomentum[i]->IsZombie()) {
			sieveholemomentum[i]->Delete();
		}
		if (!sieveholemomentumGausFit[i]->IsZombie()) {
			sieveholemomentumGausFit[i]->Delete();
		}

	}
	if(!SavesieveCheck->IsZombie()){
		SavesieveCheck->Delete();
	}
	if(!h->IsZombie()){
		h->Delete();
	}

	f1->Write();
	f1->ls();
	f1->Close();
	cutdesc.close();
}


/// used for add the first excited states in the fitting
//Recognize the save patter
// save the rec hole to the folder
// first excited state
void SavePatternHole_P1(double momentumSigmaCut=3.0){
	//search all the holes in this col and save in the folder
	std::cout<<std::endl<<std::endl;
	std::cout<<"*******Save process start ........*******"<<std::endl;
	std::cout<<"	Searching for holes in col ("<<col<<")"<<std::endl;
	std::cout<<"	CAUTION :: USING First Excited States !!!!!"<<std::endl;

	TCanvas *SaveCheckCanvas=(TCanvas *) gROOT->GetListOfCanvases()->FindObject("SieveSaveCheck");
	if(!SaveCheckCanvas){
		SaveCheckCanvas =new TCanvas("SieveSaveCheck","SieveSaveCheck",1600,1080);
	}else{
		SaveCheckCanvas->Clear();
	}
	SaveCheckCanvas->Divide(1,2);
	SaveCheckCanvas->cd(1)->Divide(2,1);
	SaveCheckCanvas->cd(2)->Divide(row_count,1);


	TString workdir_temp=WorkDir+"/FirstMomCut/";

	// check the existance of the folder, if not create the folder
	if(!boost::filesystem::is_directory(workdir_temp.Data())){
		std::cout<<"Folder :"<< workdir_temp.Data()<<"  Does Not EXIST\n Trying to Create the folder"<<std::endl;
		if(!boost::filesystem::create_directories(workdir_temp.Data())){
			std::cout<<"ERROR:: cannot creat folde, please check the permission!!!!"<<std::endl;
			exit(-1);
		}
	}

	TString CutFileName = workdir_temp + RootFileName + CutSuf;
	TString TempString(Form(CutDescFileSufSieve.Data(), FoilID, col));
	TString PlotDir(RootFileName + Form(".hcut_R_%d_%d/", FoilID, col));
	TString CutDescName = workdir_temp + RootFileName + TempString;
	// prepare the filename and folder
	//TString CutDescName = WorkDir + RootFileName + CutDescFileSufDp;
	//TString CutFileName = WorkDir + RootFileName + CutSuf; // used for save the root file and cut file
	TFile *f1=new TFile(CutFileName,"UPDATE");
	assert(f1);

	std::fstream cutdesc(CutDescName, std::ios_base::out);
	assert(cutdesc.is_open());

	SaveCheckCanvas->cd(1)->cd(1);
	// attached the root file and plot the canvas used for check
	TChain *chain = (TChain *) gROOT->FindObject("T");
	if(!chain) std::cout<<"[ERROR] CAN NOT FIND TREE"<<std::endl;
	TH2F *h=(TH2F*)gROOT->FindObject("th_vs_ph");
	if(!h)
		h=new TH2F("th_vs_ph1","th_vs_ph1",1000,-0.03,0.03,1000,-0.045,0.045);

	TH2F *SavesieveCheck = new TH2F("SaveSieveCheck", "SaveSieveCheck",
			h->GetXaxis()->GetNbins(), h->GetXaxis()->GetXmin(),
			h->GetXaxis()->GetXmax(), h->GetYaxis()->GetNbins(),
			h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax());

	TString HRS("R");
	TString filename(chain->GetFile()->GetName());
	if (filename.Contains("RHRS")) {
	} else if (filename.Contains("LHRS")) {
		HRS = "L";
	}
	chain->Project(SavesieveCheck->GetName(),Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()),Form("%s",generalcut.Data()));
	//chain->Draw(Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()),Form("%s",generalcut.Data()),"zcol");
	SavesieveCheck->Draw("zcol");
	// search how many cut file exist in the root buffer
	for (int i = 0 ; i < row_min; i ++)cutdesc<<"fEvtHdr.fRun==0"<< std::endl;
	std::vector<int> sieveIDList;
	TH2F *sievehole[row_count];
	TH1F *sieveholemomentum[row_count];    // ground states momentum
	TF1  *sieveholemomentumGausFit[row_count];
	TF1  *sieveholemomentumGausFit_p1[row_count];
	TLatex * momentumInfor1[row_count];
	TCutG *cutg;
	for (int row_iter = row_min; row_iter<row_min+row_count;row_iter++){

		cutg=(TCutG *) gROOT->FindObject(Form("hcut_R_%d_%d_%d",FoilID,col,row_iter));
		if(cutg){
			sieveIDList.push_back(row_iter-row_min);
			sievehole[row_iter-row_min] = new TH2F(
					Form("hcut_R_%d_%d_%d_hh", FoilID, col, row_iter),
					Form("hcut_R_%d_%d_%d_hh", FoilID, col, row_iter),
					h->GetXaxis()->GetNbins(), h->GetXaxis()->GetXmin(),
					h->GetXaxis()->GetXmax(), h->GetYaxis()->GetNbins(),
					h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax());
			SaveCheckCanvas->cd(2)->cd(row_iter - row_min + 1);
//			SaveCheckCanvas->cd(2)->cd(row_iter - row_min + 1)->SetLogy();
			chain->Project(sievehole[row_iter-row_min]->GetName(), Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()),Form("%s && %s",cutg->GetName(),generalcut.Data()));

			sieveholemomentum[row_iter-row_min]=new TH1F(Form("hcut_R_%d_%d_%d_h_momentum", FoilID, col, row_iter),Form("hcut_R_%d_%d_%d_momentum", FoilID, col, row_iter),600,2.1,2.25);
			chain->Project(sieveholemomentum[row_iter-row_min]->GetName(),Form("%s.gold.p",HRS.Data()),Form("%s && %s",cutg->GetName(),generalcut.Data()));
			sieveholemomentumGausFit[row_iter-row_min]=new TF1(Form("1ststatesDpgaushcut_R_%d_%d_%d", FoilID, col, row_iter),"gaus",
					sieveholemomentum[row_iter-row_min]->GetXaxis()->GetBinCenter(sieveholemomentum[row_iter-row_min]->GetMaximumBin())-0.002,
					sieveholemomentum[row_iter-row_min]->GetXaxis()->GetBinCenter(sieveholemomentum[row_iter-row_min]->GetMaximumBin())+0.002);
			sieveholemomentumGausFit[row_iter-row_min]->SetParameter(1,sieveholemomentum[row_iter-row_min]->GetXaxis()->GetBinCenter(sieveholemomentum[row_iter-row_min]->GetMaximumBin()));
			sieveholemomentum[row_iter - row_min]->Fit(
					Form("1ststatesDpgaushcut_R_%d_%d_%d", FoilID, col,
							row_iter), "R", "ep",
					sieveholemomentumGausFit[row_iter-row_min]->GetXmin(),
					sieveholemomentumGausFit[row_iter-row_min]->GetXmax());

			sieveholemomentum[row_iter - row_min]->GetXaxis()->SetRangeUser(
					sieveholemomentum[row_iter - row_min]->GetXaxis()->GetBinCenter(
							sieveholemomentum[row_iter - row_min]->GetMaximumBin())
							- 0.009,
					sieveholemomentum[row_iter - row_min]->GetXaxis()->GetBinCenter(
							sieveholemomentum[row_iter - row_min]->GetMaximumBin())
							+ 0.004);
			sieveholemomentum[row_iter - row_min]->Draw();
			sieveholemomentumGausFit[row_iter-row_min]->Draw("same");

			auto groudpcenter=sieveholemomentumGausFit[row_iter-row_min]->GetParameter(1);
			auto groudpsigma=sieveholemomentumGausFit[row_iter-row_min]->GetParameter(2);

			momentumInfor1[row_iter - row_min] =
					new TLatex(groudpcenter + 2 * groudpsigma-0.007,
							sieveholemomentumGausFit[row_iter - row_min]->GetParameter(0),
							Form("P_{0}= %2.5f",
									groudpcenter));//, groudpsigma
			momentumInfor1[row_iter - row_min]->SetTextSize(0.055);
			momentumInfor1[row_iter - row_min]->SetTextAlign(12);
			momentumInfor1[row_iter - row_min]->SetTextColor(2);
			momentumInfor1[row_iter - row_min]->Draw("same");

			SaveCheckCanvas->cd(2)->cd(row_iter - row_min + 1)->Update();
			if(groudpsigma>0.0008)groudpsigma=0.0008;
			// plot the boundary of the cut
			TLine *leftboundary=new TLine(groudpcenter-momentumSigmaCut*groudpsigma,0,groudpcenter-momentumSigmaCut*groudpsigma,(0.9*SaveCheckCanvas->cd(2)->cd(row_iter - row_min + 1)->GetUymax()));
			leftboundary->SetLineColor(3);
			leftboundary->SetLineWidth(2);
			leftboundary->Draw("same");

			TLine *rightboundary=new TLine(groudpcenter+momentumSigmaCut*groudpsigma,0,groudpcenter+momentumSigmaCut*groudpsigma,(0.9*SaveCheckCanvas->cd(2)->cd(row_iter - row_min + 1)->GetUymax()));
			rightboundary->SetLineColor(3);
			rightboundary->SetLineWidth(2);
			rightboundary->Draw("same");

			// add the second excited state gaus fit
			TH1F *hSieveP1h=(TH1F *)sieveholemomentum[row_iter - row_min]->Clone("C12.p1");
//			hSieveP1h->GetXaxis()->SetRangeUser(sieveholemomentum[row_iter - row_min]->GetXaxis()->GetXmin()
//							,groudpcenter-groudpsigma*4);

			hSieveP1h->GetXaxis()->SetRangeUser(groudpcenter-0.006
							,groudpcenter-0.004);


			// get the bin center, and used this bin center as the fit center for the P1
			double_t p1_mean=hSieveP1h->GetXaxis()->GetBinCenter(hSieveP1h->GetMaximumBin());
			std::cout<<"===>"<<p1_mean<<std::endl;
			hSieveP1h->Delete();
			sieveholemomentumGausFit_p1[row_iter-row_min]=new TF1(Form("1ststatesDpgaushcut_R_p1_%d_%d_%d", FoilID, col, row_iter),"gaus",p1_mean-0.0009,
					p1_mean+0.0009);
			sieveholemomentum[row_iter - row_min]->Fit(
								Form("1ststatesDpgaushcut_R_p1_%d_%d_%d", FoilID, col,
										row_iter), "R", "ep",
										sieveholemomentumGausFit_p1[row_iter-row_min]->GetXmin(),
										sieveholemomentumGausFit_p1[row_iter-row_min]->GetXmax());
			sieveholemomentumGausFit_p1[row_iter-row_min]->Draw("same");
			double_t p1guasMeam=sieveholemomentumGausFit_p1[row_iter-row_min]->GetParameter(1);
			double_t p1guasSigma=sieveholemomentumGausFit_p1[row_iter-row_min]->GetParameter(2);


			TLine *p1leftboundary=new TLine(p1guasMeam-(momentumSigmaCut-1)*p1guasSigma,0,p1guasMeam-(momentumSigmaCut-1)*p1guasSigma,(0.9*SaveCheckCanvas->cd(2)->cd(row_iter - row_min + 1)->GetUymax()));
			p1leftboundary->SetLineColor(5);
			p1leftboundary->SetLineWidth(2);
			p1leftboundary->Draw("same");

			TLine *p1rightboundary=new TLine(p1guasMeam+(momentumSigmaCut-1)*p1guasSigma,0,p1guasMeam+(momentumSigmaCut-1)*p1guasSigma,(0.9*SaveCheckCanvas->cd(2)->cd(row_iter - row_min + 1)->GetUymax()));
			p1rightboundary->SetLineColor(5);
			p1rightboundary->SetLineWidth(2);
			p1rightboundary->Draw("same");


			SaveCheckCanvas->cd(1)->cd(1);
			cutg->SetName(Form("hcut_R_%d_%d_%d", FoilID, col, row_iter));
			cutg->SetVarX(Form("%s.gold.ph",HRS.Data()));
			cutg->SetVarY(Form("%s.gold.th",HRS.Data()));
			cutg->SetLineColor(kMagenta);
			cutg->SetLineWidth(2);
			cutg->Draw("PL same");
			SaveCheckCanvas->cd(1)->cd(1)->SetGridx(20);
			SaveCheckCanvas->cd(1)->cd(1)->SetGridy(20);
			//plot the momentum and apply cut on the momentum

			cutg->Write("", TObject::kOverwrite); // Overwrite old cut
			if(groudpsigma>0.0008)groudpsigma=0.0008;
			cutdesc << Form("hcut_R_%d_%d_%d", FoilID, col, row_iter) << " && ";
			cutdesc << (const char*)generalcut <<" && "
					<<Form("abs(%s.gold.p-%f)<%f*%f",HRS.Data(),p1guasMeam,momentumSigmaCut-1,p1guasSigma)
					<< std::endl;

			SaveCheckCanvas->cd(1)->cd(2);

			TH2F *sieveholetemp= new TH2F(
								Form("hcut_R_%d_%d_%d_hhg", FoilID, col, row_iter),
								Form("hcut_R_%d_%d_%d_hhg", FoilID, col, row_iter),
								h->GetXaxis()->GetNbins(), h->GetXaxis()->GetXmin(),
								h->GetXaxis()->GetXmax(), h->GetYaxis()->GetNbins(),
								h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax());
			sieveholetemp->Draw("same");
			chain->Project(sieveholetemp->GetName(),Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()),Form("%s && %s && abs(%s.gold.p-%f)<%f*%f",generalcut.Data(),cutg->GetName(),HRS.Data(),groudpcenter,momentumSigmaCut,groudpsigma));

			TH2F *sieveholetempp1= new TH2F(
											Form("hcut_R_%d_%d_%d_hhp1", FoilID, col, row_iter),
											Form("hcut_R_%d_%d_%d_hhp1", FoilID, col, row_iter),
											h->GetXaxis()->GetNbins(), h->GetXaxis()->GetXmin(),
											h->GetXaxis()->GetXmax(), h->GetYaxis()->GetNbins(),
											h->GetYaxis()->GetXmin(), h->GetYaxis()->GetXmax());
			chain->Project(sieveholetempp1->GetName(),Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()),Form("%s && %s && abs(%s.gold.p-%f)<%f*%f",generalcut.Data(),cutg->GetName(),HRS.Data(),p1guasMeam,momentumSigmaCut-1,p1guasSigma));

			TLatex *eventCountLable=new TLatex(sievehole[row_iter-row_min]->GetMean(1) + 0.005,sievehole[row_iter-row_min]->GetMean(2), Form("Entries(%d,%d): %2.0f-> %2.0f p1 :%2.0f",col,row_iter, (sievehole[row_iter-row_min]->GetEntries()),sieveholetemp->GetEntries(),sieveholetempp1->GetEntries()));
			eventCountLable->SetTextSize(0.03);
			eventCountLable->SetTextAlign(12);
			eventCountLable->SetTextColor(2);
			eventCountLable->Draw("same");
			SaveCheckCanvas->cd(1)->cd(2)->SetGridx(20);
			SaveCheckCanvas->cd(1)->cd(2)->SetGridy(20);
			//sievehole[row_iter-row_min]->Delete();
			cutg->Draw("same");
			SaveCheckCanvas->Update();
			sieveholetemp->Delete();
			sieveholetempp1->Delete();

		}else{
			//if the cut does not exist,then write the cut
			cutdesc << "fEvtHdr.fRun==0" << std::endl;
		}

	}

	for(int i = row_min+row_count; i < NSieveRow; i++)
		cutdesc << "fEvtHdr.fRun==0" << std::endl;
	SaveCheckCanvas->SetName(Form("CutProfcut_R_%d_%d",FoilID, col));
	SaveCheckCanvas->Write("", TObject::kOverwrite);
	SaveCheckCanvas->SaveAs(Form("%s/%s.hcut_R_%d_%d.jpg",workdir_temp.Data(),RootFileName.Data(),FoilID, col));

	for(auto i : sieveIDList){
//	for (unsigned int i = 0; i < row_count; i++) {
		if (!sievehole[i]->IsZombie()) {
			sievehole[i]->Delete();
		}
		if (!sieveholemomentum[i]->IsZombie()) {
			sieveholemomentum[i]->Delete();
		}
		if (!sieveholemomentumGausFit[i]->IsZombie()) {
			sieveholemomentumGausFit[i]->Delete();
			sieveholemomentumGausFit_p1[i]->Delete();
		}

	}
	if(!SavesieveCheck->IsZombie()){
		SavesieveCheck->Delete();
	}
	if(!h->IsZombie()){
		h->Delete();
	}

	f1->Write();
	f1->ls();
	f1->Close();
	cutdesc.close();
}


void CurrentStatus(){
	std::cout<<"***  's' save the current cut to file *****"<<std::endl;
	std::cout<<"***  '-' decrease one hole ID         *****"<<std::endl;
	std::cout<<"***  '+' increase one hole ID         *****"<<std::endl;
	std::cout<<"Currently working on ::"<<std::endl;
	std::cout<<"	hole col :"<<col<<std::endl;
	std::cout<<"    hole row :"<<row<< "  minID:"<<row_min<<"  count:"<<row_count<<std::endl;
	std::cout<<">>>Please click the center of the hole ("<<col<<" ,"<<row<<" ) ...."<<std::endl;
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
	TFile *tempfile=new TFile("temp.root","RECREATE");

	// link the root tree and check which HRS we are working on
	TChain *chain = (TChain *) gROOT->FindObject("T");
	TString HRS("R");
	TString filename(chain->GetFile()->GetName());
	if (filename.Contains("RHRS")) {
	} else if (filename.Contains("LHRS")) {
		HRS = "L";
	}

	// check the input keys
	if(event == kKeyPress) {

		std::cout<<"Key Pressed::"<<(gPad->GetEventX())<<std::endl;
		switch ((int) (gPad->GetEventX())) {
		case '+':
			row++;
			CurrentStatus();
			break;
		case '-':
			row--;
			CurrentStatus();
			break;
		case 's':
			std::cout << "Save Button Clicked" << std::endl;
			SavePatternHole();
			SavePatternHole(300000); // with out the groud momentum cut
			SavePatternHole_P1();
			break;
		case 'q':
			{std::cout << "Quit Button Clicked" << std::endl;
			gApplication->Clear();
			gApplication->Terminate();
			break;}
		default:
			std::cout<<(char)(gPad->GetEventX())<<std::endl;
		}
	}else
	if (event==kButton1Down) {
		TH2 *h = (TH2*) select;
		gPad->GetCanvas()->FeedbackMode(kTRUE);


		// if the button is clicked
		//Rec the sieve pattern
		// get the mouse click position in histogram
		double_t x = (gPad->PadtoX(gPad->AbsPixeltoX(gPad->GetEventX())));
		double_t y = (gPad->PadtoY(gPad->AbsPixeltoY(gPad->GetEventY())));

		// create new canvas
		TCanvas *SieveRecCanvas = (TCanvas*) gROOT->GetListOfCanvases()->FindObject("SieveRecCanvas");
		if(SieveRecCanvas){
			SieveRecCanvas->Clear();
			delete SieveRecCanvas->GetPrimitive("Projection");
		}else
			SieveRecCanvas = new TCanvas("SieveRecCanvas","Projection Canvas", 1000,1000);

			SieveRecCanvas->Divide(1,2);
			SieveRecCanvas->cd(2)->Divide(4,1);

			//get the hsitogram and start rec
			SieveRecCanvas->cd(2)->cd(1);

			TH2F *selectedSievehh=(TH2F *)  gROOT->FindObject("Sieve_Selected_th_ph");
			if(selectedSievehh){
				selectedSievehh->Clear();
			}else{
			selectedSievehh = new TH2F("Sieve_Selected_th_ph",
					"Sieve_Selected_th_ph",
					100,
					h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(),
					100,
					h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
			}

			chain->Project(selectedSievehh->GetName(),
					Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()),
					Form("sqrt((%s.gold.th-%f)^2+ (%s.gold.ph-%f)^2)<0.003 && %s ",
							HRS.Data(), y, HRS.Data(), x,generalcut.Data()));
			selectedSievehh->SetContour(10);
			selectedSievehh->GetXaxis()->SetTitle(Form("%s.gold.ph",HRS.Data()));
			selectedSievehh->GetYaxis()->SetTitle(Form("%s.gold.th",HRS.Data()));
			selectedSievehh->Draw("CONT LIST");

			SieveRecCanvas->Update(); // update the canvas to let the pattern buffer in root

			// extract the contour
			TObjArray *conts = (TObjArray*) gROOT->GetListOfSpecials()->FindObject(
					"contours");
			if (!conts)
				return;
			TList *lcontour1 = (TList*) conts->At(1);
			if (!lcontour1)
				return;
			TGraph *gc1 = (TGraph*) lcontour1->First();
			if (!gc1)
				return;
			if (gc1->GetN() < 10)
				return;

			//TODO need to change the name of
			TCutG *cutg = new TCutG(Form("hcut_R_%d_%d_%d", FoilID, col, row), gc1->GetN(), gc1->GetX(), gc1->GetY());
			cutg->SetLineColor(kRed);
			cutg->SetName(Form("hcut_R_%d_%d_%d", FoilID, col, row));
			cutg->SetVarX(Form("%s.gold.ph",HRS.Data()));
			cutg->SetVarY(Form("%s.gold.th",HRS.Data()));
			cutg->Draw("same");

			SieveRecCanvas->cd(2)->cd(2);
			auto projectx = selectedSievehh->ProjectionX();
			projectx->Draw();
			projectx->Fit("gaus");

			SieveRecCanvas->cd(2)->cd(3);
			auto projecty = selectedSievehh->ProjectionY();
			projecty->Draw();
			projecty->Fit("gaus");

			// plot the cut on the canvas
			SieveRecCanvas->cd(1);

			TH2F *patternCheck=(TH2F *)  gROOT->FindObject("Sieve_Pattern_Check");
			if(patternCheck){
				patternCheck->Clear();
			}
			patternCheck = new TH2F("Sieve_Pattern_Check",
					"Sieve_Pattern_Check", h->GetXaxis()->GetNbins(),
					h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax(),
					h->GetYaxis()->GetNbins(), h->GetYaxis()->GetXmin(),
					h->GetYaxis()->GetXmax());
			chain->Project(patternCheck->GetName(),
					Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()),Form("%s",generalcut.Data()));
			patternCheck->Draw("zcol");
			cutg->Draw("same");

			TLatex *label=new TLatex(selectedSievehh->GetMean(1),selectedSievehh->GetMean(2),Form("(%d %d)",col,row));
			label->SetTextSize(0.04);
			label->SetTextColor(2);
			label->Draw("same");

			row++;
			SieveRecCanvas->Update();

			SieveRecCanvas->cd(2)->cd(4);
			TH1F *sieveholemomentum=new TH1F(Form("hcut_R_%d_%d_%d_h_momentum_check", FoilID, col, row),Form("hcut_R_%d_%d_%d_momentum_check", FoilID, col, row),600,2.1,2.25);
			chain->Project(sieveholemomentum->GetName(),Form("%s.gold.p",HRS.Data()),Form("%s && %s",cutg->GetName(),generalcut.Data()));
			sieveholemomentum->GetXaxis()->SetRangeUser(
					sieveholemomentum->GetXaxis()->GetBinCenter(
							sieveholemomentum->GetMaximumBin())
							- 0.009,
					sieveholemomentum->GetXaxis()->GetBinCenter(
							sieveholemomentum->GetMaximumBin())
							+ 0.004);
			sieveholemomentum->Draw();
			SieveRecCanvas->Update();
			SieveRecCanvas->Write();
	}
	tempfile->Write();
	tempfile->Close();
}

