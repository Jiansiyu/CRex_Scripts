/*
 * gemCluster.C
 *
 *  Created on: Jun 3, 2020
 *      Author: newdriver
 */

#include <TROOT.h>
#include <TString.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <map>
#include <TH1F.h>
#include <TCanvas.h>
#include <stdio.h>
#include <iostream>
#include <TLegend.h>

TTree *tree;
// extract the GEM cluster information
void gemCluster(){
	TString gemReplayedNamePattern="";
	TString gemPedestalNamePattern="";

	//load the information and get the plot


}

//extract the GEM pedestal informations
std::map<unsigned int,std::map<unsigned int, TH1F *> > gemPedestal(TString fname="/home/newdriver/Storage/Research/PRex_Experiment/PRex_Script/PRex_pedestal_generator/Pedestal/pedestal_LHRS_2056.root"){

	TFile *fileio= TFile::Open(fname.Data(),"READ");
	assert(fileio);
	fileio->GetObject("T",tree);
	fileio->Print();

	// check the existance of the trees and load it to the histogram
	std::map<unsigned int, TH1F *> prex_pedestalX;
	std::map<unsigned int, TH1F *> prex_pedestalY;

	for (auto gemPlane=0; gemPlane<8;  gemPlane++){
		TString pedestalPatternx(Form("prex_gems_X%d_pedestal_distribution",gemPlane));
		TString pedestalPatterny(Form("prex_gems_Y%d_pedestal_distribution",gemPlane));

		if( fileio->GetListOfKeys()->Contains(pedestalPatternx.Data())){
			//project the histogram
			prex_pedestalX[gemPlane]=(TH1F *)fileio->Get(pedestalPatternx.Data());
		}

		if (fileio->GetListOfKeys()->Contains(pedestalPatterny.Data())){
			prex_pedestalY[gemPlane]=(TH1F *) fileio->Get(pedestalPatterny.Data());
		}
	}

	//maybe need to change the name and draw type of the pedestal plot
	for (auto gemPlane=4; gemPlane<8; gemPlane++){
		// creare the plot canvas
		if(prex_pedestalX.find(gemPlane)!=prex_pedestalX.end()){
			prex_pedestalX[gemPlane]->GetXaxis()->SetRange(0,100);
			prex_pedestalX[gemPlane]->SetName(Form("PRex_X%d_Pedestal",gemPlane));
		}
		if(prex_pedestalY.find(gemPlane)!=prex_pedestalY.end()){
			prex_pedestalY[gemPlane]->GetXaxis()->SetRange(0,100);
			prex_pedestalY[gemPlane]->SetName(Form("PRex_Y%d_Pedestal",gemPlane));
		}
	}


	TCanvas *pedestalCanv= new TCanvas(Form("PRex SBS GEM Noise Level"),Form("PRex SBS GEM Noise Level"),1960,1080);
	pedestalCanv->Divide(3,2);
	pedestalCanv->Draw();
	// get the Plot
	for (auto gemPlane=prex_pedestalX.begin()->first+3; gemPlane<8; gemPlane++){
		if (prex_pedestalX.find(gemPlane)!=prex_pedestalX.end()){
			pedestalCanv->cd(gemPlane-2);
			prex_pedestalX[gemPlane]->Draw();
		}
		if(prex_pedestalY.find(gemPlane)!=prex_pedestalY.end()){
			pedestalCanv->cd(gemPlane+1);
			prex_pedestalY[gemPlane]->Draw();
		}
	}

	pedestalCanv->Update();

	std::map<unsigned int,std::map<unsigned int, TH1F *> >prex_pedestal;
	prex_pedestal[0]=prex_pedestalX;
	prex_pedestal[1]=prex_pedestalY;
	return prex_pedestal;

}

void gemSignalInfor(TString HRS="L"){
	TString rootFilename="/home/newdriver/Storage/Research/SBS/Analysis_Program/GEManalysis/prexLHRS_2139_cluster_50k_aver_highacc_v1.root";
	if (HRS.Contains("R")){
	rootFilename="/home/newdriver/Storage/Research/SBS/Analysis_Program/GEManalysis/prexRHRS_21363_cluster_50k_aver_highacc_v1.root";
	}
	std::cout<<"Working on ::"<<rootFilename.Data()<<std::endl;
	TFile *fileio= TFile::Open(rootFilename.Data(),"READ");
	assert(fileio);
	fileio->GetObject("T",tree);
	fileio->Print();

	// extact the informations
	std::map<unsigned int, TH1F *> prex_gemSignalX;
	std::map<unsigned int, TH1F *> prex_gemSignalY;

	auto gemPlane=4;
	auto gemPlaneEnd=7;

	if (HRS.Contains("R")){
	gemPlane=3;
	gemPlaneEnd=6;
	}
	for (; gemPlane<gemPlaneEnd;  gemPlane++){
		TString pedestalPatternx(Form("hClusterAdcAverageDistX_%d",gemPlane));
		TString pedestalPatterny(Form("hClusterAdcAverageDistY_%d",gemPlane));

		if( fileio->GetListOfKeys()->Contains(pedestalPatternx.Data())){
			//project the histogram
			prex_gemSignalX[gemPlane]=(TH1F *)fileio->Get(pedestalPatternx.Data());
		}

		if (fileio->GetListOfKeys()->Contains(pedestalPatterny.Data())){
			prex_gemSignalY[gemPlane]=(TH1F *) fileio->Get(pedestalPatterny.Data());
		}
	}

	TCanvas *prexGEMSignalCanv=new TCanvas("PRex GEM Signal Canv","PRex GEM Signal Canv",1960,1080);
	prexGEMSignalCanv->Divide(3,2);
	prexGEMSignalCanv->Draw();

	for (auto gemSignX_iter = prex_gemSignalX.begin(); gemSignX_iter!=prex_gemSignalX.end(); gemSignX_iter++){
		prexGEMSignalCanv->cd(gemSignX_iter->first-prex_gemSignalX.begin()->first+1);
		gemSignX_iter->second->SetTitle(Form("PRex_%sHRS_X%d",HRS.Data(),gemSignX_iter->first));
		gemSignX_iter->second->SetName(Form("PRex_%sHRS_X%d",HRS.Data(),gemSignX_iter->first));
		if (HRS.Contains("R")){
			gemSignX_iter->second->SetTitle(Form("PRex_%sHRS_X%d",HRS.Data(),gemSignX_iter->first+1));
			gemSignX_iter->second->SetName(Form("PRex_%sHRS_X%d",HRS.Data(),gemSignX_iter->first+1));
		}
		gemSignX_iter->second->GetXaxis()->SetRange(0,1500);
		gemSignX_iter->second->Draw();
	}

	for (auto gemSignY_iter = prex_gemSignalY.begin(); gemSignY_iter!=prex_gemSignalY.end(); gemSignY_iter++){
		prexGEMSignalCanv->cd(gemSignY_iter->first-prex_gemSignalY.begin()->first+4);
		gemSignY_iter->second->SetTitle(Form("PRex_%sHRS_Y%d",HRS.Data(),gemSignY_iter->first));
		gemSignY_iter->second->SetName(Form("PRex_%sHRS_Y%d",HRS.Data(),gemSignY_iter->first));
		if (HRS.Contains("R")){
			gemSignY_iter->second->SetTitle(Form("PRex_%sHRS_Y%d",HRS.Data(),gemSignY_iter->first+1));
			gemSignY_iter->second->SetName(Form("PRex_%sHRS_Y%d",HRS.Data(),gemSignY_iter->first+1));
		}
		gemSignY_iter->second->GetXaxis()->SetRange(0,1500);
		gemSignY_iter->second->Draw();
	}
	prexGEMSignalCanv->Update();
	if (rootFilename.Contains("RHRS")) {
		prexGEMSignalCanv->SaveAs("PRex_RHRS_21363.jpg");
	} else {
		prexGEMSignalCanv->SaveAs("PRex_LHRS_2139.jpg");
	}

}
void extractGEMSignalInfor(TString fname =
		"/home/newdriver/PRex/PRex_Data/GEMRootFile/prexLHRS_2141_00.root",
		TString PedestalFname =
				"/home/newdriver/Storage/Research/PRex_Experiment/PRex_Script/PRex_pedestal_generator/Pedestal/pedestal_LHRS_2056.root",
				Bool_t normalize=true) {

	// Check LHRS or RHRS
	TString HRS="LGEM.lgems";
	if (fname.Contains("LHRS")){
		std::cout<<"Working on LHRS"<<std::endl;
	}
	else
		{
		HRS="RGEM.rgems";
	}

	TFile *fileio=TFile::Open(fname.Data());
	assert(fileio);
	TTree *tree;
	fileio->GetObject("T",tree);

	if(tree->IsZombie()){
		std::cout<<"[Error]: can not find tree in the file !!!"<<std::endl;
	}else{
		std::cout<<"Total Entries in the file:"<< (tree->GetEntries())<<std::endl;
	}

	// get the HRS informations
	std::map<unsigned int, TH1F *>signalAverageADC_x;
	std::map<unsigned int, TH1F *>signalAverageADC_y;

	for (auto gemPlane=0; gemPlane<8;  gemPlane++){
		std::string checkStr(Form("%s.x%d.adc0",HRS.Data(),gemPlane));

		if(tree->GetListOfBranches()->Contains(checkStr.c_str()))
		{

			// initilized the tree
			if (signalAverageADC_x.find(gemPlane)==signalAverageADC_x.end()){
				signalAverageADC_x[gemPlane]=new TH1F(Form("GEM Signal Average ADC X%d",gemPlane),Form("GEM Signal Average ADC X%d",gemPlane),300,10,2000);
				signalAverageADC_y[gemPlane]=new TH1F(Form("GEM Signal Average ADC Y%d",gemPlane),Form("GEM Signal Average ADC Y%d",gemPlane),300,10,2000);
			}


			// fill the data in the tree
			std::string formulaPatternX(Form("Sum$(%s.x%d.strip.adc)/(%s.x%d.nhitstrips*3.0)",HRS.Data(),gemPlane,HRS.Data(),gemPlane));
			std::string formulaPatternY(Form("Sum$(%s.y%d.strip.adc)/(%s.y%d.nhitstrips*3.0)",HRS.Data(),gemPlane,HRS.Data(),gemPlane));

			// need to properly cut the data to select only the sigma cut larger than 2 ??

			tree->Project(signalAverageADC_x[gemPlane]->GetName(),formulaPatternX.c_str());
			tree->Project(signalAverageADC_y[gemPlane]->GetName(),formulaPatternY.c_str());

			std::cout<<"Project the data "<< gemPlane<<"   get::"<<signalAverageADC_x[gemPlane]->GetEntries()<<std::endl;
		}
	}
	// to plot the data, maybe need to normallized all the event to the same scale

	// The following code maybe different between LHRS and RHRS
	TCanvas *SignalAverSizeCanv=new TCanvas("GEM Signal Aver","GEM Signal Aver",1960,1080);
	SignalAverSizeCanv->Divide(3,2);
	SignalAverSizeCanv->Draw();
	for (auto gemPlane=4; gemPlane<8; gemPlane++){
		if(signalAverageADC_x.find(gemPlane)!=signalAverageADC_x.end()){
			SignalAverSizeCanv->cd(gemPlane-3);
			signalAverageADC_x[gemPlane]->Draw();
		}

		if (signalAverageADC_y.find(gemPlane)!=signalAverageADC_y.end()){
			SignalAverSizeCanv->cd(gemPlane);
			signalAverageADC_y[gemPlane]->Draw();
			//signalAverageADC_y[gemPlane]->Fit("landau");
		}
	}
	SignalAverSizeCanv->Update();

/*	TCanvas *gemSignelNoiseCanv=new TCanvas("GEM signal Noise Diagnose","GEM signal Noise Diagnose",1960,1080);
	gemSignelNoiseCanv->Divide(3,2);
	gemSignelNoiseCanv->Draw();

	//load the pedestal data and extract the pedetal distribution
	auto prex_pedestal=gemPedestal(PedestalFname.Data());
	auto prex_pedestalX=prex_pedestal[0];
	auto prex_pedestalY=prex_pedestal[1];

	// draw the plot
	for (auto gemPlane=4; gemPlane<7;  gemPlane++){
		gemSignelNoiseCanv->cd(gemPlane-3);
		if(signalAverageADC_x.find(gemPlane)!=signalAverageADC_x.end() && (prex_pedestalX.find(gemPlane)!=prex_pedestalX.end())){
			double norm=100;
			if (normalize){
				signalAverageADC_x[gemPlane]->Scale(norm/signalAverageADC_x[gemPlane]->Integral("width"));
				prex_pedestalX[gemPlane]->Scale(norm/prex_pedestalX[gemPlane]->Integral("width"));
			}
			signalAverageADC_x[gemPlane]->GetYaxis()->SetRangeUser(0,1);
			signalAverageADC_x[gemPlane]->Draw("hist");
			prex_pedestalX[gemPlane]->Draw("same hist");
		}


		gemSignelNoiseCanv->cd(gemPlane);
		if (signalAverageADC_y.find(gemPlane)!=signalAverageADC_y.end() && (prex_pedestalY.find(gemPlane)!=prex_pedestalY.end())){
			double norm=1;
			if (normalize){
				signalAverageADC_y[gemPlane]->Scale(norm/signalAverageADC_y[gemPlane]->Integral("width"));
				prex_pedestalY[gemPlane]->Scale(norm/prex_pedestalY[gemPlane]->Integral("width"));
			}
			signalAverageADC_y[gemPlane]->GetYaxis()->SetRangeUser(0,1);
			signalAverageADC_y[gemPlane]->Draw("hist");
			prex_pedestalY[gemPlane]->Draw("same hist");
		}

	}
	gemSignelNoiseCanv->Update();*/

}
