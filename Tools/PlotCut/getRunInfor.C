/*
 * getRunInfor.C
 *
 *  Created on: Jun 7, 2020
 *      Author: newdriver
 */

#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <sys/stat.h>
#include <string.h>
#include <map>
#include <iostream>

inline Bool_t IsFileExist (const std::string& name) {
	  struct stat buffer;
	  return (stat (name.c_str(), &buffer) == 0);
}

int getRunInfor(UInt_t runID,TString folder="/home/newdriver/Storage/Research/CRex_Experiment/RasterReplay/Replay/Result"){
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

	chain->GetTree();
	//THaEvent        *Event_Branch;
	ULong64_t fEvtHdr_fEvtTime;
	UInt_t fEvtHdr_fEvtNum;
	Int_t fEvtHdr_fEvtType;
	Int_t fEvtHdr_fEvtLen;
	Int_t fEvtHdr_fHelicity;
	Int_t fEvtHdr_fTargetPol;
	Int_t fEvtHdr_fRun;

	// project the branched
	chain->SetBranchAddress("fEvtHdr.fEvtTime", &fEvtHdr_fEvtTime);
	chain->SetBranchAddress("fEvtHdr.fEvtNum", &fEvtHdr_fEvtNum);
	chain->SetBranchAddress("fEvtHdr.fEvtType", &fEvtHdr_fEvtType);
	chain->SetBranchAddress("fEvtHdr.fEvtLen", &fEvtHdr_fEvtLen);
	chain->SetBranchAddress("fEvtHdr.fHelicity", &fEvtHdr_fHelicity);
	chain->SetBranchAddress("fEvtHdr.fTargetPol", &fEvtHdr_fTargetPol);
	chain->SetBranchAddress("fEvtHdr.fRun", &fEvtHdr_fRun);
	chain->SetMakeClass(1);

	// get the minimum and the maximum timestamp
	Long64_t nentries = chain->GetEntries();
	for (Long64_t jentry = 0; jentry < nentries; jentry++) {
		chain->GetEntry(jentry);
		std::cout << "time stamp::" << fEvtHdr_fEvtTime << std::endl;
	}



	return 0;
}
