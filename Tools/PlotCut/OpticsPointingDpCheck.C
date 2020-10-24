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
#include "TVector2.h"

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
TString generalcutR="R.tr.n==1 && R.vdc.u1.nclust==1&& R.vdc.v1.nclust==1 && R.vdc.u2.nclust==1 && R.vdc.v2.nclust==1";// && R.gold.p > 2.14 && R.gold.p < 2.2  ";
TString generalcutL="L.tr.n==1 && L.vdc.u1.nclust==1&& L.vdc.v1.nclust==1 && L.vdc.u2.nclust==1 && L.vdc.v2.nclust==1";//  && L.gold.p > 2.14 && L.gold.p < 2.2";

inline Bool_t IsFileExist (const std::string& name) {
    return !gSystem->AccessPathName(name.c_str());
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

TF1 *SpectroCrystalFitDp_C12(TH1F*momentumSpectro,int fitPeak=4){

    if(momentumSpectro->GetEntries()<10) {
        std::cout<<"[ERROR]:: "<<"Fit Data empty"<<std::endl;
        exit(-1);
    }

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
	momentumSpectro->Fit("fCrystalMomentum","RQ0","ep",fCrystalMomentum->GetXmin(),fCrystalMomentum->GetXmax());
	fCrystalMomentum->GetParameters(fCrystalMomentumPar);
	// get the Dp seperation for the first and second, and then project to the third and the fouth to fit the second the third excited states


	if(fitPeak>3){
		double c2GausFitPar[3];
		double c2_fitRange_Min=fCrystalMomentum->GetParameter(1)-(fCrystalMomentum->GetParameter(1)-fCrystalMomentum->GetParameter(6))*7.65407/4.43982-2*fCrystalMomentum->GetParameter(7);
		double c2_fitRange_Max=fCrystalMomentum->GetParameter(1)-(fCrystalMomentum->GetParameter(1)-fCrystalMomentum->GetParameter(6))*7.65407/4.43982+2*fCrystalMomentum->GetParameter(7);
		//fit the peak with gaussion
		momentumSpectro->Fit("gaus","R0Q","ep",c2_fitRange_Min,c2_fitRange_Max);
		momentumSpectro->GetFunction("gaus")->GetParameters(c2GausFitPar);

		double c3GausFitPar[3];
		double c3_fitRange_Min=fCrystalMomentum->GetParameter(1)-(fCrystalMomentum->GetParameter(1)-fCrystalMomentum->GetParameter(6))*9.641/4.43982-2*fCrystalMomentum->GetParameter(7);
		double c3_fitRange_Max=fCrystalMomentum->GetParameter(1)-(fCrystalMomentum->GetParameter(1)-fCrystalMomentum->GetParameter(6))*9.641/4.43982+2*fCrystalMomentum->GetParameter(7);
		// get th peak and start the fit
		momentumSpectro->Fit("gaus","R0Q","ep",c3_fitRange_Min,c3_fitRange_Max);
		momentumSpectro->GetFunction("gaus")->GetParameters(c3GausFitPar);

		double_t fCrystalGausMomentumPar[16];
		std::copy(fCrystalMomentumPar,fCrystalMomentumPar+10,fCrystalGausMomentumPar);
		std::copy(c2GausFitPar,c2GausFitPar+3,fCrystalGausMomentumPar+10);
		std::copy(c3GausFitPar,c3GausFitPar+3,fCrystalGausMomentumPar+13);

		TF1 *fCrystalGuasMomentum=new TF1("fCrystalGuasMomentum","crystalball(0)+crystalball(5)+gaus(10)+gaus(13)",c3_fitRange_Min,fgroundCrystalball->GetXmax());
		fCrystalGuasMomentum->SetParameters(fCrystalGausMomentumPar);
		momentumSpectro->Fit("fCrystalGuasMomentum","R0Q","ep",fCrystalGuasMomentum->GetXmin(),fCrystalGuasMomentum->GetXmax());
		return fCrystalGuasMomentum;

	}else if (fitPeak==3) {
		double c2GausFitPar[3];
		double c2_fitRange_Min=fCrystalMomentum->GetParameter(1)-(fCrystalMomentum->GetParameter(1)-fCrystalMomentum->GetParameter(6))*7.65407/4.43982-3*fCrystalMomentum->GetParameter(7);
		double c2_fitRange_Max=fCrystalMomentum->GetParameter(1)-(fCrystalMomentum->GetParameter(1)-fCrystalMomentum->GetParameter(6))*7.65407/4.43982+3*fCrystalMomentum->GetParameter(7);
		//fit the peak with gaussion
		momentumSpectro->Fit("gaus","","",c2_fitRange_Min,c2_fitRange_Max);
		momentumSpectro->GetFunction("gaus")->GetParameters(c2GausFitPar);

		double_t fCrystalGausMomentumPar[13];
		std::copy(fCrystalMomentumPar,fCrystalMomentumPar+10,fCrystalGausMomentumPar);
		std::copy(c2GausFitPar,c2GausFitPar+3,fCrystalGausMomentumPar+10);

		TF1 *fCrystalGuasMomentum=new TF1("fCrystalGuasMomentum","crystalball(0)+crystalball(5)+gaus(10)",c2_fitRange_Min,fgroundCrystalball->GetXmax());
		fCrystalGuasMomentum->SetParameters(fCrystalGausMomentumPar);
		momentumSpectro->Fit("fCrystalGuasMomentum","","",fCrystalGuasMomentum->GetXmin(),fCrystalGuasMomentum->GetXmax());
		return fCrystalGuasMomentum;
	}



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


TChain *LoadrootFile(UInt_t runID,TString folder="/home/newdriver/pyQuant/crex_replayed/run20201016"){
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

Int_t OpticsPointingDpCheck(UInt_t runID,TString folder="/home/newdriver/pyQuant/crex_replayed/run20201016") {
	// prepare the data
	TChain *chain=LoadrootFile(runID,folder);
	TString rootDir(folder.Data());
	TString HRS="R";
	if(runID<20000) HRS="L";

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



/*double getC12TheoreticalDp(int runID, int excitedState){
	// get the theoretical Dp value in given runID and the exxited Stetates
	std::map<int, std::map<int, double >> DpTheoreticalList;
	DpTheoreticalList[21642][0]=0.0139829;
	DpTheoreticalList[21641][0]=0.00389789;
	DpTheoreticalList[21626][0]=-0.00628822;
	DpTheoreticalList[21632][0]=-0.0156923;

//h2o o
//	-0.00576159,-0.0156143,0.00425469,0.00424658,
	DpTheoreticalList[21740][0]=-0.00576159;
	DpTheoreticalList[21762][0]=-0.0156143;
	DpTheoreticalList[21789][0]=0.00425469;
	DpTheoreticalList[21790][0]=0.00424658;

//-0.0132058,-0.0229847,-0.00326446,-0.00327241,
	DpTheoreticalList[21740][1]=-0.0132058;
	DpTheoreticalList[21762][1]=-0.0229847;
	DpTheoreticalList[21789][1]=-0.00326446;
	DpTheoreticalList[21790][1]=-0.00327241;


	//4.43982 MeV
	DpTheoreticalList[21642][1]=0.0119146;
	DpTheoreticalList[21641][1]=0.00185016;
	DpTheoreticalList[21626][1]=-0.0083155;
	DpTheoreticalList[21632][1]=-0.0176998;

	//7.65407 MeV
	//0.0104159,0.000366456,-0.00978438,-0.0191544,
	DpTheoreticalList[21642][2]= 0.0104159;
	DpTheoreticalList[21641][2]= 0.000366456;
	DpTheoreticalList[21626][2]=-0.00978438;
	DpTheoreticalList[21632][2]=-0.0191544;

	//9.641 MeV
	//0.00948953,-0.000550671,-0.0106923,-0.0200535,
	DpTheoreticalList[21642][3]=0.00948953;
	DpTheoreticalList[21641][3]=-0.000550671;
	DpTheoreticalList[21626][3]=-0.0106923;
	DpTheoreticalList[21632][3]=-0.0200535;

	return DpTheoreticalList[runID][excitedState];
}*/

UInt_t getRunID(TChain *chain){
    return (UInt_t) chain->GetMaximum("fEvtHdr.fRun");
}

double getCentralP(TChain *chain, Bool_t drawFlag=false){

    // get the run Number
    int runID=(int)chain->GetMaximum("fEvtHdr.fRun");

    TString HRS="R";
    if(runID<20000){
        HRS="L";
    }
    //TODO need to check whether this is prex/crex experiment

    double CentralP;

    if(HRS=="L"){
        TH1F *HallProbHH = new TH1F("HallLProb", "HallLProb", 1000, -1, 0);
        chain->Project(HallProbHH->GetName(),
                                         "HacL_D_LS450_FLD_DATA");
        if (HallProbHH->GetEntries() != 0) {
            CentralP = std::abs(
                    (HallProbHH->GetMean()) * 0.95282 / 0.33930);
            std::cout << "CentralMomentum is ::" << (CentralP) << std::endl;
        }
        HallProbHH->Delete();
    }else{
        //HacR_D1_NMR_SIG
        TH1F *HallR_NMR = new TH1F("HallR_NMR", "HallR_NMR", 1000, 0.7, 0.9);
        chain->Project(HallR_NMR->GetName(), "HacR_D1_NMR_SIG");
        if (HallR_NMR->GetEntries()) {
            double Mag = HallR_NMR->GetMean();
            CentralP = 2.702 * (Mag) - 1.6e-03 * (Mag) * (Mag) * (Mag);
            std::cout << "CentralMomentum is (RHRS) from NMR::" << CentralP
                      << std::endl;
        } else {
            std::cout << "\033[1;33m [Warning]\033[0m Missing HallR_NMR:"
                      << std::endl;
        }
    }
    return CentralP;
}


// read the beam E from the file
double_t getBeamE(int runID,TChain *chain,TString beamEfname="/home/newdriver/Learning/GeneralScripts/halog/beamE.txt"){
    std::map<int, double_t> beamE;
    beamE[21739]=2.1763077;
    beamE[21740]=2.1763047;
    beamE[21789]=2.1762745;
    beamE[21790]=2.1762517;
    beamE[21762]=2.1763254;

    beamE[2566]=2.175918588;
    beamE[2565]=2.175984498;
    beamE[2550]=2.17560073;
    beamE[2556]=2.1762867;
    beamE[2674]=2.1763062;
    beamE[2697]=2.1763254;
    beamE[2726]=2.1762729;

    //TODO check the root file, if it contains the hallp information, need to use this value
    {
        double beamERangeMin=2100;
        double beamERangeMax=2200;
        chain->GetListOfBranches()->Contains("HALLA_p");
        TH1F *HallEpicsBeamE = new TH1F("HallEpicsBeamE", "HallEpicsBeamE", 1000, beamERangeMin, beamERangeMax);
        chain->Project(HallEpicsBeamE->GetName(),"HALLA_p");
        double epicsBeamE=HallEpicsBeamE->GetMean();
        if ((epicsBeamE > 2170)&&(epicsBeamE < 2180)){
            std::cout<<"\033[1;32m [Infor]\033[0m Read in the Beam E in ROOT file(with correction): "<<epicsBeamE*953.4/951.1<<std::endl;
            return  epicsBeamE/1000.0*953.4/951.1;
        }
    }

    //read in the beamE information and parser
    if ((!beamEfname.IsNull()) && IsFileExist(beamEfname.Data())){
        std::cout<<"\033[1;32m [Infor]\033[0m Read in the Beam E file: "<<beamEfname.Data()<<std::endl;
        std::ifstream infile(beamEfname.Data());
        int runID_temp;
        float beamE_temp;
        while (infile >> runID_temp >> beamE_temp){
            //std::cout<<"runID:"<<runID_temp<<"   ->   "<<beamE_temp<<std::endl;
            beamE[runID_temp]=beamE_temp/1000.0;
        }
    }else{
        std::cout<<"\033[1;33m [Warning]\033[0m can not find file "<<beamEfname.Data()<<" Skip the beamE file!!!"<<std::endl;

    }

    if(beamE.find(runID)!=beamE.end()){
        std::cout<<"\033[1;32m [Infor]\033[0m Read in the Beam E file(with correction): "<<beamE[runID]*953.4/951.1<<std::endl;
        return beamE[runID]*953.4/951.1;
    }else{
        std::cout<<"\033[1;31m [CAUTION]\033[0m Can not find the Beam E for run"<<runID<<" Using default value!!!["<<__func__<<"("<<__LINE__<<")]"<<std::endl;
        return 2.17568;
    }
}

///
/// \param BeamE  Beam Energy in GeV
/// \param scatteredAngle Scattered Angle in Degree, set to 0 to use default value
/// \return
std::map<TString,double>  getWaterScatteredP(double BeamE,double  scatteredAngle){
    std::map<TString,double> res;
    Double_t amu = 931.494028 * 1.e-3;  // amu to GeV
    Double_t H2O = 18 * amu;
    Double_t TargetH = amu;
    double_t TargetO = 16 * amu;
    Double_t mass_tg = H2O;
    double beamE = BeamE;         // in GeV

    double HydMom=0.0;
    double OxyMom=0.0;

    double ScatteredAngleRad=scatteredAngle*TMath::Pi()/180.0;


    HydMom=beamE/(1+2*beamE*TMath::Sin(ScatteredAngleRad/2)*TMath::Sin(ScatteredAngleRad/2)/TargetH);
    OxyMom=beamE/(1+2*beamE*TMath::Sin(ScatteredAngleRad/2)*TMath::Sin(ScatteredAngleRad/2)/TargetO);

    res["H"]=HydMom;
    res["O"]=OxyMom;

    std::cout<<"Beam E:"<<beamE<<"  Angle:"<<scatteredAngle<<"   ep:"<<OxyMom<<std::endl;
    return  res;
}



TVector2 getBPM(UInt_t runID,TString csvfname="/home/newdriver/Learning/GeneralScripts/OptCheck/PRex/bpm_on_targ.csv"){
    if (gSystem->AccessPathName(csvfname.Data())){
        std::cout<<"\033[1;33m [Warning]\033[0m Missing csv file::"<<csvfname.Data()<<std::endl;
        exit(-1);
    }

    std::map<UInt_t,TVector2> bpmList;

    std::ifstream csvStream(csvfname.Data());
    std::string line;
    while (std::getline(csvStream,line)){
        std::istringstream s(line);
        std::string field;

        getline(s,field,',');
        TString title=field;
        if(title.Contains("run")) continue;

        UInt_t runs=std::stoi(field.c_str());

        getline(s,field,',');
        double_t bpmx=std::stof(field.c_str());

        getline(s,field,',');
        double_t bpmy=std::stof(field.c_str());

        TVector2 vec(bpmx,bpmy);

        bpmList[runs]=vec;
    }
    if(bpmList.find(runID) != bpmList.end()){
        return bpmList[runID];
    } else{
        std::cout<<"\033[1;33m [Warning]\033[0m Missing csv file::"<<csvfname.Data()<<std::endl;
        exit(-1);
    }
}

///
/// \param pos
/// \param chain
/// \param plotCanv
/// \return
TCutG *getSieveCut(TVector2 pos,TChain *chain, Bool_t plotCanv= false){
    // need to get the pos twice inoder to get the central sieve hole position correctly Maybe....
    double_t x=pos.X();
    double_t y=pos.Y();

    UInt_t runID=getRunID(chain);
    TString HRS="L";
    if (runID>20000) HRS="R";

    TH2F *PreselectedSievehh = (TH2F *) gROOT->FindObject("pre_Sieve_Selected_th_ph");
    if (PreselectedSievehh) {
        PreselectedSievehh->Clear();
    }
    PreselectedSievehh = new TH2F("pre_Sieve_Selected_th_ph", "pre_Sieve_Selected_th_ph",
                               200, -0.045,0.045, 200,
                               -0.045,0.045);
    chain->Project(PreselectedSievehh->GetName(),Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()),
                   Form("sqrt((%s.gold.th-%f)^2+ (%s.gold.ph-%f)^2)<0.003 ",
                        HRS.Data(), y, HRS.Data(), x));

    auto preSelectedSieveXh=PreselectedSievehh->ProjectionX();
    auto preSelectedSieveYh=PreselectedSievehh->ProjectionY();

    x=preSelectedSieveXh->GetXaxis()->GetBinCenter(preSelectedSieveXh->GetMaximumBin());
    y=preSelectedSieveYh->GetYaxis()->GetBinCenter(preSelectedSieveYh->GetMaximumBin());

    // update the central of the sieve hole and re-select the counte
    TH2F *selectedSievehh = (TH2F *) gROOT->FindObject("Sieve_Selected_th_ph");
    if (selectedSievehh) {
        selectedSievehh->Clear();
    }
    selectedSievehh = new TH2F("Sieve_Selected_th_ph", "Sieve_Selected_th_ph",
                               200, PreselectedSievehh->GetXaxis()->GetXmin(), PreselectedSievehh->GetXaxis()->GetXmax(), 200,
                               PreselectedSievehh->GetYaxis()->GetXmin(), PreselectedSievehh->GetYaxis()->GetXmax());
    chain->Project(selectedSievehh->GetName(),
                   Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()),
                   Form("sqrt((%s.gold.th-%f)^2+ (%s.gold.ph-%f)^2)<0.003",
                        HRS.Data(), y, HRS.Data(), x));
    selectedSievehh->SetContour(10);
    selectedSievehh->GetXaxis()->SetTitle(Form("%s.gold.ph", HRS.Data()));
    selectedSievehh->GetYaxis()->SetTitle(Form("%s.gold.th", HRS.Data()));

    // extract the contour
    TObjArray *conts = (TObjArray*) gROOT->GetListOfSpecials()->FindObject(
            "contours");
    TCutG *cut;
    if (!conts)
        return cut;
    TList *lcontour1 = (TList*) conts->At(0);
    if (!lcontour1)
        return cut;
    TGraph *gc1 = (TGraph*) lcontour1->First();
    if (!gc1)
        return cut;
    if (gc1->GetN() < 10)
        return cut;
    //TODO need to change the name of
    TCutG *cutg = new TCutG(Form("hcut_R_%ld",random()),
                            gc1->GetN(), gc1->GetX(), gc1->GetY());
    cutg->SetLineWidth(2);
    cutg->SetLineColor(kRed);
    cutg->SetVarX(Form("%s.gold.ph", HRS.Data()));
    cutg->SetVarY(Form("%s.gold.th", HRS.Data()));
    return  cutg;
}

void DynamicCanvas(){
	gStyle->SetOptStat(0);
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


    //Load ALl the root files
    std::map<int,TChain *> chainArray;
    std::map<int,int> ScanrunListArray;
    if (HRS == 'L') {
        std::cout<<"Working on LHRS"<<std::endl;
        chainArray[-2] = LoadrootFile(2566);
        chainArray[-1] = LoadrootFile(2565);
        chainArray[0] = LoadrootFile(2550);
        chainArray[1] = LoadrootFile(2556);
		chainArray[1000] = LoadrootFile(2671);
//        chainArray[999] = LoadrootFile(2671);
        chainArray[1001] = LoadrootFile(2709);

    }else {
		std::cout<<"Working on RHRS"<<std::endl;
		chainArray[-2] = LoadrootFile(21642);
		chainArray[-1] = LoadrootFile(21641);
		chainArray[0] = LoadrootFile(21626);
		chainArray[1] = LoadrootFile(21632);
//		chainArray[999] = LoadrootFile(22119);
		chainArray[1000] = LoadrootFile(22234);
		chainArray[1001] = LoadrootFile(22248);
		chainArray[1001] = LoadrootFile(22256);
	}

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
	SieveRecCanvas->cd(2)->Divide(chainArray.size(), 1);
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

	std::map<int,TH2F *>SieveThetaPhiList;

	int canvas_counter_temp=1;
	for (auto item = chainArray.begin(); item!= chainArray.end();item++){
		if((item->first)<10){
			SieveThetaPhiList[item->first]=new TH2F(Form("th_vs_ph_Dp%1d",item->first),Form("th_vs_ph_Dp%1d",item->first),1000,-0.03,0.03,1000,-0.045,0.045);
			item->second->Project(SieveThetaPhiList[item->first]->GetName(),Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()));
			SieveRecCanvas->cd(2)->cd(canvas_counter_temp);
			SieveThetaPhiList[item->first]->Draw("zcol");
//			getSieveCut(TVector2(0.0,y),item->second)->Draw("same");
			cutg->Draw("same");
		}else{
			SieveThetaPhiList[item->first]=new TH2F(Form("th_vs_ph_H20%1d",item->first),Form("th_vs_ph_H20%1d",item->first),1500,-0.03,0.03,1000,-0.045,0.045);
			item->second->Project(SieveThetaPhiList[item->first]->GetName(),Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()));
			SieveRecCanvas->cd(2)->cd(canvas_counter_temp);
			SieveThetaPhiList[item->first]->Draw("zcol");
//            getSieveCut(TVector2(x,y),item->second)->Draw("same");
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
			OptDpArrayH[item->first]=new TH1F(Form("Dp_hist:%d",item->first),Form("Dp_hist:%d",item->first),1000,-0.03,0.03);
			OptDpArrayH[item->first]->GetYaxis()->SetRange(0,10000);

			item->second->Project(OptDpArrayH[item->first]->GetName(),Form("%s.gold.dp",HRS.Data()),
			Form("%s && %s", generalcut.Data(),cutg->GetName()));
		}else{
			OptDpArrayH[item->first]=new TH1F(Form("Dp_hist:H2O%d",item->first),Form("Dp_hist:H_{2}O%d",item->first),1000,-0.03,0.03);
			item->second->Project(OptDpArrayH[item->first]->GetName(),Form("%s.gold.dp",HRS.Data()),Form("%s",cutg->GetName()));
//            OptDpArrayH[item->first]->Scale()
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


	//old value
    if (HRS == "L") { 	//0.0164994,0.00649027,-0.00405851,-0.00915162,
		DpTheoreticalMap[-2] = 0.0164994;
		DpTheoreticalMap[-1] = 0.00649027;
		DpTheoreticalMap[0] = -0.00405851;
		DpTheoreticalMap[1] = -0.00915162;
	}else {
		DpTheoreticalMap[-2] = 0.0164135;
		DpTheoreticalMap[-1] = 0.00630493;
		DpTheoreticalMap[0] = -0.00390594;
		DpTheoreticalMap[1] = -0.0134213;
	}

//	if (HRS == "L") { //0.0193804,0.00934287,-0.00123581,-0.00634336
//		DpTheoreticalMap[-2] = 0.0193804;
//		DpTheoreticalMap[-1] = 0.00934287;
//		DpTheoreticalMap[0] = -0.00123581;
//		DpTheoreticalMap[1] = -0.00634336;
//	}else {//0.0164135,0.00630493,-0.00390594,-0.0134213
//		DpTheoreticalMap[-2] = 0.0164135;
//		DpTheoreticalMap[-1] = 0.00630493;
//		DpTheoreticalMap[0] = -0.00390594;
//		DpTheoreticalMap[1] = -0.0134213;
//	}

	std::map<int, double >CentralPArray;
	for (auto item = OptDpArrayH.begin(); item!= OptDpArrayH.end();item++){
//		if((item->first)>10) continue;
		if((item->first)<10){
			item->second->SetLineColor(8+item->first);
			auto bpmVec=getBPM(getRunID(chainArray[item->first]));
			auto beamE=getBeamE(getRunID(chainArray[item->first]),chainArray[item->first]);
			legend->AddEntry((item->second),Form("C_{12} Dp:%2d%% scan run_%d  BeamPos(%f,%f), BeamE %fGeV",item->first, getRunID(chainArray[item->first]),bpmVec.X(),bpmVec.Y(),beamE));
		}else{
			item->second->SetLineColor(kRed);
            auto bpmVec=getBPM(getRunID(chainArray[item->first]));
            auto beamE=getBeamE(getRunID(chainArray[item->first]),chainArray[item->first]);
			legend->AddEntry((item->second),Form("H_{2}O Dp%d  run_%d  BeamPos(%f,%f), BeamE %fGeV",item->first-1000,getRunID(chainArray[item->first]),bpmVec.X(),bpmVec.Y(),beamE));
		}

		item->second->GetYaxis()->SetRangeUser(0,10000);
		item->second->SetLineWidth(2);

		if(item==OptDpArrayH.begin()){
			item->second->Draw();
		}else{
			item->second->Draw("same");
		}
		if((item->first)<10){
			fitFunctionsList[item->first]=SpectroCrystalFitDp_C12(item->second);
		}else{
			fitFunctionsList[item->first]=SpectroCrystalFitDp_C12(item->second,2);
		}
		fitFunctionsList[item->first]->SetLineColor(42);
		fitFunctionsList[item->first]->Draw("same");
		const int NfitPars= fitFunctionsList[item->first]->GetNpar();
		FitPars[item->first]=new double[NfitPars];
		fitFunctionsList[item->first]->GetParameters(FitPars[item->first]);


		// plot the central P informations
		if(HRS=="L"){
			TH1F *HallProbHH = new TH1F("HallLProb", "HallLProb", 1000, -1, 0);
			chainArray[item->first]->Project(HallProbHH->GetName(),
					"HacL_D_LS450_FLD_DATA", generalcut.Data());
			if (HallProbHH->GetEntries() != 0) {
				double CentralP = std::abs(
						(HallProbHH->GetMean()) * 0.95282 / 0.33930);
				std::cout << "CentralMomentum is ::" << (CentralP) << std::endl;
				TLatex *peakInfor = new TLatex(
						FitPars[item->first][1] - 5 * FitPars[item->first][2],
						FitPars[item->first][0],
						Form("P_{0}=%1.4fGeV#pm%1.3fMeV", CentralP,std::abs(
								(HallProbHH->GetRMS()) * 0.95282*1000.0 / 0.33930)));
				peakInfor->SetLineWidth(2);
				peakInfor->SetTextSize(0.02);
//				peakInfor->Draw("same");
				CentralPArray[item->first]=CentralP;
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
				std::cout << "CentralMomentum is (RHRS) from NMR::" << CentralP
						<< std::endl;

				std::cout << "CentralMomentum is ::" << (CentralP) << std::endl;
				Mag=HallR_NMR->GetRMS();
				double nmrperror=2.702 * (Mag) - 1.6e-03 * (Mag) * (Mag) * (Mag);
				TLatex *peakInfor = new TLatex(
						FitPars[item->first][1] - 5 * FitPars[item->first][2],
						FitPars[item->first][0],
						Form("P_{0}=%1.4fGeV#pn%1.3f", CentralP,nmrperror));
				peakInfor->SetLineWidth(2);
				peakInfor->SetTextSize(0.02);
//				peakInfor->Draw("same");
				CentralPArray[item->first]=CentralP;

			} else {
				std::cout << "\033[1;33m [Warning]\033[0m Missing HallR_NMR:"
						<< std::endl;
			}
		}
	}

	for (auto item = CentralPArray.begin(); item!=CentralPArray.end(); item++){
		std::cout<<" Central P ::"<<item->first<<"    -> "<<item->second<<std::endl;
	}

    DpCanvas->cd();
	// plot the bias at the carbon
    {
        for(auto iter = fitFunctionsList.begin();iter!=fitFunctionsList.end();iter++){
            int scanID=iter->first;
            if(scanID < 10){
                if(DpTheoreticalMap.find(scanID)!=DpTheoreticalMap.end()){
                    //plot the bias on the canv
                    std::cout<<"expected"<<DpTheoreticalMap[scanID]<<std::endl;
                    TLine *line=new TLine(DpTheoreticalMap[scanID],0,DpTheoreticalMap[scanID],fitFunctionsList[scanID]->GetParameter(0));
                    line->SetLineWidth(2);
                    line->Draw("same");
                    //print the bias
                    double theo=DpTheoreticalMap[scanID]*1000;
                    double meas=fitFunctionsList[scanID]->GetParameter(1)*1000;
                    TLatex *txt=new TLatex(fitFunctionsList[scanID]->GetParameter(1),fitFunctionsList[scanID]->GetParameter(0),Form("Bias:%1.3f x 10^{-3}",theo-meas));
                    txt->SetTextSize(0.02);
                    txt->Draw("same");
                }
            } else{
                // get the theoretical predicted value for the water runs
                UInt_t runID=getRunID(chainArray[scanID]);
                double beamE=getBeamE(runID,chainArray[scanID]);
                double centralP=getCentralP(chainArray[scanID]);
                double scatteredAngle=4.7655;
                if(runID > 20000) scatteredAngle=4.7717;
                auto scatteredEmap=getWaterScatteredP(beamE,scatteredAngle);

                double HydroE=0.0;
                double  OxyE=0.0;
                if (scatteredEmap.find("O")!=scatteredEmap.end()) OxyE=scatteredEmap["O"];
                if (scatteredEmap.find("H")!=scatteredEmap.end()) HydroE=scatteredEmap["H"];

                double  OxyEdp=(OxyE-centralP)/centralP;
                double  HydroEdp=(HydroE-centralP)/centralP;

                std::cout<<"Expect Position on Dp Plot::"<<OxyEdp<<","<<fitFunctionsList[scanID]->GetParameter(0)<<"  Central P"<<centralP<<std::endl;
                TLine *line=new TLine(OxyEdp,0,OxyEdp,fitFunctionsList[scanID]->GetParameter(0));
                line->SetLineWidth(2);
                line->Draw("same");

                double theo=OxyEdp;
                double meas=fitFunctionsList[scanID]->GetParameter(1);
                TLatex *txt=new TLatex(theo,fitFunctionsList[scanID]->GetParameter(0),Form("Bias:%f x 10^{-3}",1000*(theo-meas)));
                txt->SetTextSize(0.02);
                txt->SetTextColor(9);
                txt->Draw("same");

                TLine *lineH=new TLine(HydroEdp,0,HydroEdp,fitFunctionsList[scanID]->GetParameter(5));
                lineH->SetLineWidth(2);
                lineH->Draw("same");

                double  theoH=HydroEdp;
                double measH=fitFunctionsList[scanID]->GetParameter(6);
                TLatex *txth=new TLatex(theoH,fitFunctionsList[scanID]->GetParameter(5),Form("bias:%f x 10^{-3}",1000*(theoH-measH)));
                txth->SetTextSize(0.02);
                txth->SetTextColor(9);
                txth->Draw("same");
            }
        }
    }
    DpCanvas->Update();

//	for(auto item = OptDpArrayH.begin(); item!= OptDpArrayH.end();item++){
//		int i = item->first;
//		if(FitPars.find(i)!=FitPars.end()){
//			 TArrow *ar5 = new TArrow(FitPars[i][1],FitPars[i][0]/2,FitPars[i][6],FitPars[i][0]/2,15,"<|>");
//			ar5->SetAngle(60);
//			ar5->SetLineWidth(2);
//			ar5->SetLineColor(4);
//			ar5->SetFillStyle(3008);
//			ar5->SetFillColor(2);
//			ar5->Draw("same");
//
//			TLatex *txt=new TLatex(FitPars[i][1]*0.5+FitPars[i][6]*0.5,FitPars[i][0]*0.5,Form("%1.3fMeV#pm%1.3f",CentralPArray[i]*1000.0*(FitPars[i][1]-FitPars[i][6]),CentralPArray[i]*1000.0*TMath::Sqrt(fitFunctionsList[i]->GetParError(1)*fitFunctionsList[i]->GetParError(1)+fitFunctionsList[i]->GetParError(6)*fitFunctionsList[i]->GetParError(6))));
//			txt->SetLineWidth(1);
//			txt->SetTextSize(0.02);
//			txt->Draw("same");
//		}
//	}



//	for (int i =-2; i <1;i++)
//	{
//		   TArrow *ar5 = new TArrow(FitPars[i][1],FitPars[i][0],FitPars[i+1][1],FitPars[i][0],15,"<|>");
//		   ar5->SetAngle(60);
//		   ar5->SetLineWidth(4);
//		   ar5->SetLineColor(4);
//		   ar5->SetFillStyle(3008);
//		   ar5->SetFillColor(2);
//		   ar5->Draw("same");
//
//		   TLatex *txt=new TLatex(FitPars[i+1][1]+(FitPars[i][1]-FitPars[i+1][1])*0.2,FitPars[i][0],Form("#DeltaDp=%1.4f (Bias:%1.2f*10^{-4})",FitPars[i][1]-FitPars[i+1][1],10000.0*(FitPars[i][1]-FitPars[i+1][1]-(DpTheoreticalMap[i]-DpTheoreticalMap[i+1]))));
//		   txt->SetLineWidth(2);
//		   txt->SetTextSize(0.02);
//		   txt->Draw("same");
//
//	}
//
//
//	{
//		TArrow *ar5 = new TArrow(FitPars[999][1],FitPars[999][0],FitPars[999][6],FitPars[999][0],15,"<|>");
//		ar5->SetAngle(60);
//		ar5->SetLineWidth(4);
//		ar5->SetLineColor(4);
//		ar5->SetFillStyle(3008);
//		ar5->SetFillColor(2);
//		ar5->Draw("same");
//
//		TLatex *txt=new TLatex(FitPars[999][6]+(FitPars[999][1]-FitPars[999][6])*0.2,FitPars[999][0]+200,Form("#DeltaDp=%1.4f",FitPars[999][1]-FitPars[999][6]));
//		txt->SetLineWidth(2);
//		txt->SetTextSize(0.02);
//		txt->Draw("same");
//	}



	legend->Draw("same");


	//create


// 	 if(HRS=="L"){
// 	double DpTheoreticalList[8]={0.0134086,0.00337445,-0.00699935,-0.0122753,0.0154806,0.00542595,-0.00496891,-0.010256};
//	TLine *line[sizeof(DpTheoreticalList)/sizeof(double_t)];
//	for(int i =0; i < (sizeof(DpTheoreticalList)/sizeof(double_t)); i ++){
//		line[i]=new TLine(DpTheoreticalList[i],0,DpTheoreticalList[i],6000);
//		line[i]->SetLineWidth(2);
//		line[i]->Draw("same");
//	}
//	}else{
//	double DpTheoreticalList[8]={0.0141027,0.00387559,-0.00628469,-0.0157513,0.012034,0.00182792,-0.008312,-0.0177587};
//	TLine *line[sizeof(DpTheoreticalList)/sizeof(double_t)];
//	for(int i =0; i < (sizeof(DpTheoreticalList)/sizeof(double_t)); i ++){
//		line[i]=new TLine(DpTheoreticalList[i],0,DpTheoreticalList[i],6000);
//		line[i]->SetLineWidth(2);
//		line[i]->Draw("same");
//	}
//
//	}

	DpCanvas->Update();

//	auto legendHallProb = new TLegend(0.1,0.7,0.48,0.9);
//	legendHallProb->SetHeader("HallProb","C"); // option "C" allows to center the header
//	TCanvas *CentralPCanvas = (TCanvas*) gROOT->GetListOfCanvases()->FindObject("CentralPCanvas");
//	if (CentralPCanvas) {
//		CentralPCanvas->Clear();
//
//	} else {
//		CentralPCanvas = new TCanvas("CentralPCanvas", "CentralPCanvas", 600, 600);
//	}
//	CentralPCanvas->Draw();
//	std::map<int, TH1F*>HallProbh;
//	for (auto item = chainArray.begin(); item!= chainArray.end();item++){
//		HallProbh[item->first] = new TH1F("HallLProb", "HallLProb", 1000, -0.79, -0.75);
//		chainArray[item->first]->Project(HallProbh[item->first]->GetName(),
//				"HacL_D_LS450_FLD_DATA", generalcut.Data());
//		if (HallProbh[item->first]->GetEntries() != 0) {
//			double CentralP = std::abs(
//					(HallProbh[item->first]->GetMean()) * 0.95282 / 0.33930);
//			std::cout << "CentralMomentum is ::" << (CentralP) << std::endl;
//			if((item->first)<10){
//				HallProbh[item->first]->SetLineColor(8+item->first);
//				legendHallProb->AddEntry((HallProbh[item->first]),Form("C_{12} Dp:%2d%% scan",item->first));
//			}else{
//				HallProbh[item->first]->SetLineColor(kRed);
//				legendHallProb->AddEntry((HallProbh[item->first]),Form("H_{2}O Dp"));
//			}
//		}
//		if(item==chainArray.begin()){
//			HallProbh[item->first]->Draw();
//		}else{
//			HallProbh[item->first]->Draw("same");
//		}
//	}


//	// start the new canvas to draw the bias etc
//	TCanvas *DpBiasCanv=new TCanvas("Dp Bias Canvas","Dp Bias Canvas",1960,1080);
//	DpBiasCanv->SetLogy();
//	DpBiasCanv->Draw();
//	// calcualte the comman bias and subtract the common bias
//	double CommBiasC12=0.0;
//	double CommBiasH2O=0.0;
//	{
//		//for carbon, get the common bias to the 4th order
//		double CommBiasC12Sum=0.0;
//		int8_t CommBiasC12Counter=0;
//
//		double CommBiasH2OSum=0.0;
//		int8_t CommBiasH2OCounter=0;
//
//		for (auto item = OptDpArrayH.begin(); item != OptDpArrayH.end();item++) {
//			if (item->first < 10) {
//				const uint nPars=fitFunctionsList[item->first]->GetNpar();
//				double FitPars[nPars];
//				fitFunctionsList[item->first]->GetParameters(FitPars);
//
//				if (nPars>=10){
//					CommBiasC12Sum+=1000.0*(getC12TheoreticalDp(ScanrunListArray[item->first],0)-FitPars[1]);
//					CommBiasC12Counter+=1;
//
//					CommBiasC12Sum+=1000.0*(getC12TheoreticalDp(ScanrunListArray[item->first],1)-FitPars[6]);
//					CommBiasC12Counter+=1;
//
//				}
//				if (nPars>=13){
//					CommBiasC12Sum+=1000.0*(getC12TheoreticalDp(ScanrunListArray[item->first],2)-FitPars[11]);
//					CommBiasC12Counter+=1;
//
//				}
//				if(nPars >=16){
//					CommBiasC12Sum+=1000.0*(getC12TheoreticalDp(ScanrunListArray[item->first],3)-FitPars[14]);
//					CommBiasC12Counter+=1;
//				}
//
//			}else{
//
//				const uint nPars=fitFunctionsList[item->first]->GetNpar();
//				double FitPars[nPars];
//				fitFunctionsList[item->first]->GetParameters(FitPars);
//				CommBiasH2OSum+=1000.0*(getC12TheoreticalDp(ScanrunListArray[item->first],0)-FitPars[1]);
//				CommBiasH2OCounter++;
//				CommBiasH2OSum+=1000.0*(getC12TheoreticalDp(ScanrunListArray[item->first],1)-FitPars[6]);
//				CommBiasH2OCounter++;
//
//			}
//		}
//		CommBiasC12=CommBiasC12Sum/(double)CommBiasC12Counter;
//		CommBiasH2O=CommBiasH2OSum/(double)CommBiasH2OCounter;
//	}
//	DpBiasCanv->cd();
//	{
//		for (auto item = OptDpArrayH.begin(); item != OptDpArrayH.end();
//				item++) {
//			item->second->GetYaxis()->SetRangeUser(1, 100000);
//			if (item->first < 10) {
//				if (item == OptDpArrayH.begin()) {
//					item->second->Draw();
//				} else {
//					item->second->Draw("same");
//				}
//				fitFunctionsList[item->first]->Draw("same");
//
//
//				//get the excited states, calculat the seperation
//				const uint nPars=fitFunctionsList[item->first]->GetNpar();
//				double FitPars[nPars];
//				fitFunctionsList[item->first]->GetParameters(FitPars);
//				auto FitParErrors=fitFunctionsList[item->first]->GetParErrors();
//
//				if (nPars>=10){   //Ground and first excited states
//					// theoretical Seperation
//					double DpSepTheoretical=getC12TheoreticalDp(ScanrunListArray[item->first],0)-getC12TheoreticalDp(ScanrunListArray[item->first],1);
//					TLatex *txt0=new TLatex(FitPars[1],FitPars[0],Form("Dp_{0}:%1.3f Bias %1.3f #times 10^{3}",1000.0*FitPars[1],1000.0*(getC12TheoreticalDp(ScanrunListArray[item->first],0)-FitPars[1])-CommBiasC12));
//					txt0->SetTextSize(0.02);
//					txt0->Draw("same");
//
//					TLatex *txt1=new TLatex(FitPars[6],FitPars[5],Form("Dp_{1}:%1.3f Bias %1.3f #times 10^{3}",1000.0*FitPars[6],1000.0*(getC12TheoreticalDp(ScanrunListArray[item->first],1)-FitPars[6])-CommBiasC12));
//					txt1->SetTextSize(0.02);
//					txt1->Draw("same");
//
//				}
//
//				if (nPars>=13){   // second excited states
//					TLatex *txt1=new TLatex(FitPars[11],FitPars[10],Form("Dp_{2}:%1.3f Bias %1.3f #times 10^{3}",1000.0*FitPars[11],1000.0*(getC12TheoreticalDp(ScanrunListArray[item->first],2)-FitPars[11])-CommBiasC12));
//					txt1->SetTextSize(0.02);
//					txt1->Draw("same");
//				}
//
//				if(nPars >=16){   // third excited states
//					TLatex *txt1=new TLatex(FitPars[14],FitPars[13],Form("Dp_{3}:%1.3f Bias %1.3f #times 10^{3}",1000.0*FitPars[14],1000.0*(getC12TheoreticalDp(ScanrunListArray[item->first],3)-FitPars[14])-CommBiasC12));
//					txt1->SetTextSize(0.02);
//					txt1->Draw("same");
//				}
//			}
////			else{
////				// fit for the water target
////				item->second->Draw("same");
////				fitFunctionsList[item->first]->Draw("same");
////				const uint nPars=fitFunctionsList[item->first]->GetNpar();
////				double FitPars[nPars];
////				fitFunctionsList[item->first]->GetParameters(FitPars);
////				auto FitParErrors=fitFunctionsList[item->first]->GetParErrors();
////
////				TLatex *txt0=new TLatex(FitPars[1],FitPars[0]*1.5,Form("Dp_{0}:%1.3f Bias %1.3f #times 10^{3}",1000.0*FitPars[1],1000.0*(getC12TheoreticalDp(ScanrunListArray[item->first],0)-FitPars[1])-CommBiasH2O));
////				txt0->SetTextSize(0.02);
////				txt0->SetTextColor(kBlue);
////				txt0->Draw("same");
////				TLatex *txt1=new TLatex(FitPars[6],FitPars[5]*1.5,Form("Dp_{1}:%1.3f Bias %1.3f #times 10^{3}",1000.0*FitPars[6],1000.0*(getC12TheoreticalDp(ScanrunListArray[item->first],1)-FitPars[6])-CommBiasH2O));
////				txt1->SetTextSize(0.02);
////				txt1->SetTextColor(kBlue);
////				txt1->Draw("same");
////			}
//		}
//		legend->Draw("same");
//	}
//	DpBiasCanv->Update();
//
//
//	TCanvas *DpH2OBiasCanv=new TCanvas("Dp H20 Bias Canvas","Dp H2O Bias Canvas",1960,1080);
//	DpH2OBiasCanv->Draw();
//	DpH2OBiasCanv->SetLogy();
//	DpH2OBiasCanv->cd();
//	{
//		for (auto item = OptDpArrayH.begin(); item != OptDpArrayH.end();
//				item++) {
//			if (item->first > 10)
//			{
//				// fit for the water target
//				if(item->first==999){
//					item->second->Draw();
//				}else{
//					item->second->Draw("same");
//				}
//				fitFunctionsList[item->first]->Draw("same");
//				const uint nPars=fitFunctionsList[item->first]->GetNpar();
//				double FitPars[nPars];
//				fitFunctionsList[item->first]->GetParameters(FitPars);
//				auto FitParErrors=fitFunctionsList[item->first]->GetParErrors();
//
//				TLatex *txt0=new TLatex(FitPars[1],FitPars[0]*1.5,Form("Dp_{0}:%1.3f Bias %1.3f #times 10^{3}(%1.1fKeV) ",1000.0*FitPars[1],1000.0*(getC12TheoreticalDp(ScanrunListArray[item->first],0)-FitPars[1])-CommBiasH2O,1000.0*CentralPArray[item->first]*(1000.0*(getC12TheoreticalDp(ScanrunListArray[item->first],0)-FitPars[1])-CommBiasH2O)));
//				txt0->SetTextSize(0.02);
//				txt0->SetTextColor(kBlue);
//				txt0->Draw("same");
//				TLatex *txt1=new TLatex(FitPars[6],FitPars[5]*1.5,Form("Dp_{1}:%1.3f Bias %1.3f #times 10^{3}(%1.1fKeV) ",1000.0*FitPars[6],1000.0*(getC12TheoreticalDp(ScanrunListArray[item->first],1)-FitPars[6])-CommBiasH2O,1000.0*CentralPArray[item->first]*(1000.0*(getC12TheoreticalDp(ScanrunListArray[item->first],1)-FitPars[6])-CommBiasH2O)));
//				txt1->SetTextSize(0.02);
//				txt1->SetTextColor(kBlue);
//				txt1->Draw("same");
//			}
//		}
//
//		legend->Draw("same");
//	}
//	DpH2OBiasCanv->Update();
//
//
//	TCanvas *MomC12BiasCanv=new TCanvas("P Bias Canvas","P Bias Canvas",1960,1080);
//	MomC12BiasCanv->SetLogy();
//	MomC12BiasCanv->Draw();
//
//	{
//			for (auto item = OptDpArrayH.begin(); item != OptDpArrayH.end();
//					item++) {
//				item->second->GetYaxis()->SetRangeUser(1, 100000);
//				if (item->first < 10) {
//					if (item == OptDpArrayH.begin()) {
//						item->second->Draw();
//					} else {
//						item->second->Draw("same");
//					}
//					fitFunctionsList[item->first]->Draw("same");
//
//
//					//get the excited states, calculat the seperation
//					const uint nPars=fitFunctionsList[item->first]->GetNpar();
//					double FitPars[nPars];
//					fitFunctionsList[item->first]->GetParameters(FitPars);
//					auto FitParErrors=fitFunctionsList[item->first]->GetParErrors();
//
//					if (nPars>=10){   //Ground and first excited states
//						// theoretical Seperation
//						double DpSepTheoretical=getC12TheoreticalDp(ScanrunListArray[item->first],0)-getC12TheoreticalDp(ScanrunListArray[item->first],1);
//						TLatex *txt0=new TLatex(FitPars[1],FitPars[0],Form("P_{0}:%1.3fGeV Bias %1.3fKeV",FitPars[1]*CentralPArray[item->first]+CentralPArray[item->first],1000.0*(1000.0*(getC12TheoreticalDp(ScanrunListArray[item->first],0)-FitPars[1])-CommBiasC12)*CentralPArray[item->first]));
//						txt0->SetTextSize(0.02);
//						txt0->Draw("same");
//
//						TLatex *txt1=new TLatex(FitPars[6],FitPars[5],Form("P_{1}:%1.3fGeV Bias %1.3fKeV",FitPars[6]*CentralPArray[item->first]+CentralPArray[item->first],CentralPArray[item->first]*1000.0*(1000.0*(getC12TheoreticalDp(ScanrunListArray[item->first],1)-FitPars[6])-CommBiasC12)));
//						txt1->SetTextSize(0.02);
//						txt1->Draw("same");
//
//					}
//
//					if (nPars>=13){   // second excited states
//						TLatex *txt1=new TLatex(FitPars[11],FitPars[10],Form("P_{2}:%1.3fGeV Bias %1.3fKeV",FitPars[11]*CentralPArray[item->first]+CentralPArray[item->first],(1000.0*(getC12TheoreticalDp(ScanrunListArray[item->first],2)-FitPars[11])-CommBiasC12)*CentralPArray[item->first]*1000.0));
//						txt1->SetTextSize(0.02);
//						txt1->Draw("same");
//					}
//
//					if(nPars >=16){   // third excited states
//						TLatex *txt1=new TLatex(FitPars[14],FitPars[13],Form("P_{3}:%1.3fGeV Bias %1.3fKeV",FitPars[14]*CentralPArray[item->first]+CentralPArray[item->first],(1000.0*(getC12TheoreticalDp(ScanrunListArray[item->first],3)-FitPars[14])-CommBiasC12)*CentralPArray[item->first]*1000.0));
//						txt1->SetTextSize(0.02);
//						txt1->Draw("same");
//					}
//				}
////				else{
////					// fit for the water target
////					item->second->Draw("same");
////					fitFunctionsList[item->first]->Draw("same");
////					const uint nPars=fitFunctionsList[item->first]->GetNpar();
////					double FitPars[nPars];
////					fitFunctionsList[item->first]->GetParameters(FitPars);
////					auto FitParErrors=fitFunctionsList[item->first]->GetParErrors();
////
////					TLatex *txt0=new TLatex(FitPars[1],FitPars[0]*1.5,Form("Dp_{0}:%1.3f Bias %1.3f #times 10^{3}",1000.0*FitPars[1],1000.0*(getC12TheoreticalDp(ScanrunListArray[item->first],0)-FitPars[1])-CommBiasH2O));
////					txt0->SetTextSize(0.02);
////					txt0->SetTextColor(kBlue);
////					txt0->Draw("same");
////					TLatex *txt1=new TLatex(FitPars[6],FitPars[5]*1.5,Form("Dp_{1}:%1.3f Bias %1.3f #times 10^{3}",1000.0*FitPars[6],1000.0*(getC12TheoreticalDp(ScanrunListArray[item->first],1)-FitPars[6])-CommBiasH2O));
////					txt1->SetTextSize(0.02);
////					txt1->SetTextColor(kBlue);
////					txt1->Draw("same");
////				}
//			}
//			legend->Draw("same");
//		}
//	MomC12BiasCanv->Update();
//
//	// get the Momentum seperation
//	TCanvas *MomSepBiasCanv=new TCanvas("Momentum Sep Bias Canvas","Momentum Sep Bias Canvas",1960,1080);
//	MomSepBiasCanv->Draw("same");
//	MomSepBiasCanv->SetLogy();
//	{
//		for (auto item = OptDpArrayH.begin(); item != OptDpArrayH.end();
//				item++) {
//			if (item->first < 10) {
//				item->second->GetYaxis()->SetRangeUser(1, 100000);
//				if (item == OptDpArrayH.begin()) {
//					item->second->Draw();
//				} else {
//					item->second->Draw("same");
//				}
//				fitFunctionsList[item->first]->Draw("same");
//
//
//				//get the excited states, calculat the seperation
//				const uint nPars=fitFunctionsList[item->first]->GetNpar();
//				double FitPars[nPars];
//				fitFunctionsList[item->first]->GetParameters(FitPars);
//				auto FitParErrors=fitFunctionsList[item->first]->GetParErrors();
//
//				if (nPars>=10){   //Ground and first excited states
//					// theoretical Seperation
//					double MomSep=CentralPArray[item->first]*(FitPars[1]-FitPars[6]);
//					double MomBias=CentralPArray[item->first]*((getC12TheoreticalDp(ScanrunListArray[item->first],0)-getC12TheoreticalDp(ScanrunListArray[item->first],1))-(FitPars[1]-FitPars[6]));
//
//					TLatex *txt0=new TLatex(FitPars[1],FitPars[0],Form("#DeltaP_{1}:%1.3fMeV Bias %1.1fKeV",MomSep*1000.0,MomBias*1000000.0));
//					txt0->SetTextSize(0.02);
//					txt0->Draw("same");
//
//				}
//
//				if (nPars>=13){   // second excited states
//					double MomSep=CentralPArray[item->first]*(FitPars[1]-FitPars[11]);
//					double MomBias=CentralPArray[item->first]*((getC12TheoreticalDp(ScanrunListArray[item->first],0)-getC12TheoreticalDp(ScanrunListArray[item->first],2))-(FitPars[1]-FitPars[11]));
//
//					TLatex *txt0=new TLatex(FitPars[11],FitPars[10]*0.5+FitPars[5]*0.5,Form("#DeltaP_{2}:%1.3fMeV Bias %1.1fKeV",MomSep*1000.0,MomBias*1000000.0));
//					txt0->SetTextSize(0.02);
//					txt0->Draw("same");				}
//
//				if(nPars >=16){   // third excited states
//					double MomSep=CentralPArray[item->first]*(FitPars[1]-FitPars[14]);
//					double MomBias=CentralPArray[item->first]*((getC12TheoreticalDp(ScanrunListArray[item->first],0)-getC12TheoreticalDp(ScanrunListArray[item->first],3))-(FitPars[1]-FitPars[14]));
//
//					TLatex *txt0=new TLatex(FitPars[14],FitPars[13]*0.5+FitPars[10]*0.5,Form("#DeltaP_{3}:%1.3fMeV Bias %1.1fKeV",MomSep*1000.0,MomBias*1000000.0));
//					txt0->SetTextSize(0.02);
//					txt0->Draw("same");
//				}
//			}else{
//				// fit for the water target
//			}
//		}
//		legend->Draw("same");
//
//	}


}

