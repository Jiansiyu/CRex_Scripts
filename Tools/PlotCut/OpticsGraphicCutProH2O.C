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
#include <TMath.h>
#include <TH2F.h>
#include <TH1.h>
#include <TF1.h>
#include <TPaveText.h>
#include <map>
#include <vector>
#include <random>
#include <iostream>

#include <TVirtualPad.h>

#include <TObject.h>
#include <TFile.h>
#include <fstream>
#include <TLatex.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TVector3.h>
#include <TRotation.h>
#include "TImage.h"
#include <Math/Functor.h>
#include "Math/RootFinder.h"
#include "TGraphErrors.h"

int FoilID=0;

int col=3;
int row=1;
int row_min=0;
int row_count=10;
const UInt_t NSieveCol = 13;
const UInt_t NSieveRow = 7;
double CentralP;

//////////////////////////////////////////////////////////////////////////////
// Work Directory
// cut options
// Need to change
//////////////////////////////////////////////////////////////////////////////
TString prepcut;
TString generalcut;
//TString generalcutR="R.tr.n==1 && R.vdc.u1.nclust==1&& R.vdc.v1.nclust==1 && R.vdc.u2.nclust==1 && R.vdc.v2.nclust==1 && fEvtHdr.fEvtType==1 && R.gold.p > 2.14 && R.gold.p < 2.2";
TString generalcutR="R.tr.n==1";// && R.vdc.u1.nclust==1&& R.vdc.v1.nclust==1 && R.vdc.u2.nclust==1 && R.vdc.v2.nclust==1 && R.gold.p > 2.14 && R.gold.p < 2.19";
//TString generalcutL="L.tr.n==1 && L.vdc.u1.nclust==1&& L.vdc.v1.nclust==1 && L.vdc.u2.nclust==1 && L.vdc.v2.nclust==1 && fEvtHdr.fEvtType==1 && L.gold.p > 2.14 && L.gold.p < 2.19";
TString generalcutL="L.tr.n==1";// && L.vdc.u1.nclust==1&& L.vdc.v1.nclust==1 && L.vdc.u2.nclust==1 && L.vdc.v2.nclust==1  && L.gold.p > 2.14 && L.gold.p < 2.19";


///
/// \param DeltaE
/// \param BeamE
/// \return
inline double GetPointingAngle(double DeltaE, double BeamE=2.17568){
	//In pointing measurement,
	// calculate the HRS angle according to H and O seperation

	Double_t amu = 931.494028 * 1.e-3;  // amu to GeV
	Double_t H2O = 18 * amu;
	Double_t TargetH = amu;
	double_t TargetO = 16 * amu;
	Double_t mass_tg = H2O;
	double beamE = BeamE;         // in GeV
	// get the pointing measurement
	double deltaE = DeltaE; // in GeV

	double a = 4.0 * beamE * deltaE;
	double b = 2.0 * deltaE * (TargetH + TargetO)+ 2.0 * beamE * (TargetH - TargetO);
	double c = deltaE * TargetH * TargetO / beamE;
	double sin2theta = ((0. - b) - TMath::Sqrt(b * b - 4.0 * a * c)) / (2.0 * a);
	double HRSAngleReal=TMath::ASin(TMath::Sqrt(sin2theta)) * 2.0 * 180.0 / TMath::Pi();
	return HRSAngleReal;
}

///
/// \param ScatterAng
/// \param Targ_theta
/// \param Targ_phi
/// \return
inline double getP0(double ScatterAng, double Targ_theta, double Targ_phi){

    ScatterAng=ScatterAng*TMath::Pi()/180.0;
    double A=2.0*TMath::Sqrt(1+Targ_phi*Targ_phi+Targ_theta*Targ_theta)*TMath::Cos(ScatterAng);
    double B=TMath::Sqrt(2)*TMath::Sqrt(Targ_phi*Targ_phi - Targ_phi*Targ_phi*Targ_theta*Targ_theta + Targ_phi*Targ_phi*Targ_phi*Targ_phi-
            Targ_phi*Targ_phi*TMath::Cos(2*ScatterAng)-Targ_phi*Targ_phi*Targ_theta*Targ_theta*TMath::Cos(2*ScatterAng)-TMath::Power(Targ_phi,4)*TMath::Cos(2*ScatterAng));

    double C=1.0/(2*(1+TMath::Power(Targ_phi,2)));

    return TMath::ACos(C*(A-B))*180/TMath::Pi();
}

///
/// \param name
/// \return
inline Bool_t IsFileExist (const std::string& name) {
    return !gSystem->AccessPathName(name.c_str());
}

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

//        TH1F *HallEpicsBeamE = new TH1F("HallEpicsBeamE", "HallEpicsBeamE", 1000, beamERangeMin, beamERangeMax);

        TH1F *HallEpicsBeamE=(TH1F *)gROOT->FindObject("HallEpicsBeamE");
        if(HallEpicsBeamE) HallEpicsBeamE->Delete();
        HallEpicsBeamE=new TH1F("HallEpicsBeamE", "HallEpicsBeamE", 1000, beamERangeMin, beamERangeMax);

        chain->Project(HallEpicsBeamE->GetName(),"HALLA_p");
        double epicsBeamE=HallEpicsBeamE->GetMean();
        if ((epicsBeamE > 2170)&&(epicsBeamE < 2180)){
            std::cout<<"\033[1;32m [Infor]\033[0m Read in the Beam E in ROOT file(with correction): "<<epicsBeamE*953.4/951.1<<std::endl;
            return  epicsBeamE/1000.0*2182.2152/2176;
        }
    }

    //read in the beamE information and parser
    if ((!beamEfname.IsNull()) && IsFileExist(beamEfname.Data())){
        std::cout<<"\033[1;32m [Infor]\033[0m Read in the Beam E file: "<<beamEfname.Data()<<std::endl;
        std::ifstream infile(beamEfname.Data());
        int runID_temp;
        float beamE_temp;
        while (infile >> runID_temp >> beamE_temp){
            beamE[runID_temp]=beamE_temp/1000.0;
        }
    }else{
        std::cout<<"\033[1;33m [Warning]\033[0m can not find file "<<beamEfname.Data()<<" Skip the beamE file!!!"<<std::endl;

    }
    if(beamE.find(runID)!=beamE.end()){
        std::cout<<"\033[1;32m [Infor]\033[0m Read in the Beam E file(with correction): "<<beamE[runID]*953.4/951.1<<std::endl;
        return beamE[runID]*2182.2152/2176;
    }else{
        std::cout<<"\033[1;31m [CAUTION]\033[0m Can not find the Beam E for run"<<runID<<" Using default value!!!["<<__func__<<"("<<__LINE__<<")]"<<std::endl;
        exit(-1);
    }
}

/// Read the NMR/Hall Probe and calculate the cental P
/// \param chain
/// \return
double_t getCentralP(TChain *chain){
    TString HRS("R");
    TString filename(chain->GetFile()->GetName());
    if (filename.Contains("RHRS")) {
    } else if (filename.Contains("LHRS")) {
        HRS = "L";
    }

    double CentralP;
    if (HRS == "L") {
        TH1F *HallProbHH = new TH1F("HallLProb", "HallLProb", 1000, -1, 0);
        chain->Project(HallProbHH->GetName(), "HacL_D_LS450_FLD_DATA",
                       generalcut.Data());
        CentralP = std::abs((HallProbHH->GetMean()) * 0.95282 / 0.33930);
        std::cout << "CentralMomentum is (LHRS) for Hall Probe::" << (CentralP)
                  << std::endl;
    } else {
        //HacR_D1_NMR_SIG
        TH1F *HallR_NMR = new TH1F("HallR_NMR", "HallR_NMR", 1000, 0.7, 0.9);
        chain->Project(HallR_NMR->GetName(), "HacR_D1_NMR_SIG",
                       generalcut.Data());
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

double GLobalSovler_scatteredAngle;
double GLobalSovler_theta_tg;
double GLobalSovler_phi_tg;

///
/// \param HRSAngle
/// \return
double GetHRSAngle(double HRSAngle){

    TVector3 TCSX(0, -1, 0);
    TVector3 TCSZ(TMath::Sin(HRSAngle), 0, TMath::Cos(HRSAngle));
    TVector3 TCSY = TCSZ.Cross(TCSX);
    TRotation fTCSInHCS;
    fTCSInHCS.RotateAxes(TCSX, TCSY, TCSZ);
    TVector3 MomDirectionTCS(GLobalSovler_theta_tg,GLobalSovler_phi_tg,1); // target variables
    TVector3 MomDirectionHCS = fTCSInHCS*MomDirectionTCS;
    TVector3 BeamDirection(0, 0, 1);
    const Double_t ScatteringAngle = BeamDirection.Angle(MomDirectionHCS);
    return ScatteringAngle - GLobalSovler_scatteredAngle;
}

///
/// \param momentumSpectro
/// \return
TF1 *SpectroCrystalFit_H2O(TH1F*momentumSpectro){
    // locate the H peak
    auto CGroundp = momentumSpectro->GetXaxis()->GetBinCenter(momentumSpectro->GetMaximumBin());
    momentumSpectro->GetXaxis()->SetRangeUser(CGroundp - 0.03,CGroundp + 0.0044*1.7);
    momentumSpectro->GetXaxis()->SetTitle("Mom");
    momentumSpectro->GetYaxis()->SetTitle("#");

    double_t fgroudGausPar[3];
    TF1 *fgroundGaus = new TF1("groundStateGaus","gaus",CGroundp - 0.0005,CGroundp + 0.0005);
    momentumSpectro->Fit(fgroundGaus->GetName(),"RQ0","ep",fgroundGaus->GetXmin(),fgroundGaus->GetXmax());
    fgroundGaus->GetParameters(fgroudGausPar);

    TH1F *test  = (TH1F *) momentumSpectro->Clone("fitTest");
    test ->GetXaxis()->SetRangeUser(momentumSpectro->GetXaxis()->GetXmin(),fgroudGausPar[1]-10*fgroudGausPar[2]);

    auto C1stp = test ->GetXaxis()->GetBinCenter(test->GetMaximumBin());

    double_t ffirstGuasPar[3];
    TF1 *ffirstGuas=new TF1 ("firststatesgaus","gaus",C1stp-0.0015,C1stp+0.00155);
    momentumSpectro->Fit("firststatesgaus","RQ0","ep",ffirstGuas->GetXmin(),ffirstGuas->GetXmax());
    ffirstGuas->GetParameters(ffirstGuasPar);

    // change the gause fit to cristal ball
    double_t fgroundCrystalballPar[5];
    TF1 *fgroundCrystalball=new TF1("fgroundCrystal","crystalball",fgroudGausPar[1]-0.0030,fgroundGaus->GetXmax()+0.0003);
    fgroundCrystalball->SetParameters(fgroudGausPar[0],fgroudGausPar[1],fgroudGausPar[2],1.64,1.1615);
    momentumSpectro->Fit("fgroundCrystal","RQ0","same",fgroundCrystalball->GetXmin(),fgroundCrystalball->GetXmax());
    fgroundCrystalball->GetParameters(fgroundCrystalballPar);

    double_t ffirstCrystalPar[5];
    TF1 *ffirstCrystal=new TF1("ffirstCrystal","crystalball",ffirstGuasPar[1]-0.0025,ffirstGuas->GetXmax());
    ffirstCrystal->SetParameters(ffirstGuasPar[0],ffirstGuasPar[1],ffirstGuasPar[2],1.64,1.1615);
    momentumSpectro->Fit("ffirstCrystal","RQ0","ep",ffirstCrystal->GetXmin(),ffirstCrystal->GetXmax());
    ffirstCrystal->GetParameters(ffirstCrystalPar);

    double_t fCrystalMomentumPar[10];
    TF1 *fCrystalMomentum=new TF1("fCrystalMomentum","crystalball(0)+crystalball(5)",ffirstCrystal->GetXmin(),fgroundCrystalball->GetXmax());
    std::copy(fgroundCrystalballPar,fgroundCrystalballPar+5,fCrystalMomentumPar);
    std::copy(ffirstCrystalPar,ffirstCrystalPar+5,fCrystalMomentumPar+5);
    fCrystalMomentum->SetParameters(fCrystalMomentumPar);
    momentumSpectro->Fit("fCrystalMomentum","RQ0","",fCrystalMomentum->GetXmin(),fCrystalMomentum->GetXmax());
    return fCrystalMomentum;
}


///
/// \param runID
/// \param folder
/// \return
TChain *LoadrootFile(UInt_t runID,TString folder="/home/newdriver/Storage/Research/CRex_Experiment/RasterReplay/Replay/Result"){
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

///
/// \param runID
/// \param csvfname
/// \return TVector2 with the beam Position on target, Unit/meter
TVector2 getBPM(UInt_t runID,TString csvfname="bpm_on_targ.csv"){
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

        TVector2 vec(bpmx/1000.,bpmy/1000.);

        bpmList[runs]=vec;
    }
    if(bpmList.find(runID) != bpmList.end()){
        return bpmList[runID];
    } else{
        std::cout<<"\033[1;33m [Warning]\033[0m Missing csv file::"<<csvfname.Data()<<std::endl;
        exit(-1);
    }
}


/// Use get the Sieve Positions on the target cooridnation system
/// \param HRS  L/R
/// \param col  col Index
/// \param row  row index
/// \return     TVecter3 with each position
TVector3 _getSievePos(TString HRS, UInt_t col, UInt_t row){

    if (HRS == "L"){
        auto ZPos = (97.9) * 1e-2;
        // return the Sieve hole position
        assert(col < NSieveCol);
        assert(row < NSieveRow);

        // for prex experiment,
        const Double_t PositionIndex_y[]={-20.52, -18.12, -15.72, -10.92,
                                          -8.52,  -6.12,      0.0,	3.06,
                                          6.12,    12.24, 15.30,   18.30,
                                          21.42};

        const Double_t XOffset[]      = {0.0,	    6.65,	   0.0,	    0.0,    // 0 - 3
                                         6.65,	   0.0,	   0.0,    6.65,      // 4 - 7
                                         0.0,	   0,	6.65,      0,
                                         6.65};

        double_t rowID_temp=(double_t)row -3.0;

        Double_t PRex_x=rowID_temp*13.3*(1e-3)+XOffset[col]*(1e-3);
        Double_t PRex_y=PositionIndex_y[col]*(1e-3);

        assert(PRex_x<0.2 && PRex_x<-0.2);
        assert(PRex_y<0.2 && PRex_y<-0.2);

        return TVector3(PRex_x,PRex_y,ZPos);
    } else{
        auto ZPos = (97.923) * 1e-2;
        assert(col < NSieveCol);
        assert(row < NSieveRow);

        // temporary uasage
        // for prex experiment,
        const Double_t PositionIndex_y[]={-20.52, -18.12, -15.72, -10.92,
                                          -8.52,  -6.12,      0.0,	3.06,
                                          6.12,    12.24, 15.30,   18.30,
                                          21.42
        };
        const Double_t XOffset[]      = {0.0,	    6.65,	   0.0,	    0.0,    // 0 - 3
                                         6.65,	   0.0,	   0.0,    6.65,      // 4 - 7
                                         0.0,	   0,	6.65,      0,
                                         6.65};
        //  the gap of each two hole 13.3mm
        double_t rowID_temp=(double_t)row -3.0;

        Double_t PRex_x=rowID_temp*13.3*(1e-3)+XOffset[col]*(1e-3);
        Double_t PRex_y=PositionIndex_y[col]*(-1e-3);

        assert(PRex_x<0.2 && PRex_x<-0.2);
        assert(PRex_y<0.2 && PRex_y<-0.2);

        return TVector3(PRex_x,PRex_y,ZPos);
    }


}
TVector2 _getSieveAngleOnTargCoor(TVector2 beamPos,TString HRS, UInt_t col,UInt_t row,double HRSAngle = 0){
    if (HRSAngle == 0){
        if (HRS == "R"){
            HRSAngle = -4.822 * TMath::Pi()/180.0;
        }else{
            HRSAngle = 4.818 * TMath::Pi() / 180.;
        }
    }
    // canclulate the theoratical position of the sieve hole
    auto SieveHoleTCS = _getSievePos(HRS,col,row);
    const TVector3 BeamSpotHCS(beamPos.X(),beamPos.Y(),0.0);

    const Int_t a = (HRSAngle > 0) ? 1 : -1;

    // MissPoint* are in HCS
    const Double_t MissPointZ =0.0;//
    const Double_t MissPointY = 0.0;//

    TRotation fTCSInHCS; // transformations vector from TCS to HCS
    TVector3 TCSX(0, -1, 0);
    TVector3 TCSZ(TMath::Sin(HRSAngle), 0, TMath::Cos(HRSAngle));
    TVector3 TCSY = TCSZ.Cross(TCSX);
    fTCSInHCS.RotateAxes(TCSX, TCSY, TCSZ);

    TVector3 fPointingOffset(-a*MissPointZ*TMath::Cos(HRSAngle), MissPointY, -MissPointZ * TMath::Sin(HRSAngle));
    const TVector3 BeamSpotTCS = fTCSInHCS.Inverse()*(BeamSpotHCS - fPointingOffset);
    const TVector3 MomDirectionTCS = SieveHoleTCS - BeamSpotTCS;   //

    double kRealThVal = MomDirectionTCS.X() / MomDirectionTCS.Z();
    double kRealPhVal = MomDirectionTCS.Y() / MomDirectionTCS.Z();

    const Double_t x_tg = BeamSpotTCS.X() - BeamSpotTCS.Z() * kRealThVal;
    const Double_t y_tg = BeamSpotTCS.Y() - BeamSpotTCS.Z() * kRealPhVal;

    const Double_t ExtTarCor_ThetaCorr = 0.61;//0.00;//
    const Double_t ExtTarCor_DeltaCorr = 5.18;//1e36;//

    double kRealThMatrixVal = kRealThVal - x_tg * ExtTarCor_ThetaCorr;

    TVector2 SieveAnglePos(kRealThMatrixVal,kRealPhVal);
    return  SieveAnglePos;
}

std::map<UInt_t,TH1F *> getOptRange(TString HRS, UInt_t col, UInt_t row){
    // color for each run
    //run List of the Carbon Scan
    std::map<UInt_t,TH1F *> resMomh;

    return  resMomh;
}


///
/// \param runID
/// \param cutFile
/// \param folder
void getAllSievePointing(UInt_t runID=2671,
		TString cutFile =
				"/home/newdriver/Storage/Research/PRex_Workspace/PREX-MPDGEM/PRexScripts/Tools/PlotCut/Result/Water/cut_20201206/WithOutMomCut/prexLHRS_2671_-1.root.FullCut.root",
		TString folder =
				"/home/newdriver/Storage/Research/CRex_Experiment/RasterReplay/Replay/Result") {

	TChain *chain=LoadrootFile(runID,folder);
	// if the folder itself is and root file
	TString HRS="R";
	if(runID<20000){HRS="L";};

	//load the cut and load the canvas
	//plot the theta and phi, and load the cut file
	if(HRS=="L"){
		generalcut=generalcutL;
	}else{
		generalcut=generalcutR;
	}

	double CentralP=getCentralP(chain);


	TCanvas *mainPatternCanvas=(TCanvas *)gROOT->GetListOfCanvases()->FindObject("cutPro");

	if(!mainPatternCanvas){
		mainPatternCanvas=new TCanvas("cutPro","cutPro",600,600);
	}else{
		mainPatternCanvas->Clear();
	}
	mainPatternCanvas->Divide(2,2);
	mainPatternCanvas->Draw();

	mainPatternCanvas->cd(1);

	// check initial data set without cut and check it with cut
	TH2F *TargetThPh_SieveNoCutHH=(TH2F *)gROOT->FindObject("th_vs_ph_No_cut");
	if(TargetThPh_SieveNoCutHH) TargetThPh_SieveNoCutHH->Delete();
	TargetThPh_SieveNoCutHH=new TH2F("th_vs_ph_No_cut","th_vs_ph with No Cut",1000,-0.03,0.03,1000,-0.045,0.045);
	chain->Project(TargetThPh_SieveNoCutHH->GetName(),Form("%s.gold.th:%s.gold.ph",HRS.Data(),HRS.Data()));   // draw all the data with no cut
	TargetThPh_SieveNoCutHH->Draw("zcol");

	if(!(TargetThPh_SieveNoCutHH->GetEntries())){
		TPaveText *text=new TPaveText(0.2,0.2,0.8,0.8,"NDC");
		text->AddText("No Event Found !!!");
		text->SetTextColor(kRed);
		text->Draw("same");
	}

	mainPatternCanvas->cd(2);
	TH2F *TargetThPhGeneralCutHH=(TH2F *)gROOT->FindObject("th_vs_ph_GeneralCut");
	if(TargetThPhGeneralCutHH) TargetThPhGeneralCutHH->Delete();
	TargetThPhGeneralCutHH=new TH2F("th_vs_ph_GeneralCut","th_vs_ph General Cut",1000,-0.03,0.03,1000,-0.045,0.045);
	chain->Project(TargetThPhGeneralCutHH->GetName(),Form("%s.gold.th:%s.gold.ph",HRS.Data(),HRS.Data()),generalcut.Data());
	TargetThPhGeneralCutHH->Draw("zcol");
	if(!(TargetThPhGeneralCutHH->GetEntries())){
			TPaveText *text=new TPaveText(0.2,0.2,0.8,0.8,"NDC");
			text->AddText("No Event Found !!!");
			text->SetTextColor(kRed);
			text->Draw("same");
		}

	mainPatternCanvas->cd(3);
	TH2F *TargetThPhHH=(TH2F *)gROOT->FindObject("th_vs_ph");
	if(TargetThPhHH) TargetThPhHH->Delete();
	TargetThPhHH=new TH2F("th_vs_ph","th_vs_ph",1000,-0.03,0.03,1000,-0.045,0.045);

	chain->Project(TargetThPhHH->GetName(),Form("%s.gold.th:%s.gold.ph",HRS.Data(),HRS.Data()),generalcut.Data());
	TargetThPhHH->Draw("zcol");

	// read the cut file and load the cut, create the new data set with the cut
	// load the cut file
	TFile *cutFileIO=new TFile(cutFile.Data(),"READ");
	if(cutFileIO->IsZombie()){
		std::cout<<"[ERROR]:: CAN NOT FIND CUT FILE \" "<<cutFile.Data()<<"\""<<std::endl;
	}

	//loop on the files in the cut and find all the sieve hole cuts
	TCutG *sieveCut[NSieveCol][NSieveRow];
	TCut sieveAllHoleCut;
	for (int16_t col = 0; col < NSieveCol; col++){
		for (int16_t row = 0; row < NSieveRow; row++){
			auto cutg=(TCutG*)gROOT->FindObject(Form("hcut_R_%d_%d_%d", FoilID, col, row));
			if(cutg){
				sieveCut[col][row]=cutg;
				sieveCut[col][row]->SetLineWidth(2);
				sieveCut[col][row]->SetLineColor(kRed);
				sieveCut[col][row]->Draw("same");
				sieveAllHoleCut=sieveAllHoleCut||TCut(Form("hcut_R_%d_%d_%d", FoilID, col, row));

				//get the data for this canvas
				TH2F *selectedSievehh=(TH2F *)  gROOT->FindObject("Sieve_Selected_th_ph");
				if (selectedSievehh) {
					selectedSievehh->Clear();
				} else {
					selectedSievehh = new TH2F("Sieve_Selected_th_ph",
							"Sieve_Selected_th_ph", 1000,
							TargetThPhHH->GetXaxis()->GetXmin(),
							TargetThPhHH->GetXaxis()->GetXmax(), 1000,
							TargetThPhHH->GetYaxis()->GetXmin(),
							TargetThPhHH->GetYaxis()->GetXmax());
				}
				chain->Project(selectedSievehh->GetName(),Form("%s.gold.th:%s.gold.ph",HRS.Data(),HRS.Data()),Form("%s&&%s",sieveCut[col][row]->GetName(),generalcut.Data()));
				TLatex *label=new TLatex(selectedSievehh->GetMean(1),selectedSievehh->GetMean(2),Form("(%d %d)",col,row));
				label->SetTextSize(0.04);
				label->SetTextColor(2);
				label->Draw("same");
				selectedSievehh->Delete();
			}
		}
	}

	sieveAllHoleCut=sieveAllHoleCut+TCut(generalcut.Data());

	mainPatternCanvas->cd(4);
	TH2F *TargetThPh_SieveCutHH=(TH2F *)gROOT->FindObject("th_vs_ph_cut");
	if(TargetThPh_SieveCutHH) TargetThPh_SieveCutHH->Delete();
	TargetThPh_SieveCutHH=new TH2F("th_vs_ph_cut","th_vs_ph_cut",1000,-0.03,0.03,1000,-0.045,0.045);
	chain->Project(TargetThPh_SieveCutHH->GetName(),Form("%s.gold.th:%s.gold.ph",HRS.Data(),HRS.Data()),sieveAllHoleCut);
	TargetThPh_SieveCutHH->Draw("zcol");
	mainPatternCanvas->Update();

//    mainPatternCanvas->Print(Form("Pointing_run%d.pdf(",runID),"pdf");

    TCanvas * inforCanv = new TCanvas("InforCanv","InforCanv",1960,1080);
    inforCanv->Draw();
    inforCanv->cd();

    TargetThPh_SieveCutHH->Draw("zcol");
    inforCanv->Update();
    inforCanv->Print(Form("Pointing_run%d.pdf(",runID),"pdf");


    // fit each individual momentum
	TCanvas *SievePointingCanv=new TCanvas(Form("SievePointing"),Form("SievePointing"),1960,1080);
    SievePointingCanv->Draw();

	std::map<int, std::map<int, double *>> MomentumFitParArray;
	std::map<int, std::map<int, TH1F *>> SieveMomentumArray;

	// plot all the Momentum, and get  the momentum difference
    std::map<int, std::map<int, double>> scatteredAngleArray; // initial scatteredAngle
    std::map<int, std::map<int, double>> HRSAngleArray; // HRS pointing angle
    std::map<int, std::map<int, double>> scatteredAngleEnergyDiff;


    // theta and phi value
    std::map<int, std::map<int, double>> thetaValueArray;
    std::map<int, std::map<int, double>> phiValueArray;

    for (int16_t col = 0; col < NSieveCol; col++) {
		for (int16_t row = 0; row < NSieveRow; row++) {
			auto cutg = (TCutG*) gROOT->FindObject(Form("hcut_R_%d_%d_%d", FoilID, col, row));

			if (cutg) {
                SievePointingCanv->Clear();
                SievePointingCanv->Divide(1,2);
                SievePointingCanv->cd(2)->Divide(4,1);

			    SievePointingCanv->SetLogy();

				std::cout<<"Fitting col"<<col<<"row"<<row<<std::endl;
				SieveMomentumArray[col][row]=new TH1F(Form("Sieve_Col%d_Row%d_Momentum",col,row),Form("Sieve_Col%d_Row%d_Momentum",col,row),1000,2.1,2.2);
				TCut sieveCut=TCut(generalcut.Data())+TCut(Form("hcut_R_%d_%d_%d", FoilID, col, row));
				chain->Project(SieveMomentumArray[col][row]->GetName(),Form("%s.gold.dp*%f+%f",HRS.Data(),CentralP,CentralP),sieveCut);


				if(SieveMomentumArray[col][row]->GetEntries())
				{
                    //get the target theta and phi angle from the projected sieve holesm use those information together with the angle to get the central Sieve
                    TH2F *TargetThPhSieve=(TH2F *)gROOT->FindObject(Form("th_vs_ph_cut_col%d_row%d",col,row));
                    if(TargetThPhSieve) TargetThPhSieve->Delete();
                    TargetThPhSieve=new TH2F(Form("th_vs_ph_cut_col%d_row%d",col,row),Form("th_vs_ph_cut_col%d_row%d",col,row),1000,-0.03,0.03,1000,-0.045,0.045);
                    TCut sieveCut=TCut(generalcut.Data())+TCut(Form("hcut_R_%d_%d_%d", FoilID, col, row));
                    chain->Project(TargetThPhSieve->GetName(),Form("%s.gold.th:%s.gold.ph",HRS.Data(),HRS.Data()),sieveCut);

                    SievePointingCanv->cd(2)->cd(1);
                    TargetThPhSieve->Draw("zcol");
                    cutg->Draw("same");

                    auto projectx=TargetThPhSieve->ProjectionX();
                    auto projecty=TargetThPhSieve->ProjectionY();

                    projectx->GetXaxis()->SetRangeUser(projectx->GetXaxis()->GetBinCenter(projectx->GetMaximumBin())-0.003,projectx->GetXaxis()->GetBinCenter(projectx->GetMaximumBin())+0.003);
                    projectx->SetTitle(Form("Sieve_Phi_Col%d_Row%d",col,row));
                    projecty->GetXaxis()->SetRangeUser(projecty->GetXaxis()->GetBinCenter(projecty->GetMaximumBin())-0.004,projecty->GetXaxis()->GetBinCenter(projecty->GetMaximumBin())+0.004);
                    projecty->SetTitle(Form("Sieve_Theta_Col%d_Row%d",col,row));


                    SievePointingCanv->cd(2)->cd(2);
                    projectx->Fit("gaus","RQ");
                    thetaValueArray[col][row]=projectx->GetFunction("gaus")->GetParameter(1);
                    TLine *thetaValLine =new TLine(thetaValueArray[col][row],0,thetaValueArray[col][row],projectx->GetFunction("gaus")->GetParameter(0));
                    TLatex *thetaValtxt =new TLatex(thetaValueArray[col][row],projectx->GetFunction("gaus")->GetParameter(0),Form("theta %1.5f",thetaValueArray[col][row]));
                    projectx->Draw();
                    thetaValLine->Draw("same");
                    thetaValtxt->Draw("same");

                    SievePointingCanv->cd(2)->cd(3);
                    projecty->Fit("gaus","RQ");
                    phiValueArray[col][row] = projecty->GetFunction("gaus")->GetParameter(1);

                    TLine *phiValLine = new TLine(phiValueArray[col][row],0, phiValueArray[col][row], projecty->GetFunction("gaus")->GetParameter(0));
                    TLatex *phiValTxt = new TLatex(phiValueArray[col][row], projecty->GetFunction("gaus")->GetParameter(0),Form("Phi: %1.5f",phiValueArray[col][row]));

                    projecty->Draw();
                    phiValLine->Draw("same");
                    phiValTxt->Draw("same");

                    SievePointingCanv->cd(1);
					auto fitFunction=SpectroCrystalFit_H2O(SieveMomentumArray[col][row]);
					MomentumFitParArray[col][row]=new double[10];
					fitFunction->GetParameters(MomentumFitParArray[col][row]);


					//calculate the sieve angles and also the sieve theta phi
                    double  DeltaE=MomentumFitParArray[col][row][1]-MomentumFitParArray[col][row][6];
                    double scatteredAngle=GetPointingAngle(DeltaE,getBeamE(runID,chain));
                    // take the angle on the target coordination and get the central sieve position
                    scatteredAngleArray[col][row]=scatteredAngle;
                    scatteredAngleEnergyDiff[col][row]=DeltaE;


                    SieveMomentumArray[col][row]->Draw();
					fitFunction->Draw("same");
					//plot the values on the canvas, write the data
					TLatex *groundValtxt = new TLatex(MomentumFitParArray[col][row][1] + 2*MomentumFitParArray[col][row][2],MomentumFitParArray[col][row][0],Form("P=%2.5fGeV #pm %3.1fMeV",MomentumFitParArray[col][row][1],1000*fitFunction->GetParError(1)));
					groundValtxt->SetTextSize(0.05);
                    groundValtxt->SetTextAlign(12);
                    groundValtxt->SetTextColor(2);
                    groundValtxt->Draw("same");

                    TLine *groudposLine=new TLine(MomentumFitParArray[col][row][1],0,MomentumFitParArray[col][row][1],MomentumFitParArray[col][row][0]*1.1);
                    groudposLine->SetLineColor(3);
                    groudposLine->SetLineWidth(2);
                    groudposLine->Draw("same");

                    TLatex *firstValtxt = new TLatex(MomentumFitParArray[col][row][6] + 2*MomentumFitParArray[col][row][7],MomentumFitParArray[col][row][5],Form("P=%2.5fGeV #pm %3.1fMeV",MomentumFitParArray[col][row][6],1000*fitFunction->GetParError(6)));
                    firstValtxt->SetTextSize(0.05);
                    firstValtxt->SetTextAlign(12);
                    firstValtxt->SetTextColor(2);
                    firstValtxt->Draw("same");

                    TLine *firstposLine=new TLine(MomentumFitParArray[col][row][6],0,MomentumFitParArray[col][row][6],MomentumFitParArray[col][row][5]*2.0);
                    firstposLine->SetLineColor(3);
                    firstposLine->SetLineWidth(2);
                    firstposLine->Draw("same");

                    auto deltaE = MomentumFitParArray[col][row][1]-MomentumFitParArray[col][row][6];
                    auto deltaErr = TMath::Sqrt(fitFunction->GetParError(1)*fitFunction->GetParError(1)+fitFunction->GetParError(6)*fitFunction->GetParError(6));

                    TLatex *midValuetxt = new TLatex(0.5*MomentumFitParArray[col][row][6]+0.5*MomentumFitParArray[col][row][1],0.5*(MomentumFitParArray[col][row][0]+MomentumFitParArray[col][row][5]),Form("#DeltaP=%2.3f #pm %1.3fMeV",deltaE*1000.0,1000.0*deltaErr));
                    midValuetxt->SetTextSize(0.05);
                    midValuetxt->SetTextAlign(12);
                    midValuetxt->SetTextColor(2);
                    midValuetxt->Draw("same");

                    // calculate the final result of the data
                    TPaveText *pt = new TPaveText(0.1,0.8,0.4,0.9,"NDC");

                    std::cout<<"Get the HRS angle from ("<<col<<","<<row<<std::endl;
                    double a = scatteredAngleArray[col][row];
                    double theta_tg=thetaValueArray[col][row];
                    double phi_tg=phiValueArray[col][row];

                    // TODO fill the theta-phi on target
                    Bool_t UseMeasuredTargThPh= false;
                    if (UseMeasuredTargThPh){
                        // use the measured average for the theta and phi
                        GLobalSovler_phi_tg=theta_tg;
                        GLobalSovler_theta_tg=phi_tg;
                    }else{
                        // use the theoretical Value and the survey for the theta-phi
                        auto sievePos = _getSieveAngleOnTargCoor(getBPM(runID),HRS,col,row);
                        GLobalSovler_theta_tg=sievePos.X();
                        GLobalSovler_phi_tg=sievePos.Y();
                    }

                    GLobalSovler_scatteredAngle=a*TMath::Pi()/180.0;

                    // sovle the function and get the angle
                    ROOT::Math::Functor1D f1D(&GetHRSAngle);
                    ROOT::Math::RootFinder rfb(ROOT::Math::RootFinder::kBRENT);
                    rfb.SetFunction(f1D,0.0,6.0*TMath::Pi()/180.0);
                    rfb.Solve();
                    std::cout << rfb.Root()*180.0/TMath::Pi() << std::endl;
                    double hrsPointingAngle = rfb.Root()*180.0/TMath::Pi();
                    HRSAngleArray[col][row] = hrsPointingAngle;
                    pt->AddText(Form("#DeltaDp :%1.3f MeV Scatter Angle:%1.3f, HRS Angle:%1.3f",1000.0*deltaE,GetPointingAngle(deltaE,getBeamE(runID,chain)),hrsPointingAngle));
                    pt->Draw("same");

                    SievePointingCanv->Update();
					SievePointingCanv->Print(Form("Pointing_run%d.pdf",runID),"pdf");
				}
			}
		}
	}

    // draw the Momentum of the each Sieve holes
    TH2F *TargetGroundPSieve=(TH2F *)gROOT->FindObject(Form("Sieve_Ground_P"));
    if(TargetGroundPSieve) TargetGroundPSieve->Delete();
    TargetGroundPSieve=new TH2F(Form("Sieve_Ground_P"),Form("Sieve_Ground_P"),1000,-0.03,0.03,1000,-0.045,0.045);

    TCanvas *targetGroundPCanv = new TCanvas(Form("Sieve_ground_p"),Form("Sieve_ground_p"),1960,1080);
    targetGroundPCanv->Draw();
    TargetGroundPSieve->Draw("same");

    for (int16_t col = 0; col < NSieveCol; col++) {
        for (int16_t row = 0; row < NSieveRow; row++) {
            auto cutg = (TCutG *) gROOT->FindObject(Form("hcut_R_%d_%d_%d", FoilID, col, row));
            if (cutg) {
                cutg->Draw("same");
                TH2F *TargetThPhSieve=(TH2F *)gROOT->FindObject(Form("th_vs_ph_cut_col%d_row%d",col,row));
				if(TargetThPhSieve) TargetThPhSieve->Delete();
				TargetThPhSieve=new TH2F(Form("th_vs_ph_cut_col%d_row%d",col,row),Form("th_vs_ph_cut_col%d_row%d",col,row),1000,-0.03,0.03,1000,-0.045,0.045);
				TCut sieveCut=TCut(generalcut.Data())+TCut(Form("hcut_R_%d_%d_%d", FoilID, col, row));
				//project the sieve
				chain->Project(TargetThPhSieve->GetName(),Form("%s.gold.th:%s.gold.ph",HRS.Data(),HRS.Data()),sieveCut);

				auto projectx=TargetThPhSieve->ProjectionX();
				auto projecty=TargetThPhSieve->ProjectionY();

				// get central sieve position
				auto maximumX=projectx->GetXaxis()->GetBinCenter(projectx->GetMaximumBin());
				auto maximumY=projecty->GetXaxis()->GetBinCenter(projecty->GetMaximumBin());

				TLatex *txt =new TLatex(maximumX,maximumY,Form("%1.5f",MomentumFitParArray[col][row][1]));
				txt->SetTextSize(0.03);
				txt ->Draw("same");
            }
        }
    }
    targetGroundPCanv->Update();
    targetGroundPCanv->Print(Form("Pointing_run%d.pdf",runID),"pdf");

    TH2F *TargetFirstPSieve=(TH2F *)gROOT->FindObject(Form("Sieve_First_P_h"));
    if(TargetFirstPSieve) TargetFirstPSieve->Delete();
    TargetFirstPSieve=new TH2F(Form("Sieve_First_P_h"),Form("Sieve_First_P_h"),1000,-0.03,0.03,1000,-0.045,0.045);
    TCanvas *targetFirstPCanv = new TCanvas(Form("Sieve_first_p_canv"),Form("Sieve_first_p_canv"),1960,1080);
    targetFirstPCanv->Draw();
    TargetFirstPSieve->Draw();

    for (int16_t col = 0; col < NSieveCol; col++) {
        for (int16_t row = 0; row < NSieveRow; row++) {
            auto cutg = (TCutG *) gROOT->FindObject(Form("hcut_R_%d_%d_%d", FoilID, col, row));
            if (cutg) {
                cutg->Draw("same");
                TH2F *TargetThPhSieve=(TH2F *)gROOT->FindObject(Form("th_vs_ph_cut_col%d_row%d",col,row));
                if(TargetThPhSieve) TargetThPhSieve->Delete();
                TargetThPhSieve=new TH2F(Form("th_vs_ph_cut_col%d_row%d",col,row),Form("th_vs_ph_cut_col%d_row%d",col,row),1000,-0.03,0.03,1000,-0.045,0.045);
                TCut sieveCut=TCut(generalcut.Data())+TCut(Form("hcut_R_%d_%d_%d", FoilID, col, row));
                //project the sieve
                chain->Project(TargetThPhSieve->GetName(),Form("%s.gold.th:%s.gold.ph",HRS.Data(),HRS.Data()),sieveCut);

                auto projectx=TargetThPhSieve->ProjectionX();
                auto projecty=TargetThPhSieve->ProjectionY();

                // get central sieve position
                auto maximumX=projectx->GetXaxis()->GetBinCenter(projectx->GetMaximumBin());
                auto maximumY=projecty->GetXaxis()->GetBinCenter(projecty->GetMaximumBin());

                TLatex *txt =new TLatex(maximumX,maximumY,Form("%1.5f",MomentumFitParArray[col][row][6]));
                txt->SetTextSize(0.03);
                txt ->Draw("same");
            }
        }
    }
    targetFirstPCanv->Update();
    targetFirstPCanv->Print(Form("Pointing_run%d.pdf",runID),"pdf");

    //plot all the HRS angle on canvas


    // plot the angle on the canvas
    std::map<int,TH1F *>sievePointColh;
    std::map<int,TH1F *>sievePointRowh;

    // buffer all the data in the canvas
    for (int16_t col = 0; col < NSieveCol; col++) {
        for (int16_t row = 0; row < NSieveRow; row++) {
            // write the data into the histograme


            if ((HRSAngleArray.find(col)!=HRSAngleArray.end())&&(HRSAngleArray[col].find(row)!=HRSAngleArray[col].end())){
                if (!(sievePointColh.find(col)!=sievePointColh.end())){
                    sievePointColh[col] = new TH1F(Form("Sieve Col %d HRS Angle",col),Form("Sieve Col %d HRS Angle",col),10,0,10);
                    sievePointColh[col]->GetYaxis()->SetRangeUser(4.5,5.1);
                    sievePointColh[col]->SetMarkerColor(46);
                    sievePointColh[col]->SetLineWidth(2);
                    sievePointColh[col]->SetLineColor(46);
                    sievePointColh[col]->SetMarkerStyle(20);
                }

                sievePointColh[col]->SetBinContent(row+1,HRSAngleArray[col][row]);
                sievePointColh[col]->SetBinError(row+1,0.02);

                if (!(sievePointRowh.find(row) != sievePointRowh.end())){
                    sievePointRowh[row] =new TH1F(Form("Sieve Row %d HRS Angle",row),Form("Sieve Row %d HRS Angle",col),15,0,15);
                    sievePointRowh[row]->GetYaxis()->SetRangeUser(4.5,5.1);
                    sievePointRowh[row]->SetMarkerColor(46);
                    sievePointRowh[row]->SetLineWidth(2);
                    sievePointRowh[row]->SetLineColor(46);
                    sievePointRowh[row]->SetMarkerStyle(20);
                }
                sievePointRowh[row]->SetBinContent(col+1,HRSAngleArray[col][row]);
                sievePointRowh[row]->SetBinError(col+1,0.02);
            }
        }
    }
    //plot the HRS angle on the canvas
    TCanvas *pointAngColPart1Canv=new TCanvas(Form("pointAngColPart1"),Form("pointAngColPart1"),1960,1080);
    pointAngColPart1Canv->Divide(4,3);
    pointAngColPart1Canv->Draw();
    for (auto colIter = 3; colIter <=12; colIter ++){
        if(sievePointColh.find(colIter)!=sievePointColh.end()){
            pointAngColPart1Canv->cd(colIter-2);
            gStyle->SetEndErrorSize(5);
            gStyle->SetErrorX(0.);
            sievePointColh[colIter]->Draw("E1");

            // plot the reference line
            if(1) {
                double survey = 4.818;
                double final_pointVal = 4.765;
                if (HRS == "R"){
                    survey=4.822;
                    final_pointVal = 4.747;
                }

                // draw
                TLine *line=new TLine(0.5,survey,9.5,survey);
                line->SetLineWidth(2);
                line->SetLineColor(3);
                line->Draw("same");


                TLatex *text1= new TLatex(6,survey+0.07,Form("Survey: %1.3f + 0.06#circ",survey));
                text1->SetTextSize(0.04);
                text1->SetTextColorAlpha(93,0.476);

                TLatex *text2= new TLatex(6,survey-0.09,Form("Survey: %1.3f - 0.06#circ",survey));
                text2->SetTextSize(0.04);
                text2->SetTextColorAlpha(93,0.476);

                TLatex *surveyValtext= new TLatex(6,survey+0.001,Form("Survey: %1.3f #pm 0.06#circ",survey));
                surveyValtext->SetTextSize(0.04);
                surveyValtext->SetTextColorAlpha(3,0.476);


                text1->Draw("same");
                text2->Draw("same");
                surveyValtext->Draw("same");

                TLine *line1=new TLine(0.5,survey+0.06,9.5,survey+0.06);
                line1->SetLineWidth(2);
                line1->SetLineColor(93);
                line1->Draw("same");

                TLine *line2=new TLine(0.5,survey-0.06,9.5,survey-0.06);
                line2->SetLineWidth(2);
                line2->SetLineColor(93);
                line2->Draw("same");


                // line with the measurement value
                TLine *fvalueLine=new TLine(0.5,final_pointVal,9.5,final_pointVal);
                fvalueLine->SetLineWidth(2);
                fvalueLine->SetLineColor(4);
                fvalueLine->Draw("same");

                TLatex *fvaluetxt= new TLatex(6,final_pointVal+0.001,Form("Measure: %1.3f #pm 0.02#circ",final_pointVal));
                fvaluetxt->SetTextSize(0.04);
                fvaluetxt->SetTextColor(4);
                fvaluetxt->Draw("same");

            }
        }
    }
    pointAngColPart1Canv->Update();
    pointAngColPart1Canv->Print(Form("Pointing_run%d.pdf",runID),"pdf");

    TCanvas *canv_temp=new TCanvas(Form("canv_temp"),Form("canv_temp"),1960,1080);
    for (auto colIter = 3; colIter <=12; colIter ++){
        if(sievePointColh.find(colIter)!=sievePointColh.end()){

            canv_temp->Clear();
            canv_temp->cd();

            gStyle->SetEndErrorSize(5);
            gStyle->SetErrorX(0.);
            sievePointColh[colIter]->Draw("E1");

            // plot the reference line
            if(1) {
                double survey = 4.818;
                double final_pointVal = 4.765;
                if (HRS == "R") {
                    survey = 4.822;
                    final_pointVal = 4.747;
                }

                // draw
                TLine *line = new TLine(0.5, survey, 9.5, survey);
                line->SetLineWidth(2);
                line->SetLineColor(3);
                line->Draw("same");


                TLatex *text1 = new TLatex(6, survey + 0.07, Form("Survey: %1.3f + 0.06#circ", survey));
                text1->SetTextSize(0.04);
                text1->SetTextColorAlpha(93, 0.476);

                TLatex *text2 = new TLatex(6, survey - 0.09, Form("Survey: %1.3f - 0.06#circ", survey));
                text2->SetTextSize(0.04);
                text2->SetTextColorAlpha(93, 0.476);

                TLatex *surveyValtext = new TLatex(6, survey + 0.001, Form("Survey: %1.3f #pm 0.06#circ", survey));
                surveyValtext->SetTextSize(0.04);
                surveyValtext->SetTextColorAlpha(3, 0.476);


                text1->Draw("same");
                text2->Draw("same");
                surveyValtext->Draw("same");

                TLine *line1 = new TLine(0.5, survey + 0.06, 9.5, survey + 0.06);
                line1->SetLineWidth(2);
                line1->SetLineColor(93);
                line1->Draw("same");

                TLine *line2 = new TLine(0.5, survey - 0.06, 9.5, survey - 0.06);
                line2->SetLineWidth(2);
                line2->SetLineColor(93);
                line2->Draw("same");


                // line with the measurement value
                TLine *fvalueLine = new TLine(0.5, final_pointVal, 9.5, final_pointVal);
                fvalueLine->SetLineWidth(2);
                fvalueLine->SetLineColor(4);
                fvalueLine->Draw("same");

                TLatex *fvaluetxt = new TLatex(6, final_pointVal + 0.001,
                                               Form("Measure: %1.3f #pm 0.02#circ", final_pointVal));
                fvaluetxt->SetTextSize(0.04);
                fvaluetxt->SetTextColor(4);
                fvaluetxt->Draw("same");
            }
            canv_temp->Print(Form("Pointing_run%d.pdf",runID),"pdf");
        }
    }
    // write the canvas seperately


    TCanvas *pointAngRowPart1Canv=new TCanvas(Form("pointAngRowPart1Canv"),Form("pointAngRowPart1Canv"),1960,1080);
    pointAngRowPart1Canv->Divide(3,3);
    pointAngRowPart1Canv->Draw();
    for (auto rowIter = 0; rowIter<=8; rowIter++){
        if (sievePointRowh.find(rowIter)!=sievePointRowh.end()){
            pointAngRowPart1Canv->cd(rowIter+1);
            sievePointRowh[rowIter]->Draw("E1");
        }
    }
    pointAngRowPart1Canv->Update();
    pointAngRowPart1Canv->Print(Form("Pointing_run%d.pdf",runID),"pdf");



    TH2F *momentumP1=new TH2F("momentumP1","MomentmP0",100,0,100,1000,2.1,2.2);
	TH2F *momentumP2=new TH2F("momentumP2","MomentmP1",100,0,100,1000,2.1,2.2);
	for (auto iter_col=MomentumFitParArray.begin();iter_col!=MomentumFitParArray.end();iter_col++){
		for (auto iter_row=iter_col->second.begin();iter_row!=iter_col->second.end();iter_row++){
			momentumP1->Fill(iter_col->first*NSieveRow+iter_row->first,iter_row->second[1]);
			momentumP2->Fill(iter_col->first*NSieveRow+iter_row->first,iter_row->second[6]);
		}
	}

	TCanvas *canvastest2=new TCanvas("canvs2","sasas2",1000,1000);
	canvastest2->cd();

	momentumP1->SetLineWidth(2);
	momentumP1->SetMarkerStyle(20);
	momentumP1->Draw();

	momentumP2->SetLineWidth(2);
    momentumP2->SetMarkerStyle(20);
    momentumP2->SetMarkerColor(kRed);
    momentumP2->SetLineColor(kRed);
    momentumP2->Draw("same");
	canvastest2->Update();

	canvastest2->Print(Form("Pointing_run%d.pdf)",runID),"pdf");
}


Int_t OpticsGraphicCutProH20(UInt_t runID,UInt_t maximumFileas=1,TString folder="/home/newdriver/Storage/Research/CRex_Experiment/RasterReplay/Replay/Result") {
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
                if(split>maximumFileas) break;
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
                if(split > maximumFileas) break;
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

	// link the root tree and check which HRS we are working on
	TChain *chain = (TChain *) gROOT->FindObject("T");
	TString HRS("R");
	TString filename(chain->GetFile()->GetName());
	if (filename.Contains("RHRS")) {
	} else if (filename.Contains("LHRS")) {
		HRS = "L";
	}

    TH1F *eventIDhh=new TH1F("eventID","eventID",800000,0,800000);
    chain->Project(eventIDhh->GetName(),"fEvtHdr.fRun");
    int eventID=int(eventIDhh->GetMean());
    TFile *f1=new TFile(Form("./PointingCheck/result/Pointing_%d.root",eventID),"RECREATE");
    assert(f1);
	// try to extract the hall prob if this is LHRS
	double CentralP;
	if (HRS == "L") {
		TH1F *HallProbHH = new TH1F("HallLProb", "HallLProb", 1000, -1, 0);
		chain->Project(HallProbHH->GetName(), "HacL_D_LS450_FLD_DATA",
				generalcut.Data());
		CentralP = std::abs((HallProbHH->GetMean()) * 0.95282 / 0.33930);
		std::cout << "CentralMomentum is (LHRS) for Hall Probe::" << (CentralP)
				<< std::endl;
	} else {
		//HacR_D1_NMR_SIG
		TH1F *HallR_NMR = new TH1F("HallR_NMR", "HallR_NMR", 1000, 0.7, 0.9);
		chain->Project(HallR_NMR->GetName(), "HacR_D1_NMR_SIG",
				generalcut.Data());
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
	} else
		SieveRecCanvas = new TCanvas("SieveRecCanvas", "Projection Canvas",
				1000, 1000);

	SieveRecCanvas->Divide(1, 3);   // on the third line, will display the bpm correction term

	SieveRecCanvas->cd(2)->Divide(4, 1);
    SieveRecCanvas->cd(3)->Divide(4,1);
//    SieveRecCanvas->cd(4)->Divide(4,1);
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
	TH1F *momentum=new TH1F(Form("H2O gold.p run%d",eventID),Form("H2O gold.p run%d",eventID),500,2.1,2.2);
	chain->Project(momentum->GetName(),Form("%s.gold.dp*%f+%f",HRS.Data(),CentralP,CentralP),Form("%s && %s",generalcut.Data(),cutg->GetName()));
	// get the maximum bin, this should be the first excited states
	auto CGroundp=momentum->GetXaxis()->GetBinCenter(momentum->GetMaximumBin());

	momentum->GetXaxis()->SetRangeUser(CGroundp-0.02,CGroundp+0.0044*2);
	momentum->GetXaxis()->SetTitle(Form("%s.gold.dp*%f+%f",HRS.Data(),CentralP,CentralP));
	momentum->GetYaxis()->SetTitle("#");

	momentum->Draw();
	double_t fgroudGausPar[3];
	double_t ffirstGuasPar[3];
	TF1 *fgroudGaus=new TF1("groudstatesgaus","gaus",CGroundp-0.0005,CGroundp+0.0005);
	momentum->Fit("groudstatesgaus","RQ0","ep",fgroudGaus->GetXmin(),fgroudGaus->GetXmax());
	fgroudGaus->Draw("same");
	fgroudGaus->GetParameters(fgroudGausPar);

	TH1F *test=(TH1F *)momentum->Clone("fitTest");
	test->GetXaxis()->SetRangeUser(momentum->GetXaxis()->GetXmin(),fgroudGausPar[1]-10*fgroudGausPar[2]);

	auto C1stp=test->GetXaxis()->GetBinCenter(test->GetMaximumBin());

	TF1 *ffirstGuas=new TF1 ("firststatesgaus","gaus",C1stp-0.0015,C1stp+0.00155);
	momentum->Fit("firststatesgaus","RQ0","ep",ffirstGuas->GetXmin(),ffirstGuas->GetXmax());
	ffirstGuas->Draw("same");
	ffirstGuas->GetParameters(ffirstGuasPar);


	// change the gause fit to cristal ball
	double_t fgroundCrystalballPar[5];
	TF1 *fgroundCrystalball=new TF1("fgroundCrystal","crystalball",fgroudGausPar[1]-0.0030,fgroudGaus->GetXmax()+0.0003);
	fgroundCrystalball->SetParameters(fgroudGausPar[0],fgroudGausPar[1],fgroudGausPar[2],1.64,1.1615);
	momentum->Fit("fgroundCrystal","RQ0","same",fgroundCrystalball->GetXmin(),fgroundCrystalball->GetXmax());
	fgroundCrystalball->GetParameters(fgroundCrystalballPar);

	//fgroundCrystalball->Draw("same");

	double_t ffirstCrystalPar[5];
	TF1 *ffirstCrystal=new TF1("ffirstCrystal","crystalball",ffirstGuasPar[1]-0.0025,ffirstGuas->GetXmax());
	ffirstCrystal->SetParameters(ffirstGuasPar[0],ffirstGuasPar[1],ffirstGuasPar[2],1.64,1.1615);
	momentum->Fit("ffirstCrystal","RQ0","ep",ffirstCrystal->GetXmin(),ffirstCrystal->GetXmax());
	ffirstCrystal->GetParameters(ffirstCrystalPar);
	//	ffirstCrystal->Draw("same");
	// fit together
	double_t fCrystalMomentumPar[10];
	TF1 *fCrystalMomentum=new TF1("fCrystalMomentum","crystalball(0)+crystalball(5)",ffirstCrystal->GetXmin(),fgroundCrystalball->GetXmax());
	std::copy(fgroundCrystalballPar,fgroundCrystalballPar+5,fCrystalMomentumPar);
	std::copy(ffirstCrystalPar,ffirstCrystalPar+5,fCrystalMomentumPar+5);
	fCrystalMomentum->SetParameters(fCrystalMomentumPar);
	momentum->Fit("fCrystalMomentum","R","",fCrystalMomentum->GetXmin(),fCrystalMomentum->GetXmax());
	fCrystalMomentum->Draw("same");
	fCrystalMomentum->GetParameters(fCrystalMomentumPar);

	SieveRecCanvas->Update();
	// plot the reconstrcution peak
	TLine *groudposLine=new TLine(fCrystalMomentumPar[1],0,fCrystalMomentumPar[1],fgroudGausPar[0]*1.1);
	groudposLine->SetLineColor(3);
	groudposLine->SetLineWidth(2);
	groudposLine->Draw("same");

	TLine *firstposLine=new TLine(fCrystalMomentumPar[6],0,fCrystalMomentumPar[6],ffirstGuasPar[0]*2.0);
	firstposLine->SetLineColor(3);
	firstposLine->SetLineWidth(2);
	firstposLine->Draw("same");

	TPaveText *pt = new TPaveText(0.1,0.1,0.4,0.4,"NDC");
	double_t deltaE=fCrystalMomentumPar[1]-fCrystalMomentumPar[6];
	double_t deltaErr=TMath::Sqrt( (fCrystalMomentum->GetParError(1))*(fCrystalMomentum->GetParError(1))+(fCrystalMomentum->GetParError(6))*(fCrystalMomentum->GetParError(6)));

	double_t HRSAngle_Final=0.0;
//	if(HallProbHH->GetEntries()!=0)
	{
	    //Step-1 calculate the correction
        double HRSBPMCorrection=0.0;

        SieveRecCanvas->cd(3)->cd(3);
        TString imgfname=Form("/home/newdriver/Learning/GeneralScripts/halog/result/BeamE%d.jpg",eventID);
        if(!gSystem->AccessPathName(imgfname.Data())){
            TImage *img=TImage::Open(imgfname.Data());
            img->Draw();
        } else{
            std::cout << "\033[1;33m [Warning]\033[0m Missing Beam E image:"<<imgfname.Data()<< std::endl;
        }
        SieveRecCanvas->cd(1);
		// get the error
		double DeltaEMax=deltaE+deltaErr;
		double DeltaE200=deltaE+200.0/1000000.0;
		double DeltaEMin=deltaE-deltaErr;
		double AngleMax=GetPointingAngle(DeltaEMax,getBeamE(eventID,chain));
		double AngleMin=GetPointingAngle(DeltaEMin,getBeamE(eventID,chain));
		double Angletemp=GetPointingAngle(deltaE,getBeamE(eventID,chain));
		double errorAngle=0;
		if(abs(AngleMax-Angletemp) >abs(Angletemp-AngleMin) ){
			errorAngle=abs(AngleMax-Angletemp);
		}else{
			errorAngle=abs(AngleMin-Angletemp);
		}

		double errorAngle1=abs(GetPointingAngle(DeltaE200,getBeamE(eventID,chain)))-Angletemp;
		//add the correction
		double HRSAngle=GetPointingAngle(deltaE,getBeamE(eventID,chain));
		HRSAngle_Final=HRSAngle;
		pt->AddText(Form("#DeltaDp :%1.3f MeV (%1.3f#pm%1.3f#pm%1.3f Degree)",1000.0*deltaE,HRSAngle,errorAngle,errorAngle1));
	}

	pt->Draw("same");

	TLatex *t1 = new TLatex(fgroudGausPar[1] + 2 * fgroudGausPar[2],fgroudGausPar[0], Form("P=%2.5fGeV #pm %1.3fMeV", fCrystalMomentumPar[1],1000.0*(fCrystalMomentum->GetParError(1))));
	t1->SetTextSize(0.1);
	t1->SetTextAlign(12);
	t1->SetTextColor(2);
	t1->Draw("same");

	TLatex *t2 = new TLatex(ffirstGuasPar[1]+ffirstGuasPar[2]*0.5,fCrystalMomentumPar[5], Form("P=%2.5fGeV #pm %1.3fMeV", fCrystalMomentumPar[6],1000.0*(fCrystalMomentum->GetParError(6))));
	t2->SetTextSize(0.1);
	t2->SetTextAlign(12);
	t2->SetTextColor(2);
	t2->Draw("same");

	TLatex *t3 = new TLatex((ffirstGuasPar[1]+fgroudGausPar[1])/2.0,(fCrystalMomentumPar[5]+fgroudGausPar[0])*0.5, Form("#DeltaP=%2.3f #pm %1.3fMeV", 1000.0*deltaE,1000.0*deltaErr));
	t3->SetTextSize(0.1);
	t3->SetTextAlign(12);
	t3->SetTextColor(2);
	t3->Draw("same");

	//plot the bigger plot for first excited states
	SieveRecCanvas->cd(2)->cd(3);
	TH1F *groundStats=(TH1F *)momentum->Clone("H2O p.ground");
	groundStats->GetXaxis()->SetRangeUser(fCrystalMomentumPar[1]-0.002,fCrystalMomentumPar[1]+0.002);
	groundStats->Draw();
	groudposLine->Draw("same");

	SieveRecCanvas->cd(2)->cd(4);
	TH1F *firstStats=(TH1F *)momentum->Clone("H2O p.first");
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

	//ceate the data for the target focal plane variables and also save it to canvas
    {
        SieveRecCanvas->cd(3)->cd(1);
        TH1F *tgThetah=new TH1F(Form("tg_th_%d",eventID),Form("tg_th_%d",eventID),1000,-0.045,0.045);
        chain->Project(tgThetah->GetName(),Form("%s.gold.th",HRS.Data()),Form("%s && %s",generalcut.Data(),cutg->GetName()));
        tgThetah->GetXaxis()->SetRangeUser(tgThetah->GetMaximum()-0.005,tgThetah->GetMaximum()+0.005);
        tgThetah->Fit("gaus");
        tgThetah->Draw();
        auto thetaFunc=tgThetah->GetFunction("gaus");
        TLatex *txTheta=new TLatex(thetaFunc->GetParameter(1),thetaFunc->GetParameter(0),Form("theta:%f",thetaFunc->GetParameter(1)));
        txTheta->Draw("same");

        SieveRecCanvas->cd(3)->cd(2);
        TH1F *tgPhih=new TH1F(Form("tg_ph_%d",eventID),Form("tg_ph_%d",eventID),1000,-0.045,0.045);
        chain->Project(tgPhih->GetName(),Form("%s.gold.ph",HRS.Data()),Form("%s && %s",generalcut.Data(),cutg->GetName()));
        tgPhih->GetXaxis()->SetRangeUser(tgPhih->GetMaximum()-0.005,tgPhih->GetMaximum()+0.005);
        tgPhih->Fit("gaus");
        tgPhih->Draw();
        auto phiFunc=tgPhih->GetFunction("gaus");
        TLatex *txPhi=new TLatex(phiFunc->GetParameter(1),phiFunc->GetParameter(0),Form("phi:%f",phiFunc->GetParameter(1)));
        txPhi->Draw("same");

        // create file and write the data into it
        std::ofstream txtfileio("./FinalData/TargetVar/tg_variableList.txt",std::ofstream::app);
        auto writeString=Form("%5d  %1.5f   %1.5f   %1.5f   %1.5f   %1.5f   %1.5f",eventID,thetaFunc->GetParameter(1),phiFunc->GetParameter(1),thetaFunc->GetParameter(1)*180.0/3.141592654,phiFunc->GetParameter(1)*180.0/3.141592654,HRSAngle_Final,getP0(HRSAngle_Final,thetaFunc->GetParameter(1),phiFunc->GetParameter(1)));
        txtfileio << writeString<<std::endl;
        txtfileio.close();
    }

    SieveRecCanvas->SaveAs(Form("./PointingCheck/result/Pointing_%d.jpg",eventID));
	SieveRecCanvas->Write();
	hSieveHole->Delete();
	cutg->Delete();
}

