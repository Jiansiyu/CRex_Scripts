/*
 *
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
#include "TVector.h"
// used for create the folder if does not exist in the destintion folder
#include <boost/filesystem.hpp>
R__LOAD_LIBRARY(/usr/lib/x86_64-linux-gnu/libboost_filesystem.so)


int FoilID=0;

int col=3;
int row=1;
int row_min=0;
int row_count=10;
const UInt_t NSieveCol = 13;
const UInt_t NSieveRow = 7;


double clickedPosX=9999.99;
double clickedPosY=9999.99;

//////////////////////////////////////////////////////////////////////////////
// Work Directory
// cut options
// Need to change
//////////////////////////////////////////////////////////////////////////////
TString prepcut;
TString generalcut;

//CRex C-12
TString generalcutR="R.tr.n==1 && R.vdc.u1.nclust==1&& R.vdc.v1.nclust==1 && R.vdc.u2.nclust==1 && R.vdc.v2.nclust==1 && fEvtHdr.fEvtType==1";
///TString generalcutL="L.tr.n==1 && L.vdc.u1.nclust==1&& L.vdc.v1.nclust==1 && L.vdc.u2.nclust==1 && L.vdc.v2.nclust==1 && fEvtHdr.fEvtType==1";
TString generalcutL="";
// CRex Water
//TString generalcutR="R.tr.n==1 && R.vdc.u1.nclust==1&& R.vdc.v1.nclust==1 && R.vdc.u2.nclust==1 && R.vdc.v2.nclust==1 && fEvtHdr.fEvtType==1 && R.gold.p > 2.14 && R.gold.p < 2.2";
//TString generalcutL="L.tr.n==1 && L.vdc.u1.nclust==1&& L.vdc.v1.nclust==1 && L.vdc.u2.nclust==1 && L.vdc.v2.nclust==1 && L.gold.p > 2.14 && L.gold.p < 2.19";

//////////////////////////////////////////////////////////////////////////////
// Work Directory
//////////////////////////////////////////////////////////////////////////////
TString WorkDir = "/home/newdriver/Storage/Research/PRex_Workspace/PREX-MPDGEM/PRexScripts/Tools/PlotCut/Result/Cut20200701/water";


TString CutSuf = ".FullCut.root";
TString CutDescFileSufVertex = ".VertexCut.cut";
TString CutDescFileSufDp = ".DpCut.cut";
TString CutDescFileSufSieve = ".SieveCut.%d_%d.cut";
TString RootFileName;


//LHRS
int numberofSieveHoles[13]={0,0,0,5,6,5,5,6,5,5,4,3,2};
int minSieveHoles[13]=     {0,0,0,1,0,1,1,0,1,1,1,2,2};


//RHRS
//int numberofSieveHoles[13]={0,0,0,6,6,5,5,6,5,5,4,3,2};
//int minSieveHoles[13]=     {0,0,0,0,0,1,1,0,1,1,1,2,2};


///
/// \param name
/// \return
inline Bool_t IsFileExist (const std::string& name) {
    struct stat buffer;
    return (stat (name.c_str(), &buffer) == 0);
}

///
/// \param runID  the runID
/// \param maxFiles    the maximum number of files want to load
/// \param folder      the folder  that contain the data file, if the folder end with .root file,
///                    it will neglect the runID and maxFiles and just load this file
/// \return            root TChain
TChain *LoadRootFiles(UInt_t runID,UInt_t maxFiles=999,TString folder="/home/newdriver/Storage/Research/CRex_Experiment/RasterReplay/Replay/Result"){

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
/// \param chain
/// \return runID from the root file
UInt_t getRunID(TChain *chain){
    auto runID=int(chain->GetMaximum("fEvtHdr.fRun"));
    return runID;
}

TString getHRS(TChain *chain){
    if (getRunID(chain) > 20000){
        return  "R";
    }else{
        return  "L";
    }
}

/// Check the target variables on the target coordination system
/// \param runID
/// \param folder
void checkThetaPhi(UInt_t runID, TString folder="/home/newdriver/Storage/Research/PRex_Experiment/PRex_Replay/replay/Result") {

    auto chain = LoadRootFiles(runID,999,folder);

    TString HRS(getHRS(chain));
    //load the cut profile
    if (HRS == "L") {
        generalcut = generalcutL;
    } else {
        generalcut = generalcutR;
    }


    TCanvas *canv = new TCanvas("ThetaPhi_Diag", "ThetaPhi_Diag", 1960, 1080);
    canv->Divide(2, 1);
    canv->cd(2)->Divide(1, 2);
    canv->Draw();
    canv->cd(1);
    canv->cd(1)->SetGridx();
    canv->cd(1)->SetGridy();


    TH2F *TargetThPhHH = (TH2F *) gROOT->FindObject("th_vs_ph");
    if (TargetThPhHH) TargetThPhHH->Delete();
    TargetThPhHH = new TH2F("th_vs_ph", "th_vs_ph", 1000, -0.045, 0.045, 1000, -0.045, 0.045);

    chain->Project(TargetThPhHH->GetName(), Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()), generalcut.Data());
    TargetThPhHH->SetTitle(Form("Targ_Theta_Phi_run%d",runID));
    TargetThPhHH->Draw("zcol");

    //create plot
   auto  tgThetah= TargetThPhHH->ProjectionY("targ_theta_h");
   auto  tgPhih  = TargetThPhHH->ProjectionX("targ_phi_h");

   canv->cd(2)->cd(1);
   canv->cd(2)->cd(1)->SetGridx();
   canv->cd(2)->cd(1)->SetGridy();
   tgThetah->Draw();
   canv->cd(2)->cd(2);
   canv->cd(2)->cd(2)->SetGridx();
   canv->cd(2)->cd(2)->SetGridy();
   tgPhih->Draw();
   canv->Update();


}