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

//////////////////////////////////////////////////////////////////////////////
// Work Directory
// cut options
// Need to change
//////////////////////////////////////////////////////////////////////////////
TString prepcut;
TString generalcut;

//CRex C-12
TString generalcutR="R.tr.n==1 && R.vdc.u1.nclust==1&& R.vdc.v1.nclust==1 && R.vdc.u2.nclust==1 && R.vdc.v2.nclust==1 && R.gold.dp<1 && R.gold.dp > -0.1 && fEvtHdr.fEvtType==1";
TString generalcutL="L.tr.n==1 && L.vdc.u1.nclust==1&& L.vdc.v1.nclust==1 && L.vdc.u2.nclust==1 && L.vdc.v2.nclust==1 && L.gold.dp<1 && L.gold.dp > -0.1 && fEvtHdr.fEvtType==1";

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

    if (event==kButton1Down) {
        TH2 *h = (TH2*) select;

        gPad->GetCanvas()->FeedbackMode(kTRUE);

        // if the button is clicked
        //Rec the sieve pattern
        // get the mouse click position in histogram
        double_t x = (gPad->PadtoX(gPad->AbsPixeltoX(gPad->GetEventX())));
        double_t y = (gPad->PadtoY(gPad->AbsPixeltoY(gPad->GetEventY())));

        // before load to contour algorithm, get the more accurate center first
        // create new canvas
        TCanvas *SieveRecCanvas = (TCanvas*) gROOT->GetListOfCanvases()->FindObject("SieveRecCanvas");
        if(SieveRecCanvas){
            SieveRecCanvas->Clear();
            delete SieveRecCanvas->GetPrimitive("Projection");
        }else
            SieveRecCanvas = new TCanvas("SieveRecCanvas","Projection Canvas", 1000,1000);

        SieveRecCanvas->Divide(1,2);
        SieveRecCanvas->cd(1)->Divide(2,1);
        SieveRecCanvas->cd(2)->Divide(4,1);

        SieveRecCanvas->cd(1)->cd(2)->Divide(1,3);

        SieveRecCanvas->cd(1)->cd(2)->cd(1);
        //preCut
        TH2F *selectedSievePreCuthh = (TH2F *) gROOT->FindObject(
                "Sieve_Selected_th_ph_PreCut");
        if (selectedSievePreCuthh) {
            selectedSievePreCuthh->Clear();
        } else {
            selectedSievePreCuthh= new TH2F("Sieve_Selected_th_ph_PreCut",
                                            "Sieve_Selected_th_ph_PreCut", 100, h->GetXaxis()->GetXmin(),
                                            h->GetXaxis()->GetXmax(), 100, h->GetYaxis()->GetXmin(),
                                            h->GetYaxis()->GetXmax());
        }

        chain->Project(selectedSievePreCuthh->GetName(),
                       Form("%s.gold.th:%s.gold.ph", HRS.Data(), HRS.Data()),
                       Form("sqrt((%s.gold.th-%f)^2+ (%s.gold.ph-%f)^2)<0.003 && %s ",
                            HRS.Data(), y, HRS.Data(), x, generalcut.Data()));
        selectedSievePreCuthh->GetXaxis()->SetTitle(Form("%s.gold.ph",HRS.Data()));
        selectedSievePreCuthh->GetYaxis()->SetTitle(Form("%s.gold.th",HRS.Data()));
        //project to theta and phi, and start fit, get  more accurate position before pass to the counter
        auto projectxPreCut = selectedSievePreCuthh->ProjectionX();
        auto projectyPreCut = selectedSievePreCuthh->ProjectionY();
        selectedSievePreCuthh->Draw("zcol");

        SieveRecCanvas->cd(1)->cd(2)->cd(2);
        projectxPreCut->Draw();
        SieveRecCanvas->cd(1)->cd(2)->cd(3);
        projectyPreCut->Draw();
        //get the fit and update the position information
        projectxPreCut->Fit("gaus","","");
        projectyPreCut->Fit("gaus","","");

        // get the updated informations
        x=projectxPreCut->GetFunction("gaus")->GetParameter(1); //phi
        y=projectyPreCut->GetFunction("gaus")->GetParameter(1);  //theta


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
        SieveRecCanvas->cd(1)->cd(1);

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
        SieveRecCanvas->Update(); // update the canvas to let the pattern buffer in root

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

///
/// \param chain
/// \return
TVector getSieveThetaPhi(TChain *chain){

    gStyle->SetOptStat(0);

    TString HRS=getHRS(chain);

    if(HRS=="L"){
        generalcut=generalcutL;
    }else{
        generalcut=generalcutR;
    }

    TCanvas *mainPatternCanvas=(TCanvas *)gROOT->GetListOfCanvases()->FindObject("cutPro");
    if(!mainPatternCanvas){
        mainPatternCanvas=new TCanvas("cutProGetSieve","cutProGetSieve",1000,1000);
    }else{
        mainPatternCanvas->Clear();
    }

    //	TCanvas *mainPatternCanvas=new TCanvas("cut","cut",600,600);
    mainPatternCanvas->Draw();
    TH2F *TargetThPhHH=(TH2F *)gROOT->FindObject("th_vs_ph");
    if(TargetThPhHH) TargetThPhHH->Delete();
    TargetThPhHH=new TH2F("th_vs_ph","th_vs_ph",1000,-0.025,0.025,1000,-0.047,0.05);

    chain->Project(TargetThPhHH->GetName(),Form("%s.gold.th:%s.gold.ph",HRS.Data(),HRS.Data()),generalcut.Data());
    TargetThPhHH->Draw("zcol");
    mainPatternCanvas->SetGridx(10);
    mainPatternCanvas->SetGridy(10);


    mainPatternCanvas->Update();
    mainPatternCanvas->ToggleEventStatus();
    mainPatternCanvas->AddExec("ex", "DynamicCanvas()");

    //Draw
    TVector a;
    return  a;
}


///
/// \param runID
/// \param cutFile
/// \param folder
void cutProTemplate(UInt_t runID=22114,
                    TString cutFile ="/home/newdriver/Storage/Research/PRex_Workspace/PREX-MPDGEM/PRexScripts/Tools/PlotCut/Result/Final_Cut/RHRS_Cut20200518/RHRS/GroundMomCut/prexRHRS_21632_-1.root.FullCut.root",
                    TString folder ="/home/newdriver/Storage/Research/CRex_Experiment/RasterReplay/Replay/Result"){

    // check which HRS we are working on
    TString HRS="R";
    if (runID < 20000)HRS= "L";
    //load the cut profile
    if(HRS=="L"){
        generalcut=generalcutL;
    }else{
        generalcut=generalcutR;
    }

    if(!cutFile.EndsWith(".root")){
        std::cout<<"[warning]:: cut file is should end with .root, input-> "<<cutFile.Data()<<std::endl;
        exit(-1);
    }
    TFile *cutFileIO=new TFile(cutFile.Data(),"READ");
    if(cutFileIO->IsZombie()){
        std::cout<<"[ERROR]:: CAN NOT FIND CUT FILE \" "<<cutFile.Data()<<"\""<<std::endl;
        exit(-1);
    }

    TChain *chain=LoadRootFiles(runID,999,folder.Data());

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
    TargetThPhHH=new TH2F("th_vs_ph","th_vs_ph",1000,-0.025,0.025,1000,-0.047,0.05);

    chain->Project(TargetThPhHH->GetName(),Form("%s.gold.th:%s.gold.ph",HRS.Data(),HRS.Data()),generalcut.Data());
    TargetThPhHH->Draw("zcol");
    mainPatternCanvas->SetGridx(10);
    mainPatternCanvas->SetGridy(10);

    mainPatternCanvas->Update();
    mainPatternCanvas->ToggleEventStatus();
    mainPatternCanvas->AddExec("ex", "DynamicCanvas()");

}


void test(){
    auto chain = LoadRootFiles(22114,999,"/home/newdriver/Storage/Research/CRex_Experiment/RasterReplay/Replay/Result");
    getSieveThetaPhi(chain);
}