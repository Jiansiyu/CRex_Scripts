//
// Created by Siyu Jian on 7/8/20.
// Used for Create Target Burn out plot
// Usage ::
//

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TChain.h>
#include <TCut.h>
#include <TH2F.h>
#include <TH1.h>
#include <TF1.h>
#include <TPaveText.h>
#include <map>
#include <vector>
#include <random>
#include <iostream>

#include <TObject.h>
#include "TMinuit.h"
#include <TFile.h>
#include <fstream>
#include <TSystem.h>
#include <TApplication.h>
#include <TLatex.h>
#include "TGraphErrors.h"
#include "TColor.h"
#include "TLegend.h"
TString prepcut;
TString generalcut;
//TString generalcutR="R.tr.n==1 && R.vdc.u1.nclust==1&& R.vdc.v1.nclust==1 && R.vdc.u2.nclust==1 && R.vdc.v2.nclust==1 && fEvtHdr.fEvtType==1 ";
//TString generalcutL="L.tr.n==1 && L.vdc.u1.nclust==1&& L.vdc.v1.nclust==1 && L.vdc.u2.nclust==1 && L.vdc.v2.nclust==1 && fEvtHdr.fEvtType==1";
TString generalcutR="fEvtHdr.fEvtType==1 ";
TString generalcutL="fEvtHdr.fEvtType==1";

std::string DataSavePath="carbonCheck/";


double MeshXMin=15300;
double MeshXMax=76000;
double MeshXNBin=10;

double MeshYMin=29000;
double MeshYMax=68000;
double MeshYNBin=10;


inline Bool_t IsFileExist (const std::string& name) {
//	  struct stat buffer;
//	  return (stat (name.c_str(), &buffer) == 0);
    if(gSystem->AccessPathName(name.c_str())){
        return false;
    }else{
        return true;
    }
}


TChain *LoadrootFile(UInt_t runID,TString folder="/home/newdriver/Storage/Research/CRex_Experiment/RasterReplay/Replay/Result", int logLevel=0){
	TChain *chain=new TChain("T");
	TString HRS="R";
		if(runID<20000){HRS="L";};

		if(folder.EndsWith(".root")){
			chain->Add(folder.Data());
		}else{
			TString rootDir(folder.Data());
			if(runID>20000){ //RHRS
				if(IsFileExist(Form("%s/prexRHRS_%d_-1.root",rootDir.Data(),runID))){
					if(logLevel > 2)
				    std::cout<<"Add File::"<<Form("%s/prexRHRS_%d_-1.root",rootDir.Data(),runID)<<std::endl;
					chain->Add(Form("%s/prexRHRS_%d_-1.root",rootDir.Data(),runID));

					TString filename;
					int16_t split=1;
					filename=Form("%s/prexRHRS_%d_-1_%d.root",rootDir.Data(),runID,split);
					while (IsFileExist(filename.Data())){
                        if(logLevel > 2)
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
                    if(logLevel > 2)std::cout<<"Add File::"<<Form("%s/prexLHRS_%d_-1.root",rootDir.Data(),runID)<<std::endl;
					chain->Add(Form("%s/prexLHRS_%d_-1.root",rootDir.Data(),runID));

					TString filename;
					int16_t split=1;
					filename=Form("%s/prexLHRS_%d_-1_%d.root",rootDir.Data(),runID,split);
					while (IsFileExist(filename.Data())){
                        if(logLevel > 2)std::cout<<"Add File::"<<filename.Data()<<std::endl;
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
/// \param DrawMesh
/// \return
std::map<int,std::map<int,int>> MeshRasterCurrent(TChain *chain,bool DrawMesh= true,TString resultFolder="./carbonCheck"){

    //Quartz cut
    const double loadc_cutL = 600;
    const double loadc_cutR = 565;

    const double upadc_cutL = 485;
    const double upadc_cutR = 505;


    std::map<int,std::map<int,int>> MeshCellEntryList;

    double loadc_cut=loadc_cutL;
    double upadc_cut=upadc_cutL;

    int runID = (int)chain->GetMaximum("fEvtHdr.fRun");
//    std::cout<<runID<<std::endl;

    int RangeYmin=(int)chain->GetMinimum("frawcurx2");
    int RangeYmax=(int)chain->GetMaximum("frawcurx2");
    int RangeXmin=(int)chain->GetMinimum("frawcury2");
    int RangeXmax=(int)chain->GetMaximum("frawcury2");

    // generate the range for the histogram
    int PlotRangeXmin=RangeXmin-(RangeXmax-RangeXmin)/10;
    int PlotRangeXmax=RangeXmax+(RangeXmax-RangeXmin)/10;
    int  PlotRangeYmin=RangeYmin-(RangeYmax-RangeYmin)/10;
    int  PlotRangeYmax=RangeYmax+(RangeYmax-RangeYmin)/10;


    TString HRS="L";
    if(runID>20000){ //RHRS
        HRS="R";
    }
    if(HRS=="R"){
        loadc_cut=loadc_cutR;
        upadc_cut=upadc_cutR;
    }

    TH2F *currenthh=(TH2F *)gROOT->FindObject("frawcurx2 vs. frawcury2");
    if(currenthh) currenthh->Delete();
    currenthh=new TH2F("frawcurx2 vs. frawcury2","frawcurx2 vs. frawcury2",100,PlotRangeXmin,PlotRangeXmax,100,PlotRangeYmin,PlotRangeYmax);

    TH2F *currentCuthh=(TH2F *)gROOT->FindObject("frawcurx2 vs. frawcury2 cut");
    if(currentCuthh)currentCuthh->Delete();
    currentCuthh=new TH2F("frawcurx2 vs. frawcury2 cut","frawcurx2 vs. frawcury2 cut",100,PlotRangeXmin,PlotRangeXmax,100,PlotRangeYmin,PlotRangeYmax);

    TH1F *upQadc=(TH1F *)gROOT->FindObject("upQadc");
    if(upQadc) upQadc->Delete();
    upQadc=new TH1F("upQadc","upQadc",1000,400,900);

    TH1F *loQadc=(TH1F *)gROOT->FindObject("loQadc");
    if (loQadc)loQadc->Delete();
    loQadc=new TH1F("loQadc","loQadc",1000,400,900);


    TString QadcCut(Form("P.upQadc%s> %f && P.loQadc%s> %f",HRS.Data(), upadc_cut,HRS.Data(),loadc_cut));
    chain->Project(currenthh->GetName(),"frawcurx2:frawcury2");
    chain->Project(currentCuthh->GetName(),"frawcurx2:frawcury2",QadcCut.Data());
    chain->Project(upQadc->GetName(),Form("P.upQadc%s",HRS.Data()));
    chain->Project(loQadc->GetName(),Form("P.loQadc%s",HRS.Data()));

    std::cout<<"histo  "<<currentCuthh->GetEntries()<<"   vs."<<(chain->GetEntries(QadcCut.Data()))<<std::endl;

    TCanvas *canv=(TCanvas *)gROOT->GetListOfCanvases()->FindObject(Form("Canv_runID%d",runID));
    if(!canv){
        canv=new TCanvas(Form("Canv_runID%d",runID),Form("Canv_runID%d",runID),1960,1080);
    }else{
        canv->Clear();
    }

    canv->Divide(2,2);
    canv->cd(1);
    currenthh->Draw("zcol");
//    currenthh->GetXaxis()->SetRangeUser(22000,70000);
//    currenthh->GetYaxis()->SetRangeUser(2000,95000);


    canv->cd(2);
    currentCuthh->Draw("zcol");

    //draw the messh and get the number of entries
    if(DrawMesh)
    {
        for(int xIter=0; xIter<=MeshXNBin; xIter++){
            double bin=(MeshXMax-MeshXMin)/MeshXNBin;
            TLine *line=new TLine(MeshXMin+xIter*bin,MeshYMin,MeshXMin+xIter*bin,MeshYMax);
            line->SetLineColor(kRed);
            line->Draw("same");
        }
        // draw the mess on the Y dimension
        for (int yIter=0; yIter<=MeshXNBin;yIter++){
            double bin=(MeshYMax-MeshYMin)/MeshYNBin;
            TLine *line=new TLine(MeshXMin,MeshYMin+yIter*bin,MeshXMax,MeshYMin+yIter*bin);
            line->SetLineColor(kRed);
            line->Draw("same");
        }
        // get the entries with cut
        for(int xIter=0; xIter<MeshXNBin; xIter++){
            for (int yIter=0; yIter<MeshXNBin;yIter++) {
                double xbin=(MeshXMax-MeshXMin)/MeshXNBin;
                double ybin=(MeshYMax-MeshYMin)/MeshYNBin;

                double boundaryXmin=MeshXMin+xIter*xbin;
                double boundaryXmax=MeshXMin+xIter*xbin+xbin;

                double boundaryYmin=MeshYMin+yIter*ybin;
                double boundaryYmax=MeshYMin+yIter*ybin+ybin;

                TString MeshCellCut(Form("frawcury2 > %f && frawcury2 <= %f && frawcurx2 > %f && frawcurx2 <= %f",boundaryXmin,boundaryXmax,boundaryYmin, boundaryYmax));
                MeshCellEntryList[xIter][yIter]=chain->GetEntries(Form("%s && %s",QadcCut.Data(),MeshCellCut.Data()));

                // draw the number of Entries on the Canvas
                if (true){
                    TLatex *text=new TLatex(boundaryXmin+0.2*xbin, boundaryYmin+0.2*ybin,Form("%d",MeshCellEntryList[xIter][yIter]));
                    text->SetTextSize(0.03);
                    text->SetTextColor(96);
                    text->Draw("same");
                }
                canv->Update();

            }
        }
    }

    canv->cd(3);
    upQadc->Draw();
    TLine *upQadcref=new TLine(upadc_cut,0,upadc_cut,upQadc->GetMaximum());
    upQadcref->SetLineWidth(2);
    upQadcref->SetLineColor(42);
    upQadcref->Draw("same");

    canv->cd(4);
    loQadc->Draw();
    TLine *loQadcref=new TLine(loadc_cut,0,loadc_cut,loQadc->GetMaximum());
    loQadcref->SetLineWidth(2);
    loQadcref->SetLineColor(42);
    loQadcref->Draw("same");
    canv->Update();
    canv->SaveAs(Form("%s/target_check_run%d.png",resultFolder.Data(),runID));

    return MeshCellEntryList;
}


///
/// \param chain
/// \param DrawMesh
/// \return
std::map<int,std::map<int,int>> RasterCurrent(TChain *chain,bool DrawMesh= true,TString resultFolder="./carbonCheck"){

    //Quartz cut
    const double loadc_cutL = 600;
    const double loadc_cutR = 565;

    const double upadc_cutL = 485;
    const double upadc_cutR = 505;


    std::map<int,std::map<int,int>> MeshCellEntryList;

    double loadc_cut=loadc_cutL;
    double upadc_cut=upadc_cutL;

    int runID = (int)chain->GetMaximum("fEvtHdr.fRun");
    TString HRS="L";
    generalcut=generalcutL;
    if(runID>20000){ //RHRS
        HRS="R";
        generalcut=generalcutR;
    }


    if(HRS=="R"){
        loadc_cut=loadc_cutR;
        upadc_cut=upadc_cutR;
    }

    int RangeYmin=(int)chain->GetMinimum(Form("%srb.Raster2.rawcur.y",HRS.Data()));
    int RangeYmax=(int)chain->GetMaximum(Form("%srb.Raster2.rawcur.y",HRS.Data()));
    int RangeXmin=(int)chain->GetMinimum(Form("%srb.Raster2.rawcur.x",HRS.Data()));
    int RangeXmax=(int)chain->GetMaximum(Form("%srb.Raster2.rawcur.x",HRS.Data()));

    std::cout<<"X dimension:: "<<RangeXmin<<"   "<<RangeXmax<<std::endl;
    std::cout<<"Y dimension:: "<<RangeYmin<<"   "<<RangeYmax<<std::endl;

    // generate the range for the histogram
    int PlotRangeXmin=RangeXmin-(RangeXmax-RangeXmin)/50;
    int PlotRangeXmax=RangeXmax+(RangeXmax-RangeXmin)/50;
    int  PlotRangeYmin=RangeYmin-(RangeYmax-RangeYmin)/50;
    int  PlotRangeYmax=RangeYmax+(RangeYmax-RangeYmin)/50;


    int EdgeCutXmin=RangeXmin+(RangeXmax-RangeXmin)/50;
    int EdgeCutXmax=RangeXmax-(RangeXmax-RangeXmin)/50;
    int EdgeCutYmin=RangeYmin+(RangeYmax-RangeYmin)/50;
    int EdgeCutYmax=RangeYmax-(RangeYmax-RangeYmin)/50;

    std::cout<<"Xmin:"<<EdgeCutXmin<<"  xMax:"<<EdgeCutXmax<<"   Ymin:"<<EdgeCutYmin<<"   Ymax:"<<EdgeCutYmax<<std::endl;


    // Initial plot without cut
    TH2F *currenthh=(TH2F *)gROOT->FindObject("frawcurx2 vs. frawcury2");
    if(currenthh) currenthh->Delete();
    currenthh=new TH2F("frawcurx2 vs. frawcury2","frawcurx2 vs. frawcury2",100,PlotRangeXmin,PlotRangeXmax,100,PlotRangeYmin,PlotRangeYmax);
    currenthh->GetXaxis()->SetTitle(Form("%srb.Raster2.rawcur.x",HRS.Data()));
    currenthh->GetYaxis()->SetTitle(Form("%srb.Raster2.rawcur.y",HRS.Data()));

    //Cut with only the adge cut
    TH2F *currentEdgeCuthh=(TH2F *)gROOT->FindObject(Form("rawCurX2 vs. rawCury2 edge Cut %d",runID));
    if(currentEdgeCuthh) currentEdgeCuthh->Delete();
    currentEdgeCuthh=new TH2F(Form("rawCurX2 vs. rawCury2 edge Cut %d",runID),Form("rawCurX2 vs. rawCury2 edge Cut %d",runID),100,PlotRangeXmin,PlotRangeXmax,100,PlotRangeYmin,PlotRangeYmax);
    currentEdgeCuthh->GetXaxis()->SetTitle(Form("%srb.Raster2.rawcur.x",HRS.Data()));
    currentEdgeCuthh->GetYaxis()->SetTitle(Form("%srb.Raster2.rawcur.y",HRS.Data()));
    currentEdgeCuthh->GetXaxis()->SetRangeUser(EdgeCutXmin,EdgeCutXmax);
    currentEdgeCuthh->GetYaxis()->SetRangeUser(EdgeCutYmin,EdgeCutYmax);

    TH2F *currentEdgePedhh=(TH2F *)gROOT->FindObject(Form("rawCurX2 vs. rawCury2 edge Ped %d",runID));
    if(currentEdgePedhh) currentEdgePedhh->Delete();
    currentEdgePedhh=new TH2F(Form("rawCurX2 vs. rawCury2 edge Ped %d",runID),Form("rawCurX2 vs. rawCury2 edge Ped %d",runID),100,PlotRangeXmin,PlotRangeXmax,100,PlotRangeYmin,PlotRangeYmax);
    currentEdgePedhh->GetXaxis()->SetTitle(Form("%srb.Raster2.rawcur.x",HRS.Data()));
    currentEdgePedhh->GetYaxis()->SetTitle(Form("%srb.Raster2.rawcur.y",HRS.Data()));
    currentEdgePedhh->GetXaxis()->SetRangeUser(EdgeCutXmin,EdgeCutXmax);
    currentEdgePedhh->GetYaxis()->SetRangeUser(EdgeCutYmin,EdgeCutYmax);

    // Start the Edge cut and the lower the Qadc cut
    TH2F *currentEdgePedCuthh=(TH2F *)gROOT->FindObject(Form("rawCurX2 vs. rawCury2 edge Ped Cut %d",runID));
    if(currentEdgePedCuthh) currentEdgePedCuthh->Delete();
    currentEdgePedCuthh=new TH2F(Form("rawCurX2 vs. rawCury2 edge Ped Cut %d",runID),Form("rawCurX2 vs. rawCury2 edge Ped Cut %d",runID),100,PlotRangeXmin,PlotRangeXmax,100,PlotRangeYmin,PlotRangeYmax);
    currentEdgePedCuthh->GetXaxis()->SetTitle(Form("%srb.Raster2.rawcur.x",HRS.Data()));
    currentEdgePedCuthh->GetYaxis()->SetTitle(Form("%srb.Raster2.rawcur.y",HRS.Data()));
    currentEdgePedCuthh->GetXaxis()->SetRangeUser(EdgeCutXmin,EdgeCutXmax);
    currentEdgePedCuthh->GetYaxis()->SetRangeUser(EdgeCutYmin,EdgeCutYmax);


    TH2F *meshcellhh=(TH2F *) gROOT->FindObject("MeshCell");
    if(meshcellhh) meshcellhh->Delete();
    meshcellhh=new TH2F("MeshCell","MeshCell",100,PlotRangeXmin,PlotRangeXmax,100,PlotRangeYmin,PlotRangeYmax);
    meshcellhh->GetXaxis()->SetTitle(Form("%srb.Raster2.rawcur.x",HRS.Data()));
    meshcellhh->GetYaxis()->SetTitle(Form("%srb.Raster2.rawcur.y",HRS.Data()));
    meshcellhh->GetXaxis()->SetRangeUser(EdgeCutXmin,EdgeCutXmax);
    meshcellhh->GetYaxis()->SetRangeUser(EdgeCutYmin,EdgeCutYmax);

    TH2F *targetPoshh=(TH2F *)gROOT->FindObject("targX_vs_targY");
    if (targetPoshh) targetPoshh->Delete();
    targetPoshh=new TH2F("targX_vs_targY","targX_vs_targY",500,-8,8,500,-8,8);
    targetPoshh->GetYaxis()->SetTitle("TargetY");
    targetPoshh->GetXaxis()->SetTitle("TargetX");

    TH1F *targetXPoshh=(TH1F *)gROOT->FindObject("targX");
    if(targetXPoshh) targetXPoshh->Delete();
    targetXPoshh=new TH1F("targX","targX",500,-8,8);


    TH1F *targetYPoshh=(TH1F *)gROOT->FindObject("targY");
    if(targetYPoshh) targetYPoshh->Delete();
    targetYPoshh=new TH1F("targY","targY",500,-8,8);



    TH1F *upQadc=(TH1F *)gROOT->FindObject("upQadc");
    if(upQadc) upQadc->Delete();
    upQadc=new TH1F("upQadc","upQadc",1000,400,900);

    TH1F *loQadc=(TH1F *)gROOT->FindObject("loQadc");
    if (loQadc)loQadc->Delete();
    loQadc=new TH1F("loQadc","loQadc",1000,400,900);

    TString edgeCut(Form("%srb.Raster2.rawcur.x > %d && %srb.Raster2.rawcur.x < %d && %srb.Raster2.rawcur.y > %d && %srb.Raster2.rawcur.y < %d",HRS.Data(),EdgeCutXmin,HRS.Data(),EdgeCutXmax,HRS.Data(),EdgeCutYmin,HRS.Data(),EdgeCutYmax));
    TString QadcPedCut(Form("P.upQadc%s < %f && P.loQadc%s < %f",HRS.Data(), upadc_cut,HRS.Data(),loadc_cut));
    TString QadcCut(Form("P.upQadc%s> %f && P.loQadc%s> %f",HRS.Data(), upadc_cut,HRS.Data(),loadc_cut));


    chain->Project(currenthh->GetName(),Form("%srb.Raster2.rawcur.y:%srb.Raster2.rawcur.x",HRS.Data(),HRS.Data()),Form("%s",generalcut.Data()));
    chain->Project(currentEdgeCuthh->GetName(),Form("%srb.Raster2.rawcur.y:%srb.Raster2.rawcur.x",HRS.Data(),HRS.Data()),Form("%s",generalcut.Data()));
    chain->Project(currentEdgePedCuthh->GetName(),Form("%srb.Raster2.rawcur.y:%srb.Raster2.rawcur.x",HRS.Data(),HRS.Data()),Form("%s && %s",QadcCut.Data(),generalcut.Data()));
    chain->Project(currentEdgePedhh->GetName(),Form("%srb.Raster2.rawcur.y:%srb.Raster2.rawcur.x",HRS.Data(),HRS.Data()),Form("%s && %s",QadcPedCut.Data(),generalcut.Data()));
    chain->Project(upQadc->GetName(),Form("P.upQadc%s",HRS.Data()),Form("%s",generalcut.Data()));
    chain->Project(loQadc->GetName(),Form("P.loQadc%s",HRS.Data()),Form("%s",generalcut.Data()));
    chain->Project(targetPoshh->GetName(),"targy:targx",Form("%s",generalcut.Data()));
    chain->Project(targetXPoshh->GetName(),"targx",Form("%s",generalcut.Data()));
    chain->Project(targetYPoshh->GetName(),"targy",Form("%s",generalcut.Data()));

    TCanvas *canv=(TCanvas *)gROOT->GetListOfCanvases()->FindObject(Form("Canv_runID%d",runID));
    if(!canv){
        canv=new TCanvas(Form("Canv_runID%d",runID),Form("Canv_runID%d",runID),1960,1080);
    }else{
        canv->Clear();
    }
    canv->Divide(4,2);
    canv->cd(1);
    currenthh->Draw("zcol");
    {
        TLine *line0=new TLine(EdgeCutXmin,EdgeCutYmin,EdgeCutXmax,EdgeCutYmin);
        line0->SetLineColor(2);
        line0->Draw("same");

        TLine *line1=new TLine(EdgeCutXmin,EdgeCutYmax,EdgeCutXmax,EdgeCutYmax);
        line1->SetLineColor(2);
        line1->Draw("same");

        TLine *line2=new TLine(EdgeCutXmin,EdgeCutYmin,EdgeCutXmin,EdgeCutYmax);
        line2->SetLineColor(2);
        line2->Draw("same");

        TLine *line3=new TLine(EdgeCutXmax,EdgeCutYmin,EdgeCutXmax,EdgeCutYmax);
        line3->SetLineColor(2);
        line3->Draw("same");
    }

    canv->cd(2);
    currentEdgeCuthh->Draw("zcol");

    canv->cd(3);
    currentEdgePedCuthh->Draw("zcol");

    if(DrawMesh)
    {

        std::cout<<"Using Boundary: X"<<MeshXMin << " -> "<<MeshXMax<<"  Y :"<<MeshYMin<<" --> "<< MeshYMax<<std::endl;

        for(int xIter=0; xIter<=MeshXNBin; xIter++){
            double bin=(MeshXMax-MeshXMin)/MeshXNBin;
            TLine *line=new TLine(MeshXMin+xIter*bin,MeshYMin,MeshXMin+xIter*bin,MeshYMax);
            line->SetLineColor(kRed);
            line->Draw("same");
        }
        // draw the mess on the Y dimension
        for (int yIter=0; yIter<=MeshXNBin;yIter++){
            double bin=(MeshYMax-MeshYMin)/MeshYNBin;
            TLine *line=new TLine(MeshXMin,MeshYMin+yIter*bin,MeshXMax,MeshYMin+yIter*bin);
            line->SetLineColor(kRed);
            line->Draw("same");
        }
        // get the entries with cut
        for(int xIter=0; xIter<MeshXNBin; xIter++){
            for (int yIter=0; yIter<MeshXNBin;yIter++) {
                double xbin=(MeshXMax-MeshXMin)/MeshXNBin;
                double ybin=(MeshYMax-MeshYMin)/MeshYNBin;

                double boundaryXmin=MeshXMin+xIter*xbin;
                double boundaryXmax=MeshXMin+xIter*xbin+xbin;

                double boundaryYmin=MeshYMin+yIter*ybin;
                double boundaryYmax=MeshYMin+yIter*ybin+ybin;

                TString MeshCellCut(Form("%srb.Raster2.rawcur.y > %f && %srb.Raster2.rawcur.y < %f && %srb.Raster2.rawcur.x >%f && %srb.Raster2.rawcur.x < %f",HRS.Data(),boundaryYmin,HRS.Data(),boundaryYmax,HRS.Data(),boundaryXmin,HRS.Data(),boundaryXmax));
                MeshCellEntryList[xIter][yIter]=chain->GetEntries(Form("%s && %s",QadcCut.Data(),MeshCellCut.Data()));

                chain->Project(meshcellhh->GetName(),Form("%srb.Raster2.rawcur.y:%srb.Raster2.rawcur.x",HRS.Data(),HRS.Data()),Form("%s && %s",QadcCut.Data(),MeshCellCut.Data()));
                meshcellhh->Draw("same");

                // draw the number of Entries on the Canvas
                if (true){
                    TLatex *text=new TLatex(boundaryXmin+0.2*xbin, boundaryYmin+0.2*ybin,Form("%d",MeshCellEntryList[xIter][yIter]));
                    text->SetTextSize(0.03);
                    text->SetTextColor(96);
                    text->Draw("same");
                }
                canv->Update();

            }
        }
    }

    canv->cd(4);
    currentEdgePedhh->Draw("zcol");


    canv->cd(5);
    upQadc->Draw();
    TLine *upQadcref=new TLine(upadc_cut,0,upadc_cut,upQadc->GetMaximum());
    upQadcref->SetLineWidth(2);
    upQadcref->SetLineColor(42);
    upQadcref->Draw("same");

    canv->cd(6);
    loQadc->Draw();
    TLine *loQadcref=new TLine(loadc_cut,0,loadc_cut,loQadc->GetMaximum());
    loQadcref->SetLineWidth(2);
    loQadcref->SetLineColor(42);
    loQadcref->Draw("same");


    canv->cd(7);
    canv->cd(7)->SetGridx();
    canv->cd(7)->SetGridy();
    targetPoshh->Draw("zcol");

    canv->cd(8)->Divide(1,2);
    canv->cd(8)->cd(1);
    targetXPoshh->Draw();
    targetXPoshh->Fit("gaus");
    TLatex *text1=new TLatex(targetXPoshh->GetFunction("gaus")->GetParameter(1),targetXPoshh->GetFunction("gaus")->GetParameter(0),Form("%1.3f",targetXPoshh->GetFunction("gaus")->GetParameter(1)));
    text1->Draw("same");

    canv->cd(8)->cd(2);
    targetYPoshh->Draw();
    targetYPoshh->Fit("gaus");
    TLatex *text2=new TLatex(targetYPoshh->GetFunction("gaus")->GetParameter(1),targetYPoshh->GetFunction("gaus")->GetParameter(0),Form("%1.3f",targetYPoshh->GetFunction("gaus")->GetParameter(1)));
    text2->Draw("same");

    canv->SaveAs(Form("%s/target_%d.jpg",resultFolder.Data(),runID));
    return MeshCellEntryList;
}

double GetCharge(int runID){
    double TotalCharge=0.0;

    return TotalCharge;
}

/// \param referenceMesh, the Initial Mesh(The target is not damanged)
/// \param targetMesh ,
/// \return    thickness ratio
double NormalizeRatio(std::map<int,std::map<int,int>> referenceMesh,std::map<int,std::map<int,int>> targetMesh,std::map<int,std::map<int,double>> & meshedRatio){

    int ref_EdgeSum=0;
    int tag_EdgeSum=0;

    int ref_CenterSum=0;
    int tag_CenterSum=0;

    for(auto xIter=referenceMesh.begin();xIter!=referenceMesh.end();xIter++){
        for(auto yItter=(xIter->second).begin(); yItter!=(xIter->second).end();yItter++){
            int xMeshID=xIter->first;
            int yMeshID=yItter->first;

            //used for calculate the ratio
//            if((xIter->first==0)||(xIter->first==9)||(yItter->first==0)||(yItter->first==9))
//            if((xIter->first==0)&&(yItter->first==9))
            {
                ref_EdgeSum+=yItter->second;
                tag_EdgeSum+=targetMesh[xIter->first][yItter->first];
            }

            if((xMeshID>2 && xMeshID <6)&&(yMeshID>2 && yMeshID <6)){
                //calculate the ratio
                ref_CenterSum+=yItter->second;
                tag_CenterSum+=targetMesh[xIter->first][yItter->first];
            }
        }
    }

    double referenceFactor=(double ) tag_EdgeSum/ref_EdgeSum;
    double ratio=(double ) tag_CenterSum/((double )ref_CenterSum*referenceFactor);
    // used for calculate the ratio
    for(auto xIter=referenceMesh.begin();xIter!=referenceMesh.end();xIter++) {
        for (auto yItter = xIter->second.begin(); yItter != xIter->second.end(); yItter++) {
            //calculate the ratio
            int refEntries= yItter->second;
            int targEntries=targetMesh[xIter->first][yItter->first];
            // fill the value ratio
            meshedRatio[xIter->first][yItter->first]=(double ) targEntries/(refEntries*referenceFactor);
        }
    }
    return  ratio;
}


void  GetBoundary(std::map<int,int> runList, TString folder="/home/newdriver/pyQuant/prex_replayed/rootfile",TString resultFolder="./carbonCheck"){

    std::vector<int> runIDs;
    for (auto iter= runList.begin();iter!=runList.end();iter++){
        runIDs.push_back(iter->second);
    }
//    UInt_t runIDs[]={2140,2141,2182,2181,2204,2202,2199,2291,2300};

    int boundaryXmin=0;
    int boundaryXmax=1000000;

    int boundaryYmin=0;
    int boundaryYmax=1000000;

//    for (int i =0; i < sizeof(runIDs)/sizeof(UInt_t); i++)
    for(auto item : runIDs)
    {
//        std::cout<<runIDs[i]<<std::endl;
//        int runID=runIDs[i];
        std::cout<< item<<std::endl;
        int runID=item;


        TString HRS="L";
        if(runID>20000){ //RHRS
            HRS="R";
        }

        auto chain=LoadrootFile(runID,folder);
        int RangeYmin=(int)chain->GetMinimum(Form("%srb.Raster2.rawcur.y",HRS.Data()));
        int RangeYmax=(int)chain->GetMaximum(Form("%srb.Raster2.rawcur.y",HRS.Data()));
        int RangeXmin=(int)chain->GetMinimum(Form("%srb.Raster2.rawcur.x",HRS.Data()));
        int RangeXmax=(int)chain->GetMaximum(Form("%srb.Raster2.rawcur.x",HRS.Data()));

        std::cout<<"X dimension:: "<<RangeXmin<<"   "<<RangeXmax<<std::endl;
        std::cout<<"Y dimension:: "<<RangeYmin<<"   "<<RangeYmax<<std::endl;

        // generate the range for the histogram
        int  PlotRangeXmin=RangeXmin-(RangeXmax-RangeXmin)/10;
        int  PlotRangeXmax=RangeXmax+(RangeXmax-RangeXmin)/10;
        int  PlotRangeYmin=RangeYmin-(RangeYmax-RangeYmin)/10;
        int  PlotRangeYmax=RangeYmax+(RangeYmax-RangeYmin)/10;

        int EdgeCutXmin=RangeXmin+(RangeXmax-RangeXmin)/20;
        int EdgeCutXmax=RangeXmax-(RangeXmax-RangeXmin)/20;
        int EdgeCutYmin=RangeYmin+(RangeYmax-RangeYmin)/20;
        int EdgeCutYmax=RangeYmax-(RangeYmax-RangeYmin)/20;

        std::cout<<"Xmin:"<<EdgeCutXmin<<"  xMax:"<<EdgeCutXmax<<"   Ymin:"<<EdgeCutYmin<<"   Ymax:"<<EdgeCutYmax<<std::endl;

        if (EdgeCutXmin > boundaryXmin) boundaryXmin = EdgeCutXmin;
        if (EdgeCutXmax < boundaryXmax) boundaryXmax = EdgeCutXmax;
        if (EdgeCutYmin > boundaryYmin) boundaryYmin = EdgeCutYmin;
        if (EdgeCutYmax < boundaryYmax) boundaryYmax = EdgeCutYmax;

    }

    MeshXMin=boundaryXmin;
    MeshXMax=boundaryXmax;
    MeshYMin=boundaryYmin;
    MeshYMax=boundaryYmax;


    std::cout<<"Boundary: X"<<MeshXMin << " -> "<<MeshXMax<<"  Y :"<<MeshYMin<<" --> "<< MeshYMax<<std::endl;
}

///
/// \param runFile
/// \param folder
void TargetThicknessCal(TString runFile="",TString folder="/home/newdriver/pyQuant/prex_replayed/rootfile"){


    std::map<int,int> runList;
    if(!runFile.IsNull()){
        // if the input the run list file
        std::string  line;
        std::ifstream fin;
        fin.open(runFile.Data());
        if (fin.is_open()){
            while (getline(fin,line)){
                line=line.substr(0,line.find("#",0));// remove the commend line
                if(!line.empty()){   // decode the line
                    std::stringstream  ss(line);
                    std::vector<int> value;
                    while(getline(ss,line,',')){
                        value.push_back(std::stoi(line));
                    }
                    if (value.size()==2){
                        runList[value[0]]=value[1];
                    }
                    value.clear();
                }
            }
        }else{
            std::cout << "\033[1;31m [ERROR]\033[0m"<<"  Can NOT load file:: "<<runFile.Data();
        }
        fin.close();
        if(runList.find(0)==runList.end()){
            std::cout<<"\033[1;31m [ERROR]\033[0m"<< "Input file format error\n"<<
            "\t ===> Accept Format: \n"<<
            " \t\t\tID, runID \n \t\t\tID=0 will be used as reference run!!"<<std::endl;
        }
    } else{

//        runList[0]=2140;
//        runList[1]=2182;
//        runList[2]=2204;
//        runList[3]=2291;
//        runList[4]=2300;

        runList[0]=21179;
        runList[1]=21181;
        runList[2]=21290;
        runList[3]=21415;



//        runList[0]=21014;
//        runList[1]=21180;
//        runList[2]=21195;
//        runList[3]=21290;
//        runList[0]=2140;
//        runList[1]=2140;
//        runList[2]=2141;
//        runList[3]=2291;
//        runList[4]=2299;
//        runList[5]=2300;
    }

    GetBoundary(runList);

    std::map<int, std::map<int, std::map<int, double>>> TargThicknessRationList;

    std::cout<<"Run Lis"<<std::endl;
    for (auto i = runList.begin(); i!=runList.end();i++){
        std::cout<<'\t'<<i->first<<"   "<<i->second<<std::endl;
    }

    std::map<int,double> reletiveThickness;
    reletiveThickness[0]=1.000;

    auto referenceMesh=RasterCurrent(LoadrootFile(runList[0],folder));
    auto runIter=(runList.begin());
    for (runIter++;runIter!=runList.end();runIter++){
        std::cout<<runIter->second<<std::endl;
        int runID=runIter->second;
        auto meshInfor=RasterCurrent(LoadrootFile(runIter->second,folder));

        std::map<int,std::map<int,double>> meshedCellRatio;
        double thickness=NormalizeRatio(referenceMesh,meshInfor, meshedCellRatio);
        std::cout<<"runID"<<runIter->second<<"   "<<thickness<<std::endl;
        reletiveThickness[runIter->first]=thickness;
        TargThicknessRationList[runIter->first]=meshedCellRatio;

        //TODO, need to get the plot and write to PDF

        //plot the Data row by row
        {
            int NcellX=meshedCellRatio.size();
            int NcellY=0;
            for (auto xIter=meshedCellRatio.begin();xIter!=meshedCellRatio.end();xIter++) {
                if (NcellY < xIter->second.size()) NcellY=xIter->second.size();
            }
            // write the data into the
            // Y iter    Histo on X dimension
            std::map<int,TH1F *> thicknesshh;
            //std::map<int,TH2F *> thicknessMaphh;

            //initialize the plot
            int colorTemp=17;
            for (auto xIter=meshedCellRatio.begin();xIter!=meshedCellRatio.end();xIter++) {
                for(auto yItter=xIter->second.begin();yItter!=xIter->second.end();yItter++){
                    //
                    if(thicknesshh.find(yItter->first)==thicknesshh.end()){

                        TString plotName(Form("run%d_Y_%c_hh", runID, char(yItter->first+65)));
                        std::cout<<yItter->first<<"    "<<plotName.Data()<<std::endl;

                        thicknesshh[yItter->first]=(TH1F *) gROOT->FindObject(plotName.Data());
                        if(thicknesshh[yItter->first]) thicknesshh[yItter->first]->Delete();
                        thicknesshh[yItter->first]=new TH1F(plotName.Data(),plotName.Data(),15,-3,12);
                        thicknesshh[yItter->first]->SetLineColor(colorTemp+10);
                        thicknesshh[yItter->first]->SetLineWidth(2);
                        thicknesshh[yItter->first]->SetMarkerStyle(20);
                        thicknesshh[yItter->first]->SetMarkerColor(colorTemp+10);
                        thicknesshh[yItter->first]->GetYaxis()->SetRangeUser(0.7,1.2);
                        colorTemp=colorTemp+10;
                    }// if not initilized, initialLize the plot
                    thicknesshh[yItter->first]->Fill(xIter->first,meshedCellRatio[xIter->first][yItter->first]);
                    std::cout<<"y"<<yItter->first<<"  x "<<xIter->first<<"  value"<<meshedCellRatio[xIter->first][yItter->first]<<std::endl;
                    thicknesshh[yItter->first]->SetBinError(xIter->first,meshedCellRatio[xIter->first][yItter->first]*0.015);
                }
            }
            // draw the plot in canvas
            {
                TCanvas *meshedCanv = (TCanvas *) gROOT->GetListOfCanvases()->FindObject(
                        Form("MeshedCellCanv_runID%d", runID));
                if (!meshedCanv) {
                    meshedCanv = new TCanvas(Form("MeshedCellCanv_runID%d", runID),
                                             Form("MeshedCellCanv_runID%d", runID), 1960, 1080);
                } else {
                    meshedCanv->Clear();
                }
                meshedCanv->Draw();
                meshedCanv->cd();
                TLegend *lgend = new TLegend(0.1, 0.7, 0.48, 0.9);
                for (auto iter = thicknesshh.begin(); iter != thicknesshh.end(); iter++) {
                    lgend->AddEntry(iter->second, iter->second->GetName());
                    if (iter == thicknesshh.begin()) {
                        iter->second->Draw("HIST P");
                    } else {
                        iter->second->Draw("same HIST P");
                    }
                }
                lgend->Draw("same");
                meshedCanv->Update();
                meshedCanv->SaveAs(Form("carbonCheck/%s.jpg", meshedCanv->GetName()));
            }

            // get the relative thickness colored map
            {
                TCanvas *thicknessRatiohhCanv=(TCanvas *) gROOT->GetListOfCanvases()->FindObject(Form("Thickness_Ratio_Map_run_%d",runID));
                if(thicknessRatiohhCanv)thicknessRatiohhCanv->Delete();
                thicknessRatiohhCanv=new TCanvas(Form("Thickness_Ratio_Map_run_%d",runID),Form("Thickness_Ratio_Map_run_%d",runID),1960,1080);
                thicknessRatiohhCanv->Draw();

                TH2F *thicknessMaphh=(TH2F *) gROOT->FindObject(Form("Relative Thickness Map %d",runID));
                if(thicknessMaphh)thicknessMaphh->Delete();
                thicknessMaphh=new TH2F(Form("Relative Thickness Map %d",runID),Form("Relative Thickness Map %d",runID),MeshXNBin,MeshXMin,MeshXMax,MeshYNBin,MeshYMin,MeshYMax);

                // get the mesh
                for(int xIter=0; xIter<MeshXNBin; xIter++) {
                    for (int yIter = 0; yIter < MeshXNBin; yIter++) {
                        double xbin = (MeshXMax - MeshXMin) / MeshXNBin;
                        double ybin = (MeshYMax - MeshYMin) / MeshYNBin;

                        double boundaryXmin = MeshXMin + xIter * xbin;
                        double boundaryXmax = MeshXMin + xIter * xbin + xbin;

                        double boundaryYmin = MeshYMin + yIter * ybin;
                        double boundaryYmax = MeshYMin + yIter * ybin + ybin;

                        double  binCenterX=(boundaryXmax+boundaryXmin)/2;
                        double  binCenterY=(boundaryYmax+boundaryYmin)/2;

                        // get the cell center, and fill it with the color
                        if((meshedCellRatio.find(xIter)!=meshedCellRatio.end())&&(meshedCellRatio[xIter].find(yIter)!=meshedCellRatio[xIter].end())){
                            thicknessMaphh->Fill(binCenterX,binCenterY,meshedCellRatio[xIter][yIter]);
                        }
                    }
                }

                TH2F *currentEdgePedhh=(TH2F *)gROOT->FindObject(Form("rawCurX2 vs. rawCury2 edge Ped %d",runID));
                if(currentEdgePedhh){
                    double  xMax = currentEdgePedhh->GetXaxis()->GetXmax();
                    double  xMin = currentEdgePedhh->GetXaxis()->GetXmin();

                    double  yMax = currentEdgePedhh->GetYaxis()->GetXmax();
                    double  yMin = currentEdgePedhh->GetYaxis()->GetXmin();

                //    thicknessMaphh->GetXaxis()->SetRange(xMin,xMax);
                //    thicknessMaphh->GetYaxis()->SetRange(yMin,yMax);
                }

                thicknessMaphh->Draw("zcol");
                thicknessRatiohhCanv->Update();
                thicknessRatiohhCanv->SaveAs(Form("carbonCheck/%s.jpg", thicknessRatiohhCanv->GetName()));
            }


        }

    }

    // plot the Error canvas
    {
        TCanvas *targetThicknessCanv = new TCanvas(Form("Run%d", runList[0]), Form("Run%d", runList[0]), 2080, 1960);


        double thicknessX[runList.size()];
        double thicknessY[runList.size()];
        double thicknessXErr[runList.size()];
        double thicknessYErr[runList.size()];

        for (auto item = reletiveThickness.begin(); item != reletiveThickness.end(); item++) {
            thicknessX[item->first] = item->first;
            thicknessY[item->first] = item->second;
            thicknessXErr[item->first] = 0;
            thicknessYErr[item->first] = item->second * 0.015;
        }

        auto geprex = new TGraphErrors(runList.size() + 1, thicknessX, thicknessY, thicknessXErr, thicknessYErr);
        geprex->GetYaxis()->SetRangeUser(0.8, 1.2);
        geprex->GetXaxis()->SetRangeUser(-2, 7);
        geprex->GetXaxis()->SetLimits(-2, 7);
        geprex->SetTitle("Target D9-208Pb7-D10");
        geprex->GetXaxis()->SetTitle("RunID");
        geprex->GetYaxis()->SetTitle("Thickness");
        geprex->SetLineWidth(2);
        geprex->SetLineColor(6);
        geprex->SetMarkerStyle(20);
        geprex->SetMarkerColor(6);
        geprex->Draw("ap");
        {
            if (true) {

                for (auto item = reletiveThickness.begin(); item != reletiveThickness.end(); item++) {
                    TLatex *text = new TLatex(item->first, item->second + 0.05, Form("%1.4f", item->second));
                    text->SetTextColor(94);
                    text->Draw("same");
                }
            }

            if (true) {
                TLine *line1 = new TLine(geprex->GetXaxis()->GetXmin(), 1.0, geprex->GetXaxis()->GetXmax(), 1.0);
                line1->SetLineColor(88);

                TLine *line2 = new TLine(geprex->GetXaxis()->GetXmin(), 1.05, geprex->GetXaxis()->GetXmax(), 1.05);
                line2->SetLineColor(92);

                TLine *line3 = new TLine(geprex->GetXaxis()->GetXmin(), 0.95, geprex->GetXaxis()->GetXmax(), 0.95);
                line3->SetLineColor(92);

                line1->Draw("same");
                line2->Draw("same");
                line3->Draw("same");
            }
        }
    }

    // plot all calv and save it into pdf files
/*    {

       std::map<int, double> reletiveCharge;
        reletiveCharge[0] =0.0;
        reletiveCharge[1] =2.97;
        reletiveCharge[2] =5.09;
        reletiveCharge[3] =15.54;


        // plot the charge thickness canvas
        TCanvas *targetThicknessCanv=new TCanvas(Form("Charge Run%d",runList[0]),Form("CharefRun%d",runList[0]),2080,1960);
        double thicknessX[runList.size()];
        double thicknessY[runList.size()];
        double thicknessXErr[runList.size()];
        double thicknessYErr[runList.size()];

        for (auto item = reletiveThickness.begin();item!=reletiveThickness.end();item++){
            if (reletiveCharge.find(item->first)==reletiveCharge.end()) continue;
            thicknessX[item->first]=reletiveCharge[item->first];
            thicknessY[item->first]=item->second;
            thicknessXErr[item->first]=0;
            thicknessYErr[item->first]=item->second*0.03;
            std::cout<<"Write "<< item->first<<" Charge"<<reletiveCharge[item->first]<<std::endl;
        }
        auto geprex=new TGraphErrors(21,thicknessX,thicknessY,thicknessXErr,thicknessYErr);
        geprex->GetYaxis()->SetRangeUser(0.8,1.2);
        geprex->GetXaxis()->SetRangeUser(-1,20);
        geprex->GetXaxis()->SetLimits(-1,20);
        geprex->SetTitle("Target D9-208Pb7-D10");
        geprex->GetXaxis()->SetTitle("RunID");
        geprex->GetYaxis()->SetTitle("Thickness");
        geprex->SetLineWidth(2);
        geprex->SetLineColor(6);
        geprex->SetMarkerStyle(20);
        geprex->SetMarkerColor(6);
        geprex->Draw("ap");
        {
            if (true){

                for (auto item = reletiveThickness.begin();item!=reletiveThickness.end();item++){
                    if(reletiveCharge.find(item->first)!=reletiveCharge.end()){
                    TLatex *text=new TLatex(reletiveCharge[item->first],item->second+0.05,Form("%1.4f",item->second));
                    text->SetTextColor(94);
                    text->Draw("same");
                    }
                }
            }

            if(true){
                TLine *line1=new TLine(geprex->GetXaxis()->GetXmin(),1.0,geprex->GetXaxis()->GetXmax(),1.0);
                line1->SetLineColor(88);

                TLine *line2=new TLine(geprex->GetXaxis()->GetXmin(),1.05,geprex->GetXaxis()->GetXmax(),1.05);
                line2->SetLineColor(92);

                TLine *line3=new TLine(geprex->GetXaxis()->GetXmin(),0.95,geprex->GetXaxis()->GetXmax(),0.95);
                line3->SetLineColor(92);

                line1->Draw("same");
                line2->Draw("same");
                line3->Draw("same");
            }
        }

    }*/

}



///
/// \param runID, the run number that want to check file. It will look for the root file that matches the runID
/// \param folde, root file location
/// \param resultFolder: folders that used for save the generated plot
void TargetCheck(UInt_t runID, TString folder="/home/newdriver/pyQuant/prex_replayed/rootfile",TString resultFolder="./carbonCheck"){
    auto chain=LoadrootFile(runID,folder);
    RasterCurrent(chain, false,resultFolder.Data());
}

