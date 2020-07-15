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
#include "iostream"
#include "TSystem.h"
#include "TGraphErrors.h"
#include "iostream"
#include "fstream"
#include "sstream"


TString prepcut;
TString generalcut;
TString generalcutR="R.tr.n==1 && R.vdc.u1.nclust==1&& R.vdc.v1.nclust==1 && R.vdc.u2.nclust==1 && R.vdc.v2.nclust==1 && R.gold.p > 2.14 && R.gold.p < 2.2 && fEvtHdr.fEvtType==1 ";
TString generalcutL="L.tr.n==1 && L.vdc.u1.nclust==1&& L.vdc.v1.nclust==1 && L.vdc.u2.nclust==1 && L.vdc.v2.nclust==1  && L.gold.p > 2.14 && L.gold.p < 2.2 && fEvtHdr.fEvtType==1";

std::string DataSavePath="carbonCheck/";


inline Bool_t IsFileExist (const std::string& name) {
//	  struct stat buffer;
//	  return (stat (name.c_str(), &buffer) == 0);
    if(gSystem->AccessPathName(name.c_str())){
        return false;
    }else{
        return true;
    }
}


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



struct MeshStruc{
    int Xid=0;
    int Yid=0;

    double xmin=0;
    double xmax=0;
    double ymin=0;
    double ymax=0;

    int Entries=0;
};

///
/// \param chain
/// \param DrawMesh
/// \return
std::map<int,std::map<int,int>> MeshRasterCurrent(TChain *chain,bool DrawMesh= true,TString resultFolder="./carbonCheck"){

    const double MeshXMin=25000;
    const double MeshXMax=70000;
    const double MeshXNBin=10;

    const double MeshYMin=13000;
    const double MeshYMax=80000;
    const double MeshYNBin=10;

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


/// \param referenceMesh
/// \param targetMesh
/// \return
std::map<int,std::map<int,int>> NormalizeCell(std::map<int,std::map<int,int>> referenceMesh,std::map<int,std::map<int,int>> targetMesh){

    std::map<int,std::map<int,int>> a ;
    return  a;

}


/// \param referenceMesh, the Initial Mesh(The target is not damanged)
/// \param targetMesh ,
/// \return    thickness ratio
double NormalizeRatio(std::map<int,std::map<int,int>> referenceMesh,std::map<int,std::map<int,int>> targetMesh){

    int ref_EdgeSum=0;
    int tag_EdgeSum=0;

    int ref_CenterSum=0;
    int tag_CenterSum=0;

    for(auto xIter=referenceMesh.begin();xIter!=referenceMesh.end();xIter++){
        for(auto yItter=xIter->second.begin(); yItter!=xIter->second.end();yItter++){

            int xMeshID=xIter->first;
            int yMeshID=yItter->first;

            if((xIter->first==0)||(xIter->first==9)||(yItter->first==0)||(yItter->first==9)){
                ref_EdgeSum+=yItter->second;
                tag_EdgeSum+=targetMesh[xIter->first][yItter->first];
            }
            if((xMeshID>2 && xMeshID <6)&&(yMeshID>2 && yMeshID <6)){
             ref_CenterSum+=yItter->second;
             tag_CenterSum+=targetMesh[xIter->first][yItter->first];
            }
        }
    }

    //calculate the ratio
    double NomalizeFact=(double )ref_EdgeSum/tag_EdgeSum;

    double ratio=(double)ref_CenterSum/(tag_CenterSum*NomalizeFact);
    std::cout<<ratio<<std::endl;
    return  ratio;
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
        runList[0]=2140;
        runList[1]=2140;
        runList[2]=2141;
        runList[3]=2291;
        runList[4]=2299;
        runList[5]=2300;
    }

    std::cout<<"Run Lis"<<std::endl;
    for (auto i = runList.begin(); i!=runList.end();i++){
        std::cout<<'\t'<<i->first<<"   "<<i->second<<std::endl;
    }

    std::map<int,double> reletiveThickness;
    reletiveThickness[0]=1.000;

    auto referenceMesh=MeshRasterCurrent(LoadrootFile(runList[0],folder));
    auto runIter=(runList.begin());
    for (runIter++;runIter!=runList.end();runIter++){
        std::cout<<runIter->second<<std::endl;
        auto meshInfor=MeshRasterCurrent(LoadrootFile(runIter->second,folder));

        double thickness=NormalizeRatio(referenceMesh,meshInfor);
        std::cout<<"runID"<<runIter->second<<"   "<<thickness<<std::endl;
        reletiveThickness[runIter->first]=thickness;
    }

    // plot the Error canvas

    TCanvas *targetThicknessCanv=new TCanvas(Form("Run%d",runList[0]),Form("Run%d",runList[0]),2080,1960);
    double thicknessX[runList.size()];
    double thicknessY[runList.size()];
    double thicknessXErr[runList.size()];
    double thicknessYErr[runList.size()];

    for (auto item = reletiveThickness.begin();item!=reletiveThickness.end();item++){
        thicknessX[item->first]=item->first;
        thicknessY[item->first]=item->second;
        thicknessXErr[item->first]=0;
        thicknessYErr[item->first]=item->second*0.015;
    }

    auto geprex=new TGraphErrors(runList.size()+1,thicknessX,thicknessY,thicknessXErr,thicknessYErr);
    geprex->GetYaxis()->SetRangeUser(0.8,1.2);
    geprex->GetXaxis()->SetRangeUser(-2,7);
    geprex->GetXaxis()->SetLimits(-2,7);
    geprex->SetTitle("Target D9-208Pb10-D10");
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
             TLatex *text=new TLatex(item->first,item->second+0.05,Form("%1.4f",item->second));
             text->SetTextColor(94);
             text->Draw("same");
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
}



///
/// \param runID, the run number that want to check file. It will look for the root file that matches the runID
/// \param folde, root file location
/// \param resultFolder: folders that used for save the generated plot
void TargetCheck(UInt_t runID, TString folder="/home/newdriver/pyQuant/prex_replayed/rootfile",TString resultFolder="./carbonCheck"){
    auto chain=LoadrootFile(runID,folder);
    MeshRasterCurrent(chain, false,resultFolder.Data());
}