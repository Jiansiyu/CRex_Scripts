#include <TROOT.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TF1.h>
#include <TStyle.h>
#include <TColor.h>
#include <TGaxis.h>
#include <TCanvas.h>
#include <TPaveLabel.h>
#include <TPaveText.h>
#include <TGraphErrors.h>
#include <fstream>
#include <iostream>
#include <TCut.h>

void Qsq(TString targ){ 

  gROOT->Reset();
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1111);
  gStyle->SetTitleYOffset(1.3);
//  gStyle->SetStatH(0.3);
//  gStyle->SetStatW(0.3);
//  gStyle->SetTitleH(0.1);
//  gStyle->SetTitleW(0.3);
  gROOT->ForceStyle();
  TGaxis::SetMaxDigits(3);
  

  char* prex_path = getenv("PREXROOT");

 
  const double z_scale = 100;
  // Cut for FBUS
  double upadc_cutL = 485;
  double upadc_cutR = 505;
  double USLPMT = 5370;
  double USLgain = 7.5;//7.5e5
  double USRPMT = 5401;
  double USRgain = 6.46;//6.46e5
  // RHRS = R.*
  // LHRS = L.*
  int run_num = (int)T->GetMaximum("fEvtHdr.fRun");
//    int run_num = 1983;
      TChain *T = new TChain("T");

  if(run_num < 10000){
      //LHRS
//      T->Add(Form("prex_counting/prexLHRS_%d_50000.root",run_num));
       T->Add(Form("%s/prexLHRS_%d_200000_test.root",prex_path,run_num));
      TCut trig_cut = "";
//      TCut trig_cut = "fEvtHdr.fEvtType==2";
//      TCut trig_cut = "(P.evtypebits&2)==2";
      // only one cluster in each VDC plane
      TCut vdc_cut = "L.vdc.u1.nclust==1&&L.vdc.v1.nclust==1&&L.vdc.u2.nclust==1&&L.vdc.v2.nclust==1";
      //track cut on theta and phi
      TCut tr_cut = "(L.tr.th[0]<0.05&&L.tr.th[0]>-0.2&&L.tr.ph[0]<0.1&&L.tr.ph[0]>-0.1)";
      //track cut on target
      TCut tg_cut = "(L.tr.tg_th[0]<0.055&&L.tr.tg_th[0]>-0.055&&L.tr.tg_ph[0]>-0.018&&L.tr.tg_ph[0]<0.026)";

      // FBUS adccuts
      TCut adc_cut_up = Form("P.upQadcL> %f", upadc_cutL);
      //FBUS <adccuts
      TCut adc_cut_up_ped = Form("P.upQadcL< %f", upadc_cutL);
      //cut on radiative tail
//      TCut x_cut = "(L.tr.x[0]+0.9*L.tr.th[0]) > -0.068";

      TCut x_cut = "";
      TCut cut_basic = trig_cut + vdc_cut;
      TCut cut = cut_basic + tr_cut + tg_cut;
      TCut cut_wadc = cut + adc_cut_up;
      TCut cut_wadc_ped = cut + adc_cut_up_ped;
 

      TCanvas *c2 = new TCanvas("c2","c2",1000,1000);
      TH1F* qsq = new TH1F("qsq", Form("LHRS %s Q^{2} for Run %d",targ.Data(),run_num), 100, 0, 0.015);


      T->Project(qsq->GetName(), "EK_L.Q2[0]", cut_wadc); 
      qsq->SetXTitle("Q^{2} (GeV/c)^{2}");
      qsq->SetLineColor(kBlack);
      qsq->Draw();
      c2->SaveAs(Form("QsqL_%s_run%d.jpg",targ.Data(),run_num));


    }else{
      // RHRS
      T->Add(Form("%s/prexRHRS_%d_200000_test.root",prex_path,run_num));
      TCut trig_cut = "";
      TCut vdc_cut = "R.vdc.u1.nclust==1&&R.vdc.v1.nclust==1&&R.vdc.u2.nclust==1&&R.vdc.v2.nclust==1";
      TCut tr_cut = "(R.tr.th[0]<0.05&&R.tr.th[0]>-0.2&&R.tr.ph[0]<0.1&&R.tr.ph[0]>-0.1)";
      TCut tg_cut = "(R.tr.tg_th[0]<0.05&&R.tr.tg_th[0]>-0.055&&R.tr.tg_ph[0]>-0.026&&R.tr.tg_ph[0]<0.016)";

      // FBUS adccuts
      TCut adc_cut_up = Form("P.upQadcR> %f", upadc_cutR);
      //FBUS <adccuts
      TCut adc_cut_up_ped = Form("P.upQadcR< %f", upadc_cutR);

      TCut x_cut = "";
      TCut cut_basic = trig_cut + vdc_cut;
      TCut cut = cut_basic + tr_cut + tg_cut;

      TCut cut_wadc = cut + adc_cut_up;
      TCut cut_wadc_ped = cut + adc_cut_up_ped;

      TCanvas *c2 = new TCanvas("c2","c2",1000,1000);
      TH1F* qsq = new TH1F("qsq", Form("RHRS %s Q^{2} for Run %d",targ.Data(),run_num), 100, 0, 0.015);

      T->Project(qsq->GetName(), "EK_R.Q2[0]", cut_wadc);                 
      qsq->SetXTitle("Q^{2} (GeV/c)^{2}");
      qsq->SetLineColor(kBlack);
      qsq->Draw();
      c2->SaveAs(Form("QsqR_%s_run%d.jpg",targ.Data(),run_num));


    }

}

