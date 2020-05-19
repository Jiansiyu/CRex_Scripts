/*
 * DrawText.C
 *
 *  Created on: Apr 3, 2020
 *      Author: newdriver
 */
#include <unistd.h>
#include <Riostream.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TMath.h>
#include <TProfile.h>
#include <TString.h>
#include <TText.h>
#include <TCut.h>
#include <TLine.h>
#include <map>
#include <string>
#include <TMathBase.h>
#include <TF1.h>
#include <TPaveText.h>

void DrawText(TString infile="/home/newdriver/Storage/Research/PRex_Workspace/PREX-MPDGEM/PRexScripts/Tools/PlotCut/Result/Cut20200322/Test/LHRS_EventNewNewRun/LargeDataSetVersion/WithOutMomCut/Sieve.Full.test", int NKineID=4){
	gStyle->SetOptStat(0);
	TCanvas *a=new TCanvas("a","a",1000,1000);
	a->SetLogy();

	a->cd();
	TTree *thm = new TTree("thm","MPD Histo Mode data");
	TString rform="kineID/I:kx/D:kth/D:ky/D:kph/D:kbeamE/D:kbeamX/D:kbeamY/D:EvtID/I:kDp/D:kelectronP/D";
	Int_t nrow = thm->ReadFile(infile.Data(),rform);
	thm->Print();
	TH1F *hRealMomentumCentralSieve = new TH1F(Form("hMomentumKi_centralSieve"),
					Form("hMomentumKin_centralSieve"), 500,
					2.16, 2.178);

	int Row=3;
	int Col=6;
	std::string sieveholeCut=Form("((kineID%98)%7==%d)&&((kineID%98)/7==%d)",Row,Col);

	thm->Project(hRealMomentumCentralSieve->GetName(),"kelectronP","kineID==45+98*3");

	TH1F *momentum=(TH1F *)hRealMomentumCentralSieve->Clone("momentum");
	hRealMomentumCentralSieve->Draw("hist");
			 // start the fit functions
	 //if (hRealMomentumCentralSieve[i]->GetEntries()>1000)
	 {
		 auto CGroundp=momentum->GetXaxis()->GetBinCenter(momentum->GetMaximumBin());
		 auto C1stp=CGroundp-0.00443891;
		 hRealMomentumCentralSieve->GetXaxis()->SetRangeUser(CGroundp-0.0044*3,CGroundp+0.0044*2);

		 double_t fgroudGausPar[3];
		 double_t ffirstGuasPar[3];
		 TF1 *fgroudGaus=new TF1("groudstatesgaus","gaus",CGroundp-0.0005,CGroundp+0.0005);
		 momentum->Fit("groudstatesgaus","R","ep",fgroudGaus->GetXmin(),fgroudGaus->GetXmax());
		 fgroudGaus->GetParameters(fgroudGausPar);

		 TF1 *ffirstGuas=new TF1 ("firststatesgaus","gaus",C1stp-0.0006,C1stp+0.0004);
		 momentum->Fit("firststatesgaus","R","ep",ffirstGuas->GetXmin(),ffirstGuas->GetXmax());
		 //ffirstGuas->Draw("same");
		 ffirstGuas->GetParameters(ffirstGuasPar);
		 // change the gause fit to cristal ball
		double_t fgroundCrystalballPar[5];
		TF1 *fgroundCrystalball=new TF1("fgroundCrystal","crystalball",fgroudGausPar[1]-0.0030,fgroudGaus->GetXmax()+0.0003);
		fgroundCrystalball->SetParameters(fgroudGausPar[0],fgroudGausPar[1],fgroudGausPar[2],1.64,1.1615);
		momentum->Fit("fgroundCrystal","R","same",fgroundCrystalball->GetXmin(),fgroundCrystalball->GetXmax());
		fgroundCrystalball->GetParameters(fgroundCrystalballPar);
		//fgroundCrystalball->Draw("same");
		double_t ffirstCrystalPar[5];
		TF1 *ffirstCrystal=new TF1("ffirstCrystal","crystalball",ffirstGuasPar[1]-0.0025,ffirstGuas->GetXmax());
		ffirstCrystal->SetParameters(ffirstGuasPar[0],ffirstGuasPar[1],ffirstGuasPar[2],1.64,1.1615);
		momentum->Fit("ffirstCrystal","R","ep",ffirstCrystal->GetXmin(),ffirstCrystal->GetXmax());
		ffirstCrystal->GetParameters(ffirstCrystalPar);

		// fit together
		double_t fCrystalMomentumPar[10];
		TF1 *fCrystalMomentum=new TF1("fCrystalMomentum","crystalball(0)+crystalball(5)",ffirstCrystal->GetXmin(),fgroundCrystalball->GetXmax());
		std::copy(fgroundCrystalballPar,fgroundCrystalballPar+5,fCrystalMomentumPar);
		std::copy(ffirstCrystalPar,ffirstCrystalPar+5,fCrystalMomentumPar+5);
		fCrystalMomentum->SetParameters(fCrystalMomentumPar);
		momentum->Fit("fCrystalMomentum","","",fCrystalMomentum->GetXmin(),fCrystalMomentum->GetXmax());
		fCrystalMomentum->Draw("same");
		fCrystalMomentum->GetParameters(fCrystalMomentumPar);

		TPaveText *pt = new TPaveText(0.1,0.8,0.3,0.9,"NDC");
		pt->AddText(Form("%1.3f MeV (%2.2f\%%)",1000.0*(fCrystalMomentumPar[1]-fCrystalMomentumPar[6]),100.0*abs(abs(fCrystalMomentumPar[1]-fCrystalMomentumPar[6])-0.00443891)/0.00443891));
		pt->Draw("same");

	 }


}


