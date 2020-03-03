/*
 * GEM_rawView.C
 *
 *  Created on: Feb 23, 2020
 *      Author: newdriver
 */
#include <string>
#include <TChain.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <map>
#include <TLatex.h>
#include <iostream>
int GEM_rawView(const char* filename, int evtID=0){
	TChain *chain = new TChain("GEMHit");
	chain->Add(filename);
	std::map<int, TH1F*> striphX;
	std::map<int, TH1F*> striphY;


	TCanvas * c1=new TCanvas("c1","c1",1000,1200);
	c1->Divide(2,6);
	c1->Draw();
	// 6 GEM detectors
	// project to the histogram and plot to canvas
	TLatex *t2[6], *t1[6];
	while(evtID++<chain->GetEntries()){
		std::cout<<"Working on Event "<<evtID<<std::endl;
		c1->SetTitle(Form("GEMRawDisplay_evtID%d",evtID));
		for(int i =0; i < 3; i ++){
				striphX[i]=new TH1F(Form("striph_moudle%d_x",i),Form("striph_moudle%d_x",i),600,0,600);
				striphY[i]=new TH1F(Form("striph_moudle%d_y",i),Form("striph_moudle%d_y",i),600,0,600);
		}
		for(int i =3; i < 6; i ++){
				striphX[i]=new TH1F(Form("striph_moudle%d_x",i),Form("striph_moudle%d_x",i),1600,0,1600);
				striphY[i]=new TH1F(Form("striph_moudle%d_y",i),Form("striph_moudle%d_y",i),1600,0,1600);
		}
		for(int i =0; i<6; i ++){
			// project X
			std::cout<<"\tProject X ... "<<std::endl;
			chain->Project(striphX[i]->GetName(),Form("strip"),Form("detID==%d && planeID==0 && evtID==%d",i,evtID));
			c1->cd(2*i+1);
			striphX[i]->Draw("HIST");
			c1->cd(2*i+1)->Update();
			// search for the fired strips
			std::string stripsXstr;
			std::cout<<"\t\t";
			for (int element =0; element < striphX[i]->GetNbinsX();element++){
				if(striphX[i]->GetBinContent(element)!=0){
				std::cout<<(element)<<", ";
				stripsXstr+=Form("%d  ",element);}
			}
			std::cout<<std::endl;
			// plot the seached first strips
			t1[i] = new TLatex(c1->cd(2*i+1)->GetTopMargin(),c1->cd(2*i+1)->GetLeftMargin(), Form("%s",stripsXstr.c_str()));
					t1[i]->SetTextSize(0.1);
					t1[i]->SetTextAlign(12);
					t1[i]->SetTextColor(2);
					t1[i]->Draw("same");

			c1->cd(2*i+1)->Update();
			std::cout<<"\tProject X  Done"<<striphX[i]->GetEntries()<<std::endl;

			// project Y
			std::cout<<"\tProject Y  "<<std::endl;
			chain->Project(striphY[i]->GetName(),Form("strip"),Form("detID==%d && planeID==1 && evtID==%d",i,evtID));
			c1->cd(2*i+2);
			striphY[i]->Draw("HIST");
			c1->cd(2*i+2)->Update();
			std::string stripsYstr;
			std::cout << "\t\t";
			for (int element = 0; element < striphY[i]->GetNbinsX();
					element++) {
				if (striphY[i]->GetBinContent(element) != 0){
					std::cout << (element) << ", ";
				stripsYstr+=Form("%d  ",element);}
			}
			std::cout << std::endl;
			// plot the seached first strips
			t2[i] = new TLatex(c1->cd(2*i+2)->GetLeftMargin(),c1->cd(2*i+2)->GetTopMargin(), Form("%s",stripsYstr.c_str()));
					t2[i]->SetTextSize(0.1);
					t2[i]->SetTextAlign(12);
					t2[i]->SetTextColor(2);
					t2[i]->Draw("same");
			c1->cd(2*i+2)->Update();
			std::cout<<"\tProject X  Done"<<std::endl;


		}

		c1->Update();
		getchar();

		for(int i =0; i < 6; i ++){
			striphX[i]->Delete();
			striphY[i]->Delete();
			t1[i]->Delete();
			t2[i]->Delete();
		}

	}

	return 1;
}



