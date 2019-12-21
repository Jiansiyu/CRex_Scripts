/*
 * clumAlign.cxx
 *
 *  Created on: Sep 22, 2019
 *      Author: newdriver
 *
 *  Used for  do the GEM position correction
 *
 *  Ideas:
 *
 *  	Analysis the .dat file
 *  	generate the align result
 *  	get the residual, and change the alignment parameter accordingly
 *  	'SED' the new alignment result to the database
 *  	repeat those steps
 *  	plot the  convergence curve
 *
 */

#include <string>
#include <iostream>
#include <TTree.h>
#include <TBranch.h>
#include <TFile.h>
#include <TString.h>
#include <vector>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <map>
#include <TPaveText.h>
#include <TStyle.h>
#include <ctime>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <TChain.h>
#include <TProfile.h>
class gemCorrectionlog {
public:
	gemCorrectionlog(){};

	gemCorrectionlog operator << (std::string a){
		std::ofstream myfile;
		myfile.open("log.log", std::ios::ate | std::ios::app);
		std::string logStri(Form("[%s] Log: %s",getDataTime().c_str(),a.c_str()));
		std::cout<<logStri.c_str()<<std::endl;
		myfile<<logStri.c_str()<<"\n";

		myfile.close();
		return *this;
	}
private:
	std::string getDataTime(){
	    std::time_t t = std::time(0);   // get time now
	    std::tm* now = std::localtime(&t);
	    return Form("%d-%d-%d %d:%d:%d",now->tm_year+1900,now->tm_mon+1,now->tm_mday,now->tm_hour,now->tm_min,now->tm_sec);
	}

};

gemCorrectionlog gemLog;


struct gemCorrection{

public:
	gemCorrection(){};
	gemCorrection(Int_t chamberID, double_t x, double_t y, double_t z) :
			X(x), Y(y), Z(z), ChamberID(chamberID) {
	};

	double_t GetX(){return X;};
	double_t GetY(){return Y;};
	double_t GetZ(){return Z;};
	Int_t GetChamberID(){return ChamberID;};

	inline void Print(){
		gemLog<<Form("Chamber :: %d,   correction -> %f;    %f;    %f", ChamberID,X,Y,Z);
	}

/*	inline gemCorrection operator - (gemCorrection const & a,gemCorrection const & b){
		if(a.ChamberID!=b.ChamberID){
			std::cout<<"[WORNING]:: The two chamber does not match, this is not expected"<< std::endl;
		}
		return(gemCorrection(a.ChamberID,a.X-b.X,a.Y-b.Y,a.Z-b.Z));
	}*/
protected:
	double_t X;
	double_t Y;
	double_t Z;
	Int_t ChamberID;

};

inline bool IsfileExist(const std::string &fname){
	struct stat buffer;
	return (stat (fname.c_str(), &buffer) == 0);
}

//call the X-term to analysis the data
bool RawAnalysis(Int_t runNumber){
	//call the bash script to analysis the data
	gemLog <<Form("Start Analysis run : %d", runNumber);
	system(Form("./scripts/analysisRaw.sh %d ",runNumber));

	gemLog<<"Finish analysis";
	return true;
}



inline std::vector<Int_t> GetChamberList(){
	std::vector<Int_t> gemChamberList;
	gemChamberList.push_back(1);
	gemChamberList.push_back(2);
	gemChamberList.push_back(3);
	gemChamberList.push_back(4);
	gemChamberList.push_back(5);
	gemChamberList.push_back(6);
	return gemChamberList;
}

inline std::map<Int_t,std::string> GetDBPositionNameList(Int_t runNumber, std::string xy="x"){
	std::map<Int_t,std::string> DBPositionNameList;
	std::string hrsprofix;
	if(runNumber>20000){
		hrsprofix="RGEM.rgems";
	}else{
		hrsprofix="LGEM.lgems";
	}
	for (auto ChamberID : GetChamberList()){
		DBPositionNameList[ChamberID]=Form("%s.%s%d.position",hrsprofix.c_str(),xy.c_str(),ChamberID);
		gemLog<< DBPositionNameList[ChamberID].c_str();
	}
	return DBPositionNameList;
}


std::map<Int_t, gemCorrection > ReadDatabase(std::string gemdbFname,Int_t runNumber){

	//Get the database
	std::map<Int_t, gemCorrection > result;
	std::ifstream filestream;
	std::map<std::string, std::string> DBBuffer;
	filestream.open( gemdbFname.c_str(), std::ifstream::in);
	std::string line;
	while( std::getline(filestream, line)){
		std::istringstream is_line(line);
		std::string key;
		  if( std::getline(is_line, key, '=') )
		  {
		    std::string value;
		    if( std::getline(is_line, value) )
		      DBBuffer[key.c_str()]=value.c_str();
		  }
	}

	// search for the keys
	auto DBPositionNameList=GetDBPositionNameList(runNumber);
	for (auto DBChamberParaIter=DBPositionNameList.begin(); DBChamberParaIter!=DBPositionNameList.end(); DBChamberParaIter++ ){
		for (auto bufferiter = DBBuffer.begin(); bufferiter!=DBBuffer.end();bufferiter++){
			auto bufferkey=bufferiter->first;
			std::size_t found=bufferkey.find(DBChamberParaIter->second);
			if(found!=std::string::npos){
				gemLog<<Form("%d   key ::%s   value:: %s ", DBChamberParaIter->first,DBChamberParaIter->second.c_str(),DBBuffer[bufferkey].c_str());
				std::istringstream iss(DBBuffer[bufferkey].c_str());
				std::vector<double_t> correction;
				correction.clear();
				for (std::string s; iss>>s;){
					correction.push_back(atof(s.c_str()));
				}
				if(correction.size()==3){
					result[DBChamberParaIter->first]= gemCorrection(DBChamberParaIter->first,correction.at(0),correction.at(1),correction.at(2));
				}
			}
		}
	}

	return result;

}

// analysis the ROOT file that generated, and get the residual. Used for correct the X-Y alignment
std::map<Int_t, gemCorrection > ResidualEvaluation(std::string rootFileName,std::string db_fname,Int_t runID, std::string workdir, std::string evalueTerm="transition"){

	std::cout<<"Analysis file "<< rootFileName.c_str()<<std::endl;
	std::map<Int_t, gemCorrection > gemResidual;

	// load the root files
	TChain *gemChain=new TChain("T");
	gemChain->Add(rootFileName.c_str());

	// need to load the database to get the z-dimension etc.
	std::map<Int_t, gemCorrection > gemDBValue=ReadDatabase(db_fname, runID);

	std::map<Int_t, TH1F*> gemCoordResidualX;
	std::map<Int_t, TH1F*> gemCoordResidualY;

	// for translation correction
	if(evalueTerm=="transition"){
		// get the GEM run list and project to histo
		std::vector<Int_t> ChamberList=GetChamberList();

		gemCoordResidualX.clear();
		gemCoordResidualY.clear();
		for(auto chamberID : ChamberList){
			gemCoordResidualX[chamberID]=new TH1F(Form("chamber%d_gemx_coordresidual",chamberID),Form("chamber%d_gemx_coordresidual",chamberID),1000,-0.05,0.05);
			gemCoordResidualY[chamberID]=new TH1F(Form("chamber%d_gemy_coordresidual",chamberID),Form("chamber%d_gemy_coordresidual",chamberID),1000,-0.05,0.05);
		}

		// project the data to histo
		for(auto chamberID : ChamberList){
			if(runID>20000){
				gemChain->Project(gemCoordResidualX[chamberID]->GetName(),Form("R.tr.x + R.tr.th * %f - RGEM.rgems.x%d.coord.pos",gemDBValue[chamberID].GetZ(),chamberID));
				gemChain->Project(gemCoordResidualY[chamberID]->GetName(),Form("R.tr.y + R.tr.ph * %f - RGEM.rgems.y%d.coord.pos",gemDBValue[chamberID].GetZ(),chamberID));
			}else{
				gemChain->Project(gemCoordResidualX[chamberID]->GetName(),Form("L.tr.x + L.tr.th * %f - LGEM.lgems.x%d.coord.pos",gemDBValue[chamberID].GetZ(),chamberID));
				gemChain->Project(gemCoordResidualY[chamberID]->GetName(),Form("L.tr.y + L.tr.ph * %f - LGEM.lgems.y%d.coord.pos",gemDBValue[chamberID].GetZ(),chamberID));
			}
			gemCoordResidualX[chamberID]->Fit("gaus");
			gemCoordResidualY[chamberID]->Fit("gaus");


			// get the fit parameters
			double_t MeanX=gemCoordResidualX[chamberID]->GetFunction("gaus")->GetParameter(1);// get the mean of the
			double_t MeanErrX=gemCoordResidualX[chamberID]->GetFunction("gaus")->GetParError(1);// get the mean of the

			// get the fit parameters
			double_t MeanY=gemCoordResidualY[chamberID]->GetFunction("gaus")->GetParameter(1); // get the mean of the
			double_t MeanErrY=gemCoordResidualY[chamberID]->GetFunction("gaus")->GetParError(1); // get the mean of the

			std::cout<<" X mean "<<MeanX<< "   err:"<<MeanErrX<<"   Y mean: "<<MeanY<<"  err:"<<MeanErrY<<std::endl;

			gemResidual[chamberID]=gemCorrection(chamberID,MeanX,MeanY,0);
		}
		// plot and return the evaluation result
		TCanvas *x=new TCanvas("x","x",1000,1000);
		TCanvas *y=new TCanvas("y","y",1000,1000);
		x->Divide(3,2);
		y->Divide(3,2);
		for(auto chamberID : ChamberList){
			x->cd(chamberID);
			gemCoordResidualX[chamberID]->Draw();
			y->cd(chamberID);
			gemCoordResidualY[chamberID]->Draw();
		}
		x->Update();
		y->Update();
		x->Draw();
		y->Draw();
		x->SaveAs(Form("%s/x.jpg",workdir.c_str()));
		y->SaveAs(Form("%s/y.jpg",workdir.c_str()));
	}else{ // for GEM rotation correction

	}
	return gemResidual;
}

// used for evaluation the Z alignment
// check the note book for detail
// Attention: need to align the x-y peer to align the z
//Question : How to evaluate the the boundary for the plot. Maybe need a good way to
std::map<Int_t, gemCorrection > ProfileZEvaluation(std::string rootFileName,std::string db_fname,Int_t runID, std::string workdir,std::map<Int_t, gemCorrection > gemDB, std::string evalueTerm="transition"){

	gemLog<<"Start Evaluate the Z dimension";
	gemLog<<Form("Analysis file ::: %s", rootFileName.c_str());

	std::map<Int_t, gemCorrection > gemResidual;

	// load the root files
	TChain *gemChain=new TChain("T");
	gemChain->Add(rootFileName.c_str());

	//need to load the z-dimension to do the projection, maybe latter will need to add to the root file
	std::string HRSarm = "R";
	if(runID < 20000){
		HRSarm="L";
	}

	std::map<Int_t, TH2F *> vdcGEMResidualhh_xth; // plot the
	std::map<Int_t, TH2F *> vdcGEMResidualhh_yph; // plot the

	// general cut on the GEMs and VDCs

	for (auto chamberID : GetChamberList()){
		vdcGEMResidualhh_xth[chamberID]=new TH2F(Form("chamber%d_residual_x_th",chamberID),Form("chamber%d_residual_x_th",chamberID),400,-0.02,0.02,400,-0.015,0.015);
		vdcGEMResidualhh_yph[chamberID]=new TH2F(Form("chamber%d_residual_y_ph",chamberID),Form("chamber%d_residual_y_ph",chamberID),400,-0.02,0.02,400,-0.015,0.015);

		std::string projectresidualX;
		std::string projectresidualY;
		if(HRSarm == "R"){
			projectresidualX=Form("R.tr.x + R.tr.th * %f - RGEM.rgems.x%d.coord.pos : R.tr.th ",gemDB[chamberID].GetZ(),chamberID);
			projectresidualY=Form("R.tr.y + R.tr.ph * %f - RGEM.rgems.y%d.coord.pos : R.tr.ph ",gemDB[chamberID].GetZ(),chamberID);

		}else{
			projectresidualX=Form("L.tr.x + L.tr.th * %f - LGEM.lgems.x%d.coord.pos : L.tr.th ",gemDB[chamberID].GetZ(),chamberID);
			projectresidualY=Form("L.tr.y + L.tr.ph * %f - LGEM.lgems.y%d.coord.pos : L.tr.ph ",gemDB[chamberID].GetZ(),chamberID);
		}

		gemLog<<"Plotting";
		gemLog<<projectresidualX.c_str();
		gemLog<<projectresidualY.c_str();

		// PROJECT THE DATA TO THE HISTOGRAM
		gemChain->Project(vdcGEMResidualhh_xth[chamberID]->GetName(),projectresidualX.c_str());
		gemChain->Project(vdcGEMResidualhh_yph[chamberID]->GetName(),projectresidualY.c_str());

	}

	TCanvas *Canvasprofxth=new TCanvas("profileX","profileX",600,800);
	TCanvas *Canvasprofyph=new TCanvas("profileY","profileY",600,800);
	Canvasprofxth->Divide(3,2);
	Canvasprofyph->Divide(3,2);
	Canvasprofxth->Draw();
	Canvasprofyph->Draw();

	for(auto chamberID : GetChamberList()){
		Canvasprofxth->cd(chamberID);
		vdcGEMResidualhh_xth[chamberID]->ProfileX()->Draw();
		Canvasprofyph->cd(chamberID);
		vdcGEMResidualhh_yph[chamberID]->ProfileX()->Draw();
	}

	//fit the functions and get the result
	for(auto chamberID : GetChamberList()){
		auto residualprofx=vdcGEMResidualhh_xth[chamberID]->ProfileX();
		residualprofx->Fit("pol1","","",-0.015,0.01);
		auto residualprofy=vdcGEMResidualhh_yph[chamberID]->ProfileX();
		residualprofy->Fit("pol1","","",-0.015,0.01);
		gemLog<<Form("Chamber %d :::   x fit   %f    y fit:  %f", chamberID,residualprofx->GetFunction("pol1")->GetParameter("p1"), residualprofy->GetFunction("pol1")->GetParameter("p1"));

		double correction=0.5*(residualprofx->GetFunction("pol1")->GetParameter("p1") + residualprofy->GetFunction("pol1")->GetParameter("p1") );
		gemResidual[chamberID]=gemCorrection(chamberID,gemDB[chamberID].GetX(),gemDB[chamberID].GetY(),gemDB[chamberID].GetZ()-correction);
	}

	// check the alignment result
	std::map<Int_t, TH2F *> vdcGEMResidualhh_xth_corrected; // plot the
	std::map<Int_t, TH2F *> vdcGEMResidualhh_yph_corrected; // plot the
	for (auto chamberID : GetChamberList()){
		vdcGEMResidualhh_xth_corrected[chamberID]=new TH2F(Form("chamber%d_residual_x_th_corrected",chamberID),Form("chamber%d_residual_x_th_corrected",chamberID),400,-0.02,0.015,400,-0.015,0.015);
		vdcGEMResidualhh_yph_corrected[chamberID]=new TH2F(Form("chamber%d_residual_y_ph_corrected",chamberID),Form("chamber%d_residual_y_ph_corrected",chamberID),400,-0.02,0.015,400,-0.015,0.015);

		std::string projectresidualX;
		std::string projectresidualY;
		if(HRSarm == "R"){
			projectresidualX=Form("R.tr.x + R.tr.th * %f - RGEM.rgems.x%d.coord.pos : R.tr.th ",gemResidual[chamberID].GetZ(),chamberID);
			projectresidualY=Form("R.tr.y + R.tr.ph * %f - RGEM.rgems.y%d.coord.pos : R.tr.ph ",gemResidual[chamberID].GetZ(),chamberID);

		}else{
			projectresidualX=Form("L.tr.x + L.tr.th * %f - LGEM.lgems.x%d.coord.pos : L.tr.th ",gemResidual[chamberID].GetZ(),chamberID);
			projectresidualY=Form("L.tr.y + L.tr.ph * %f - LGEM.lgems.y%d.coord.pos : L.tr.ph ",gemResidual[chamberID].GetZ(),chamberID);
		}

		gemLog<<"Plotting";
		gemLog<<projectresidualX.c_str();
		gemLog<<projectresidualY.c_str();

		// PROJECT THE DATA TO THE HISTOGRAM
		gemChain->Project(vdcGEMResidualhh_xth_corrected[chamberID]->GetName(),projectresidualX.c_str());
		gemChain->Project(vdcGEMResidualhh_yph_corrected[chamberID]->GetName(),projectresidualY.c_str());

		Canvasprofxth->cd(chamberID);
		auto profX=vdcGEMResidualhh_xth_corrected[chamberID]->ProfileX();
		profX->SetLineColor(3);
		//profX->Draw("same");
		Canvasprofyph->cd(chamberID);
		auto profY=vdcGEMResidualhh_yph_corrected[chamberID]->ProfileX();
		profY->SetLineColor(3);
		//profY->Draw("same");
	}
	Canvasprofxth->SaveAs(Form("run%d_xth.jpg",runID));
	Canvasprofyph->SaveAs(Form("run%d_yph.jpg",runID));

	return gemResidual;
}


// used for evaluation the Z alignment
// check the note book for detail
// Attention: need to align the x-y peer to align the z
//Question : How to evaluate the the boundary for the plot. Maybe need a good way to
std::map<Int_t, gemCorrection > ProfileZEvaluation1(std::string rootFileName,std::string db_fname,Int_t runID, std::string workdir,std::map<Int_t, gemCorrection > gemDB, std::string evalueTerm="transition"){

	gemLog<<"Start Evaluate the Z dimension";
	gemLog<<Form("Analysis file ::: %s", rootFileName.c_str());

	std::map<Int_t, gemCorrection > gemResidual;

	// load the root files
	TChain *gemChain=new TChain("T");
	gemChain->Add(rootFileName.c_str());

	//need to load the z-dimension to do the projection, maybe latter will need to add to the root file
	std::string HRSarm = "R";
	if(runID < 20000){
		HRSarm="L";
	}

	std::map<Int_t, TH2F *> vdcGEMResidualhh_xth; // plot the
	std::map<Int_t, TH2F *> vdcGEMResidualhh_yph; // plot the

	// general cut on the GEMs and VDCs

	for (auto chamberID : GetChamberList()){
		vdcGEMResidualhh_xth[chamberID]=new TH2F(Form("chamber%d_residual_x_th",chamberID),Form("chamber%d_residual_x_th",chamberID),400,-0.02,0.02,400,-0.015,0.02);
		vdcGEMResidualhh_yph[chamberID]=new TH2F(Form("chamber%d_residual_y_ph",chamberID),Form("chamber%d_residual_y_ph",chamberID),400,-0.02,0.02,400,-0.015,0.02);

		std::string projectresidualX;
		std::string projectresidualY;
		if(HRSarm == "R"){
			projectresidualX=Form("R.tr.x + R.tr.th * %f - RGEM.rgems.x%d.coord.pos : R.tr.th ",gemDB[chamberID].GetZ(),chamberID);
			projectresidualY=Form("R.tr.y + R.tr.ph * %f - RGEM.rgems.y%d.coord.pos : R.tr.ph ",gemDB[chamberID].GetZ(),chamberID);

		}else{
			projectresidualX=Form("L.tr.x + L.tr.th * %f - LGEM.lgems.x%d.coord.pos : L.tr.th ",gemDB[chamberID].GetZ(),chamberID);
			projectresidualY=Form("L.tr.y + L.tr.ph * %f - LGEM.lgems.y%d.coord.pos : L.tr.ph ",gemDB[chamberID].GetZ(),chamberID);
		}

		gemLog<<"Plotting";
		gemLog<<projectresidualX.c_str();
		gemLog<<projectresidualY.c_str();

		// PROJECT THE DATA TO THE HISTOGRAM
		gemChain->Project(vdcGEMResidualhh_xth[chamberID]->GetName(),projectresidualX.c_str());
		gemChain->Project(vdcGEMResidualhh_yph[chamberID]->GetName(),projectresidualY.c_str());

	}

	TCanvas *Canvasprofxth=new TCanvas("profileX1","profileX1",600,800);
	TCanvas *Canvasprofyph=new TCanvas("profileY1","profileY1",600,800);
	Canvasprofxth->Divide(3,2);
	Canvasprofyph->Divide(3,2);
	Canvasprofxth->Draw();
	Canvasprofyph->Draw();

	for(auto chamberID : GetChamberList()){
		Canvasprofxth->cd(chamberID);
		vdcGEMResidualhh_xth[chamberID]->Draw("zcol");
		Canvasprofyph->cd(chamberID);
		vdcGEMResidualhh_yph[chamberID]->Draw("zcol");
	}
/*
	//fit the functions and get the result
	for(auto chamberID : GetChamberList()){
		auto residualprofx=vdcGEMResidualhh_xth[chamberID]->ProfileX();
		residualprofx->Fit("pol1","","",-0.015,0.01);
		auto residualprofy=vdcGEMResidualhh_yph[chamberID]->ProfileX();
		residualprofy->Fit("pol1","","",-0.015,0.01);
		gemLog<<Form("Chamber %d :::   x fit   %f    y fit:  %f", chamberID,residualprofx->GetFunction("pol1")->GetParameter("p1"), residualprofy->GetFunction("pol1")->GetParameter("p1"));

		double correction=0.5*(residualprofx->GetFunction("pol1")->GetParameter("p1") + residualprofy->GetFunction("pol1")->GetParameter("p1") );
		gemResidual[chamberID]=gemCorrection(chamberID,gemDB[chamberID].GetX(),gemDB[chamberID].GetY(),gemDB[chamberID].GetZ()-correction);
	}

	// check the alignment result
	std::map<Int_t, TH2F *> vdcGEMResidualhh_xth_corrected; // plot the
	std::map<Int_t, TH2F *> vdcGEMResidualhh_yph_corrected; // plot the
	for (auto chamberID : GetChamberList()){
		vdcGEMResidualhh_xth_corrected[chamberID]=new TH2F(Form("chamber%d_residual_x_th_corrected",chamberID),Form("chamber%d_residual_x_th_corrected",chamberID),400,-0.02,0.015,400,-0.015,0.015);
		vdcGEMResidualhh_yph_corrected[chamberID]=new TH2F(Form("chamber%d_residual_y_ph_corrected",chamberID),Form("chamber%d_residual_y_ph_corrected",chamberID),400,-0.02,0.015,400,-0.015,0.015);

		std::string projectresidualX;
		std::string projectresidualY;
		if(HRSarm == "R"){
			projectresidualX=Form("R.tr.x + R.tr.th * %f - RGEM.rgems.x%d.coord.pos : R.tr.th ",gemResidual[chamberID].GetZ(),chamberID);
			projectresidualY=Form("R.tr.y + R.tr.ph * %f - RGEM.rgems.y%d.coord.pos : R.tr.ph ",gemResidual[chamberID].GetZ(),chamberID);

		}else{
			projectresidualX=Form("L.tr.x + L.tr.th * %f - LGEM.lgems.x%d.coord.pos : L.tr.th ",gemResidual[chamberID].GetZ(),chamberID);
			projectresidualY=Form("L.tr.y + L.tr.ph * %f - LGEM.lgems.y%d.coord.pos : L.tr.ph ",gemResidual[chamberID].GetZ(),chamberID);
		}

		gemLog<<"Plotting";
		gemLog<<projectresidualX.c_str();
		gemLog<<projectresidualY.c_str();

		// PROJECT THE DATA TO THE HISTOGRAM
		gemChain->Project(vdcGEMResidualhh_xth_corrected[chamberID]->GetName(),projectresidualX.c_str());
		gemChain->Project(vdcGEMResidualhh_yph_corrected[chamberID]->GetName(),projectresidualY.c_str());

		Canvasprofxth->cd(chamberID);
		auto profX=vdcGEMResidualhh_xth_corrected[chamberID]->ProfileX();
		profX->SetLineColor(3);
		profX->Draw("same");
		Canvasprofyph->cd(chamberID);
		auto profY=vdcGEMResidualhh_yph_corrected[chamberID]->ProfileX();
		profY->SetLineColor(3);
		profY->Draw("same");
	}
	Canvasprofxth->SaveAs(Form("run%d_xth.jpg",runID));
	Canvasprofyph->SaveAs(Form("run%d_yph.jpg",runID));*/

	return gemResidual;
}

struct prexTrackEvent{

	prexTrackEvent(int16_t ID, double_t x, double_t y, double_t z){
		DetectorID=ID;
		this->x=x;
		this->y=y;
		this->z=z;
	}
private:
	int16_t DetectorID;
	double_t x;
	double_t y;
	double_t z;

};
struct prexGEMevent{

};
struct prextrack{
public:
	prextrack(int entry, double_t vdcx, double_t vdcy, double_t vdcz,
			double_t vdcth, double_t vdcph, double_t gemx, double_t gemy,
			double_t gemz, double_t gemth, double_t gemph) {

	};
	~prextrack(){};

	void clone(prextrack track){

		this->entry=track.entry;

		this->vdc_x=track.vdc_x;
		this->vdc_y=track.vdc_y;
		this->vdc_z=track.vdc_z;
		this->vdc_ph=track.vdc_ph;
		this->vdc_th=track.vdc_th;

		this->gem_x=track.gem_x;
		this->gem_y=track.gem_y;
		this->gem_z=track.gem_z;
		this->gem_ph=track.gem_ph;
		this->gem_th=track.gem_th;
	}
private:
	int64_t entry;
	double_t vdc_x;
	double_t vdc_y;
	double_t vdc_z;
	double_t vdc_th;
	double_t vdc_ph;

	double_t gem_x;
	double_t gem_y;
	double_t gem_z;
	double_t gem_th;
	double_t gem_ph;
};
//New Z dimension align method
//Maybe need a lot of statistics to get plot
//
std::map<Int_t, gemCorrection > ProfileZEvaluation2(std::string rootFileName,std::string db_fname,Int_t runID, std::string workdir,std::map<Int_t, gemCorrection > gemDB, std::string evalueTerm="transition"){

	gemLog<<"Start Evaluate the Z dimension";
	gemLog<<Form("Analysis file ::: %s", rootFileName.c_str());

	std::map<Int_t, gemCorrection > gemResidual;

	// load the root files
	TChain *gemChain=new TChain("T");
	gemChain->Add("/home/newdriver/Storage/Research/PRex_Experiment/prex_analyzer/replay/Result/prexRHRS_21363_-1.root");
	gemChain->Add("/home/newdriver/Storage/Research/PRex_Experiment/prex_analyzer/replay/Result/prexRHRS_21363_-1_1.root");
	gemChain->Add("/home/newdriver/Storage/Research/PRex_Experiment/prex_analyzer/replay/Result/prexRHRS_21363_-1_2.root");

	//gemChain->Add(rootFileName.c_str());

	//need to load the z-dimension to do the projection, maybe latter will need to add to the root file
	std::string HRSarm = "R";
	if(runID < 20000){
		HRSarm="L";
	}

	auto Entries=gemChain->GetEntries();

	//



	return gemResidual;
}

//correct the database
void SED2Database(std::string gemdbFname, std::string key, std::string value, std::string backupName){
	gemLog<<"Start Modify the database...."<<Form("Modify the database :::: ./scripts/db_modify.sh %s  %s  %s  %s", gemdbFname.c_str(), key.c_str(), value.c_str(), backupName.c_str());
	system(Form("./scripts/db_modify.sh %s  %s  %s  %s", gemdbFname.c_str(), key.c_str(), value.c_str(), backupName.c_str()));
}

inline std::string GenerateResultFolder(Int_t runID){
	std::time_t t = std::time(0);   // get time now
	std::tm* now = std::localtime(&t);
	std::string saveFolder= Form("./Result/Run%0d_%04d_%02d_%02d_%02d_%02d_%02d",runID, now->tm_year + 1900, now->tm_mon + 1,
			now->tm_mday, now->tm_hour, now->tm_min, now->tm_sec);
	mkdir(saveFolder.c_str(),0777);
	return saveFolder.c_str();
}

inline Bool_t DatabaseModify(std::string DB_name,Int_t runID, std::map<Int_t, gemCorrection > gemDBValue){

	for(auto dbValueIter=gemDBValue.begin(); dbValueIter!=gemDBValue.end(); dbValueIter++){
		SED2Database(DB_name.c_str(),Form("RGEM.rgems.x%d.position ",dbValueIter->first), Form("\" %5.8f %5.8f %5.8f\"", dbValueIter->second.GetX(),dbValueIter->second.GetY(),dbValueIter->second.GetZ()),Form("%s_bk",DB_name.c_str()));
		SED2Database(DB_name.c_str(),Form("RGEM.rgems.y%d.position ",dbValueIter->first), Form("\" %5.8f %5.8f %5.8f\"", dbValueIter->second.GetX(),dbValueIter->second.GetY(),dbValueIter->second.GetZ()),Form("%s_bk",DB_name.c_str()));
	}
	system(Form("diff %s %s_bk",DB_name.c_str(), DB_name.c_str()));
	return true;
}

void clumAlign(Int_t runNumber){
	//
	//Generate the filename and folders that used
	std::string RawPath="/home/newdriver/PRex/PRex_Data/Raw";
	std::string RawfName;
	std::string gemDBName;
	std::string backupFolder;

	// analyze the raw data
	//RawAnalysis(runNumber);
	if(const char *DB_path=std::getenv("DB_DIR")){
		gemLog<<Form("DB PATH : %s",DB_path);
	}

	// create the folder to store the result
	if(runNumber>20000){
		gemDBName=Form("%s/db_RGEM.rgems.dat",std::getenv("DB_DIR"));
	}else{
		gemDBName=Form("%s/db_LGEM.lgems.dat",std::getenv("DB_DIR"));
	}
	gemLog<<Form("DB location :: %s", gemDBName.c_str());

int16_t temp=0;
while(temp++<=10)
{
	//Generate the folder used for save the optimized result
	std::string WorkFolder=GenerateResultFolder(runNumber);

	//read the database and used for generate the new adjustable parameters
	std::map<Int_t, gemCorrection > gemDBOriginalValue=ReadDatabase(gemDBName.c_str(),runNumber);
	// copy the current database to the new direction
	system(Form("cp -r %s %s/database_save.dat",gemDBName.c_str(), WorkFolder.c_str()));


	// call the script to analysis the data

	RawAnalysis(runNumber);

	std::string rootfilePath="./replay/Result";

	system(Form("cp %s/prexRHRS_%d.root %s",rootfilePath.c_str(),runNumber,WorkFolder.c_str()));

	//call the script to evaluate the current database
	std::map<Int_t, gemCorrection > gemResidual=
			ResidualEvaluation(Form("%s/prexRHRS_%d.root",rootfilePath.c_str(),runNumber),gemDBName.c_str(),runNumber,WorkFolder.c_str());


	std::map<Int_t, gemCorrection > gemNewDB;
	for(auto chamberID : GetChamberList()){

		gemNewDB[chamberID]= gemCorrection(gemDBOriginalValue[chamberID].GetChamberID(),
				gemDBOriginalValue[chamberID].GetX()+gemResidual[chamberID].GetX(),
				gemDBOriginalValue[chamberID].GetY()+gemResidual[chamberID].GetY(),
				gemDBOriginalValue[chamberID].GetZ());
	}

	DatabaseModify(gemDBName.c_str(),runNumber,gemNewDB);

}

}


// iteration scripts to get all the parameter aligned
void clumAlignIter(Int_t runNumber) {
	//Generate the filename and folders that used
	std::string RawPath = "/home/newdriver/PRex/PRex_Data/Raw";
	std::string RawfName;
	std::string gemDBName;
	std::string backupFolder;

	// analyze the raw data
	//RawAnalysis(runNumber);
	if (const char *DB_path = std::getenv("DB_DIR")) {
		gemLog << Form("DB PATH : %s", DB_path);
	}

	// create the folder to store the result
	if (runNumber > 20000) {
		gemDBName = Form("%s/db_RGEM.rgems.dat", std::getenv("DB_DIR"));
	} else {
		gemDBName = Form("%s/db_LGEM.lgems.dat", std::getenv("DB_DIR"));
	}
	gemLog << Form("DB location :: %s", gemDBName.c_str());

	int16_t temp = 0;

	while (temp++ <= 1000)
	{
		//Generate the folder used for save the optimized result
		std::string WorkFolder = GenerateResultFolder(runNumber);

		//read the database and used for generate the new adjustable parameters
		std::map<Int_t, gemCorrection> gemDBOriginalValue = ReadDatabase(
				gemDBName.c_str(), runNumber);
		// copy the current database to the new direction
		system(
				Form("cp -r %s %s/database_original.dat", gemDBName.c_str(),
						WorkFolder.c_str()));

		std::string rootfilePath = "./replay/Result";

		// call the script to analysis the data
		RawAnalysis(runNumber);


		system(
				Form("cp %s/prexRHRS_%d.root %s", rootfilePath.c_str(),
						runNumber, WorkFolder.c_str()));

		//call the script to evaluate the current database
		std::map<Int_t, gemCorrection> gemResidual = ResidualEvaluation(
				Form("%s/prexRHRS_%d.root", rootfilePath.c_str(), runNumber),
				gemDBName.c_str(), runNumber, WorkFolder.c_str());

		std::map<Int_t, gemCorrection> gemNewDB;
		for (auto chamberID : GetChamberList()) {

			gemNewDB[chamberID] = gemCorrection(
					gemDBOriginalValue[chamberID].GetChamberID(),
					gemDBOriginalValue[chamberID].GetX()
							+ gemResidual[chamberID].GetX(),
					gemDBOriginalValue[chamberID].GetY()
							+ gemResidual[chamberID].GetY(),
					gemDBOriginalValue[chamberID].GetZ());
		}

		DatabaseModify(gemDBName.c_str(), runNumber, gemNewDB);

		// copy the current database to the new direction
		system(Form("cp -r %s %s/database_new.dat", gemDBName.c_str(),
						WorkFolder.c_str()));

		// after align the x-y, re-analyze the data, and start evaluate the z-dimension
		RawAnalysis(runNumber);

		// copy the file and rename
		system(Form("cp %s/prexRHRS_%d.root %s/prexRHRS_%d_xy_corrected.root", rootfilePath.c_str(),
								runNumber, WorkFolder.c_str(),runNumber));

		// call the scripts to evaluate the current database Z-dimension
		auto gemZcorrectedDB=ProfileZEvaluation(Form("%s/prexRHRS_%d.root", rootfilePath.c_str(), runNumber),
				gemDBName.c_str(),runNumber,WorkFolder.c_str(), ReadDatabase(gemDBName.c_str(), runNumber)
				);

		DatabaseModify(gemDBName.c_str(), runNumber, gemZcorrectedDB);
		system(Form("cp -r %s %s/database_zcorrected_new.dat", gemDBName.c_str(),
						WorkFolder.c_str()));
		system(Form("cp  /home/newdriver/Storage/Research/PRex_Workspace/PREX-MPDGEM/PRexScripts/Tools/ClumAlign/*.jpg %s/",WorkFolder.c_str()));



	}

}

void alignmentcheck(TString rootFname, TString DBname, int16_t runNumber){

	// check the alignment result for x-y
	std::map<Int_t, gemCorrection> gemResidual = ResidualEvaluation(
			rootFname.Data(),
			DBname.Data(), runNumber, "./");


	// check the alignment result for z
	auto gemZcorrectedDB=ProfileZEvaluation(rootFname.Data(),
			DBname.Data(),runNumber,"./", ReadDatabase(DBname.Data(), runNumber)
			);
	auto gemZcorrectedDB1=ProfileZEvaluation1(rootFname.Data(),
			DBname.Data(),runNumber,"./", ReadDatabase(DBname.Data(), runNumber)
			);
}

