/*
 * gemvdcReplay.C
 *
 *  Created on: Oct 6, 2019
 *      Author: newdriver
 */
#include <stdlib.h>
#include <iostream>
#include <TROOT.h>

R__LOAD_LIBRARY(../bobParityData/libParity.so)
R__LOAD_LIBRARY(libTreeSearch-GEM.so)
R__LOAD_LIBRARY(/home/newdriver/Storage/Research/PRex_Workspace/PREX-MPDGEM/PRexAnalyzer/libprexCounting.so)


using namespace std;

Bool_t IsFileExist(const Char_t * fname);


void gemvdcReplay(Int_t runNo=0, Int_t lastevt=-1){
	char infile[300];
	  //static const char* replay_dir_prefix = "./%s";
	  THaApparatus* HRSR = new THaHRS("R","Right arm HRS");
	  gHaApps->Add( HRSR );

	  THaApparatus* HRSL = new THaHRS("L","Left arm HRS");
	  gHaApps->Add( HRSL );
	 //   HRSL->AddDetector( new TriFadcXscin("s0","s0 scintillator",kFALSE) );
	 //   HRSL->AddDetector( new THaVDC("vdc", "Vertical Drift Chamber"));

	  // add detectors that are not in the default config
	  //  HRSL->AddDetector( new THaScintillator("s0","Scintillator S0"));


	  // done in rootlogon.C:  gSystem->Load("libParity.so");
	  gHaApps->Add(new ParityData("P","HAPPEX Data"));

	  if(runNo>20000){
	  PREXStand *rprex= new PREXStand("RGEM", "GEM stand Right Arm");
	  MPDGEMTracker *rgems = new MPDGEMTracker("rgems", "Collection of GEMs in Right arm");
	  rprex->AddDetector(rgems);
	  gHaApps->Add(rprex);
	  gHaPhysics->Add( new THaGoldenTrack( "RGEM.gold", "RGEM Golden Track", "RGEM" ));
	  }
	  else
	  {
	  PREXStand *lprex= new PREXStand("LGEM", "GEM stand Right Arm");
	  MPDGEMTracker *lgems = new MPDGEMTracker("lgems", "Collection of GEMs in Right arm");
	  lprex->AddDetector(lgems);
	  gHaApps->Add(lprex);
	  gHaPhysics->Add( new THaGoldenTrack( "LGEM.gold", "LGEM Golden Track", "LGEM" ));
	  }

	  // e-p scattering
	  Double_t amu = 931.494*1.e-3;  // amu to GeV
	  Double_t he4 = 4*amu;
	  Double_t pb208 = 208*amu;
	  double_t c12   = 12 *amu;
	  //Double_t mass_tg  = he4; // Helium
	  Double_t mass_tg = c12;

	  // Electron kinematics
	  THaPhysicsModule* EK_r = new THaElectronKine("EK_R",
						       "Electron kinematics in HRS-R",
						       "R",mass_tg);
	  THaPhysicsModule* EK_l = new THaElectronKine("EK_L",
						       "Electron kinematics in HRS-L",
						       "L",mass_tg);

	  gHaPhysics->Add( EK_r );
	  gHaPhysics->Add( EK_l );

	  gHaPhysics->Add( new THaGoldenTrack( "R.gold", "HRS-R Golden Track", "R" ));
	  gHaPhysics->Add( new THaGoldenTrack( "L.gold", "HRS-L Golden Track", "L" ));


	  THaScalerEvtHandler *lscaler = new THaScalerEvtHandler("Left","HA scaler event type 140");
	  lscaler->SetDebugFile("scaler/LeftScaler.txt");
	  gHaEvtHandlers->Add (lscaler);

	  THaScalerEvtHandler *rscaler = new THaScalerEvtHandler("Right","HA scaler event type 140");
	  rscaler->SetDebugFile("scaler/RightScaler.txt");
	  gHaEvtHandlers->Add (rscaler);


	  THaAnalyzer* analyzer = new THaAnalyzer;

	  char outname[300];
	  if(runNo>20000)
	    sprintf(outname,"Result/prexRHRS_%d_%d.root",runNo, lastevt);
	  else
	    sprintf(outname,"Result/prexLHRS_%d_%d.root",runNo, lastevt);
	  analyzer->SetOutFile( outname );
	  analyzer->SetCutFile("./cuts/onlana.cuts");
	  if(runNo>20000)
	    analyzer->SetOdefFile("./outputDef/output_R.def");
	  else  analyzer->SetOdefFile("./outputDef/output_L.def");
	  analyzer->SetSummaryFile("./log/summary_example.log"); // optional

	  TString oldfilename="";
	  int found;
	  Int_t FirstEventNum = 0;
	  THaRun *oldrun=0, *run;
	  Bool_t exit = false;

	  if(lastevt>=0) lastevt+= FirstEventNum;

	  THaEvent* event = new THaEvent;

	  for (Int_t nsplit=0;!exit; nsplit++){
	    found = 0;
	    if(runNo>20000)
	      sprintf(infile,"/home/newdriver/PRex/PRex_Data/Raw/prexRHRS_%d.dat.%d",runNo,nsplit);
	    else
	      sprintf(infile,"/home/newdriver/PRex/PRex_Data/Raw/prexLHRS_%d.dat.%d",runNo,nsplit);

	    if( IsFileExist(infile) )
	      found = 1;

	    if( !found || oldfilename== infile)
	      {
		exit = true;
	      }
	    else
	      {
		oldfilename = infile;

		if (oldrun)
		  {
		    run = new THaRun(*oldrun);
		    run->SetFilename(infile);
		  }
		else
		  {
		    run = new THaRun(infile);
		  }

		if(lastevt>=0) run->SetLastEvent(lastevt);
		run->SetFirstEvent(FirstEventNum);

		try{
		  analyzer->Process(run);
		}
		catch( exception& e)
		  {
		    cerr << "Unhandled exception during replay: " << e.what() << endl;
		    cerr << "Exiting." << endl;
		    run->Close();
		    break;
		  }

		run->Close();
		if( !oldrun ) oldrun = run;

	      }//while file exist
	  }//nsplit loop
}

Bool_t IsFileExist(const Char_t * fname)
{
  fstream testfile;

  testfile.open(fname, ios_base::in);
  Bool_t isopen = testfile.is_open();
  testfile.close();

  return isopen;
}

