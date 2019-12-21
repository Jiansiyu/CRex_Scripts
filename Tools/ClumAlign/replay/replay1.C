/*
 * A more flexiable replay script for PRexII experiment
 * Author: Siyu
 * sj9va@virginia.edu
 */


//#include <TROOT.h>
//#include <THaApparatus.h>
//#include <PREXStand.h>
//#include <MPDGEMTracker.h>
//#include <THaHRS.h>
//#include <ParityData.h>
//#include <THaGoldenTrack.h>
//#include <THaAnalyzer.h>
//#include <THaTrackingDetector.h>
//#include <THaRun.h>
//#include <THaScalerEvtHandler.h>
//#include <THaElectronKine.h>
//#include <TString.h>
//#include <fstream>
//#include <vector>
//#include <THaIdealBeam.h>
//#include <THaDetector.h>
//#include <THaReactionPoint.h>

R__LOAD_LIBRARY(/home/newdriver/Storage/Research/PRex_Experiment/prex_analyzer/bobParityData/libParity.so)
R__LOAD_LIBRARY(libTreeSearch-GEM.so)
R__LOAD_LIBRARY(/home/newdriver/Storage/Research/PRex_Workspace/PREX-MPDGEM/PRexAnalyzer/libprexCounting.so)

using namespace std;


inline  Bool_t IsFileExist(TString fname){
	ifstream f(fname.Data());
	return f.good();
}

inline Bool_t CreateFolder(TString folder){
return 1;
}


void replay1(Int_t runNo, Int_t lastevt=-1,Bool_t EnableGEM=kTRUE, TString defaultDataPath="/home/newdriver/PRex/PRex_Data/Raw",TString defaultSavePath="./Result", TString fileNamePatter="prex%sHRS_%d.dat.%d", Int_t MaxFilesplit=100){

	if (lastevt <0){
		std::cout<<"Will Replay all the events in files"<<std::endl;
	}{
		std::cout<<"Will replay "<<lastevt<<"  events in file"<<std::endl;
	}

	if(system("ls $DB_DIR")!=0){
		std::cout<<" [ERROR]:: CANNOT FIND database in folder "<<getenv("DB_DIR")<<std::endl;
		exit(-1);
	}

	TString HRSarm="R";
	if(runNo<20000){
		HRSarm="L";
	}

	// search the file in folder and skip files that did not exist
	std::vector<TString> RawFileList;
	for(Int_t split=0; split<MaxFilesplit; split++){
		TString filename=defaultDataPath+"/"+Form(fileNamePatter.Data(),HRSarm.Data(),runNo,split);
		if(IsFileExist(filename.Data())){
			RawFileList.push_back(filename);
		}else{
			break;
		}
	}

	std::cout<<"oooooooOOOOOOOOO00000000OOOOOOOOOooooooooo"<<std::endl;
	std::cout<<" Start Analysis The Following files:"<<std::endl;
	for (auto file :RawFileList ){
		std::cout<<file.Data()<<std::endl;
	}


	//Need to check whether folder exist, if not need to create one
	THaApparatus *HRSAppa;
	if(HRSarm=="R"){

		HRSAppa=new THaHRS("R","Right arm HRS");
	}else if (HRSarm=="L") {
		HRSAppa=new THaHRS("L","Left arm HRS");
	}
	gHaApps->Add(HRSAppa);
	gHaApps->Add(new ParityData("P","HAPPEX Data"));

	// add the GEM detector to the apparatus
	if(EnableGEM && (HRSarm=="R")){

		PREXStand  *rprex= new PREXStand("RGEM", "GEM stand Right Arm");
		MPDGEMTracker *rgems = new MPDGEMTracker("rgems", "Collection of GEMs in Right arm");
		rprex->AddDetector(rgems);
		gHaApps->Add(rprex);
		gHaPhysics->Add( new THaGoldenTrack( "RGEM.gold", "RGEM Golden Track", "RGEM" ));

	}else if (EnableGEM && (HRSarm=="L")) {
		PREXStand *lprex= new PREXStand("LGEM", "GEM stand Right Arm");
		MPDGEMTracker *lgems = new MPDGEMTracker("lgems", "Collection of GEMs in Right arm");
		lprex->AddDetector(lgems);
		gHaApps->Add(lprex);
		gHaPhysics->Add( new THaGoldenTrack( "LGEM.gold", "LGEM Golden Track", "LGEM" ));
	}

	// start the Ideal Beam configure
    // Ideal beam (perfect normal incidence and centering)
	THaIdealBeam* ib = new THaIdealBeam("IB", "Ideal beam");
    gHaApps->Add(ib);

    // e-p scattering
    Double_t amu = 931.494*1.e-3;  // amu to GeV
    Double_t he4 = 4*amu;
    Double_t pb208 = 208*amu;
    double_t c12   = 12 *amu;

    Double_t mass_tg = c12;

    // Single-arm electron kinematics for the one spectrometer we have set up.
    // We assume a carbon-12 target (12 AMU)
    gHaPhysics->Add( new THaPrimaryKine( Form("%s.ekine",HRSarm.Data()),
                                       Form("%sHRS electron kinematics",HRSarm.Data()),
									   HRSarm.Data(), 0.511e-3, mass_tg ));

    // Vertex position calculated from RHRS golden track and ideal beam
    // (will have poor resolution if raster is on)
    gHaPhysics->Add( new THaReactionPoint( Form("%s.vx",HRSarm.Data()),
                                         Form("Vertex %s",HRSarm.Data()),
                                         HRSarm.Data(), "IB" ));

    // Electron kinematics
	if (HRSarm == "R") {
		THaPhysicsModule* EK_r = new THaElectronKine("EK_R",
				"Electron kinematics in HRS-R", "R", mass_tg);
		gHaPhysics->Add(EK_r);
	    gHaPhysics->Add( new THaGoldenTrack( "R.gold", "HRS-R Golden Track", "R" ));

	} else if (HRSarm == "L") {

		THaPhysicsModule* EK_l = new THaElectronKine("EK_L",
				"Electron kinematics in HRS-L", "L", mass_tg);
		gHaPhysics->Add(EK_l);
	    gHaPhysics->Add( new THaGoldenTrack( "L.gold", "HRS-L Golden Track", "L" ));
	}

  if(HRSarm == "L"){
    THaScalerEvtHandler *lscaler = new THaScalerEvtHandler("Left","HA scaler event type 140");
    lscaler->SetDebugFile("./log/LeftScaler.txt");
    gHaEvtHandlers->Add (lscaler);
  }else{
    THaScalerEvtHandler *rscaler = new THaScalerEvtHandler("Right","HA scaler event type 140");
    rscaler->SetDebugFile("./log/RightScaler.txt");
    gHaEvtHandlers->Add (rscaler);
  }


  THaAnalyzer* analyzer = new THaAnalyzer;

  char outname[300];
  if(HRSarm == "R")
    sprintf(outname,"Result/prexRHRS_%d_%d.root",runNo, lastevt);
  else
    sprintf(outname,"Result/prexLHRS_%d_%d.root",runNo, lastevt);

  analyzer->SetOutFile( outname );
  analyzer->SetCutFile("./cuts/onlana.cuts");

  if(HRSarm == "R")
    analyzer->SetOdefFile("./outputDef/output_R.def");
  else  analyzer->SetOdefFile("./outputDef/output_L.def");

  analyzer->SetSummaryFile("./log/summary_example.log"); // optional




  TString oldfilename="";
  int found;
  Int_t FirstEventNum = 0;
  THaRun *oldrun=0, *run;
  Bool_t exit = false;

  THaEvent* event = new THaEvent;
  for (auto infile : RawFileList){

	if(!oldrun){
		run = new THaRun(infile.Data());
		oldrun = run;
	}else{
	    run = new THaRun(*oldrun);
	    run->SetFilename(infile.Data());
	}

	if(lastevt>=0) run->SetLastEvent(lastevt);

	run->SetFirstEvent(FirstEventNum);

	try {
		analyzer->Process(run);

	 } catch (exception& e) {
		cerr << "Unhandled exception during replay: " << e.what() << endl;
		cerr << "Exiting." << endl;
		run->Close();
		break;
	  }

	run->Close();

  }

}


void replay(Int_t runNo=0, Int_t lastevt=-1, Int_t splitID=0){
  //  R. Michaels, May 2014
  //  Steering script for Hall A analyzer
 /* int runNo, hrs;/prex_analyzer
  char infile[300];
  cout<<"Run number please (-1 to exit):";
  cin>>runNo;
  cout<<"HRS? (1 Left, 2 Right) ";
  cin>>hrs;*/
  
  /*if(runNo<20000){
    gSystem->Load("../ParityData/libParity_LHRS.so"); 
  }else{
    gSystem->Load("../ParityData/libParity_RHRS.so");
   }*/
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
  Double_t pb208 = 208*amu;
  
  Double_t mass_tg  = pb208; // Helium

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
  
  THaEvent* event = new THaEvent;
  for (Int_t nsplit=0;nsplit<1;nsplit++){  
  if(runNo>20000)
  sprintf(infile,"/home/newdriver/PRex/PRex_Data/Raw/prexRHRS_%d.dat.%d",runNo,splitID);
  //sprintf(infile,"/adaqfs/home/a-onl/siyu/RawData/prexRHRS_%d.dat.%d",runNo,nsplit);
  else
  sprintf(infile,"/home/newdriver/PRex/PRex_Data/Raw/prexLHRS_%d.dat.%d",runNo,splitID);
  cout<<"replay: Try file "<<infile<<endl;
  THaRun *run;
  run = new THaRun(infile);
 // run->SetDataRequired( THaRunBase::kDate );
  TDatime now;
  run->SetDate( now );
  run->SetDataRequired( 0 );  // it was: ( THaRunBase::kDate );

  analyzer->EnableBenchmarks();
  analyzer->SetEvent( event );
  char outname[300];
  if(runNo>20000)
  sprintf(outname,"Result/prexRHRS_%d.root",runNo, splitID);
  else
  sprintf(outname,"Result/prexLHRS_%d.root",runNo, splitID);
  analyzer->SetOutFile( outname );
  analyzer->SetCutFile("cuts/onlana.cuts");
 if(runNo>20000)
  analyzer->SetOdefFile("outputDef/output_R.def");
 else  analyzer->SetOdefFile("outputDef/output_L.def");
   analyzer->SetSummaryFile("log/summary_example.log"); // optional
//
  run->SetLastEvent(lastevt);   // Number of events to process
  analyzer->Process(*run);
  system("touch finish.txt");
  }
}
