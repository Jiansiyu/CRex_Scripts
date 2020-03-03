/*
 * standlone version of replay scripts
 *
 */

#include <THaApparatus.h>
#include <THaHRS.h>
#include <ParityData.h>
#include <THaPhysicsModule.h>
#include <THaElectronKine.h>
#include <THaGoldenTrack.h>
#include <THaScalerEvtHandler.h>
#include <THaAnalyzer.h>
#include <THaRun.h>

#include <TROOT.h>
#include <fstream>
#include <string.h>
#include <string>
#include <iostream>
#include <limits.h>
#include <unistd.h>

struct HRSrunInfor{
	//list of target
	  // e-p scattering
	  Double_t amu = 931.494028*1.e-3;  // amu to GeV
	  Double_t he4 = 4*amu;
	  Double_t pb208 = 208*amu;
	  Double_t c12 = 12*amu-.511e-3*6;
	  Double_t ca40= 40*amu;
	  Double_t ca48= 48*amu;

	  // config the path used in the decoding
	  // config file > auto config
	  std::string currentPath;
	  std::string runPath;
	  std::string configPath;
	  std::string RawDataPath;
	  std::string LogPath;
	  std::string ResultSavePath;

public:

	virtual HRSrunInfor(std::string="runConfig.txt"){};
	virtual ~HRSrunInfor(){};
	double getTargetMass(){
		return c12;
	};

	void ReadConfig(void){

	};

	std::string GetRunConfigPath(){
		return configPath.c_str();
	}
	std::string GetRawDataPath(){
		return RawDataPath.c_str();
	}
	std::string GetLogPath(){
		return LogPath.c_str();
	}
	std::string GetResultPath(){
		return ResultSavePath.c_str();
	}

};

HRSrunInfor *runInfor;


Bool_t IsFileExist(const Char_t * fname){
	std::fstream testfile;
	testfile.open(fname, std::ios_base::in);
	Bool_t isopen = testfile.is_open();
	testfile.close();
	return isopen;
}


std::string getexepath()
{
  char result[ PATH_MAX ];
  ssize_t count = readlink( "/proc/self/exe", result, PATH_MAX );
  return std::string( result, (count > 0) ? count : 0 );
}

//try to find the run path if needed
void EnvCheck(){
	//check whether run on the ifarm or run locally
	std::string uname(system("whoami"));

	// if run on a local computer
	if(uname == "newdriver"){  // run locally,
		std::string CurrentPath=getexepath();

	}else{  // if run on JLab IFARM


	}
}

void main(Int_t runNo=0, std::string Config="Config.txt",std::string workDir="", Int_t lastevt=-1){
	std::string HRS="L";
	if(runNo>20000){
		HRS="R";
	}
	runInfor=new HRSrunInfor();

	char infile[300];
	THaApparatus *HRSApparatus;
	if(HRS=="R"){
		HRSApparatus=new THaHRS("R","Right arm HRS");
	}else{
		HRSApparatus=new THaHRS("L","Left arm HRS");
	}
	gHaApps->Add(HRSApparatus);
	gHaApps->Add(new ParityData("P","HAPPEX Data"));
	//electron kinematics
	THaPhysicsModule *EK_HRS;
	if(HRS=="R"){
		EK_HRS= new THaElectronKine("EK_R",
			       "Electron kinematics in HRS-R",
			       "R",runInfor->getTargetMass());
	}else{

		EK_HRS=new THaElectronKine("EK_L",
			       "Electron kinematics in HRS-L",
			       "L",runInfor->getTargetMass());
	}
	gHaPhysics->Add(EK_HRS);
	if(HRS=="R"){
		gHaPhysics->Add( new THaGoldenTrack( "R.gold", "HRS-R Golden Track", "R" ));
	}else{
		gHaPhysics->Add( new THaGoldenTrack( "L.gold", "HRS-L Golden Track", "L" ));
	}

	// work DIR
	if(HRS=="L"){
	  THaScalerEvtHandler *lscaler = new THaScalerEvtHandler("Left","HA scaler event type 140");
	  lscaler->SetDebugFile(Form("%s/LeftScaler.txt",runInfor->GetLogPath().c_str()));
	  gHaEvtHandlers->Add (lscaler);
	}else{
	  THaScalerEvtHandler *rscaler = new THaScalerEvtHandler("Right","HA scaler event type 140");
	  rscaler->SetDebugFile(Form("%s/RightScaler.txt",runInfor->GetLogPath().c_str()));
	  gHaEvtHandlers->Add (rscaler);
	}

	THaAnalyzer *analyzer=new THaAnalyzer;
	char outname[300];
	if (HRS=="R")
		sprintf(outname, "Result/prexRHRS_%d_%d.root", runNo, lastevt);
	else
		sprintf(outname, "Result/prexLHRS_%d_%d.root", runNo, lastevt);
	analyzer->SetOutFile(outname);
	analyzer->SetCutFile(Form("%s/cuts/onlana.cuts",runInfor->GetRunConfigPath().c_str()));
	if (runNo > 20000)
		analyzer->SetOdefFile(Form("%s/outputDef/output_R.def",runInfor->GetRunConfigPath().c_str()));
	else
		analyzer->SetOdefFile(Form("%s/outputDef/output_L.def",runInfor->GetRunConfigPath().c_str()));
	analyzer->SetSummaryFile(Form("%s/log/summary_example.log",runInfor->GetRunConfigPath().c_str())); // optional

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
	      sprintf(infile,"%s/prexRHRS_%d.dat.%d",runInfor->GetRawDataPath().c_str(),runNo,nsplit);
	    else
	      sprintf(infile,"%s/prexLHRS_%d.dat.%d",runInfor->GetRawDataPath().c_str(),runNo,nsplit);

	    if( IsFileExist(infile) )
	      found = 1;

	    if( !found || oldfilename== infile) // did not find the file exit
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
		else // if this is the first run
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

