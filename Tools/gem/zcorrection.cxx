{

    TCanvas canvas;
    canvas.Draw();
	TH1F *residual;
	
	TH1F *residualRMS=new TH1F ("residualRMS","residual",2000,-0.5,0.3);
	residualRMS->GetYaxis()->SetRangeUser(0,0.01);
	double deltaZmin=-0.5;
	double deltaZmax=0.3;
	double deltaZbin=0.001;
	
	int counter=0;
	for(double deltaZ=deltaZmin;deltaZ<deltaZmax;deltaZ+=deltaZbin){
	//
	residual=new TH1F ("residual","residual",1000,-0.05,0.05);
	std::string command(Form("(R.tr.x-RGEM.tr.x-%f*RGEM.tr.th)",deltaZ));
	T->Project("residual",command.c_str());
	residual->Fit("gaus","QM","",-0.02,0.02);
	double a[3];
	residual->GetFunction("gaus")->GetParameters(a);
	residualRMS->Fill(deltaZ,a[2]);
	if((counter++)%10==1){
	residualRMS->Draw("hist");
	canvas.Update();
	}
	residual->Delete();
	
	}
		
	residualRMS->Draw("hist");
	canvas.Update();


}
