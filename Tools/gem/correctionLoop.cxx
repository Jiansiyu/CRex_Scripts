{
 TCanvas canvas;
 canvas.Divide(3,2);
 

 
 canvas.cd(1);
 T->Draw("RGEM.rgems.x1.coord.3Dresid>>gemResidX1(1000,-0.01,0.01)");
 gemResidX1->Fit("gaus","","",-0.002,0.002);


 canvas.cd(2);
 T->Draw("RGEM.rgems.x2.coord.3Dresid>>gemResidX2(1000,-0.01,0.01)");
 gemResidX2->Fit("gaus","","",-0.002,0.002);

 canvas.cd(3);
 T->Draw("RGEM.rgems.x3.coord.3Dresid>>gemResidX3(1000,-0.01,0.01)");
 gemResidX3->Fit("gaus","","",-0.002,0.002);

 canvas.cd(4);
 T->Draw("RGEM.rgems.x4.coord.3Dresid>>gemResidX4(1000,-0.01,0.01)");
 gemResidX4->Fit("gaus","","",-0.001,0.001);

 canvas.cd(5);
 T->Draw("RGEM.rgems.x5.coord.3Dresid>>gemResidX5(1000,-0.01,0.01)");
 gemResidX5->Fit("gaus","","",-0.002,0.002);

 canvas.cd(6);
 T->Draw("RGEM.rgems.x6.coord.3Dresid>>gemResidX6(1000,-0.01,0.01)");
 gemResidX6->Fit("gaus","","",-0.002,0.002);
 
/* 
 // initiallize all the Histogram
 for(int i = 0 ; i < 6 ; i ++){
 canvas.cd(i+1);
 gemResidX[i]=new TH1F(Form("Resid_x%d",i+1),Form("Resid_x%d",i+1),1000,-0.02,0.02);
 gemResidY[i]=new TH1F(Form("Resid_y%d",i+1),Form("Resid_y%d",i+1),1000,-0.02,0.02);
 
 T->Draw(Form("RGEM.rgems.x%d.coord.resid >> gemResidX[%d]",i+1,i));
 //T->Draw(Form("RGEM.rgems.y%d.coord.resid>>gemResidY[%d]",i+1,i));
 gemResidX[i]->GetXaxis()->SetRangeUser(-0.02,0.02);
 }
 */
 canvas.Update();
 
 
 TCanvas canvasy;
 canvasy.Divide(3,2);
 

 
 canvasy.cd(1);
 T->Draw("RGEM.rgems.y1.coord.3Dresid>>gemResidY1(1000,-0.01,0.01)");
 gemResidY1->Fit("gaus","","",-0.002,0.002);


 canvasy.cd(2);
 T->Draw("RGEM.rgems.y2.coord.3Dresid>>gemResidY2(1000,-0.01,0.01)");
 gemResidY2->Fit("gaus","","",-0.002,0.002);

 canvasy.cd(3);
 T->Draw("RGEM.rgems.y3.coord.3Dresid>>gemResidY3(1000,-0.01,0.01)");
 gemResidY3->Fit("gaus","","",-0.002,0.002);

 canvasy.cd(4);
 T->Draw("RGEM.rgems.y4.coord.3Dresid>>gemResidY4(1000,-0.01,0.01)");
 gemResidY4->Fit("gaus","","",-0.002,0.002);

 canvasy.cd(5);
 T->Draw("RGEM.rgems.y5.coord.3Dresid>>gemResidY5(1000,-0.01,0.01)");
 gemResidY5->Fit("gaus","","",-0.002,0.002);

 canvasy.cd(6);
 T->Draw("RGEM.rgems.y6.coord.3Dresid>>gemResidY6(1000,-0.01,0.01)");
 //T->Draw("RGEM.rgems.y6.coord.resid>>gemResidY6(1000,-0.01,0.01)","","Y+");
 gemResidY6->Fit("gaus","","",-0.002,0.002);
 
 canvasy.Update();
 
};
