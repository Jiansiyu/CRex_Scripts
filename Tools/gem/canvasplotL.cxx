{
 TCanvas canvas;
 canvas.Divide(3,2);
 

 
 canvas.cd(1);
 T->Draw("LGEM.lgems.x1.coord.resid>>gemResidX1(1000,-0.01,0.01)");
 gemResidX1->Fit("gaus","","",-0.004,0.004);


 canvas.cd(2);
 T->Draw("LGEM.lgems.x2.coord.resid>>gemResidX2(1000,-0.01,0.01)");
 gemResidX2->Fit("gaus","","",-0.004,0.004);

 canvas.cd(3);
 T->Draw("LGEM.lgems.x3.coord.resid>>gemResidX3(1000,-0.01,0.01)");
 gemResidX3->Fit("gaus","","",-0.004,0.004);

 canvas.cd(4);
 T->Draw("LGEM.lgems.x4.coord.resid>>gemResidX4(1000,-0.01,0.01)");
 gemResidX4->Fit("gaus","","",-0.002,0.002);

 canvas.cd(5);
 T->Draw("LGEM.lgems.x5.coord.resid>>gemResidX5(1000,-0.01,0.01)");
 gemResidX5->Fit("gaus","","",-0.002,0.002);

 canvas.cd(6);
 T->Draw("LGEM.lgems.x6.coord.resid>>gemResidX6(1000,-0.01,0.01)");
 gemResidX6->Fit("gaus","","",-0.002,0.002);
 
/* 
 // initiallize all the Histogram
 for(int i = 0 ; i < 6 ; i ++){
 canvas.cd(i+1);
 gemResidX[i]=new TH1F(Form("Resid_x%d",i+1),Form("Resid_x%d",i+1),1000,-0.02,0.02);
 gemResidY[i]=new TH1F(Form("Resid_y%d",i+1),Form("Resid_y%d",i+1),1000,-0.02,0.02);
 
 T->Draw(Form("LGEM.lgems.x%d.coord.resid >> gemResidX[%d]",i+1,i));
 //T->Draw(Form("LGEM.lgems.y%d.coord.resid>>gemResidY[%d]",i+1,i));
 gemResidX[i]->GetXaxis()->SetRangeUser(-0.02,0.02);
 }
 */
 canvas.Update();
 
 
 TCanvas canvasy;
 canvasy.Divide(3,2);
 

 
 canvasy.cd(1);
 T->Draw("LGEM.lgems.y1.coord.resid>>gemResidY1(1000,-0.01,0.01)");
 gemResidY1->Fit("gaus","","",-0.004,0.004);


 canvasy.cd(2);
 T->Draw("LGEM.lgems.y2.coord.resid>>gemResidY2(1000,-0.01,0.01)");
 gemResidY2->Fit("gaus","","",-0.004,0.004);

 canvasy.cd(3);
 T->Draw("LGEM.lgems.y3.coord.resid>>gemResidY3(1000,-0.01,0.01)");
 gemResidY3->Fit("gaus","","",-0.004,0.004);

 canvasy.cd(4);
 T->Draw("LGEM.lgems.y4.coord.resid>>gemResidY4(1000,-0.01,0.01)");
 gemResidY4->Fit("gaus","","",-0.002,0.002);

 canvasy.cd(5);
 T->Draw("LGEM.lgems.y5.coord.resid>>gemResidY5(1000,-0.01,0.01)");
 gemResidY5->Fit("gaus","","",-0.002,0.002);

 canvasy.cd(6);
 T->Draw("LGEM.lgems.y6.coord.resid>>gemResidY6(1000,-0.01,0.01)");
 T->Draw("LGEM.lgems.y6.coord.resid>>gemResidY6(1000,-0.01,0.01)","","Y+");
 gemResidY6->Fit("gaus","","",-0.002,0.002);
 
 canvasy.Update();
 
 TCanvas canvasz;
 canvasz.Divide(3,2);
 canvasz.cd(1);
 
 T->Draw("LGEM.lgems.x1.coord.resid:LGEM.tr.th>>h1(1000,-0.03,0.03,1000,-0.03,0.03)");
 h1->Fit("pol1","","",-0.02,0.02);
 TLine a1(-0.03,0.0,0.03,0.0);
 a1.SetLineWidth(2);
 a1.SetLineColor(3);
 a1.Draw("same");
 canvasz.cd(2);
 T->Draw("LGEM.lgems.x2.coord.resid:LGEM.tr.th>>h2(1000,-0.03,0.03,1000,-0.03,0.03)");
 h2->Fit("pol1","","",-0.02,0.02);
 a1.Draw("same");

 canvasz.cd(3);
 T->Draw("LGEM.lgems.x3.coord.resid:LGEM.tr.th>>h3(1000,-0.03,0.03,1000,-0.03,0.03)");
 h3->Fit("pol1","","",-0.02,0.02);
 a1.Draw("same");

 canvasz.cd(4);
 T->Draw("LGEM.lgems.x4.coord.resid:LGEM.tr.th>>h4(1000,-0.03,0.03,1000,-0.03,0.03)");
 h4->Fit("pol1","","",-0.02,0.02);
 a1.Draw("same");

 canvasz.cd(5);
 T->Draw("LGEM.lgems.x5.coord.resid:LGEM.tr.th>>h5(1000,-0.03,0.03,1000,-0.03,0.03)");
 h5->Fit("pol1","","",-0.02,0.02);
 a1.Draw("same");

 canvasz.cd(6);
 T->Draw("LGEM.lgems.x6.coord.resid:LGEM.tr.th>>h6(1000,-0.03,0.03,1000,-0.03,0.03)");
 h6->Fit("pol1","","",-0.02,0.02);
 a1.Draw("same");
 
 TCanvas canvaresidY;
 canvaresidY.Divide(3,2);
 canvaresidY.cd(1);
 T->Draw("LGEM.lgems.y1.coord.resid:LGEM.tr.ph>>hr1(1000,-0.03,0.03,1000,-0.03,0.03)");
 hr1->Fit("pol1","","",-0.02,0.02);
 a1.SetLineWidth(2);
 a1.SetLineColor(3);
 a1.Draw("same");
 
 canvaresidY.cd(2);
 T->Draw("LGEM.lgems.y2.coord.resid:LGEM.tr.ph>>hr2(1000,-0.03,0.03,1000,-0.03,0.03)");
 hr2->Fit("pol1","","",-0.02,0.02);
 a1.SetLineWidth(2);
 a1.SetLineColor(3);
 a1.Draw("same");

 canvaresidY.cd(3);
 T->Draw("LGEM.lgems.y3.coord.resid:LGEM.tr.ph>>hr3(1000,-0.03,0.03,1000,-0.03,0.03)");
 hr3->Fit("pol1","","",-0.02,0.02);
 a1.SetLineWidth(2);
 a1.SetLineColor(3);
 a1.Draw("same");

 canvaresidY.cd(4);
 T->Draw("LGEM.lgems.y4.coord.resid:LGEM.tr.ph>>hr4(1000,-0.03,0.03,1000,-0.03,0.03)");
 hr4->Fit("pol1","","",-0.02,0.02);
 a1.SetLineWidth(2);
 a1.SetLineColor(3);
 a1.Draw("same");

 canvaresidY.cd(5);
 T->Draw("LGEM.lgems.y5.coord.resid:LGEM.tr.ph>>hr5(1000,-0.03,0.03,1000,-0.03,0.03)");
 hr5->Fit("pol1","","",-0.02,0.02);
 a1.SetLineWidth(2);
 a1.SetLineColor(3);
 a1.Draw("same");

 canvaresidY.cd(6);
 T->Draw("LGEM.lgems.y6.coord.resid:LGEM.tr.ph>>hr6(1000,-0.03,0.03,1000,-0.03,0.03)");
 hr6->Fit("pol1","","",-0.02,0.02);
 a1.SetLineWidth(2);
 a1.SetLineColor(3);
 a1.Draw("same");


 
 TCanvas canvasxy;
 canvasxy.Divide(3,2);
 canvasxy.cd(1);
  T->Draw("LGEM.lgems.y1.coord.trkpos:LGEM.lgems.x1.coord.trkpos>>hh1(1000,-0.3,0.3,1000,-0.3,0.3)","","zcol");
  canvasxy.cd(2);
  T->Draw("LGEM.lgems.y2.coord.trkpos:LGEM.lgems.x2.coord.trkpos>>hh2(1000,-0.3,0.3,1000,-0.3,0.3)","","zcol");
  canvasxy.cd(3);
  T->Draw("LGEM.lgems.y3.coord.trkpos:LGEM.lgems.x3.coord.trkpos>>hh3(1000,-0.3,0.3,1000,-0.3,0.3)","","zcol");
  canvasxy.cd(4);
  T->Draw("LGEM.lgems.y4.coord.trkpos:LGEM.lgems.x4.coord.trkpos>>hh4(1000,-0.3,0.3,1000,-0.3,0.3)","","zcol");
  canvasxy.cd(5);
  T->Draw("LGEM.lgems.y5.coord.trkpos:LGEM.lgems.x5.coord.trkpos>>hh5(1000,-0.3,0.3,1000,-0.3,0.3)","","zcol");
  canvasxy.cd(6);
  T->Draw("LGEM.lgems.y6.coord.trkpos:LGEM.lgems.x6.coord.trkpos>>hh6(1000,-0.3,0.3,1000,-0.3,0.3)","","zcol");
 
 std::cout<<"-----------------Start Fit resid -------"<<std::endl;
 TCanvas canvas5;
 canvas5.Divide(3,2);
 
 canvas5.cd(1);
 T->Draw("LGEM.lgems.y2.hit.pos-LGEM.lgems.y1.hit.pos>>hhh1(1000,1,1)");
 std::cout<<">>>> y2"<<std::endl;
 hhh1->Fit("gaus","","",-0.01,0.01);


 canvas5.cd(2);
 T->Draw("LGEM.lgems.y3.hit.pos-LGEM.lgems.y1.hit.pos >>hhh2(1000,1,1)");
 std::cout<<">>>> y3"<<std::endl;
 hhh2->Fit("gaus","","",-0.015,0.015);
 
 canvas5.cd(3);
 T->Draw("LGEM.lgems.y4.hit.pos-LGEM.lgems.y1.hit.pos>>hhh3(1000,1,1)"); 
 std::cout<<">>>> y4"<<std::endl;
 hhh3->Fit("gaus","","",-0.02,0.02);
 
 canvas5.cd(4);
 T->Draw("LGEM.lgems.y5.hit.pos-LGEM.lgems.y1.hit.pos>>hhh4(1000,1,1)");
 std::cout<<">>>> y5"<<std::endl;
 hhh4->Fit("gaus","","",-0.02,0.02);
 
 
 canvas5.cd(5);
 T->Draw("LGEM.lgems.y6.hit.pos-LGEM.lgems.y1.hit.pos>>hhh5(1000,1,1)");
 std::cout<<">>>> y6"<<std::endl;
 hhh5->Fit("gaus","","",-0.02,0.02);
 
 canvas5.Draw();
 TCanvas canvas6;
 canvas6.Divide(3,2);
 canvas6.cd(1);
 std::cout<<"-----------------Start Fit resid XX-------"<<std::endl;
 T->Draw("LGEM.lgems.x2.hit.pos-LGEM.lgems.x1.hit.pos>>hhh6(1000,1,1)");
 std::cout<<">>>> x2"<<std::endl;
 hhh6->Fit("gaus","","",-0.015,0.015);

 canvas6.cd(2);
 T->Draw("LGEM.lgems.x3.hit.pos-LGEM.lgems.x1.hit.pos>>hhh7(1000,1,1)");
 std::cout<<">>>> x3"<<std::endl;
 hhh7->Fit("gaus","","",-0.02,0.02);
 
 canvas6.cd(3);
 T->Draw("LGEM.lgems.x4.hit.pos-LGEM.lgems.x1.hit.pos>>hhh8(1000,1,1)");
 std::cout<<">>>> x4"<<std::endl;
 hhh8->Fit("gaus","","",-0.02,0.02);
 
 canvas6.cd(4);
 T->Draw("LGEM.lgems.x5.hit.pos-LGEM.lgems.x1.hit.pos>>hhh9(1000,1,1)");
  std::cout<<">>>> x5"<<std::endl;
 hhh9->Fit("gaus","","",-0.025,0.025);
 
 canvas6.cd(5);
 T->Draw("LGEM.lgems.x6.hit.pos-LGEM.lgems.x1.hit.pos>>hhh10(1000,1,1)");
 std::cout<<">>>> x6"<<std::endl;
 hhh10->Fit("gaus","","",-0.025,0.025);
 
 canvas6.Draw();
 
};
