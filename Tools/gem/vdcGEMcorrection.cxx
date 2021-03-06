{
 TPaveText *pt;
 TCanvas canvas;
 canvas.Divide(3,2);
 

 std::cout<<"----------------- X Dimension--------------"<<std::endl;
 std::cout<<"====> x1"<<std::endl;
 canvas.cd(1);
 T->Draw("R.tr.x + R.tr.th * 1.161 - RGEM.rgems.x1.coord.pos >>gemResidX1(1000,-0.05,0.05)");
 gemResidX1->Fit("gaus","","",-0.03,0.05);
std::cout<<std::endl;

 std::cout<<"====> x2"<<std::endl;
 canvas.cd(2);
 T->Draw("R.tr.x + R.tr.th * 1.7979800 - RGEM.rgems.x2.coord.pos>>gemResidX2(1000,-0.05,0.05)");
 gemResidX2->Fit("gaus","","",-0.03,0.015);
 gemResidX2->GetFunction("gaus")->GetParameter(2);
 //pt = new TPaveText(0.2,0.7,0.4,0.85,"NDC");
 //pt->AddText(Form("Mean:%f,\n RMS:%f",gemResidX2->GetFunction("gaus")->GetParameter(2),gemResidX2->GetFunction("gaus")->GetParameter(2)));
 //pt->Draw();

 std::cout<<std::endl;

 std::cout<<"====> x3"<<std::endl;
 canvas.cd(3);
 T->Draw("R.tr.x + R.tr.th * 2.0902131 - RGEM.rgems.x3.coord.pos>>gemResidX3(1000,-0.05,0.05)");
 gemResidX3->Fit("gaus","","",-0.03,0.015);
 std::cout<<std::endl;
 
 std::cout<<"====> x4"<<std::endl;
 canvas.cd(4);
 T->Draw("R.tr.x + R.tr.th * 2.7165651 - RGEM.rgems.x4.coord.pos>>gemResidX4(1000,-0.05,0.05)");
 gemResidX4->Fit("gaus","","",-0.03,0.02);
 
 pt = new TPaveText(0.2,0.7,0.4,0.85,"NDC");
//pt->AddText(Form("Fit Sigma: %5.2f",((TF1 *)(gemResidX4-->GetListOfFunctions()->FindObject("gaus")))->GetParameter(2)));
 pt->SetLineColor(0);
 pt->SetFillColor(0);
 pt->SetShadowColor(0);
 pt->SetTextSize(0.04);
// pt->Draw();
 std::cout<<std::endl;

 std::cout<<"====> x5"<<std::endl;
 canvas.cd(5);
 T->Draw("R.tr.x + R.tr.th * 2.8749137 - RGEM.rgems.x5.coord.pos>>gemResidX5(1000,-0.05,0.05)");
 gemResidX5->Fit("gaus","","",-0.03,0.02);
 std::cout<<std::endl;
 
 std::cout<<"====> x6"<<std::endl;
 canvas.cd(6);
 T->Draw("R.tr.x + R.tr.th * 2.9976041 - RGEM.rgems.x6.coord.pos>>gemResidX6(1000,-0.05,0.05)");
 gemResidX6->Fit("gaus","","",-0.03,0.02);
 std::cout<<std::endl;
 std::cout<<std::endl;
 
 std::cout<<"----------------- Y Dimension--------------"<<std::endl;
 TCanvas canvas1;
 canvas1.Divide(3,2);
 
 std::cout<<"====> 71"<<std::endl;
 canvas1.cd(1);
 T->Draw("R.tr.y + R.tr.ph * 1.161 - RGEM.rgems.y1.coord.pos >>gemResidY1(1000,-0.05,0.05)");
 gemResidY1->Fit("gaus","","",-0.02,0.01);
 std::cout<<std::endl;
 
 std::cout<<"====> y2"<<std::endl;
 canvas1.cd(2);
 T->Draw("R.tr.y + R.tr.ph * 1.7979800 - RGEM.rgems.y2.coord.pos>>gemResidY2(1000,-0.05,0.05)");
 gemResidY2->Fit("gaus","","",-0.02,0.02);
 std::cout<<std::endl;
 
 std::cout<<"====> y3"<<std::endl;
 canvas1.cd(3);
 T->Draw("R.tr.y + R.tr.ph * 2.0902131 - RGEM.rgems.y3.coord.pos>>gemResidY3(1000,-0.05,0.05)");
 gemResidY3->Fit("gaus","","",-0.02,0.02);
 std::cout<<std::endl;
 
 std::cout<<"====> y4"<<std::endl;
 canvas1.cd(4);
 T->Draw("R.tr.y + R.tr.ph * 2.7165651 - RGEM.rgems.y4.coord.pos>>gemResidY4(1000,-0.05,0.05)");
 gemResidY4->Fit("gaus","","",-0.02,0.02);
 std::cout<<std::endl;
 
 std::cout<<"====> y5"<<std::endl;
 canvas1.cd(5);
 T->Draw("R.tr.y + R.tr.ph * 2.8749137 - RGEM.rgems.y5.coord.pos>>gemResidY5(1000,-0.05,0.05)");
 gemResidY5->Fit("gaus","","",-0.02,0.02);
 std::cout<<std::endl;
 
 std::cout<<"====> y6"<<std::endl;
 canvas1.cd(6);
 T->Draw("R.tr.y + R.tr.ph * 2.9976041 - RGEM.rgems.y6.coord.pos>>gemResidY6(1000,-0.05,0.05)");
 gemResidY6->Fit("gaus","","",-0.02,0.02);
 std::cout<<std::endl;
 
 
 
};
