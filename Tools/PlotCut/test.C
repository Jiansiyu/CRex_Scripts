#include "TCanvas.h"
#include "TImage.h"
void test(){
    TCanvas *canv=new TCanvas("a","a",1960,1080);
    canv->Draw();
    TImage *img=TImage::Open("/home/newdriver/Learning/GeneralScripts/halog/result/BeamE3105.jpg");
    canv->cd();
    img->Draw();
    canv->Update();
}
