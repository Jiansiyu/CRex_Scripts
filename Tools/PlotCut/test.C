#include "TCanvas.h"
#include "TImage.h"
void test(){
    TCanvas *canv=new TCanvas("a","a",1960,1080);
    canv->Draw();
    TImage *img=TImage::Open("https://raw.githubusercontent.com/Jiansiyu/GeneralScripts/master/halog/result/BeamE1696.jpg");
    canv->cd();
    img->Draw();
    canv->Update();
}
