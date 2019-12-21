//R__LOAD_LIBRARY(../ParityData/libParity.so)
//using namespace std;

#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

Bool_t IsFileExist(const Char_t * fname);

void setup(Int_t runNo=0, Int_t lastevt=-1){

	char infile[300];





}

Bool_t IsFileExist(const Char_t * fname)
{
  fstream testfile;

  testfile.open(fname, ios_base::in);
  Bool_t isopen = testfile.is_open();
  testfile.close();

  return isopen;
}
