#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <Riostream.h>
#include "auxiliaryPt.h"

#endif
using namespace std;

void test()// where the output figures will be
{
   //samples:
    int nEntry=0;
    double x[20];
    string inputFile =Form("dataSource/corryield_pt_Bs.txt");
    ifstream in;
    in.open(inputFile.c_str());

    string tmp;
    getline(in,tmp);//ignore first line/ the header
    //  while(in.good()){
    //      in >> x[0] >> x[1] >> x[2] >> x[3] >> x[4] >> x[5] >> x[6] >> x[7];
    while(in >> x[0] >> x[1] >> x[2] >> x[3] >> x[4] >> x[5] >> x[6] >> x[7])
      {
        if(nEntry==0){// low-pt bin, in different rapidity window
	  binLow[nEntry] = x[0];
	  bs_low[nEntry] = x[2]; //central value
          cout<< "GATE low: "<<bs_low[0]<<"\t x2= "<<x[2]<<endl;
        }

	if(nEntry>0){
	  binHigh[nEntry-1] = x[0];
	  bs_high[nEntry-1] = x[2];
	  cout<< "GATE high: "<<bs_high[nEntry-1]<<"\t x2= "<<x[2]<<endl;
	}

		
        nEntry++;
    }//while loop
  in.close();//close input file

   
  cout<< "GATE 1: "<<binLow[0]<<"\t central= "<<bs_low[0]<<endl;
  for(int ie=0;ie<3;ie++)
    cout<< "GATE 2: "<<binHigh[ie]<<"\t central= "<<bs_high[ie]<<endl;
  
}
