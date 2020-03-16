/*
Macro to plot xsec and ratio vs <Npart> for Bs and Bp)

Input: txt files in inputDir, with 7 columns: Npart central val, stat, systUp, systDown, glbUp,glbDown

Output: xsec vs pt, ratio vs pt.

 */
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <Riostream.h>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"

#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TMath.h"

#include "TPaveStats.h"
#include "TROOT.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"

#include "CMS_lumi.C"
#include "tdrstyle.C"

#include "auxiliaryCent.h"
#include "auxiliaryRef.h"

#endif
using namespace std;

void plotCent(bool bSavePlots     = 1,
            bool bDoDebug         = 1, //  figure out if things are read properly
            bool whichPlot        = 1, //0 is x-sec, 1 is for ratio
	      bool drawRef        = 1, //draw Ref (for ratio only
            const char* inputDir  = "dataSource", //inptu txt files
            const char* outputDir = "figs")// where the output figures will be
{
  gSystem->mkdir(Form("./%s/png",outputDir), kTRUE);
  gSystem->mkdir(Form("./%s/pdf",outputDir), kTRUE);

  //set the style
  setTDRStyle();
 
  //samples:
  const unsigned int nMes      = 2;
  const char* inputFileType_cent[3] = {"corryield_cent_30_90", "corryield_cent_0_30",                                                     "corryield_cent_0_90"};
  const char* inputFileType_ratio[3] = {"ratio_cent_30_90","ratio_cent_0_30",                                                              "ratio_cent_0_90"};
  
  const char* mesonName[nMes]  = {"Bs", "Bp"};
  Int_t endMes = nMes;
  
  for (Int_t ib=0; ib<endMes; ib++){
    int nEntry=0;
    for(Int_t icent=0; icent<3; icent++){
      ifstream in;
      string inputFileName = Form("%s/%s_%s.txt",inputDir,inputFileType_cent[icent],mesonName[ib]);
      if(whichPlot==1) inputFileName = Form("%s/%s.txt",inputDir,inputFileType_ratio[icent]);
      cout << "########## Input file name: " << inputFileName << endl;
      
      in.open(inputFileName.c_str());
      if (!in.is_open()) {
	cout << "input file " << inputFileName << " cannot be open" << endl;
	continue;
      }
      
      //get first line, the header, and discard it
      string tmpstrg;
      getline(in,tmpstrg);//ignore first line/ the header
      double x[20]={0};

      
      while(in >> x[0] >> x[1] >> x[2] >> x[3] >> x[4] >> x[5] >> x[6]){
	glbSystDown = x[5]*100;
	glbSystUp   = x[6]*100;
	cout<<"x= "<<x[0]<<endl;
	double tmpStat = x[2]*x[1];
	double tmpSysL = x[3]*x[1];
	double tmpSysH = x[4]*x[1];
	if(nEntry<2){// centrality bins
	  binLow[nEntry] = x[0];
	  if(ib==0){//bs
	    bs_low[nEntry] = x[1]; //central value
	    //stat uncert
	   
	    bs_low_yStatL[nEntry] = tmpStat;
	    bs_low_yStatH[nEntry] = bs_low_yStatL[nEntry];
	    //bin width
	    bs_low_xErrL[nEntry] = boxWidth;
	    bs_low_xErrH[nEntry] = bs_low_xErrL[nEntry];
	    //systm. uncert
      	    bs_low_ySystL[nEntry] = tmpSysL;
	    bs_low_ySystH[nEntry] = tmpSysH;

	    cout<<"Gate 1a: "<<nEntry<<"\t"<<bs_low_yStatL[nEntry]<<"\t"<<bs_low_yStatH[nEntry]<<endl;
	    cout<<"Gate 1b: "<<nEntry<<"\t"<<bs_low_ySystL[nEntry]<<"\t"<<bs_low_ySystH[nEntry]<<endl;
	  }
	  else{//bp
	    //central value
	    bpl_low[nEntry] = x[1];
	    //stat uncert
	    bpl_low_yStatL[nEntry] = tmpStat;
	    bpl_low_yStatH[nEntry] = bpl_low_yStatL[nEntry];
	    //bin width
	    bpl_low_xErrL[nEntry] = boxWidth;
	    bpl_low_xErrH[nEntry] = bpl_low_xErrL[nEntry];
	    //systm. uncert
	    bpl_low_ySystL[nEntry] = tmpSysL;
	    bpl_low_ySystH[nEntry] = tmpSysH;
            
	    cout<<"Gate 2a: "<<nEntry<<"\t"<<bpl_low_yStatL[nEntry]<<"\t"<<bpl_low_yStatH[nEntry]<<endl;
	    cout<<"Gate 2b: "<<nEntry<<"\t"<<bpl_low_ySystL[nEntry]<<"\t"<<bpl_low_ySystH[nEntry]<<endl;
	  }
	}
	if(nEntry==2){// minbias bin
	  binHigh[0] = x[0];
	  if(ib==0){//bs
	    //central value
	    bs_high[nEntry-2] = x[1];
	    //stat uncert
	    bs_high_yStatL[nEntry-2] = x[2]*x[1];
	    bs_high_yStatH[nEntry-2] = bs_high_yStatL[nEntry-2];
	    //bin width
	    bs_high_xErrL[nEntry-2] = boxWidth;
	    bs_high_xErrH[nEntry-2] = bs_high_xErrL[nEntry-2];
	    //systm. uncert
	    bs_high_ySystL[nEntry-2] = x[3]*x[1];
	    bs_high_ySystH[nEntry-2] = x[4]*x[1];

	    cout<<"Gate 3a: "<<nEntry<<"\t"<<bs_high_yStatL[nEntry-2]<<"\t"<<bs_high_yStatH[nEntry-2]<<endl;
	    cout<<"Gate 3b: "<<nEntry<<"\t"<<bs_high_ySystL[nEntry-2]<<"\t"<<bs_high_ySystH[nEntry-2]<<endl;
	    
	  }else{//bp
	    //central value
	    bpl_high[nEntry-2] = x[1];
	    //stat uncert
	    bpl_high_yStatL[nEntry-2] = x[2]*x[1];
	    bpl_high_yStatH[nEntry-2] = bpl_high_yStatL[nEntry-2];
	    //bin width
	    bpl_high_xErrL[nEntry-2] = boxWidth;
	    bpl_high_xErrH[nEntry-2] = bpl_high_xErrL[nEntry-2];
	    //systm. uncert
	    bpl_high_ySystL[nEntry-2] = x[3]*x[1];
	    bpl_high_ySystH[nEntry-2] = x[4]*x[1];

	    cout<<"Gate 4a: "<<nEntry<<"\t"<<bpl_high_yStatL[nEntry-2]<<"\t"<<bpl_high_yStatH[nEntry-2]<<endl;
	    cout<<"Gate 4b: "<<nEntry<<"\t"<<bpl_high_ySystL[nEntry-2]<<"\t"<<bpl_high_ySystH[nEntry-2]<<endl;
	  }
	}//high-pt bins
	nEntry++;
      }//reading input file line by line
      in.close();//close input file
    }// for icent
    if(bDoDebug){
   
      for(int i=0; i<endMes; i++){
	if(ib==0)
	  cout<<"Element_cent " << i << "\t binLow[i] = "<< binLow[i] << "\t statUncertL = "<<bs_low_yStatL[i]<< "\t statUncertH = "<<bs_low_yStatH[i]<<"\t systUncertL = "<<bs_low_ySystL[i]<<"\t systUncertH = "<<bs_low_ySystH[i]<< endl;
		
	if(ib==1)
	  cout<<"Element_cent " << i << "\t binLow[i] = "<< binLow[i] << "\t statUncertL = "<<bpl_low_yStatL[i]<< "\t statUncertH = "<<bpl_low_yStatH[i]<<"\t systUncertL = "<<bpl_low_ySystL[i]<<"\t systUncertH = "<<bpl_low_ySystH[i]<< endl;
	  }
	
      
    }//bDebug
    cout<<"@@@@@@@@@@@ Finished meson "<<mesonName[ib]<<endl;
  }//for each meson,ib
  
  //----------------------------------------------------------------
  // gr = new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);
    // Bs
    TGraphAsymmErrors *pgBs_low = new TGraphAsymmErrors(nBinsLow, binLow, bs_low,
                                                        bs_low_xErrL, bs_low_xErrH,
                                                        bs_low_yStatL,bs_low_yStatH);
    TGraphAsymmErrors *pgBs_high= new TGraphAsymmErrors(nBinsHigh,binHigh,bs_high,
                                                        bs_high_xErrL, bs_high_xErrH,
                                                        bs_high_yStatL,bs_high_yStatH);
    // Bplus
    TGraphAsymmErrors *pgBpl_low = new TGraphAsymmErrors(nBinsLow, binLow, bpl_low,
                                                        bpl_low_xErrL, bpl_low_xErrH,
                                                        bpl_low_yStatL,bpl_low_yStatH);
    TGraphAsymmErrors *pgBpl_high= new TGraphAsymmErrors(nBinsHigh,binHigh,bpl_high,
                                                         bpl_high_xErrL, bpl_high_xErrH,
                                                         bpl_high_yStatL,bpl_high_yStatH);
      

    //Systmeatic uncertainty
    // Bs
    TGraphAsymmErrors *pgBs_syst_low    = new TGraphAsymmErrors(nBinsLow, binLow, bs_low,
                                                                bs_low_xErrL, bs_low_xErrH,
                                                                bs_low_ySystL,bs_low_ySystH);
    TGraphAsymmErrors *pgBs_syst_high   = new TGraphAsymmErrors(nBinsHigh,binHigh,bs_high,
                                                                bs_high_xErrL, bs_high_xErrH,
                                                                bs_high_ySystL,bs_high_ySystH);

    // Bp
    TGraphAsymmErrors *pgBpl_syst_low    = new TGraphAsymmErrors(nBinsLow, binLow,  bpl_low,
                                                                 bpl_low_xErrL, bpl_low_xErrH,
                                                                 bpl_low_ySystL, bpl_low_ySystH);
    TGraphAsymmErrors *pgBpl_syst_high   = new TGraphAsymmErrors(nBinsHigh,binHigh, bpl_high,
                                                                 bpl_high_xErrL,bpl_high_xErrH,
                                                                 bpl_high_ySystL,bpl_high_ySystH);
   
    //==========================================
    //------------------------------------------
    TGraphAsymmErrors *pgRatio_low    = new TGraphAsymmErrors(nBinsLow, binLow, bs_low,
                                                              bs_low_xErrL, bs_low_xErrH,
                                                              bs_low_yStatL,bs_low_yStatH);
    TGraphAsymmErrors *pgRatio_high   = new TGraphAsymmErrors(nBinsHigh,binHigh,bs_high,
                                                              bs_high_xErrL, bs_high_xErrH,
                                                              bs_high_yStatL,bs_high_yStatH);

    TGraphAsymmErrors *pgRatio_syst_low = new TGraphAsymmErrors(nBinsLow, binLow, bs_low,
                                                                bs_low_xErrL, bs_low_xErrH,
                                                                bs_low_ySystL,bs_low_ySystH);
    TGraphAsymmErrors *pgRatio_syst_high= new TGraphAsymmErrors(nBinsHigh,binHigh,bs_high,
                                                                bs_high_xErrL, bs_high_xErrH,
                                                                bs_high_ySystL,bs_high_ySystH);

    //=========== reference
     //================== reference
    TGraphAsymmErrors * FragBand = new TGraphAsymmErrors(BandBin,BandXCent,BandY,BandXErrCent,BandXErrCent,BandYErr,BandYErr);
    FragBand->SetName("BandErr");
    FragBand->SetMarkerStyle(20);
    FragBand->SetMarkerSize(0.8);
    FragBand->SetFillColorAlpha(kGreen,0.5);
    FragBand->SetFillStyle(3004);
    FragBand->SetLineWidth(2);
    FragBand->SetLineColor(kGreen);

    
  // // **************** marker setup
    // Bs
    // marker style
    pgBs_low->SetMarkerStyle(markerLow[0]);
    pgBs_high->SetMarkerStyle(markerHigh[0]);
    
    pgBpl_low->SetMarkerStyle(markerLow[1]);
    pgBpl_high->SetMarkerStyle(markerHigh[1]);
      
    pgRatio_low->SetMarkerStyle(markerRatio[1]);
    pgRatio_high->SetMarkerStyle(markerRatio[0]);
      
    // marker size
    pgBs_low->SetMarkerSize(markerSizeLow[0]);
    pgBs_high->SetMarkerSize(markerSizeHigh[0]);
      
    pgBpl_low->SetMarkerSize(markerSizeLow[1]);
    pgBpl_high->SetMarkerSize(markerSizeHigh[1]);
        
    pgRatio_low->SetMarkerSize(markerSizeRatio[0]);
    pgRatio_high->SetMarkerSize(markerSizeRatio[1]);

    // marker color
    pgBs_low->SetMarkerColor(colorLow[0]);
    pgBs_high->SetMarkerColor(colorHigh[0]);
      
    pgBpl_low->SetMarkerColor(colorLow[1]);
    pgBpl_high->SetMarkerColor(colorHigh[1]);
        
    pgRatio_low->SetMarkerColor(colorRatio[0]);
    pgRatio_high->SetMarkerColor(colorRatio[1]);
    
    // systematic boxes
    pgBs_syst_low->SetFillColorAlpha(kBlue-9,0.5);
    pgBs_syst_high->SetFillColorAlpha(kBlue-9,0.5);
    pgBpl_syst_low->SetFillColorAlpha(kGreen-9,0.5);
    pgBpl_syst_high->SetFillColorAlpha(kGreen-9,0.5);
 
    pgRatio_syst_low->SetFillColorAlpha(colorRatio[0],0.2);
    pgRatio_syst_high->SetFillColorAlpha(colorRatio[1],0.2);
   
  
    //-------------------------------------------
    TF1 *f4 = new TF1("f4","1",0,420);
    f4->SetLineWidth(1);
    f4->SetLineColor(1);
    f4->SetLineStyle(1);
    f4->GetYaxis()->SetTitle(yAxName[whichPlot]);
    f4->GetXaxis()->SetTitle(xAxName[0]);
    f4->GetXaxis()->CenterTitle(kTRUE);
    f4->GetYaxis()->SetRangeUser(5e2,5e6);
    if(whichPlot==1) f4->GetYaxis()->SetRangeUser(0.0,1.5);
    //f4->GetXaxis()->SetNdivisions(-6);

   //---------------- general stuff
   TLatex *lat = new TLatex();
   lat->SetNDC();

  // // ##################################################### x-sec canvas
    TCanvas *pc1 = new TCanvas("pc1","pc1");
    f4->Draw();// axis
  
    CMS_lumi(pc1,19011,0);
    
    if(whichPlot==0)// x-section
    {
      gPad->SetLogy();
  
      pgBs_syst_low->Draw("2");
      pgBpl_syst_low->Draw("2");
      
      pgBs_low->Draw("P");
      pgBpl_low->Draw("P");
        
//        pgBs_syst_high->Draw("2");
//        pgBs_high->Draw("P");
//
//        pgBpl_syst_high->Draw("2");
//        pgBpl_high->Draw("P");
    }
    else{//ratio plot
          if(drawRef){
            FragBand->Draw("5same");
	    

        }
            pgRatio_syst_low->Draw("2");
            pgRatio_low->Draw("P");
            
	    //  pgRatio_syst_high->Draw("2");
	    // pgRatio_high->Draw("P");
	 
    }
    
//supplemental info on plot:
    if(whichPlot==0){
        lat->SetTextFont(42);
        lat->SetTextSize(ltxSetTextSize2);
     lat->DrawLatex(xsec_ltxText1_xStart,xsec_ltxText1_yStart,"#splitline{|y|<2.4}{p_{T}>10 GeV/c}");

        lat->SetTextFont(42);
        lat->SetTextSize(ltxSetTextSize2);
        lat->DrawLatex(xsec_ltxText1_xStart,xsec_ltxText1_yStart-0.6,Form("Global uncert.:#pm %.1f %%",glbSystDown));

        // legend
        TLegend *legXSec = new TLegend(legXsec_xLowStart,legXsec_y,legXsec_xLowEnd,legXsec_y+0.1,"B_{s}^{0}  B^{+}","brNDC");
        legXSec->SetBorderSize(0);
       
        legXSec->SetTextSize(ltxSetTextSize2);
        legXSec->SetLineColor(1);
        legXSec->SetLineStyle(1);
        legXSec->SetLineWidth(1);
        legXSec->SetFillColor(19);
        legXSec->SetFillStyle(0);
        legXSec->SetTextFont(42);
	// legXSec->SetHeader("","L");
        legXSec->SetNColumns(2);
        legXSec->SetColumnSeparation(0.05);
        legXSec->AddEntry(pgBs_low," ","p");
        legXSec->AddEntry(pgBpl_low," ","p");

	legXSec->Draw();
	
    }else{
        lat->SetTextSize(ltxSetTextSize1);
        lat->SetTextFont(42);
        lat->SetTextFont(132);
        lat->DrawLatex(ratio_ltxText1_xStart,ratio_ltxText1_yStart,"#frac{B_{s}^{0}}{B^{+}}");
        
        lat->SetTextFont(42);
        lat->SetTextSize(ltxSetTextSize2);
        lat->DrawLatex(ratio_ltxText2_xStart,ratio_ltxText2_yStart,"#splitline{|y|<2.4}{p_{T}>10 GeV/c}");

        lat->SetTextFont(42);
        lat->SetTextSize(ltxSetTextSize2);
        lat->DrawLatex(legRatio_xLowStart-0.05,legRatio_y-0.08,Form("Global uncert.:#pm %.1f %%",glbSystDown));

	if(drawRef)
	  {
	    TLegend *legRatioRef = new TLegend(legRatioRef_xLowStart,legRatioRef_y,legRatioRef_xLowEnd,legRatio_y+0.1,"f_{s}/f_{u} reference","brNDC");
	    legRatioRef->SetBorderSize(0);
	    legRatioRef->SetTextFont(132);
	    
	    legRatioRef->SetTextSize(ltxSetTextSize2);
	    legRatioRef->SetLineColor(1);
	    legRatioRef->SetLineStyle(1);
	    legRatioRef->SetLineWidth(1);
	    legRatioRef->SetFillColor(19);
	    legRatioRef->SetFillStyle(0);
	    TLegendEntry *entry1Ref = legRatioRef->AddEntry("FragBand","PDG","P");
	    entry1Ref->SetTextFont(42);
	    entry1Ref->SetFillStyle(1001);
	    entry1Ref->SetMarkerStyle(25);
	    entry1Ref->SetMarkerSize(1.4);
	    entry1Ref->SetMarkerColor(kGreen);
	    entry1Ref->SetLineWidth(5);
	    // TLegendEntry *entry2Ref = legRatioRef->AddEntry("pgRatio_high","|y| < 2.4","p");
	    // entry2Ref->SetTextFont(42);
	    // entry2Ref->SetMarkerStyle(markerRatio[1]);
	    // entry2Ref->SetMarkerSize(1.7);
	    // entry2Ref->SetFillStyle(1001);
    
	    legRatioRef->Draw();

	  }

    }
    
    
  // gPad->RedrawAxis();
   pc1->Update();

   if(bSavePlots)
   {
       if (whichPlot==0)
       {
           pc1->SaveAs(Form("%s/pdf/xsec_vsCent.pdf",outputDir));
           pc1->SaveAs(Form("%s/png/xsec_vsCent.png",outputDir));
       }else{
	 pc1->SaveAs(Form("%s/pdf/ratio_vsCent_ref%d.pdf",outputDir,drawRef));
	 pc1->SaveAs(Form("%s/png/ratio_vsCent_ref%d.png",outputDir,drawRef));
        }
   }

}
