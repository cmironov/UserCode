/*
Macro to plot xsec and ratio vs pt for Bs and Bp)

Input: txt files in inputDir, with 8 columns: ptmin, ptmax, central val, stat, systUp, systDown, glbUp,glbDown

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

#include "auxiliaryPt.h"
#include "auxiliaryRef.h"

#endif
using namespace std;

void plotPt(bool bSavePlots       = 1,
            bool bDoDebug         = 0, //  figure out if things are read properly
            bool whichPlot        = 1, //0 is x-sec, 1 is for ratio
            bool drawRef          = 1, //draw Ref (for ratio only
            const char* inputDir  = "dataSource", //inptu txt files
            const char* outputDir = "figs")// where the output figures will be
{
    gSystem->mkdir(Form("./%s/png",outputDir), kTRUE);
    gSystem->mkdir(Form("./%s/pdf",outputDir), kTRUE);

  //set the style
    setTDRStyle();
 
  //samples:
    const unsigned int nMes      = 2;
    const char* inputFileType[2] = {"corryield_pt", "ratio_pt"};
    const char* mesonName[nMes]  = {"Bs", "Bp"};
    
    Int_t endMes = nMes;
  
    for (Int_t ib=0; ib<endMes; ib++){
        ifstream in;
        string inputFileName = Form("%s/%s_%s.txt",inputDir,inputFileType[whichPlot],mesonName[ib]);
        if(whichPlot==1) inputFileName = Form("%s/%s.txt",inputDir,inputFileType[whichPlot]);
        cout << "########## Input file name: " << inputFileName << endl;
            
        in.open(inputFileName.c_str());
        if (!in.is_open()) {
            cout << "input file " << inputFileName << " cannot be open" << endl;
            continue;
        }
        double x[20];
        
        //get first line, the header, and discard it
        string tmpstrg;
        getline(in,tmpstrg);//ignore first line/ the header
        
        int nEntry=0;
        while(in >> x[0] >> x[1] >> x[2] >> x[3] >> x[4] >> x[5] >> x[6] >> x[7])
	  {
	    glbSystDown = x[6]*100;
            glbSystUp   = x[7]*100;
            
            if(nEntry==0){// low-pt bin, in different rapidity window
	      binLow[nEntry]    = (x[1]-x[0])/2 + x[0];
              if(ib==0){//bs
		bs_low[nEntry] = x[2]; //central value
                //stat uncert
        bs_low_yStatL[nEntry] = x[3]*x[2];
        bs_low_yStatH[nEntry] = bs_low_yStatL[nEntry];
                //bin width
                  bs_low_xErrL[nEntry] = (x[1]-x[0])/2;
        bs_low_xErrH[nEntry] = bs_low_xErrL[nEntry];
                //systm. uncert
		bs_low_ySystL[nEntry] = x[4]*x[2];
		bs_low_ySystH[nEntry] = x[5]*x[2];
                                     
		if(bDoDebug){
		  cout << "Central  \t statUncertL \t systUncerL " << endl;
		  cout<< bs_low[nEntry] << "\t" << bs_low_yStatL[nEntry] << "\t" <<   bs_low_ySystL[nEntry] << endl;
		  cout<<"==========================================="<<endl;
		}
	      }
	      else{//bp
                   //central value
		bpl_low[nEntry] = x[2];
		//stat uncert
		bpl_low_yStatL[nEntry] = x[3]*x[2];
		bpl_low_yStatH[nEntry] = bpl_low_yStatL[nEntry];
		//bin width
        bpl_low_xErrL[nEntry] = (x[1]-x[0])/2;
		bpl_low_xErrH[nEntry] = bpl_low_xErrL[nEntry];
		//systm. uncert
		bpl_low_ySystL[nEntry] = x[4]*x[2];
		bpl_low_ySystH[nEntry] = x[5]*x[2];
		
		if(bDoDebug){
		  cout << "Central  \t statUncertL  \t systUncertL " << endl;
		  cout<< bpl_low[nEntry] << "\t" << bpl_low_yStatL[nEntry] << "\t" <<   bpl_low_ySystL[nEntry] << endl;
		  cout<<"==========================================="<<endl;
		}
	      }
            }
	    if(nEntry>0){// high-pt bins, |y|<2.4
	      binHigh[nEntry-1] = (x[1]-x[0])/2 + x[0];
	      if(ib==0){//bs
		//central value
		bs_high[nEntry-1] = x[2];
		//stat uncert
		bs_high_yStatL[nEntry-1] = x[3]*x[2];
		bs_high_yStatH[nEntry-1] = bs_high_yStatL[nEntry-1];
		//bin width
        bs_high_xErrL[nEntry-1] = (x[1]-x[0])/2;
		bs_high_xErrH[nEntry-1] = bs_high_xErrL[nEntry-1];
		//systm. uncert
		bs_high_ySystL[nEntry-1] = x[4]*x[2];
		bs_high_ySystH[nEntry-1] = x[5]*x[2];
		
		if(bDoDebug){
		  if(nEntry==1) cout << "Central  \t statUncertL  \t systUncertL "<< endl;
		  cout<< bs_high[nEntry-1] << "\t" << bs_high_yStatL[nEntry-1] << "\t" << bs_high_ySystL[nEntry-1] << endl;
		}
	      }else{//bp
		//central value
		bpl_high[nEntry-1] = x[2];
		//stat uncert
		bpl_high_yStatL[nEntry-1] = x[3]*x[2];
		bpl_high_yStatH[nEntry-1] = bpl_high_yStatL[nEntry-1];
		//bin width
              bpl_high_xErrL[nEntry-1] = (x[1]-x[0])/2;
		bpl_high_xErrH[nEntry-1] = bpl_high_xErrL[nEntry-1];
		//systm. uncert
		bpl_high_ySystL[nEntry-1] = x[4]*x[2];
		bpl_high_ySystH[nEntry-1] = x[5]*x[2];
	      }
	    }//high-pt bins
	    nEntry++;
	  }//reading input file line by line
      
        cout<<"@@@@@@@@@@@ Finished meson "<<mesonName[ib]<<endl;
	if(bDoDebug){
	  if(ib==0) cout<<"Element_low = "<< binLow[0] <<"\t bs_low[0]= " << bs_low[0] << "\t statUncertL = "<<bs_low_yStatL[0]<<endl;
	
	  if(ib==1)
	    {
	      for(int i=0; i<3; i++){
		cout<<"Element_high " << i << "\t binHigh[i] = "<< binHigh[i] <<"\t bs_high[i]= " << bs_high[i] << "\t statUncertL = "<<bs_high_yStatL[i]<<"\t systUncertL = "<<bs_high_ySystL[i]<< endl;
	      }
	    }
	}//bDebug
        in.close();//close input file
	
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
    // pgBs_syst_high->Draw("APL");
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

    //================== reference
    TGraphAsymmErrors * FragBand = new TGraphAsymmErrors(BandBin,BandX,BandY,BandXErr,BandXErr,BandYErr,BandYErr);
    FragBand->SetName("BandErr");
    FragBand->SetMarkerStyle(20);
    FragBand->SetMarkerSize(0.8);
    FragBand->SetFillColorAlpha(kGreen,0.5);
    FragBand->SetFillStyle(3004);
    FragBand->SetLineWidth(2);
    FragBand->SetLineColor(kGreen);
    
    //==================
  // // **************** marker setup
    // Bs
    // marker style
    pgBs_low->SetMarkerStyle(markerLow[0]);
    pgBs_high->SetMarkerStyle(markerHigh[0]);
    
    pgBpl_low->SetMarkerStyle(markerLow[1]);
    pgBpl_high->SetMarkerStyle(markerHigh[1]);
      
    pgRatio_low->SetMarkerStyle(markerRatio[0]);
    pgRatio_high->SetMarkerStyle(markerRatio[1]);
      
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
    TF1 *f4 = new TF1("f4","1",5,50);
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
	pgBs_syst_high->Draw("2");
        pgBs_low->Draw("P");
        pgBs_high->Draw("P");
        
        pgBpl_syst_low->Draw("2");
	pgBpl_syst_high->Draw("2");
        pgBpl_low->Draw("P");
        pgBpl_high->Draw("P");
    }
    else{//ratio plot
            pgRatio_syst_low->Draw("2");
            pgRatio_low->Draw("P");
            
            pgRatio_syst_high->Draw("2");
            pgRatio_high->Draw("P");
        
        if(drawRef){
            FragBand->Draw("5same");
        }
    }
    
//supplemental info on plot:
    if(whichPlot==0){
    
        
        lat->SetTextFont(42);
        lat->SetTextSize(ltxSetTextSize2);
        lat->DrawLatex(xsec_ltxText1_xStart,xsec_ltxText1_yStart,"Cent. 0-90%");

	lat->SetTextFont(42);
        lat->SetTextSize(ltxSetTextSize2);
        lat->DrawLatex(xsec_ltxText1_xStart,xsec_ltxText1_yStart-0.7,Form("Global uncert.:#pm %.1f %%",glbSystDown));

	// legend
	TLegend *legXSec = new TLegend(legXsec_xLowStart,legXsec_y,legXsec_xLowEnd,legXsec_y+0.15,"B_{s}^{0}                    B^{+}","brNDC");
        legXSec->SetBorderSize(0);
       
        legXSec->SetTextSize(ltxSetTextSize2);
        legXSec->SetLineColor(1);
        legXSec->SetLineStyle(1);
        legXSec->SetLineWidth(1);
        legXSec->SetFillColor(19);
        legXSec->SetFillStyle(0);
	legXSec->SetTextFont(42);
	legXSec->SetNColumns(2);
	legXSec->SetColumnSeparation(0.0);
	legXSec->AddEntry(pgBs_low,"1.5 < |y| < 2.4","p");
	legXSec->SetTextFont(132);
	legXSec->AddEntry(pgBpl_low," ","p");
	legXSec->AddEntry(pgBs_high,"|y| < 2.4","p");
	legXSec->SetTextFont(132);
	legXSec->AddEntry(pgBpl_high," ","p");

	legXSec->Draw();
	
    }else{
        lat->SetTextSize(ltxSetTextSize1);
        lat->SetTextFont(42);
        lat->SetTextFont(132);
        lat->DrawLatex(ratio_ltxText1_xStart,ratio_ltxText1_yStart,"#frac{B_{s}^{0}}{B^{+}}");
        
        lat->SetTextFont(42);
        lat->SetTextSize(ltxSetTextSize2);
        lat->DrawLatex(ratio_ltxText2_xStart,ratio_ltxText2_yStart,"Cent. 0-90%");
       
        lat->SetTextFont(42);
        lat->SetTextSize(ltxSetTextSize2);
        lat->DrawLatex(legRatio_xLowStart-0.05,legRatio_y-0.1,Form("Global uncert.:#pm %.1f %%",glbSystDown));
        
        // legend
        TLegend *legRatio = new TLegend(legRatio_xLowStart,legRatio_y,legRatio_xLowEnd,legRatio_y+0.1,NULL,"brNDC");
        legRatio->SetBorderSize(0);
        legRatio->SetTextFont(132);
        legRatio->SetTextSize(ltxSetTextSize2);
        legRatio->SetLineColor(1);
        legRatio->SetLineStyle(1);
        legRatio->SetLineWidth(1);
        legRatio->SetFillColor(19);
        legRatio->SetFillStyle(0);
        TLegendEntry *entry1 = legRatio->AddEntry("pgRatio_low","1.5 < |y| < 2.4","p");
        entry1->SetTextFont(42);
        entry1->SetMarkerStyle(markerRatio[0]);
        entry1->SetMarkerSize(1.7);
        entry1->SetFillStyle(1001);
        
        TLegendEntry *entry2 = legRatio->AddEntry("pgRatio_high","|y| < 2.4","p");
        entry2->SetTextFont(42);
        entry2->SetMarkerStyle(markerRatio[1]);
        entry2->SetMarkerSize(1.7);
        entry2->SetFillStyle(1001);
    
        legRatio->Draw();
	//---------------
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
           pc1->SaveAs(Form("%s/pdf/xsec_vsPt.pdf",outputDir));
           pc1->SaveAs(Form("%s/png/xsec_vsPt.png",outputDir));
       }else{
            pc1->SaveAs(Form("%s/pdf/ratio_vsPt_ref%d.pdf",outputDir,drawRef));
            pc1->SaveAs(Form("%s/png/ratio_vsPt_ref%d.png",outputDir,drawRef));
        }
   }

}
