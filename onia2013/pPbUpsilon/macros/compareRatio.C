#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
#include <fstream>

#include <sstream>
#include <string>
#include <vector>

#include <Riostream.h>

#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"

#include "./dats/data_2013.h"
#include "./dats/data_2013_ref.h"
#endif

// // what rapidity bins available
const int nEtaIntervals                   = 2;
const char* dimuYCM_legend[nEtaIntervals] = {"|y_{CM}|< 1","|y_{CM}|< 2.4"};


// data sample

const int nBin4 = 4; 
const int nBin2 = 2; 
const char* pairLegend[nBin2]   = {"TRK-TRK",
				   "GLB-GLB"};
const char* sampleLegend[nBin4] = {"pPb 5TeV",
				   "PbPb 2.76 TeV",
				   "pp 2.76 TeV",
				   "pp 7 TeV"};

std::string recoType[nBin4]   = {"pp",
				 "RegIt",
				 "HI",
				 "HI pub"};
std::string trgType[nBin2]    = {"Dimuon0","DoubleMu3"};
std::string trkType[nBin2]    = {"Glb-Glb","Trk-trk"};

double centBins_compare4[nBin4]        = { 1, 3, 5, 7};
double centBinsEdges_compare4[nBin4+1] = {0, 2, 4, 6, 8};
double centBins_compare2[nBin2]        = { 1, 3};
double centBinsEdges_compare2[nBin2+1] = {0, 2, 4};

void compareRatio(bool bSavePlots   = true, 
		  bool compTrking   = true,
		  bool compGlbTrk   = false,
		  bool compTrg7tev  = false,
		  bool doPM1        = false,
		  const char* figNamePrefix="sanityRatio")
{
  gROOT->Macro("/Users/eusmartass/Documents/cmswrite/hin-10-006/logon.C+");
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(5);
  

  TGraphErrors *pgCent_2s1s_trig7   = new TGraphErrors(2); // Dimuon0, DoubleMu3
  TGraphErrors *pgCent_3s1s_trig7   = new TGraphErrors(2); // 

  TGraphErrors *pgCent_2s1s_trking  = new TGraphErrors(2); // 50-100: trk-trk, glb-glb
  TGraphErrors *pgCent_3s1s_trking  = new TGraphErrors(2); // 

  TGraphErrors *pgCent_2s1s_aa50100_reco = new TGraphErrors(4); // 276 AA:  0: pp 1: regit 2: hi zhen tree 3: hi publish

  // fill the graphs:
  for (int iBin=0; iBin<4; iBin++)
    {
      switch (iBin)
	{
	case 0:
	  pgCent_2s1s_aa50100_reco->SetPoint(iBin,centBins_compare4[iBin], ppreco_gg_ypm24_cent50100_2s1s[0]);
	  pgCent_2s1s_aa50100_reco->SetPointError(iBin,0,ppreco_gg_ypm24_cent50100_2s1s_statError[0]);

	  if(doPM1)
	    {
	      //AA pp reco
	      pgCent_2s1s_aa50100_reco->SetPoint(iBin,centBins_compare4[iBin], ppreco_gg_ypm1_cent50100_2s1s[0]);
	      pgCent_2s1s_aa50100_reco->SetPointError(iBin,0,ppreco_gg_ypm1_cent50100_2s1s_statError[0]);
	      
	      // trigger 7TeV Dimuon0
	      pgCent_2s1s_trig7->SetPoint(iBin,centBins_compare4[iBin], pp7_yPM1_2s1s[0]);
	      pgCent_3s1s_trig7->SetPoint(iBin,centBins_compare4[iBin], pp7_yPM1_3s1s[0]);
	      
	      pgCent_2s1s_trig7->SetPointError(iBin,0, pp7_yPM1_2s1s_statError[0]);
	      pgCent_3s1s_trig7->SetPointError(iBin,0, pp7_yPM1_3s1s_statError[0]);

	      // TT
	      pgCent_2s1s_trking->SetPoint(iBin,centBins_compare4[iBin], ppreco_ypm1_cent50100_2s1s[0]);
	      pgCent_2s1s_trking->SetPointError(iBin,0, ppreco_ypm1_cent50100_2s1s_statError[0]);
	    }

	  break;
	case 1:
	  pgCent_2s1s_aa50100_reco->SetPoint(iBin,centBins_compare4[iBin], regit_gg_ypm24_cent50100_2s1s[0]);
	  pgCent_2s1s_aa50100_reco->SetPointError(iBin,0,regit_gg_ypm24_cent50100_2s1s_statError[0]);

	  if(doPM1)
	    {
	      // regit
	      pgCent_2s1s_aa50100_reco->SetPoint(iBin,centBins_compare4[iBin], regit_gg_ypm1_cent50100_2s1s[0]);
	      pgCent_2s1s_aa50100_reco->SetPointError(iBin,0,regit_gg_ypm1_cent50100_2s1s_statError[0]);

	      // trigger 7TeV: DblMu3
	      pgCent_2s1s_trig7->SetPoint(iBin,centBins_compare4[iBin], dblmu3_pp7_yPM1_2s1s[0]);
	      pgCent_3s1s_trig7->SetPoint(iBin,centBins_compare4[iBin], dblmu3_pp7_yPM1_3s1s[0]);

	      pgCent_2s1s_trig7->SetPointError(iBin,0, dblmu3_pp7_yPM1_2s1s_statError[0]);
	      pgCent_3s1s_trig7->SetPointError(iBin,0, dblmu3_pp7_yPM1_3s1s_statError[0]);

	      // GG
	      pgCent_2s1s_trking->SetPoint(iBin,centBins_compare4[iBin], ppreco_gg_ypm1_cent50100_2s1s[0]);
	      pgCent_2s1s_trking->SetPointError(iBin,0, ppreco_gg_ypm1_cent50100_2s1s_statError[0]);
	    }
	  break;
	case 2:
	  pgCent_2s1s_aa50100_reco->SetPoint(iBin,centBins_compare4[iBin], zhen_gg_ypm24_cent50100_2s1s[0]);
	  pgCent_2s1s_aa50100_reco->SetPointError(iBin,0,zhen_gg_ypm24_cent50100_2s1s_statError[0]);

	  if(doPM1)
	    {
	      pgCent_2s1s_aa50100_reco->SetPoint(iBin,centBins_compare4[iBin], zhen_gg_ypm1_cent50100_2s1s[0]);
	      pgCent_2s1s_aa50100_reco->SetPointError(iBin,0,zhen_gg_ypm1_cent50100_2s1s_statError[0]);
	    }
	  break;
	case 3:
	  pgCent_2s1s_aa50100_reco->SetPoint(iBin,centBins_compare4[iBin], pub_gg_ypm24_cent50100_2s1s[0]);
	  pgCent_2s1s_aa50100_reco->SetPointError(iBin,0,pub_gg_ypm24_cent50100_2s1s_statError[0]);

	  break;

	}
    }

  
  cout << "!!!! Got all the numebrs1"<<endl;

  // drawing 2010 #A4C639
   Color_t myColor[2] = {kBlue+2,kGreen+2};// 2s1s, 3s1s
  // Color_t myColor[5] = {"#00008B","#8DB600"};// 2s1s, 3s1s
  int myMarker[2]    = {33,27}; 

  // 2s/1s
  pgCent_2s1s_trig7->SetMarkerColor(myColor[0]);
  pgCent_2s1s_trig7->SetLineColor(myColor[0]);
  pgCent_2s1s_trig7->SetMarkerStyle(myMarker[0]);
  pgCent_2s1s_trig7->SetMarkerSize(2);

  pgCent_2s1s_trking->SetMarkerColor(myColor[0]);
  pgCent_2s1s_trking->SetLineColor(myColor[0]);
  pgCent_2s1s_trking->SetMarkerStyle(myMarker[0]);
  pgCent_2s1s_trking->SetMarkerSize(2);

  pgCent_2s1s_aa50100_reco->SetMarkerColor(myColor[0]);
  pgCent_2s1s_aa50100_reco->SetLineColor(myColor[0]);
  pgCent_2s1s_aa50100_reco->SetMarkerStyle(myMarker[0]);
  pgCent_2s1s_aa50100_reco->SetMarkerSize(2);

  // 3s/1s
  pgCent_3s1s_trig7->SetMarkerColor(myColor[1]);
  pgCent_3s1s_trig7->SetLineColor(myColor[1]);
  pgCent_3s1s_trig7->SetMarkerStyle(myMarker[1]);
  pgCent_3s1s_trig7->SetMarkerSize(2);

  pgCent_3s1s_trking->SetMarkerColor(myColor[1]);
  pgCent_3s1s_trking->SetLineColor(myColor[1]);
  pgCent_3s1s_trking->SetMarkerStyle(myMarker[1]);
  pgCent_3s1s_trking->SetMarkerSize(2);
 
  // drawing
  // define axes: 

  TH1D *phAxis4 = new TH1D("phAxis4",";;#varUpsilon(ns)/#varUpsilon(1S)",
			   nBin4,centBinsEdges_compare4);
  phAxis4->SetDirectory(0);
  for (int ib=0; ib<nBin4; ib++)
    {
      phAxis4->GetXaxis()->SetBinLabel(ib+1,Form("%s",recoType[ib].c_str()));
    }
  phAxis4->GetXaxis()->CenterLabels(true);
  phAxis4->GetXaxis()->CenterTitle(true);
  phAxis4->GetXaxis()->LabelsOption("h");
  phAxis4->SetTickLength(0,"X");
  phAxis4->GetYaxis()->SetRangeUser(0.,0.5);  


  TH1D *phAxis2 = new TH1D("phAxis2",";;#varUpsilon(ns)/#varUpsilon(1S)",
			   nBin2,centBinsEdges_compare2);
  phAxis2->SetDirectory(0);
  for (int ib=0; ib<nBin2; ib++)
    {
      if(compGlbTrk)
	phAxis2->GetXaxis()->SetBinLabel(ib+1,Form("%s",trkType[ib].c_str()));
      if(compTrg7tev)
	phAxis2->GetXaxis()->SetBinLabel(ib+1,Form("%s",trgType[ib].c_str()));
    }
  phAxis2->GetXaxis()->CenterLabels(true);
  phAxis2->GetXaxis()->CenterTitle(true);
  phAxis2->GetXaxis()->LabelsOption("h");
  phAxis2->SetTickLength(0,"X");
  phAxis2->GetYaxis()->SetRangeUser(0.,0.5);  


  //------------- prompt
  TCanvas *c1 = new TCanvas("c1","c1");
 
  if(compTrking)
    phAxis4->Draw();
  else 
    phAxis2->Draw();
 
  // writing: 
  TLatex lx; 
  lx.SetNDC();
  lx.SetTextFont(42);
  lx.SetTextAlign(11);
  lx.DrawLatex(0.6, 1. - c1->GetTopMargin()*2.5,"CMS Private");

  if(compTrking || compGlbTrk)
    {
      lx.DrawLatex(0.6, 1. - c1->GetTopMargin()*4,Form("%s",sampleLegend[1]));
      lx.DrawLatex(0.6, 1. - c1->GetTopMargin()*5.5,Form("Cent. 50-100 %%"));
    }
  if(compTrg7tev)
    lx.DrawLatex(0.6, 1. - c1->GetTopMargin()*4,Form("%s",sampleLegend[3]));

  if(doPM1) lx.DrawLatex(0.2, 1. - c1->GetTopMargin()*7.5,"|y_{CM}| < 1");
  else lx.DrawLatex(0.2, 1. - c1->GetTopMargin()*7.5,"|y_{CM}| < 2.4");

 
  //------
  TLegend *leg1 = new TLegend(0.25,0.85,0.35,0.95);
  leg1->SetFillStyle(0);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->SetMargin(0.2);
  leg1->SetTextSize(0.035);

  TLine l1;
 
  if(compTrking)
    {
      pgCent_2s1s_aa50100_reco->Draw("P z");
      l1.DrawLine(0,0.15,8,0.15);
    }

  if(compGlbTrk)
    {
      pgCent_2s1s_trking->Draw("P z");
      l1.DrawLine(0,0.13,4,0.13);
    }

  if(compTrg7tev)
    {
      pgCent_2s1s_trig7->Draw("P z");
      pgCent_3s1s_trig7->Draw("P z");
      l1.DrawLine(0,0.2,4,0.2);
      l1.DrawLine(0,0.31,4,0.31);
    }

 
 
  leg1->AddEntry(pgCent_2s1s_trking,"#varUpsilon(2S)/#varUpsilon(1S)","P");
  if(compTrg7tev) leg1->AddEntry(pgCent_3s1s_trking,"#varUpsilon(3S)/#varUpsilon(1S)","P");
  
  leg1->Draw();
 
  
   gPad->RedrawAxis();

  
   if(bSavePlots)
     { 
       c1->SaveAs(Form("figs/checks/%s_compTrking%d_compGT%d_compTrg%d_ypm1%d.pdf",figNamePrefix,compTrking,compGlbTrk,compTrg7tev,doPM1));
       c1->SaveAs(Form("figs/checks/%s_compTrking%d_compGT%d_compTrg%d_ypm1%d.png",figNamePrefix,compTrking,compGlbTrk,compTrg7tev,doPM1));          
     }
   return;
}

