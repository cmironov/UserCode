#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TLegend.h"

#include "dats/data_2013.h"
#endif


void plot_Rcp_rap(bool isPaper=false, 
		  bool doCmsYearComp=false,
		  bool plotSingleRatio=true
		)
{
  gROOT->Macro("/Users/eusmartass/Documents/cmswrite/hin-10-006/logon.C+");
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(5);
  
  int nEtaBins = 2; 

  TGraphErrors *gEta_rMB50100_1s   = new TGraphErrors(nEtaBins, etaBins, ups_rMB50100_1s_EtaBins, etaBinsError, ups_rMB50100_1s_EtaBinsError);
  TGraphErrors *gEta_r02050100_1s   = new TGraphErrors(nEtaBins, etaBins, ups_r02050100_1s_EtaBins, etaBinsError, ups_r02050100_1s_EtaBinsError);
 
  TGraphErrors *gEta_rMB50100_1s2s   = new TGraphErrors(nEtaBins, etaBins, ups_rMB50100_1s2s_EtaBins, etaBinsError, ups_rMB50100_1s2s_EtaBinsError);
  TGraphErrors *gEta_r02050100_1s2s   = new TGraphErrors(nEtaBins, etaBins, ups_r02050100_1s2s_EtaBins, etaBinsError, ups_r02050100_1s2s_EtaBinsError);

  TGraphErrors *gEta_rMB50100_1s3s   = new TGraphErrors(nEtaBins, etaBins, ups_rMB50100_1s3s_EtaBins, etaBinsError, ups_rMB50100_1s3s_EtaBinsError);
  TGraphErrors *gEta_r02050100_1s3s   = new TGraphErrors(nEtaBins, etaBins, ups_r02050100_1s3s_EtaBins, etaBinsError, ups_r02050100_1s3s_EtaBinsError);
 
  // drawing 2010
  gEta_rMB50100_1s->SetMarkerColor(kRed);
  gEta_rMB50100_1s->SetMarkerStyle(24);
  gEta_r02050100_1s->SetMarkerColor(kRed);
  gEta_r02050100_1s->SetMarkerStyle(20);

  gEta_rMB50100_1s2s->SetMarkerColor(kBlue);
  gEta_rMB50100_1s2s->SetMarkerStyle(24);
  gEta_r02050100_1s2s->SetMarkerColor(kBlue);
  gEta_r02050100_1s2s->SetMarkerStyle(20);

  gEta_rMB50100_1s3s->SetMarkerColor(kGreen+2);
  gEta_rMB50100_1s3s->SetMarkerStyle(25);
  gEta_r02050100_1s3s->SetMarkerColor(kGreen+2);
  gEta_r02050100_1s3s->SetMarkerStyle(21);

  // marker size
  gEta_rMB50100_1s->SetMarkerSize(1.2);
  gEta_rMB50100_1s2s->SetMarkerSize(1.2);
  gEta_r02050100_1s->SetMarkerSize(1.2);
  gEta_r02050100_1s2s->SetMarkerSize(1.2);

  gEta_rMB50100_1s3s->SetMarkerSize(1.2);
  gEta_r02050100_1s3s->SetMarkerSize(1.2);
  // drawing
  // general stuff
  TF1 *f1;
  f1 = new TF1("f1","1",-1.93,1.93);
  f1->SetLineWidth(1);
  f1->GetXaxis()->SetTitle("y_{CM}");
  f1->GetYaxis()->SetTitle("Cen. X/(Cent. 50-100%)");
  f1->GetYaxis()->SetRangeUser(0.0,2.5);
  f1->GetXaxis()->CenterTitle(kTRUE);

  // sqrt(sig_lumi(6%)*sig_lumi(6%)+sig_taa(5.7%)*sig_taa(5.7%)) = 0.083
  TBox *lumi = new TBox(0.0,0.917,0.1,1.083);
  lumi->SetFillColor(kGray+1);

  TLatex *pre = new TLatex(0.2,2.3,"CMS Preliminary");//0.78125
  pre->SetTextSize(0.05);
  pre->Draw();

  TLatex *l1 = new TLatex(0.2,2.1,"pPb #sqrt{s_{NN}} = 5 TeV");
  l1->SetTextSize(0.05);
 
  TLatex *latex1 = new TLatex(0.2,1.9,Form("L_{int} = 5.7 nb^{-1}"));
  latex1->SetTextSize(0.05);

  //------------- prompt
  TCanvas *c1 = new TCanvas("c1","c1");
  f1->Draw();
  f1->Draw("same");
  pre->Draw();
  l1->Draw();
  latex1->Draw();
 
  TLegend *leg10 = new TLegend(0.5,0.5,0.8,0.6);
  leg10->SetFillStyle(0);
  leg10->SetFillColor(0);
  leg10->SetBorderSize(0);
  leg10->SetMargin(0.2);
  leg10->SetTextSize(0.035);
  if(plotSingleRatio)
    {
      gEta_rMB50100_1s->Draw("PZ");
      gEta_r02050100_1s->Draw("PZ");
      leg10->SetTextColor(kRed);
      leg10->AddEntry(gEta_rMB50100_1s,"#Upsilon(1S) P=50-100%","");
      leg10->AddEntry(gEta_rMB50100_1s,"C=0-100%","P");
      leg10->AddEntry(gEta_r02050100_1s,"C=0-20%","P");
     
    }
  else
    {
      gEta_rMB50100_1s3s->Draw("PZ");
      gEta_r02050100_1s3s->Draw("PZ"); 
      gEta_rMB50100_1s2s->Draw("PZ");
      gEta_r02050100_1s2s->Draw("PZ"); 

      leg10->SetTextColor(kGreen+2);
      leg10->AddEntry(gEta_rMB50100_1s3s,"#Upsilon(3S)/#Upsilon(1S) P=50-100%","");
      leg10->AddEntry(gEta_rMB50100_1s3s,"C=0-100%","P");
      leg10->AddEntry(gEta_r02050100_1s3s,"C=0-20%","P");
    }
  leg10->Draw();
  TLegend *leg1 = new TLegend(0.15,0.3,0.3,0.45);
  leg1->SetFillStyle(0);
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->SetMargin(0.2);
  leg1->SetTextSize(0.035);
  leg1->SetTextColor(kBlue);
 if(!plotSingleRatio)
   {
     leg1->AddEntry(gEta_rMB50100_1s2s,"#Upsilon(2S)/#Upsilon(1S) P=50-100%","");
     leg1->AddEntry(gEta_rMB50100_1s2s,"C=0-100%","P");
     leg1->AddEntry(gEta_r02050100_1s2s,"C=0-20%","P");
     leg1->Draw();
   }


  if(plotSingleRatio)
    {
      c1->SaveAs("ups_RcpVsY_single_slides.pdf");
      c1->SaveAs("ups_RcpVsY_single_slides.png");
    }
  else
    {
      c1->SaveAs("ups_RcpVsY_double_slides.pdf");
      c1->SaveAs("ups_RcpVsY_double_slides.png");
    }

  gPad->RedrawAxis();

  return;
}
