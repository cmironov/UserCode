#if !defined(__CINT__) || defined(__MAKECINT__)

#include <Riostream.h>
#include <TSystem.h>
#include <TProfile.h>
#include <TBrowser.h>
#include <TROOT.h>
#include <TGraph.h>
#include <TNtuple.h>
#include <TString.h>
#include <TH1D.h>
#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TInterpreter.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>

#endif


void v2SummaryPlots_cent()
{
  gROOT->Macro("/Users/eusmartass/Software/utilities/setStyle.C+");
  gStyle->SetOptFit(0);

  const char* signal[5]   = {"","NSig","NPr","NNp","NBkg"};
  const char* legend[5]   = {"","J/#psi","Prompt J/#psi","Non-prompt J/#psi","Background"};
  int choseSignal         = 1; // 1:inclusive 2:prompt 3:non-prompt
  const char* chosenSignal= signal[choseSignal];
  const char* outputName[1]  = {"nominalFit_trig"};

  // options
  bool bSavePlots      = true; 
  bool bAddBkg         = false;
  bool bAddVtx         = false;
  bool bAddNoFlat      = false;
  bool bAddTrigger     = true;

  const char* eventPlane[2] = {"","EP: etHFp & etHFm"};
  double rapIntegrated[2]   = {0.0, 2.4}; 
  double ptIntegrated[3]    = {3, 6.5, 40};  
  // inclusive numbers
  // centrlaity bins: 0-5, 5-10, 10-30, 30-60
  const int ncentbins1                    = 4;
  double ncoll1[ncentbins1]                = {381.3,   329.4,  224.3,   89.9};
  // centrlaity bins: 0-10, 10-30, 30-60
  const int ncentbins                    = 4;
  double ncoll[ncentbins]                = {355.4, 261.4178, 187.1470, 89.9};
  double ncoll_err[ncentbins]            = { 0.,     0.,       0.,      0.};

  // sgn
  double jpsi_centbins[ncentbins]        = {0.039, 0.074, 0.067, 0.067};
  double jpsi_centbins_err[ncentbins]    = {0.025, 0.022, 0.020,0.021};

  // bkg
  double bkg_centbins[ncentbins]        = {0.0999, 0.1116, 0.2048};
  double bkg_centbins_err[ncentbins]    = {0.0253, 0.0198, 0.0400};
  
   // vtx10cm signal
  double jpsiVtx10_centbins[ncentbins]        = {0.0546, 0.0815, 0.0353};
  double jpsiVtx10_centbins_err[ncentbins]    = {0.0290, 0.0176, 0.0249};

  // noFlattening
  double jpsiNoFlat_centbins[ncentbins]       = {0.0364, 0.0807, 0.0290};
  double jpsiNoFlat_centbins_err[ncentbins]   = {0.0273, 0.0166, 0.0236};

  // autocorrelations
  double jpsiAutoCorr_centbins[ncentbins]     = {0.0563, 0.0677, 0.0718};
  double jpsiAutoCorr_centbins_err[ncentbins] = {0.0272, 0.0165, 0.0231};

  // triggers
  double jpsiL1NHitTrig_centbins[ncentbins]       = {0.0361, 0.0528, 0.0560, 0.0408};
  double jpsiL1NHitTrig_centbins_err[ncentbins]   = {0.0262, 0.0215,0.0238, 0.0220};

  double jpsiL1NHitTrig_onlyCow_centbins[ncentbins]     = {0.0970, 0.1116, 0.0698, 0.0403};
  double jpsiL1NHitTrig_onlyCow_centbins_err[ncentbins] = {0.0384, 0.0313, 0.0338, 0.0327};

  double jpsiL1NHitTrig_noCow_centbins[ncentbins]       = {-0.0231, 0.0045, 0.0380, 0.0434};
  double jpsiL1NHitTrig_noCow_centbins_err[ncentbins]   = { 0.0354, 0.0292, 0.0336, 0.0297};

  double jpsiL2pt3Trig_centbins[ncentbins]       = {0.0395, 0.0610, 0.0231};
  double jpsiL2pt3Trig_centbins_err[ncentbins]   = { 0.0321, 0.0199,0.0278};

  double jpsiL3ptOpenTrig_centbins[ncentbins]    = {0.0005, 0.0237, 0.0629};
  double jpsiL3ptOpenTrig_centbins_err[ncentbins]= { 0.0574, 0.0358, 0.0570};


  TGraphErrors *pg_jpsi_centbins = new TGraphErrors(ncentbins,ncoll,jpsi_centbins,ncoll_err,jpsi_centbins_err);
  TGraphErrors *pg_bkg_centbins  = new TGraphErrors(ncentbins,ncoll,bkg_centbins,ncoll_err,bkg_centbins_err);
 
 
  TGraphErrors *pg_jpsiVtx10_centbins        = new TGraphErrors(ncentbins,ncoll,jpsiVtx10_centbins,ncoll_err,jpsiVtx10_centbins_err);
  TGraphErrors *pg_jpsiNoFlat_centbins       = new TGraphErrors(ncentbins,ncoll,jpsiNoFlat_centbins,ncoll_err,jpsiNoFlat_centbins_err);
  TGraphErrors *pg_jpsiAutoCorr_centbins     = new TGraphErrors(ncentbins,ncoll,jpsiAutoCorr_centbins,ncoll_err,jpsiAutoCorr_centbins_err);
 
  TGraphErrors *pg_jpsiL1NHitTrig_centbins           = new TGraphErrors(ncentbins,ncoll,jpsiL1NHitTrig_centbins,ncoll_err,jpsiL1NHitTrig_centbins_err);
  TGraphErrors *pg_jpsiL1NHitTrig_onlyCow_centbins   = new TGraphErrors(ncentbins,ncoll,jpsiL1NHitTrig_onlyCow_centbins,ncoll_err,jpsiL1NHitTrig_onlyCow_centbins_err);
  TGraphErrors *pg_jpsiL1NHitTrig_noCow_centbins     = new TGraphErrors(ncentbins,ncoll,jpsiL1NHitTrig_noCow_centbins,ncoll_err,jpsiL1NHitTrig_noCow_centbins_err);

  TGraphErrors *pg_jpsiL2pt3Trig_centbins    = new TGraphErrors(ncentbins,ncoll,jpsiL2pt3Trig_centbins,ncoll_err,jpsiL2pt3Trig_centbins_err);
  TGraphErrors *pg_jpsiL3ptOpenTrig_centbins = new TGraphErrors(ncentbins,ncoll,jpsiL3ptOpenTrig_centbins,ncoll_err,jpsiL3ptOpenTrig_centbins_err);
  //_____________________________________________________________________________
  // centrality dependence
  // N_{part} axis
  TCanvas *pcCentrality = new TCanvas("pcCentrality","pcCentrality");
  pcCentrality->cd();
  TH1F *phPadCentrality = new TH1F("phPadCentrality",";N_{part};v_{2};",400,0,400);
  phPadCentrality->GetXaxis()->SetLabelSize(20);
  phPadCentrality->GetXaxis()->SetLabelFont(43);
  phPadCentrality->GetXaxis()->SetTitleSize(30);
  phPadCentrality->GetXaxis()->SetTitleFont(43);
  phPadCentrality->GetXaxis()->SetTitleOffset(1.1);
  phPadCentrality->GetXaxis()->CenterTitle();

  phPadCentrality->GetYaxis()->SetLabelSize(20);
  phPadCentrality->GetYaxis()->SetLabelFont(43);
  phPadCentrality->GetYaxis()->SetTitleSize(32);
  phPadCentrality->GetYaxis()->SetTitleFont(43);
  phPadCentrality->GetYaxis()->SetTitleOffset(1.1);
  phPadCentrality->GetYaxis()->CenterTitle();

  phPadCentrality->SetMaximum(0.18);
  phPadCentrality->SetMinimum(-0.05);
  if(bAddBkg)
    {
      phPadCentrality->SetMaximum(0.4);
    }
  phPadCentrality->Draw();
  
  pg_jpsi_centbins->SetMarkerStyle(20);
  pg_jpsi_centbins->SetMarkerSize(1.8);
  pg_jpsi_centbins->SetMarkerColor(602);
  pg_jpsi_centbins->Draw("pz");

  if(bAddBkg)
    {
      pg_bkg_centbins->SetMarkerStyle(24);
      pg_bkg_centbins->SetMarkerSize(1.8);
      pg_bkg_centbins->SetMarkerColor(kBlue);
      pg_bkg_centbins->Draw("pz");
      TLegend *legCent = new TLegend(0.7,0.5,0.9,0.6);
      legCent->SetFillColor(0);
      legCent->SetBorderSize(0);
      legCent->SetTextSize(0.03);
      legCent->AddEntry(pg_jpsi_centbins,"Signal","P");
      legCent->AddEntry(pg_bkg_centbins,"Bkg","P");
      legCent->Draw("same");
    }

 if(bAddVtx)
    {
      pg_jpsiVtx10_centbins->SetMarkerStyle(24);
      pg_jpsiVtx10_centbins->SetMarkerSize(1.8);
      pg_jpsiVtx10_centbins->SetMarkerColor(kBlue);
      pg_jpsiVtx10_centbins->Draw("p");
      
      TLegend *legVtx = new TLegend(0.19,0.54,0.52,0.69);
      legVtx->SetFillColor(0);
      legVtx->SetBorderSize(0);
      legVtx->SetTextSize(0.03);
      legVtx->AddEntry(pg_jpsi_centbins,"Default: |z_{vtx}|<25 cm","P");
      legVtx->AddEntry(pg_jpsiVtx10_centbins,"|z_{vtx}|<10 cm","P");
      legVtx->Draw("same");
    }

  
  if(bAddNoFlat)
    {
      pg_jpsiNoFlat_centbins->SetMarkerStyle(24);
      pg_jpsiNoFlat_centbins->SetMarkerSize(1.8);
      pg_jpsiNoFlat_centbins->SetMarkerColor(kBlue);
      //   pg_jpsiNoFlat_centbins->Draw("p");

      pg_jpsiAutoCorr_centbins->SetMarkerStyle(27);
      pg_jpsiAutoCorr_centbins->SetMarkerSize(1.8);
      pg_jpsiAutoCorr_centbins->SetMarkerColor(kBlue);
      pg_jpsiAutoCorr_centbins->Draw("p");
      
      TLegend *legEp = new TLegend(0.44,0.31,0.77,0.46);
      legEp->SetFillColor(0);
      legEp->SetBorderSize(0);
      legEp->SetTextSize(0.03);
     //  legEp->AddEntry(pg_jpsi_centbins,"Default: w/ flattening","P");
//       legEp->AddEntry(pg_jpsiNoFlat_centbins,"No flattening","P");

      legEp->AddEntry(pg_jpsi_centbins,"+- && -+","P");
      legEp->AddEntry(pg_jpsiAutoCorr_centbins,"++ && --","P");
      
      legEp->Draw("same");

    }

  if(bAddTrigger)
    {
      pg_jpsiL1NHitTrig_centbins->SetMarkerStyle(24);
      pg_jpsiL2pt3Trig_centbins->SetMarkerStyle(27);
      pg_jpsiL3ptOpenTrig_centbins->SetMarkerStyle(28);
      
      pg_jpsiL1NHitTrig_centbins->SetMarkerColor(kBlue);
      pg_jpsiL2pt3Trig_centbins->SetMarkerColor(kBlue);
      pg_jpsiL3ptOpenTrig_centbins->SetMarkerColor(kBlue);
      
      pg_jpsiL1NHitTrig_centbins->SetMarkerSize(1.8);
      pg_jpsiL2pt3Trig_centbins->SetMarkerSize(1.8);
      pg_jpsiL3ptOpenTrig_centbins->SetMarkerSize(1.8);

      pg_jpsiL1NHitTrig_onlyCow_centbins->SetMarkerStyle(24);
      pg_jpsiL1NHitTrig_noCow_centbins->SetMarkerStyle(24);
	
      pg_jpsiL1NHitTrig_onlyCow_centbins->SetMarkerColor(kGreen+2);
      pg_jpsiL1NHitTrig_noCow_centbins->SetMarkerColor(kMagenta);

      pg_jpsiL1NHitTrig_onlyCow_centbins->SetMarkerSize(1.8);
      pg_jpsiL1NHitTrig_noCow_centbins->SetMarkerSize(1.8);

      pg_jpsiL1NHitTrig_centbins->Draw("P");
      //   pg_jpsiL2pt3Trig_centbins->Draw("p");
      // pg_jpsiL3ptOpenTrig_centbins->Draw("p");

     
      pg_jpsiL1NHitTrig_onlyCow_centbins->Draw("P");
      pg_jpsiL1NHitTrig_noCow_centbins->Draw("P");
      
      TLegend *legTrig = new TLegend(0.2,0.15,0.57,0.29);
      legTrig->SetFillColor(0);
      legTrig->SetBorderSize(0);
      legTrig->SetTextSize(0.03);
      legTrig->AddEntry(pg_jpsi_centbins,"Default: all triggers","P");
      legTrig->AddEntry(pg_jpsiL1NHitTrig_centbins,"HLT_HIL1DoubleMu0_HighQ","P");
      //   legTrig->AddEntry(pg_jpsiL2pt3Trig_centbins,"HLT_HIL2DoubleMu3","P");
      //  legTrig->AddEntry(pg_jpsiL3ptOpenTrig_centbins,"HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy","P");
      legTrig->AddEntry(pg_jpsiL1NHitTrig_noCow_centbins,"HLT_HIL1DoubleMu0_HighQ+NoCowboy","P");
      legTrig->AddEntry(pg_jpsiL1NHitTrig_onlyCow_centbins,"HLT_HIL1DoubleMu0_HighQ+OnlyCowboy","P");
      legTrig->Draw("same");
    }

  //_______ stuff to write
  TLatex *tex_c1 = new TLatex(0.53,0.92,"CMS Preliminary");
  tex_c1->SetNDC();
  tex_c1->SetTextAlign(13);
  tex_c1->SetTextFont(43);
  tex_c1->SetTextSize(25);
  tex_c1->SetLineWidth(1);
  tex_c1->Draw();

  TLatex *tex_c2 = new TLatex(0.53,0.86,"PbPb  #sqrt{s_{NN}} = 2.76 TeV");
  tex_c2->SetNDC();
  tex_c2->SetTextAlign(13);
  tex_c2->SetTextFont(43);
  tex_c2->SetTextSize(25);
  tex_c2->SetLineWidth(2);
  tex_c2->Draw();

  TLatex *tex_c3 = new TLatex(0.53,0.80,"L_{int} = 150 #mub^{-1}");
  tex_c3->SetNDC();
  tex_c3->SetTextAlign(13);
  tex_c3->SetTextFont(43);
  tex_c3->SetTextSize(25);
  tex_c3->SetLineWidth(2);
  tex_c3->Draw();

  TLatex *lt1 = new TLatex();
  lt1->SetNDC();
  lt1->SetTextSize(0.04);
  lt1->DrawLatex(0.18,0.89,Form("Inclusive %s",legend[choseSignal]));  // what signal is
  lt1->SetTextSize(0.038);
  lt1->DrawLatex(0.18,0.83,Form("|y| < %.1f",rapIntegrated[1]));       // rapidity
  lt1->DrawLatex(0.18,0.77,Form("%.1f < p_{T} < %.0f GeV/c", ptIntegrated[1], ptIntegrated[2])); 
  //  lt1->DrawLatex(0.18,0.71,Form("%s",eventPlane[1]));

  if(!bAddTrigger)
    {
      lt1->DrawLatex(0.25, 0.2, Form("30-60%%"));
      lt1->DrawLatex(0.45, 0.2, Form("20-30%%"));
      lt1->DrawLatex(0.6, 0.2, Form("10-20%%"));
      lt1->DrawLatex(0.8, 0.2, Form("0-10%%"));
    }
  else
    {
      lt1->DrawLatex(0.25, 0.36, Form("30-60%%"));
      lt1->DrawLatex(0.45, 0.36, Form("20-30%%"));
      lt1->DrawLatex(0.6, 0.36, Form("10-20%%"));
      lt1->DrawLatex(0.8, 0.36, Form("0-10%%"));
    }

  if(bSavePlots)
    {
      pcCentrality->SaveAs(Form("figs/png/%s_%s_v2_centBins.png",chosenSignal,outputName[0]));
      pcCentrality->SaveAs(Form("figs/pdf/%s_%s_v2_centBins.pdf",chosenSignal,outputName[0]));
    }
 
}



