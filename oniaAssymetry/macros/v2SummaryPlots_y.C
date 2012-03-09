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


void v2SummaryPlots_y()
{
  //gROOT->Macro("/Users/eusmartass/Software/utilities/setStyle.C+");
  gStyle->SetOptFit(0);

  const char* signal[5]      = {"","NSig","NPr","NNp","NBkg"};
  const char* legend[5]      = {"","J/#psi","Prompt J/#psi","Non-prompt J/#psi","Background"};
  int choseSignal            = 1; // 1:inclusive 2:prompt 3:non-prompt
  const char* chosenSignal   = signal[choseSignal];
  const char* outputName[1]  = {"nominalFit_etaGap"};
  // options
  bool bSavePlots      = true; 
  bool bAddBkg         = false;
  bool bAddVtx         = false;
  bool bAddNoFlat      = true;
  bool bAddTrigger     = false;

  const char* eventPlane[2] = {"","EP: etHFp & etHFm"};
  double rapIntegrated[2]   = {0.0, 2.4}; 
  double ptIntegrated[3]    = {3, 6.5, 40}; 
  // pt_rap bins, high_pt
  const int nybins_highPt = 3;
  double ybins_highPt_center[nybins_highPt]          = {0.6, 1.4, 2.0};
  double ybins_highPt_center_err[nybins_highPt]      = {0.6, 0.2, 0.4};
  double ybins_highPt_center_errGhost[nybins_highPt] = {0.,0.,0.};

  // sgn
  double jpsi_highPt_ptybins[nybins_highPt]          = {0.071, 0.080, 0.064};
  double jpsi_highPt_ptybins_err[nybins_highPt]      = {0.017, 0.026, 0.024};
  double jpsi_highPt_ptybins_errGhost[nybins_highPt] = {0.,0.,0 };

  // bkg
  double bkg_highPt_ptybins[nybins_highPt]         = {0.0418, 0.1694, 0.1378};
  double bkg_highPt_ptybins_err[nybins_highPt]     = {0.0303, 0.0304, 0.0257};

  // correlations
  double jpsiAutoCorr_highPt_ptybins[nybins_highPt]        = {0.1255,0.0593};
  double jpsiAutoCorr_highPt_ptybins_err[nybins_highPt]    = {0.0182,0.0265};

  // finer bins 
  const int nyfinebins_highPt = 5;
  double yfinebins_highPt_center[nyfinebins_highPt]          = {0.3, 0.9, 1.4, 1.8, 2.2};
  double yfinebins_highPt_center_err[nyfinebins_highPt]      = {0.3, 0.3, 0.2, 0.2, 0.2};
  double yfinebins_highPt_center_errGhost[nyfinebins_highPt] = {0.,  0.,  0.,  0.,  0.};

  double jpsiCorr_highPt_ptyfinebins[nyfinebins_highPt]      = {
    0.0703,
    0.0738,
    0.0800,
    0.0685,
    0.0590
  };

  double jpsiCorr_highPt_ptyfinebins_err[nyfinebins_highPt]  = {
    0.0255,
    0.0231,
    0.0256,
    0.0280,
    0.0486
  };

  // triggers
  double jpsiL1NHitTrig_highPt_ptybins[nybins_highPt]      = {0.0585, 0.0522, 0.0426};
  double jpsiL1NHitTrig_highPt_ptybins_err[nybins_highPt]  = {0.0180, 0.0273, 0.0264};

  double jpsiL1NHitTrig_highPt_noCow_ptybins[nybins_highPt]      = {0.0508, 0.0103, 0.0251};
  double jpsiL1NHitTrig_highPt_noCow_ptybins_err[nybins_highPt]  = {0.0251, 0.0386, 0.0335};

  double jpsiL1NHitTrig_highPt_onlyCow_ptybins[nybins_highPt]      = {0.0709, 0.0950, 0.0738};
  double jpsiL1NHitTrig_highPt_onlyCow_ptybins_err[nybins_highPt]  = {0.0258, 0.0385, 0.0426};

  // vtx10cm signal
  double jpsiVtx10_highPt_ptybins[nybins_highPt]           = {0.0710, 0.0826, 0.0747};
  double jpsiVtx10_highPt_ptybins_err[nybins_highPt]       = {0.0185, 0.0274, 0.0259};

  // pt_rap bins, high_pt
  const int nybins_lowPt = 1;
  double ybins_lowPt_center[nybins_lowPt]          = {2.0};
  double ybins_lowPt_center_err[nybins_lowPt]      = {0.0};

  //sgn
  double jpsi_lowPt_ptybins[nybins_lowPt]          = {0.028};
  double jpsi_lowPt_ptybins_err[nybins_lowPt]      = {0.022};

  // bkg
  double bkg_lowPt_ptybins[nybins_lowPt]           = {0.0972};
  double bkg_lowPt_ptybins_err[nybins_lowPt]       = {0.0091};

  // vtx10cm signal
  double jpsiVtx10_lowPt_ptybins[nybins_lowPt]           = {0.0564};
  double jpsiVtx10_lowPt_ptybins_err[nybins_lowPt]       = {0.0246};

  // noFlattening
  double jpsiNoFlat_lowPt_ptybins[nybins_lowPt]          = {0.0786};
  double jpsiNoFlat_lowPt_ptybins_err[nybins_lowPt]      = {0.0244};

  // triggers
  double jpsiL1NHitTrig_lowPt_ptybins[nybins_lowPt]      = {0.0638};
  double jpsiL1NHitTrig_lowPt_ptybins_err[nybins_lowPt]  = {0.0251};

  double jpsiL2pt3Trig_lowPt_ptybins[nybins_lowPt]       = {0.0710};
  double jpsiL2pt3Trig_lowPt_ptybins_err[nybins_lowPt]   = {0.0659};

  double jpsiL3ptOpenTrig_lowPt_ptybins[nybins_lowPt]    = {0.0216};
  double jpsiL3ptOpenTrig_lowPt_ptybins_err[nybins_lowPt]= {0.0703};

  double jpsiL1NHitTrig_lowPt_noCow_ptybins[nybins_lowPt]      = {0.0550};
  double jpsiL1NHitTrig_lowPt_noCow_ptybins_err[nybins_lowPt]  = {0.0338};

  double jpsiL1NHitTrig_lowPt_onlyCow_ptybins[nybins_lowPt]      = {0.0761};
  double jpsiL1NHitTrig_lowPt_onlyCow_ptybins_err[nybins_lowPt]  = {0.0374};

  // autocorrelations
  double jpsiAutoCorr_lowPt_ptybins[nybins_lowPt]        = {0.1418};
  double jpsiAutoCorr_lowPt_ptybins_err[nybins_lowPt]    = {0.0246};


  TGraphErrors *pg_jpsi_highPt_ptybins      = new TGraphErrors(nybins_highPt,ybins_highPt_center,jpsi_highPt_ptybins,ybins_highPt_center_errGhost,jpsi_highPt_ptybins_err); 
  TGraphErrors *pg_bkg_highPt_ptybins       = new TGraphErrors(nybins_highPt,ybins_highPt_center,bkg_highPt_ptybins,ybins_highPt_center_errGhost,bkg_highPt_ptybins_err); 
  TGraphErrors *pg_jpsi_highPt_ptybinsGhost = new TGraphErrors(nybins_highPt,ybins_highPt_center,jpsi_highPt_ptybins,ybins_highPt_center_err,jpsi_highPt_ptybins_errGhost);   

  TGraphErrors *pg_jpsiAutoCorr_highPt_ptybins    = new TGraphErrors(nybins_highPt,ybins_highPt_center,jpsiAutoCorr_highPt_ptybins,ybins_highPt_center_err,jpsiAutoCorr_highPt_ptybins_err);
  TGraphErrors *pg_jpsiAutoCorr_lowPt_ptybins     = new TGraphErrors(nybins_lowPt,ybins_lowPt_center,jpsiAutoCorr_lowPt_ptybins,ybins_lowPt_center_err,jpsiAutoCorr_lowPt_ptybins_err);

  // correlations in y:
  TGraphErrors *pg_jpsiCorr_highPt_ptyfinebins    = new TGraphErrors(nyfinebins_highPt,yfinebins_highPt_center,jpsiCorr_highPt_ptyfinebins,yfinebins_highPt_center_err,jpsiCorr_highPt_ptyfinebins_err);

  TGraphErrors *pg_jpsi_lowPt_ptybins             = new TGraphErrors(nybins_lowPt,ybins_lowPt_center,jpsi_lowPt_ptybins,ybins_lowPt_center_err,jpsi_lowPt_ptybins_err);
  TGraphErrors *pg_bkg_lowPt_ptybins              = new TGraphErrors(nybins_lowPt,ybins_lowPt_center,bkg_lowPt_ptybins,ybins_lowPt_center_err,bkg_lowPt_ptybins_err);
  TGraphErrors *pg_jpsiVtx10_lowPt_ptybins        = new TGraphErrors(nybins_lowPt,ybins_lowPt_center,jpsiVtx10_lowPt_ptybins,ybins_lowPt_center_err,jpsiVtx10_lowPt_ptybins_err);
  TGraphErrors *pg_jpsiVtx10_highPt_ptybins        = new TGraphErrors(nybins_highPt,ybins_highPt_center,jpsiVtx10_highPt_ptybins,ybins_highPt_center_err,jpsiVtx10_highPt_ptybins_err);
  TGraphErrors *pg_jpsiNoFlat_lowPt_ptybins       = new TGraphErrors(nybins_lowPt,ybins_lowPt_center,jpsiNoFlat_lowPt_ptybins,ybins_lowPt_center_err,jpsiNoFlat_lowPt_ptybins_err);
  TGraphErrors *pg_jpsiL1NHitTrig_lowPt_ptybins   = new TGraphErrors(nybins_lowPt,ybins_lowPt_center,jpsiL1NHitTrig_lowPt_ptybins,ybins_lowPt_center_err,jpsiL1NHitTrig_lowPt_ptybins_err);
  TGraphErrors *pg_jpsiL1NHitTrig_lowPt_noCow_ptybins   = new TGraphErrors(nybins_lowPt,ybins_lowPt_center,jpsiL1NHitTrig_lowPt_noCow_ptybins,ybins_lowPt_center_err,jpsiL1NHitTrig_lowPt_noCow_ptybins_err);
  TGraphErrors *pg_jpsiL1NHitTrig_lowPt_onlyCow_ptybins   = new TGraphErrors(nybins_lowPt,ybins_lowPt_center,jpsiL1NHitTrig_lowPt_onlyCow_ptybins,ybins_lowPt_center_err,jpsiL1NHitTrig_lowPt_onlyCow_ptybins_err);
  TGraphErrors *pg_jpsiL2pt3Trig_lowPt_ptybins    = new TGraphErrors(nybins_lowPt,ybins_lowPt_center,jpsiL2pt3Trig_lowPt_ptybins,ybins_lowPt_center_err,jpsiL2pt3Trig_lowPt_ptybins_err);
  TGraphErrors *pg_jpsiL3ptOpenTrig_lowPt_ptybins = new TGraphErrors(nybins_lowPt,ybins_lowPt_center,jpsiL3ptOpenTrig_lowPt_ptybins,ybins_lowPt_center_err,jpsiL3ptOpenTrig_lowPt_ptybins_err);

  TGraphErrors *pg_jpsiL1NHitTrig_highPt_ptybins   = new TGraphErrors(nybins_highPt,ybins_highPt_center,jpsiL1NHitTrig_highPt_ptybins,ybins_highPt_center_err,jpsiL1NHitTrig_highPt_ptybins_err);
  TGraphErrors *pg_jpsiL1NHitTrig_highPt_noCow_ptybins   = new TGraphErrors(nybins_highPt,ybins_highPt_center,jpsiL1NHitTrig_highPt_noCow_ptybins,ybins_highPt_center_err,jpsiL1NHitTrig_highPt_noCow_ptybins_err);
  TGraphErrors *pg_jpsiL1NHitTrig_highPt_onlyCow_ptybins   = new TGraphErrors(nybins_highPt,ybins_highPt_center,jpsiL1NHitTrig_highPt_onlyCow_ptybins,ybins_highPt_center_err,jpsiL1NHitTrig_highPt_onlyCow_ptybins_err);

  // rapidity dependence
  TCanvas *pcY = new TCanvas("pcY","pcY");
  pcY->cd();
  TH1F *pcPadY = new TH1F("pcPadY",";|y|;v_{2};",100,0,2.4);
  pcPadY->GetXaxis()->SetLabelSize(20);
  pcPadY->GetXaxis()->SetLabelFont(43);
  pcPadY->GetXaxis()->SetTitleSize(27);
  pcPadY->GetXaxis()->SetTitleFont(43);
  pcPadY->GetXaxis()->SetTitleOffset(1.2);
  pcPadY->GetXaxis()->CenterTitle();

  pcPadY->GetYaxis()->SetLabelSize(20);
  pcPadY->GetYaxis()->SetLabelFont(43);
  pcPadY->GetYaxis()->SetTitleSize(32);
  pcPadY->GetYaxis()->SetTitleFont(43);
  pcPadY->GetYaxis()->SetTitleOffset(1.1);
  pcPadY->GetYaxis()->CenterTitle();

  pcPadY->SetMaximum(0.25);
  pcPadY->SetMinimum(-0.1);
  //pcPadY->SetMaximum(0.18);
  //pcPadY->SetMinimum(-0.05);

  if(bAddBkg)
  {
    pcPadY->SetMaximum(0.4);
  }
  pcPadY->Draw();
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
  //-----------------------


  pg_jpsi_highPt_ptybins->SetMarkerStyle(20);
  pg_jpsi_highPt_ptybins->SetMarkerSize(2.0);
  pg_jpsi_highPt_ptybins->SetMarkerColor(kRed+2);
  pg_jpsi_highPt_ptybins->SetLineColor(kRed+2);
  pg_jpsi_highPt_ptybins->Draw("[P]");

  pg_jpsi_highPt_ptybinsGhost->SetMarkerStyle(21);
  pg_jpsi_highPt_ptybinsGhost->SetMarkerSize(1.8);
  pg_jpsi_highPt_ptybinsGhost->SetMarkerColor(kBlue+2);
  pg_jpsi_highPt_ptybinsGhost->SetLineColor(kBlue+2);
  //pg_jpsi_highPt_ptybinsGhost->Draw("P");

  pg_jpsi_lowPt_ptybins->SetMarkerStyle(21);
  pg_jpsi_lowPt_ptybins->SetMarkerSize(1.8);
  pg_jpsi_lowPt_ptybins->SetMarkerColor(kBlue+2);
  pg_jpsi_lowPt_ptybins->SetLineColor(kBlue+2);
  pg_jpsi_lowPt_ptybins->Draw("[P]");

  if(bAddBkg)
  {
    pg_bkg_highPt_ptybins->SetMarkerStyle(24);
    pg_bkg_highPt_ptybins->SetMarkerSize(1.8);
    pg_bkg_highPt_ptybins->SetMarkerColor(kRed+2);
    pg_bkg_highPt_ptybins->Draw("[P]");

    pg_bkg_lowPt_ptybins->SetMarkerStyle(25);
    pg_bkg_lowPt_ptybins->SetMarkerSize(1.8);
    pg_bkg_lowPt_ptybins->SetMarkerColor(kBlue+2);
    pg_bkg_lowPt_ptybins->SetLineColor(kBlue+2);
    pg_bkg_lowPt_ptybins->Draw("[P]");

    TLegend *legY = new TLegend(0.7,0.14,0.9,0.26);
    legY->SetFillColor(0);
    legY->SetBorderSize(0);
    legY->SetTextSize(0.03);
    legY->AddEntry(pg_jpsi_lowPt_ptybins,"Signal","P");
    legY->AddEntry(pg_bkg_lowPt_ptybins,"Bkg","P");
    legY->Draw("same");
  }

  if(bAddVtx)
  {
    pg_jpsiVtx10_highPt_ptybins->SetMarkerStyle(24);
    pg_jpsiVtx10_highPt_ptybins->SetMarkerSize(1.8);
    pg_jpsiVtx10_highPt_ptybins->SetMarkerColor(kBlue+2);
    pg_jpsiVtx10_highPt_ptybins->SetLineColor(kBlue+2);
    pg_jpsiVtx10_highPt_ptybins->Draw("[p]");

    pg_jpsiVtx10_lowPt_ptybins->SetMarkerStyle(25);
    pg_jpsiVtx10_lowPt_ptybins->SetMarkerSize(1.8);
    pg_jpsiVtx10_lowPt_ptybins->SetMarkerColor(kBlue+2);
    pg_jpsiVtx10_lowPt_ptybins->SetLineColor(kBlue+2);
    pg_jpsiVtx10_lowPt_ptybins->Draw("[p]");

    TLegend *legVtx = new TLegend(0.55,0.15,0.875,0.297);
    legVtx->SetFillColor(0);
    legVtx->SetBorderSize(0);
    legVtx->SetTextSize(0.03);
    legVtx->AddEntry(pg_jpsi_highPt_ptybins,"Default: |z_{vtx}|<25 cm","P");
    legVtx->AddEntry(pg_jpsiVtx10_highPt_ptybins,"|z_{vtx}|<10 cm","P");
    legVtx->Draw("same");
  }


  if(bAddNoFlat)
  {
    pg_jpsiNoFlat_lowPt_ptybins->SetMarkerStyle(24);
    pg_jpsiNoFlat_lowPt_ptybins->SetMarkerSize(1.8);
    pg_jpsiNoFlat_lowPt_ptybins->SetMarkerColor(kBlue+2);
    //  pg_jpsiNoFlat_lowPt_ptybins->Draw("p");

    pg_jpsiAutoCorr_highPt_ptybins->SetMarkerStyle(27);
    pg_jpsiAutoCorr_highPt_ptybins->SetMarkerSize(1.8);
    pg_jpsiAutoCorr_highPt_ptybins->SetMarkerColor(kRed+2);
    pg_jpsiAutoCorr_highPt_ptybins->SetLineColor(kRed+2);
    //pg_jpsiAutoCorr_highPt_ptybins->Draw("PZ");

    pg_jpsiAutoCorr_lowPt_ptybins->SetMarkerStyle(27);
    pg_jpsiAutoCorr_lowPt_ptybins->SetMarkerSize(1.8);
    pg_jpsiAutoCorr_lowPt_ptybins->SetMarkerColor(kBlue);
    //  pg_jpsiAutoCorr_lowPt_ptybins->Draw("p");

    pg_jpsiCorr_highPt_ptyfinebins->SetMarkerStyle(24);
    pg_jpsiCorr_highPt_ptyfinebins->SetMarkerSize(1.8);
    pg_jpsiCorr_highPt_ptyfinebins->SetMarkerColor(kRed+2);
    pg_jpsiCorr_highPt_ptyfinebins->SetLineColor(kRed+2);
    //pg_jpsiCorr_highPt_ptyfinebins->SetMarkerColor(kBlue+7);
    pg_jpsiCorr_highPt_ptyfinebins->Draw("[P]");

    TLegend *legEp = new TLegend(0.55,0.15,0.875,0.297);
    legEp->SetFillColor(0);
    legEp->SetBorderSize(0);
    legEp->SetTextSize(0.03);
    //legEp->AddEntry(pg_jpsi_lowPt_ptybins,"Default: w/ flattening","P");
    //legEp->AddEntry(pg_jpsiNoFlat_lowPt_ptybins,"No flattening","P");

    //legEp->AddEntry(pg_jpsi_lowPt_ptybins,"+- && -+","P");
    //legEp->AddEntry(pg_jpsiAutoCorr_lowPt_ptybins,"++ && --","P");

    legEp->Draw("same");

  }

  if(bAddTrigger)
  {
    pg_jpsiL1NHitTrig_lowPt_ptybins->SetMarkerStyle(25);
    pg_jpsiL1NHitTrig_lowPt_onlyCow_ptybins->SetMarkerStyle(25);
    pg_jpsiL1NHitTrig_lowPt_noCow_ptybins->SetMarkerStyle(25);

    pg_jpsiL2pt3Trig_lowPt_ptybins->SetMarkerStyle(27);
    pg_jpsiL3ptOpenTrig_lowPt_ptybins->SetMarkerStyle(28);

    pg_jpsiL1NHitTrig_lowPt_ptybins->SetMarkerColor(kBlue+2);
    pg_jpsiL1NHitTrig_lowPt_ptybins->SetLineColor(kBlue+2);
    pg_jpsiL1NHitTrig_highPt_ptybins->SetMarkerColor(kBlue+2);
    pg_jpsiL1NHitTrig_highPt_ptybins->SetLineColor(kBlue+2);

    pg_jpsiL2pt3Trig_lowPt_ptybins->SetMarkerColor(kBlue+2);
    pg_jpsiL3ptOpenTrig_lowPt_ptybins->SetMarkerColor(kBlue+2);

    pg_jpsiL1NHitTrig_lowPt_ptybins->SetMarkerSize(1.8);
    pg_jpsiL2pt3Trig_lowPt_ptybins->SetMarkerSize(1.8);
    pg_jpsiL3ptOpenTrig_lowPt_ptybins->SetMarkerSize(1.8);

    pg_jpsiL1NHitTrig_highPt_ptybins->SetMarkerStyle(24);
    pg_jpsiL1NHitTrig_highPt_onlyCow_ptybins->SetMarkerStyle(24);
    pg_jpsiL1NHitTrig_highPt_noCow_ptybins->SetMarkerStyle(24);


    pg_jpsiL1NHitTrig_highPt_onlyCow_ptybins->SetMarkerColor(kGreen+2);
    pg_jpsiL1NHitTrig_highPt_onlyCow_ptybins->SetLineColor(kGreen+2);
    pg_jpsiL1NHitTrig_highPt_noCow_ptybins->SetMarkerColor(kMagenta);
    pg_jpsiL1NHitTrig_highPt_noCow_ptybins->SetLineColor(kMagenta);

    pg_jpsiL1NHitTrig_lowPt_onlyCow_ptybins->SetMarkerColor(kGreen+2);
    pg_jpsiL1NHitTrig_lowPt_onlyCow_ptybins->SetLineColor(kGreen+2);
    pg_jpsiL1NHitTrig_lowPt_noCow_ptybins->SetMarkerColor(kMagenta);
    pg_jpsiL1NHitTrig_lowPt_noCow_ptybins->SetLineColor(kMagenta);

    pg_jpsiL1NHitTrig_highPt_ptybins->SetMarkerSize(1.8);
    pg_jpsiL1NHitTrig_highPt_onlyCow_ptybins->SetMarkerSize(1.8);
    pg_jpsiL1NHitTrig_highPt_noCow_ptybins->SetMarkerSize(1.8);

    pg_jpsiL1NHitTrig_lowPt_ptybins->SetMarkerSize(1.8);
    pg_jpsiL1NHitTrig_lowPt_onlyCow_ptybins->SetMarkerSize(1.8);
    pg_jpsiL1NHitTrig_lowPt_noCow_ptybins->SetMarkerSize(1.8);

    pg_jpsi_lowPt_ptybins->Draw("[P]");
    pg_jpsi_highPt_ptybins->Draw(" [P] ");
    pg_jpsiL1NHitTrig_lowPt_ptybins->Draw("[P]");
    pg_jpsiL1NHitTrig_highPt_ptybins->Draw("[P]");
    pg_jpsiL1NHitTrig_lowPt_noCow_ptybins->Draw("[P]");
    pg_jpsiL1NHitTrig_highPt_noCow_ptybins->Draw("[P]");
    pg_jpsiL1NHitTrig_lowPt_onlyCow_ptybins->Draw("[P]");
    pg_jpsiL1NHitTrig_highPt_onlyCow_ptybins->Draw("[P]");

    //pg_jpsiL2pt3Trig_lowPt_ptybins->Draw("p");
    //pg_jpsiL3ptOpenTrig_lowPt_ptybins->Draw("p");

    TLegend *legTrig = new TLegend(0.18,0.16,0.47,0.33);
    //TLegend *legTrig = new TLegend(0.18,0.16,0.47,0.30);
    legTrig->SetFillColor(0);
    legTrig->SetBorderSize(0);
    legTrig->SetTextSize(0.03);
    legTrig->AddEntry(pg_jpsi_highPt_ptybins,"Default: all triggers","P");
    legTrig->AddEntry(pg_jpsiL1NHitTrig_highPt_ptybins,"HLT_HIL1DoubleMu0_HighQ","P");
    //legTrig->AddEntry(pg_jpsiL2pt3Trig_lowPt_ptybins,"HLT_HIL2DoubleMu3","P");
    //legTrig->AddEntry(pg_jpsiL3ptOpenTrig_lowPt_ptybins,"HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy","P");
    legTrig->AddEntry(pg_jpsiL1NHitTrig_highPt_noCow_ptybins,"HLT_HIL1DoubleMu0_HighQ+Sailor","P");
    legTrig->AddEntry(pg_jpsiL1NHitTrig_highPt_onlyCow_ptybins,"HLT_HIL1DoubleMu0_HighQ+Cowboy","P");
    legTrig->Draw("same");
  }

  TLatex *lt1 = new TLatex();
  lt1->SetNDC();
  lt1->SetTextSize(0.04);
  lt1->DrawLatex(0.18,0.89,Form("Inclusive %s",legend[choseSignal]));  // what signal is
  lt1->SetTextSize(0.038);
  lt1->DrawLatex(0.18,0.83,Form("Cent. 10 - 60 %%")); 
  // lt1->DrawLatex(0.18,0.77,Form("%s",eventPlane[1]));

  if(!bAddTrigger)
  {
    TLegend *leg1 = new TLegend(0.21,0.14,0.44,0.26);
    leg1->SetFillColor(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.03);
    leg1->AddEntry(pg_jpsi_highPt_ptybins,Form("%.1f < p_{T} < %.0f GeV/c", ptIntegrated[1], ptIntegrated[2]),"P");
    leg1->AddEntry(pg_jpsi_lowPt_ptybins,Form("%.0f < p_{T} < %.1f GeV/c", ptIntegrated[0], ptIntegrated[1]),"P");
    leg1->Draw("same");
  }
  else
  {
    TLegend *leg1 = new TLegend(0.62,0.25,0.85,0.36);
    //TLegend *leg1 = new TLegend(0.44,0.35,0.67,0.46);
    leg1->SetFillColor(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.03);
    //   leg1->AddEntry(pg_jpsi_highPt_ptybins,Form("%.1f < p_{T} < %.0f GeV/c", ptIntegrated[1], ptIntegrated[2]),"P");
    leg1->AddEntry(pg_jpsi_lowPt_ptybins,Form("%.0f < p_{T} < %.1f GeV/c", ptIntegrated[0], ptIntegrated[1]),"P");
    leg1->Draw("same");
  }
  if(bSavePlots)
  {
    pcY->SaveAs(Form("figs/png/%s_%s_v2_ptyBins.png",chosenSignal,outputName[0]));
    pcY->SaveAs(Form("figs/pdf/%s_%s_v2_ptyBins.pdf",chosenSignal,outputName[0]));
  }
  //________________________________________________________________________________________

}



