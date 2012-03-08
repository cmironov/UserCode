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


void v2SummaryPlots_pt()
{
  gROOT->Macro("/Users/eusmartass/Software/utilities/setStyle.C+");
  gStyle->SetOptFit(0);

  const char* signal[5]   = {"","NSig","NPr","NNp","NBkg"};
  const char* legend[5]   = {"","J/#psi","Prompt J/#psi","Non-prompt J/#psi","Background"};
  int choseSignal         = 1; // 1:inclusive 2:prompt 3:non-prompt
  const char* chosenSignal= signal[choseSignal];
  const char* outputName[1]  = {"nominalFit"};
 
  // options
  bool bSavePlots      = true; 
  bool bAddBkg         = false;
  bool bAddVtx         = false;
  bool bAddNoFlat      = false;
  bool bAddTrigger     = false;

  const char* eventPlane[2] = {"","EP: etHFp & etHFm"};
  double rapIntegrated[2]   = {0.0, 2.4}; 
  
  // inclusive numbers
  // pt bins: 
  const int nptbins                      = 3;
  double ptbins_center[nptbins]          = {7.3, 9.0, 13.4}; // <pt> values in the AN
  double ptbins_center_err[nptbins]      = {0., 0., 0.};
  double ptbins_center_errGhost[nptbins] = {0., 0., 0.};

  // signal
  double jpsi_ptbins[nptbins]            = {0.075, 0.1, 0.043};
  double jpsi_ptbins_err[nptbins]        = {0.022, 0.022, 0.020};
  double jpsi_ptbins_errGhost[nptbins]   = {0.,     0.,     0., };

  // bkg
  double bkg_ptbins[nptbins]             = {0.1485, 0.2101, 0.0298};
  double bkg_ptbins_err[nptbins]         = {0.0277, 0.0350, 0.0604};

  // vtx10cm signal
  double jpsiVtx10_ptbins[nptbins]       = {0.0567, 0.0912, 0.0266};
  double jpsiVtx10_ptbins_err[nptbins]   = {0.0258, 0.0252, 0.0297};

  // noFlattening
  double jpsiNoFlat_ptbins[nptbins]          = {0.0313, 0.1094, 0.0432};
  double jpsiNoFlat_ptbins_err[nptbins]      = {0.0247, 0.0235, 0.0330};

  // triggers
  double jpsiL1NHitTrig_ptbins[nptbins]      = {0.0272, 0.1059, 0.0405};
  double jpsiL1NHitTrig_ptbins_err[nptbins]   = {0.0254, 0.0238, 0.0330};

  double jpsiL2pt3Trig_ptbins[nptbins]       = {-0.0214, 0.0925, 0.0514};
  double jpsiL2pt3Trig_ptbins_err[nptbins]   = { 0.0346, 0.0281, 0.0357};

  double jpsiL3ptOpenTrig_ptbins[nptbins]    = {-0.0281, 0.0966, 0.0456};
  double jpsiL3ptOpenTrig_ptbins_err[nptbins]= { 0.0647, 0.0598, 0.0579};

  //autocorrelations
  double jpsiAutoCorr_ptbins[nptbins]        = {0.0638, 0.1051, 0.1462};
  double jpsiAutoCorr_ptbins_err[nptbins]    = {0.0242, 0.0232, 0.0320};


  TGraphErrors *pg_jpsi_ptbins      = new TGraphErrors(nptbins,ptbins_center,jpsi_ptbins,ptbins_center_errGhost,jpsi_ptbins_err);
  TGraphErrors *pg_jpsi_ptbinsGhost = new TGraphErrors(nptbins,ptbins_center,jpsi_ptbins,ptbins_center_err,jpsi_ptbins_errGhost);
 
  TGraphErrors *pg_jpsiAutoCorr_ptbins     = new TGraphErrors(nptbins,ptbins_center,jpsiAutoCorr_ptbins,ptbins_center_err,jpsiAutoCorr_ptbins_err);

  TGraphErrors *pg_bkg_ptbins              = new TGraphErrors(nptbins,ptbins_center,bkg_ptbins,ptbins_center_errGhost,bkg_ptbins_err);
 
  TGraphErrors *pg_jpsiVtx10_ptbins        = new TGraphErrors(nptbins,ptbins_center,jpsiVtx10_ptbins,ptbins_center_errGhost,jpsiVtx10_ptbins_err);
  TGraphErrors *pg_jpsiNoFlat_ptbins       = new TGraphErrors(nptbins,ptbins_center,jpsiNoFlat_ptbins,ptbins_center_errGhost,jpsiNoFlat_ptbins_err);
  TGraphErrors *pg_jpsiL1NHitTrig_ptbins   = new TGraphErrors(nptbins,ptbins_center,jpsiL1NHitTrig_ptbins,ptbins_center_errGhost,jpsiL1NHitTrig_ptbins_err);
  TGraphErrors *pg_jpsiL2pt3Trig_ptbins    = new TGraphErrors(nptbins,ptbins_center,jpsiL2pt3Trig_ptbins,ptbins_center_errGhost,jpsiL2pt3Trig_ptbins_err);
  TGraphErrors *pg_jpsiL3ptOpenTrig_ptbins = new TGraphErrors(nptbins,ptbins_center,jpsiL3ptOpenTrig_ptbins,ptbins_center_errGhost,jpsiL3ptOpenTrig_ptbins_err);
  //_____________________________________________________________________________
  // pt dependence
  TCanvas *pcPt = new TCanvas("pcPt","pcPt");
  pcPt->cd();
  TH1F *pcPadPT = new TH1F("pcPadPT",";p_{T} (GeV/c);v_{2};",100,0,40);
  pcPadPT->GetXaxis()->SetLabelSize(20);
  pcPadPT->GetXaxis()->SetLabelFont(43);
  pcPadPT->GetXaxis()->SetTitleSize(27);
  pcPadPT->GetXaxis()->SetTitleFont(43);
  pcPadPT->GetXaxis()->SetTitleOffset(1.2);
  pcPadPT->GetXaxis()->CenterTitle();

  pcPadPT->GetYaxis()->SetLabelSize(20);
  pcPadPT->GetYaxis()->SetLabelFont(43);
  pcPadPT->GetYaxis()->SetTitleSize(32);
  pcPadPT->GetYaxis()->SetTitleFont(43);
  pcPadPT->GetYaxis()->SetTitleOffset(1.1);
  pcPadPT->GetYaxis()->CenterTitle();

  pcPadPT->SetMaximum(0.18);
  pcPadPT->SetMinimum(-0.05);
  if(bAddBkg)
    {
      pcPadPT->SetMaximum(0.4);
    }
  pcPadPT->Draw();

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
  //________________________________________ 
  pg_jpsi_ptbins->SetMarkerStyle(20);
  pg_jpsi_ptbins->SetMarkerSize(1.8);
  pg_jpsi_ptbins->SetMarkerColor(602);
  pg_jpsi_ptbins->Draw("pz");

  pg_jpsi_ptbinsGhost->SetMarkerStyle(20);
  pg_jpsi_ptbinsGhost->SetMarkerSize(1.8);
  pg_jpsi_ptbinsGhost->SetMarkerColor(602);
  pg_jpsi_ptbinsGhost->Draw("P");

  if(bAddBkg)
    {
      pg_bkg_ptbins->SetMarkerStyle(24);
      pg_bkg_ptbins->SetMarkerSize(1.8);
      pg_bkg_ptbins->SetMarkerColor(kBlue);
      pg_bkg_ptbins->Draw("p");
      TLegend *legPt = new TLegend(0.7,0.5,0.9,0.6);
      legPt->SetFillColor(0);
      legPt->SetBorderSize(0);
      legPt->SetTextSize(0.03);
      legPt->AddEntry(pg_jpsi_ptbins,"Signal","P");
      legPt->AddEntry(pg_bkg_ptbins,"Bkg","P");
      legPt->Draw("same");
    }

  if(bAddVtx)
    {
      pg_jpsiVtx10_ptbins->SetMarkerStyle(24);
      pg_jpsiVtx10_ptbins->SetMarkerSize(1.8);
      pg_jpsiVtx10_ptbins->SetMarkerColor(kBlue);
      pg_jpsiVtx10_ptbins->Draw("p");
      
      TLegend *legVtx = new TLegend(0.55,0.15,0.875,0.297);
      legVtx->SetFillColor(0);
      legVtx->SetBorderSize(0);
      legVtx->SetTextSize(0.03);
      legVtx->AddEntry(pg_jpsi_ptbins,"Default: |z_{vtx}|<25 cm","P");
      legVtx->AddEntry(pg_jpsiVtx10_ptbins,"|z_{vtx}|<10 cm","P");
      legVtx->Draw("same");
    }

  
  if(bAddNoFlat)
    {
      pg_jpsiNoFlat_ptbins->SetMarkerStyle(24);
      pg_jpsiNoFlat_ptbins->SetMarkerSize(1.8);
      pg_jpsiNoFlat_ptbins->SetMarkerColor(kBlue);
      //  pg_jpsiNoFlat_ptbins->Draw("p");

      pg_jpsiAutoCorr_ptbins->SetMarkerStyle(27);
      pg_jpsiAutoCorr_ptbins->SetMarkerSize(1.8);
      pg_jpsiAutoCorr_ptbins->SetMarkerColor(kBlue);
      pg_jpsiAutoCorr_ptbins->Draw("p");
      
      TLegend *legEp = new TLegend(0.55,0.15,0.875,0.297);
      legEp->SetFillColor(0);
      legEp->SetBorderSize(0);
      legEp->SetTextSize(0.03);
      // legEp->AddEntry(pg_jpsi_ptbins,"Default: w/ flattening","P");
      // legEp->AddEntry(pg_jpsiNoFlat_ptbins,"No flattening","P");
      
      legEp->AddEntry(pg_jpsi_ptbins,"+- && -+","P");
      legEp->AddEntry(pg_jpsiAutoCorr_ptbins,"++ && --","P");

      legEp->Draw("same");

    }

  if(bAddTrigger)
    {
      pg_jpsiL1NHitTrig_ptbins->SetMarkerStyle(24);
      pg_jpsiL2pt3Trig_ptbins->SetMarkerStyle(27);
      pg_jpsiL3ptOpenTrig_ptbins->SetMarkerStyle(28);
      
      pg_jpsiL1NHitTrig_ptbins->SetMarkerColor(kBlue);
      pg_jpsiL2pt3Trig_ptbins->SetMarkerColor(kBlue);
      pg_jpsiL3ptOpenTrig_ptbins->SetMarkerColor(kBlue);
      
      pg_jpsiL1NHitTrig_ptbins->SetMarkerSize(1.8);
      pg_jpsiL2pt3Trig_ptbins->SetMarkerSize(1.8);
      pg_jpsiL3ptOpenTrig_ptbins->SetMarkerSize(1.8);

      pg_jpsiL1NHitTrig_ptbins->Draw("P");
      pg_jpsiL2pt3Trig_ptbins->Draw("p");
      pg_jpsiL3ptOpenTrig_ptbins->Draw("p");
      
      TLegend *legTrig = new TLegend(0.4,0.15,0.57,0.29);
      legTrig->SetFillColor(0);
      legTrig->SetBorderSize(0);
      legTrig->SetTextSize(0.03);
      legTrig->AddEntry(pg_jpsi_ptbins,"Default: all triggers","P");
      legTrig->AddEntry(pg_jpsiL1NHitTrig_ptbins,"HLT_HIL1DoubleMu0_HighQ","P");
      legTrig->AddEntry(pg_jpsiL2pt3Trig_ptbins,"HLT_HIL2DoubleMu3","P");
      legTrig->AddEntry(pg_jpsiL3ptOpenTrig_ptbins,"HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy","P");
      legTrig->Draw("same");
    }

  TLatex *lt1 = new TLatex();
  lt1->SetNDC();
  lt1->SetTextSize(0.04);
  lt1->DrawLatex(0.18,0.89,Form("Inclusive %s",legend[choseSignal]));  // what signal is
  lt1->SetTextSize(0.038);
  lt1->DrawLatex(0.18,0.83,Form("|y| < %.1f",rapIntegrated[1]));       // rapidity
  lt1->DrawLatex(0.18,0.77,Form("Cent. 10 - 60 %%")); 
  //  lt1->DrawLatex(0.18,0.71,Form("%s",eventPlane[1]));
  
  if(bSavePlots)
    {
      pcPt->SaveAs(Form("figs/png/%s_%s_v2_ptBins.png",chosenSignal,outputName[0]));
      pcPt->SaveAs(Form("figs/pdf/%s_%s_v2_ptBins.pdf",chosenSignal,outputName[0]));
    }

}



