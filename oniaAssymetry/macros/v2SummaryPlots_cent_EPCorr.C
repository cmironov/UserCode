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
//  gROOT->Macro("/Users/eusmartass/Software/utilities/setStyle.C+");
  gROOT->Macro("./rootlogon.C");
  gStyle->SetOptFit(0);

  const char* signal[5]   = {"","NSig","NPr","Pr","NBkg"};
  const char* legend[5]   = {"","Inclusive J/#psi","Non-prompt J/#psi", "Prompt J/#psi","Background"};
  int choseSignal         = 1; // 1:inclusive 2:prompt 3:non-prompt
  const char* chosenSignal= signal[choseSignal];
//  const char* outputName[1]  = {"nominal_bit1_cowboy_sailor"};
  const char* outputName[1]  = {"nominal_bkg"};

  // options
  bool bSavePlots      = true; 
  bool bAddBkg         = true;
  bool bAddVtx         = false;
  bool bAddNoFlat      = false;
  bool bAddTrigger     = false; 
  bool bAddAutoCor     = false;

  const char* eventPlane[2] = {"","EP: etHFp & etHFm"};
  double rapIntegrated[2]   = {0.0, 2.4}; 
  double ptIntegrated[2]    = {6.5, 40};  
//  double ptIntegrated[3]    = {3, 6.5, 40};  
  // inclusive numbers
  // centrlaity bins: 0-5, 5-10, 10-30, 30-60
  const int ncentbins1                    = 4;
  double ncoll1[ncentbins1]                = {381.3,   329.4,  224.3,   89.9};
  // centrlaity bins: 0-10, 10-20, 20-30, 30-60
  const int ncentbins                    = 4;
  double ncoll[ncentbins]                = {355.4, 261.4178, 187.1470, 89.9};
  double ncoll_err[ncentbins]            = { 0.,     0.,       0.,      0.};

  char dirname[512];
  sprintf(dirname,"./a1_corrV2.root");
  TFile *f1 = new TFile(dirname);
  if(!f1->IsOpen()) { cout << "cannot open a1_corrV2.root" << endl; return;}
  
  char histname[200];
  //inclusive
  sprintf(histname,"nominal_NSig_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsi_centbins = (TGraphErrors*)f1->Get(histname);
  if (!pg_jpsi_centbins) { cout << "cannot load nominal_NSig case." << endl; return;}
  // prompt
  sprintf(histname,"nominal_NPr_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_pr_jpsi_centbins = (TGraphErrors*)f1->Get(histname);
  // non-prompt
  sprintf(histname,"nominal_NNp_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_npr_jpsi_centbins = (TGraphErrors*)f1->Get(histname);

  // systematic studies
  sprintf(histname,"nominal_NBkg_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_bkg_centbins = (TGraphErrors*)f1->Get(histname);

  sprintf(histname,"zVtxLT10_NSig_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiVtx10_centbins = (TGraphErrors*)f1->Get(histname);
  sprintf(histname,"noFlat_NSig_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiNoFlat_centbins = (TGraphErrors*)f1->Get(histname);
  sprintf(histname,"autoCorr_NSig_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiAutoCorr_centbins = (TGraphErrors*)f1->Get(histname);
 
  sprintf(histname,"bit1_NSig_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_centbins = (TGraphErrors*)f1->Get(histname);
  sprintf(histname,"cowboy_NSig_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_onlyCow_centbins = (TGraphErrors*)f1->Get(histname);
  sprintf(histname,"sailor_NSig_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL1NHitTrig_noCow_centbins = (TGraphErrors*)f1->Get(histname);

/*  sprintf(histname,"_NSig_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL2pt3Trig_centbins = (TGraphErrors*)f1->Get(histname);
  sprintf(histname,"_NSig_rap%.1f-%.1f_pT%.1f-%.1f",rapIntegrated[0],rapIntegrated[1],ptIntegrated[0],ptIntegrated[1]);
  TGraphErrors *pg_jpsiL3ptOpenTrig_centbins = (TGraphErrors*)f1->Get(histname);*/

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

  phPadCentrality->SetMaximum(0.25);
  phPadCentrality->SetMinimum(-0.1);
  //phPadCentrality->SetMaximum(0.18);
  //phPadCentrality->SetMinimum(-0.05);
  if(bAddBkg)
    {
      phPadCentrality->SetMaximum(0.4);
    }
  phPadCentrality->Draw();
  
  switch(choseSignal){
  case 1:
    pg_jpsi_centbins->SetMarkerStyle(20);
    pg_jpsi_centbins->SetMarkerSize(1.8);
    pg_jpsi_centbins->SetMarkerColor(602);
    pg_jpsi_centbins->SetLineColor(602);
    pg_jpsi_centbins->Draw("[P]");
    break;
  case 2:
    pg_pr_jpsi_centbins->SetMarkerStyle(20);
    pg_pr_jpsi_centbins->SetMarkerSize(1.8);
    pg_pr_jpsi_centbins->SetMarkerColor(602);
    pg_pr_jpsi_centbins->SetLineColor(602);
    pg_pr_jpsi_centbins->Draw("[P]");
    break;
  case 3:
    pg_npr_jpsi_centbins->SetMarkerStyle(20);
    pg_npr_jpsi_centbins->SetMarkerSize(1.8);
    pg_npr_jpsi_centbins->SetMarkerColor(602);
    pg_npr_jpsi_centbins->SetLineColor(602);
    pg_npr_jpsi_centbins->Draw("[P]");
    break;
  default:
    cout<<"Pick a valid signal!!"<<endl;
  }
  

  if(bAddBkg)
    {
      pg_bkg_centbins->SetMarkerStyle(24);
      pg_bkg_centbins->SetMarkerSize(1.8);
      pg_bkg_centbins->SetMarkerColor(kBlue);
      pg_bkg_centbins->SetLineColor(kBlue);
      pg_bkg_centbins->Draw("[P]");
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
      pg_jpsiVtx10_centbins->SetLineColor(kBlue);
      pg_jpsiVtx10_centbins->Draw("[P]");
      
      TLegend *legVtx = new TLegend(0.19,0.625,0.50,0.736);
      //TLegend *legVtx = new TLegend(0.19,0.61,0.52,0.72);
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
//      pg_jpsiL2pt3Trig_centbins->SetMarkerStyle(27);
//      pg_jpsiL3ptOpenTrig_centbins->SetMarkerStyle(28);
      
      pg_jpsiL1NHitTrig_centbins->SetMarkerColor(kBlue);
      pg_jpsiL1NHitTrig_centbins->SetLineColor(kBlue);
//      pg_jpsiL2pt3Trig_centbins->SetMarkerColor(kBlue);
//      pg_jpsiL3ptOpenTrig_centbins->SetMarkerColor(kBlue);
      
      pg_jpsiL1NHitTrig_centbins->SetMarkerSize(1.8);
//      pg_jpsiL2pt3Trig_centbins->SetMarkerSize(1.8);
//      pg_jpsiL3ptOpenTrig_centbins->SetMarkerSize(1.8);

      pg_jpsiL1NHitTrig_onlyCow_centbins->SetMarkerStyle(24);
      pg_jpsiL1NHitTrig_noCow_centbins->SetMarkerStyle(24);
	
      pg_jpsiL1NHitTrig_onlyCow_centbins->SetMarkerColor(kGreen+2);
      pg_jpsiL1NHitTrig_noCow_centbins->SetMarkerColor(kMagenta);
      pg_jpsiL1NHitTrig_onlyCow_centbins->SetLineColor(kGreen+2);
      pg_jpsiL1NHitTrig_noCow_centbins->SetLineColor(kMagenta);

      pg_jpsiL1NHitTrig_onlyCow_centbins->SetMarkerSize(1.8);
      pg_jpsiL1NHitTrig_noCow_centbins->SetMarkerSize(1.8);

      pg_jpsiL1NHitTrig_centbins->Draw("[P]");
      //pg_jpsiL1NHitTrig_centbins->Draw("P");
      //   pg_jpsiL2pt3Trig_centbins->Draw("p");
      // pg_jpsiL3ptOpenTrig_centbins->Draw("p");

     
      pg_jpsiL1NHitTrig_onlyCow_centbins->Draw("[P]");
      pg_jpsiL1NHitTrig_noCow_centbins->Draw("[P]");
      
      TLegend *legTrig = new TLegend(0.2,0.15,0.57,0.31);
      legTrig->SetFillColor(0);
      legTrig->SetBorderSize(0);
      legTrig->SetTextSize(0.03);
      legTrig->AddEntry(pg_jpsi_centbins,"Default: all triggers","P");
      legTrig->AddEntry(pg_jpsiL1NHitTrig_centbins,"HLT_HIL1DoubleMu0_HighQ","P");
      //   legTrig->AddEntry(pg_jpsiL2pt3Trig_centbins,"HLT_HIL2DoubleMu3","P");
      //  legTrig->AddEntry(pg_jpsiL3ptOpenTrig_centbins,"HLT_HIL3DoubleMuOpen_Mgt2_OS_NoCowboy","P");
      legTrig->AddEntry(pg_jpsiL1NHitTrig_noCow_centbins,"HLT_HIL1DoubleMu0_HighQ+Sailor","P");
      legTrig->AddEntry(pg_jpsiL1NHitTrig_onlyCow_centbins,"HLT_HIL1DoubleMu0_HighQ+Cowboy","P");
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
  lt1->DrawLatex(0.18,0.89,Form("%s",legend[choseSignal]));  // what signal is
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
      lt1->DrawLatex(0.25, 0.34, Form("30-60%%"));
      lt1->DrawLatex(0.45, 0.34, Form("20-30%%"));
      lt1->DrawLatex(0.6, 0.34, Form("10-20%%"));
      lt1->DrawLatex(0.8, 0.34, Form("0-10%%"));
    }

  if(bSavePlots)
    {
      gSystem->mkdir(Form("./figs/png/"),kTRUE);
      gSystem->mkdir(Form("./figs/pdf/"),kTRUE);
      pcCentrality->SaveAs(Form("figs/png/%s_%s_v2_centBins.png",chosenSignal,outputName[0]));
      pcCentrality->SaveAs(Form("figs/pdf/%s_%s_v2_centBins.pdf",chosenSignal,outputName[0]));
    }
 
}
